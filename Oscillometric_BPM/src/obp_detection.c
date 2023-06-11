#include <assert.h>
#include "obp_detection.h"
#include "obp_const.h"
#include "vector/c_vector.h"
#include "log_wrapper.h"

struct sOBPDetection {
    Vector pressure_data;
    Vector oscillation_data;
    Vector max_Amp;
    Vector max_time;
    Vector min_Amp;
    Vector min_time;
    Vector heart_rate;
    Vector omwe_data;
    Vector omwe_time; 

    double ratio_SBP;
    double ratio_DBP;
    double max_valid_HR;
    double min_valid_HR;
    double prominence;
    int min_data_size;
    int min_peak_time;
    double sampling_rate;
    int min_number_peaks;
    double cutoff_Hyst;

    double result_MAP;
    double result_SBP;
    double result_DBP;
    bool engouh_data;
};

static uint8_t ref_counter = 0;
static OBPDetection_t sInstance = {
    .ratio_SBP = 0.57f,
    .ratio_DBP = 0.70f,
    .min_number_peaks = 10,
    .max_valid_HR = kMAX_HEART_RATE,
    .min_valid_HR = kMIN_HEART_RATE,
    .sampling_rate = kSAMPLE_RATE,
    .prominence = kPROMINENCE,
    .min_data_size = 1200,
    .min_peak_time = 300,
    .cutoff_Hyst = kCUTOFF_HYST,
    .engouh_data = false
};

/* PRIVATE API*/

static double lerp(float a, float b, float t)
{
    return a + t * (b - a);
}

static double get_ratio(double lower_bound, double upper_bound, double value)
{
    return ((value - lower_bound) / (upper_bound -lower_bound));
}

static double get_average(Vector *vec)
{
    double average = 0.0f;
    if (vector_is_empty(vec))
    {
        return average;
    }

    VECTOR_FOR_EACH(vec, current)
    {
        average += ITERATOR_GET_AS(double, &current);
    }
    return (average / vec->size);
}

static double get_pressure_at(OBPDetection_t *inst, int time)
{
    double average;
    int hrSampleHalf = (inst->sampling_rate * (int)get_average(&inst->heart_rate)) / 120;

    if (inst->pressure_data.size > (time + hrSampleHalf))
    {
        log_debug("first pressure: %d, last pressure: %d, size: %d", time - hrSampleHalf, time + hrSampleHalf, inst->pressure_data.size);
        double first_pressure = *((double*)vector_get(&inst->pressure_data, time - hrSampleHalf));
        double last_pressure = *((double*)vector_get(&inst->pressure_data, time + hrSampleHalf));
        average = (first_pressure + last_pressure)/2;
    }
    else
    {
        log_info("Trying to get pressure at time %d with hrSampleHalf %d and data size %ld\n", time, 
                                                                                            hrSampleHalf, 
                                                                                            inst->pressure_data.size);
        average = *((double*)vector_get(&inst->pressure_data, time));
    }

    return average;
}

static bool is_heart_rate_valid(double heart_rate)
{
    if (heart_rate > kMAX_HEART_RATE || heart_rate < kMIN_HEART_RATE)
    {
        return false;
    }

    return true;
}
static bool is_valid_maximal(OBPDetection_t *inst)
{
    bool is_valid = false;

    assert(inst->oscillation_data.size >= 2);

    size_t second_from_last = inst->oscillation_data.size - 2;
    const double test_value = VECTOR_GET_AS(double, &inst->oscillation_data, second_from_last);
    double test_sample_number = inst->oscillation_data.size - 1;
    
    if (vector_is_empty(&inst->max_time))
    {
        vector_push_back(&inst->max_time, (void*)&test_sample_number);
        vector_push_back(&inst->max_Amp, (void*)&test_value);
    } else
    {
        assert(inst->max_time.size != 0);
        assert(inst->max_Amp.size != 0);
        double last_maxtime = *((double*)vector_back(&inst->max_time));
        if ((test_sample_number - last_maxtime) < inst->min_peak_time)
        {
            double last_max_amp = *((double*)vector_back(&inst->max_Amp));
            if (last_max_amp < test_value)
            {
                log_debug("Update last value");
                vector_pop_back(&inst->max_Amp);
                vector_pop_back(&inst->max_time);
                vector_push_back(&inst->max_Amp, (void*)&test_value);
                vector_push_back(&inst->max_time, (void*)&test_sample_number);
            } else
            {
                /*Ignore this max*/
            }
        } else
        {
            log_debug("Push new value");

            vector_push_back(&inst->max_Amp, (void*)&test_value);
            vector_push_back(&inst->max_time, (void*)&test_sample_number);
        }

        if (inst->max_time.size > 1)
        {
            size_t _index = inst->max_time.size - 2;
            double second_from_last = VECTOR_GET_AS(double, &inst->max_time, _index); 
            double last_maxtime = *((double*)vector_back(&inst->max_time));
            double new_heart_rate = (60.0 * inst->sampling_rate) / (double) (last_maxtime - second_from_last);
            log_debug("Heart rate: %0.4f", new_heart_rate);
            if (is_heart_rate_valid(new_heart_rate))
            {
                vector_push_back(&inst->heart_rate, (void*)&new_heart_rate);
                is_valid = true;
            } else
            {
                log_debug("Clear pulse");
                /*Into here is a invalid pulse entered*/
                vector_clear(&inst->max_Amp);
                vector_clear(&inst->max_time);
                vector_clear(&inst->min_Amp);
                vector_clear(&inst->min_time);
                vector_clear(&inst->heart_rate);
                vector_push_back(&inst->max_time, (void*)&test_sample_number);
                vector_push_back(&inst->max_Amp, (void*)&test_value);
                is_valid = false;
            }
        }
    }

    return is_valid;
}

static bool check_maximal(OBPDetection_t *inst)
{
    bool is_valid = false;
    if (inst->oscillation_data.size <= inst->min_data_size)
    {
        return false;
    }
    size_t start_index = inst->oscillation_data.size - 2;
    double data = VECTOR_GET_AS(double, &inst->oscillation_data, start_index);
    // log_debug("data: %0.4f", data);
    if (data <= inst->prominence)
    {
        return false;
    }

    size_t local_domain = inst->oscillation_data.size - 3;
    size_t end_index = inst->oscillation_data.size;
    size_t local_max_idx = local_domain;
    double local_max_value = VECTOR_GET_AS(double, &inst->oscillation_data, local_max_idx);
    for (size_t i = local_domain; i < end_index; i++)
    {
        double new_value = VECTOR_GET_AS(double, &inst->oscillation_data, i);
        if (new_value > local_max_value)
        {
            local_max_value = new_value;
            local_max_idx = i;
        }
    }
    if ((end_index - local_max_idx) == 2)
    {
        is_valid = is_valid_maximal(inst);
    }

    return is_valid;
}

static void find_minimal(OBPDetection_t *inst)
{
    if (inst->max_Amp.size >= 2)
    {
        size_t first_time_index = inst->max_time.size - 2;
        double first_max = VECTOR_GET_AS(double, &inst->max_time, first_time_index);
        double last_max = *((double*)vector_back(&inst->max_time)); 
        double min_pos = first_max;
        double min_value = VECTOR_GET_AS(double, &inst->oscillation_data, first_max);
        for (size_t i = first_max; i <= last_max; i++)
        {
            double new_value = VECTOR_GET_AS(double, &inst->oscillation_data, i);
            if (new_value < min_value)
            {
                min_value = new_value;
                min_pos = i;
            }
        }
        double distance = min_pos - first_max;

        if (inst->min_time.size == (inst->max_time.size - 1))
        {
            vector_pop_back(&inst->min_Amp);
            vector_push_back(&inst->min_Amp, &min_value);
            vector_pop_back(&inst->min_time);
            distance += first_max;
            vector_push_back(&inst->min_time, &distance);
        } else
        {
            distance += first_max;
            vector_push_back(&inst->min_Amp, &min_value);
            vector_push_back(&inst->min_time, &distance);
        }
    }
}

static bool is_enough_data(OBPDetection_t *inst)
{
    bool is_enough = false;
    log_debug("max AMP size %d", inst->max_Amp.size);
    if (inst->max_Amp.size <= inst->min_number_peaks)
    {
        return false;
    }
    log_debug("enter...");
    double max_el = 0.0f;
    VECTOR_FOR_EACH(&inst->max_Amp, iter)
    {
        double local_value = ITERATOR_GET_AS(double, &iter);
        if (max_el < local_value)
        {
            max_el = local_value;
        }
    }
    double last_value = *((double*)vector_back(&inst->max_Amp));
    double second_from_last = VECTOR_GET_AS(double, &inst->max_Amp, inst->max_Amp.size - 3);
    double first_from_last = VECTOR_GET_AS(double, &inst->max_Amp, inst->max_Amp.size - 2);
    log_debug("max_el %0.4f, last value: %0.4f, second_from last: %0.4f, first from last: %0.4f", max_el, last_value, second_from_last, first_from_last);
    if (max_el > 1.5 && ((last_value < second_from_last && last_value < first_from_last) ||
                            (last_value < 2 * inst->prominence)))
    {
        double cut_off = max_el * (inst->ratio_DBP - inst->cutoff_Hyst);

        if ((second_from_last < cut_off) && (first_from_last < cut_off) && (last_value < cut_off))
        {
            is_enough = true;
        }
    }

    return is_enough;
}

static void find_OWME(OBPDetection_t *inst)
{
    log_debug("find owme");
    size_t j = 0;

    for (size_t i = 0; i < inst->min_time.size - 1; ++i)
    {
        double amp_min1 = *((double*)vector_get(&inst->min_Amp, j));
        double amp_max1 = *((double*)vector_get(&inst->max_Amp, j));
        double time_max1 = *((double*)vector_get(&inst->max_time, j));

        double amp_min2 =*((double*)vector_get(&inst->min_Amp, j + 1));
        double amp_max2 = *((double*)vector_get(&inst->max_Amp, j + 1));
        double time_max2 =  *((double*)vector_get(&inst->max_time, j + 1));
        j++;
        double time_min1 = VECTOR_GET_AS(double, &inst->min_time, i);
        double time_min2 = VECTOR_GET_AS(double, &inst->min_time, i + 1); 
        assert(time_min1 > time_max1);
        assert(time_min2 > time_max2);

        double lerp_max = lerp(amp_max1, amp_max2, get_ratio(time_max1, time_max2, time_min1));
        double lerp_min = lerp(amp_min1, amp_min2, get_ratio(time_min1, time_min2, time_max2));

        lerp_max -= amp_min1;
        amp_max2 -= lerp_min;
        vector_push_back(&inst->omwe_data, &lerp_max);
        vector_push_back(&inst->omwe_time, &time_min1);
        vector_push_back(&inst->omwe_data, &amp_max2);
        vector_push_back(&inst->omwe_time, &time_max2);
    }
}

static void find_MAP(OBPDetection_t *inst)
{
    log_debug("Enter...");
    int imax_OWME = 0;
    double max_OWME = 0.0f;
    int counter = 0;
    VECTOR_FOR_EACH(&inst->omwe_data, current)
    {
        double current_value = ITERATOR_GET_AS(double, &current);
        if (current_value > max_OWME)
        {
            max_OWME = current_value;
            imax_OWME = counter;
        }

        counter ++;
    }

    int max_time = *((double*)vector_get(&inst->omwe_time, imax_OWME));

    inst->result_MAP = get_pressure_at(inst, max_time);
    double sbp_search = inst->ratio_SBP * max_OWME;
    double ubSBP = 0.0f;
    double lbSBP = 0.0f;
    int lbSTime = 0;
    int ubSTime = 0;

    for (int i = 0; i < imax_OWME; i++)
    {
        double omwerSBP = *((double*)vector_get(&inst->omwe_data, i));
        // log_debug("value: %0.4f sbp_search: %0.4f", omwerSBP, sbp_search);
        if (omwerSBP > sbp_search)
        {
            ubSBP = omwerSBP;
            log_debug("ubSBP: %0.4f, i: %d", ubSBP, i);
            ubSTime = *((double*)vector_get(&inst->omwe_time, i));
            log_debug("ubSTime: %d", ubSTime);
            i -= 1;
            lbSBP = *((double*)vector_get(&inst->omwe_data, i));
            log_debug("lbSBP: %0.4f", lbSBP);
            lbSTime = *((double*)vector_get(&inst->omwe_time, i));
            log_debug("lbSTime: %d", lbSTime);
            break;
        }
    }
    double lerp_SBP_time = lerp(lbSTime, ubSTime, get_ratio(lbSBP, ubSBP, sbp_search));
    log_debug("lert_SBP: %0.4f", lerp_SBP_time);
    inst->result_SBP = get_pressure_at(inst, lerp_SBP_time);

    double dpb_search = inst->ratio_DBP * max_OWME;
    double ubDBP = 0;
    double lbDBP = 0;
    int lbDTime = 0;
    int ubDTime = 0;

    for (int i = imax_OWME; i < inst->omwe_data.size; i++)
    {
        double omwerDBP = *((double*)vector_get(&inst->omwe_data, i));
        if (omwerDBP < dpb_search)
        {
            lbDBP = omwerDBP;
            lbDTime =  *((double*)vector_get(&inst->omwe_time, i));
            i--;
            ubDBP = *((double*)vector_get(&inst->omwe_data, i));
            ubDTime = *((double*)vector_get(&inst->omwe_time, i));
            break;
        }
    }

    if (lbDBP != 0)
    {
        int lerpDBPTime = (int)lerp(ubDTime, lbDTime, 1.0 - get_ratio(lbDBP, ubDBP, dpb_search));
        inst->result_DBP = get_pressure_at(inst, lerpDBPTime);
    } else 
    {
        log_warn("Couldn't find DBP\n");
    }
}

/* PUBLIC API */
OBPDetection_t *OBP_get_instance()
{
    ref_counter ++;
    vector_setup(&sInstance.pressure_data, 1000, sizeof(double));
    vector_setup(&sInstance.oscillation_data, 1000, sizeof(double));
    vector_setup(&sInstance.max_Amp, 1000, sizeof(double));
    vector_setup(&sInstance.max_time, 1000, sizeof(double));
    vector_setup(&sInstance.min_Amp, 1000, sizeof(double));
    vector_setup(&sInstance.min_time, 1000, sizeof(double));
    vector_setup(&sInstance.heart_rate, 1000, sizeof(double));
    vector_setup(&sInstance.omwe_data, 1000, sizeof(double));
    vector_setup(&sInstance.omwe_time, 1000, sizeof(double));
    return &sInstance;
}

void OBP_set_SBP_ratio(OBPDetection_t *instance, double val)
{
    instance->ratio_SBP = val;
}

void OBP_set_DBP_ratio(OBPDetection_t *instance, double val)
{
    instance->ratio_DBP = val;
}

void OBP_set_min_number_of_peaks(OBPDetection_t *instance, int val)
{
    instance->min_number_peaks = val;
}

void OBP_process_sample(OBPDetection_t *instance, double pressure, double oscillation)
{
    bool new_max = false;
    vector_push_back(&instance->pressure_data, (void*)&pressure);
    vector_push_back(&instance->oscillation_data, (void*)&oscillation);
    if (check_maximal(instance))
    {
        find_minimal(instance);
        if (is_enough_data(instance))
        {
            log_debug("Enter...");
            find_OWME(instance);
            find_MAP(instance);
            instance->engouh_data = true;
        }
        new_max = true;
    }
}

double OBP_get_current_heart_rate(const OBPDetection_t *instance)
{
    return *((double*)vector_back(&instance->heart_rate));
}

double OBP_get_average_heart_rate(const OBPDetection_t *instance)
{
    return get_average(&instance->heart_rate);
}

double OBP_get_MAP(const OBPDetection_t *instance)
{
    return instance->result_MAP;
}

double OBP_get_SBP(const OBPDetection_t *instance)
{
    return instance->result_SBP;
}

double OBP_get_DBP(const OBPDetection_t *instance)
{
    return instance->result_DBP;
}

bool OBP_is_enough_data(const OBPDetection_t *inst)
{
    return inst->engouh_data;
}

void OBP_reset_instance(OBPDetection_t *instance)
{
    /*Perform reset all measurement variable before next step*/
    instance->result_MAP = 0.0f;
    instance->result_DBP = 0.0f;
    instance->result_SBP = 0.0f;
}
