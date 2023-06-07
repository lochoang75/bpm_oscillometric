#include "obp_detection.h"
#include "obp_const.h"
#include "c_vector.h"

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
};

static uint8_t ref_counter = 0;
static OBPDetection_t sInstance;

/* PRIVATE API*/

static double lerp(float a, float b, float t)
{
    return a + t * (b - a);
}

static double get_ratio(double lower_bound, double upper_bound, double value)
{
    return ((value - lower_bound) / (upper_bound -lower_bound));
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

    if (inst->oscillation_data.size < 2)
    {
        return false;
    }

    size_t second_from_last = inst->oscillation_data.size - 2;
    const double test_value = VECTOR_GET_AS(double, &inst->oscillation_data, second_from_last);
    size_t test_sample_number =inst->oscillation_data.size - 1;
    
    if (vector_is_empty(&inst->max_time))
    {
        vector_push_back(&inst->max_time, &test_sample_number);
        vector_push_back(&inst->max_Amp, &test_value);
    } else
    {
        if (vector_is_empty(&inst->max_time) ||
            vector_is_empty(&inst->max_Amp))
        {
            return false;
        }

        double last_maxtime = *((double*)vector_back(&inst->max_time));
        if ((test_sample_number - last_maxtime) < kMIN_PEAK_TIME)
        {
            double last_max_amp = *((double*)vector_back(&inst->max_Amp));
            if (last_max_amp < test_value)
            {
                vector_pop_back(&inst->max_Amp);
                vector_pop_back(&inst->max_time);
                vector_push_back(&inst->max_Amp, &test_value);
                vector_push_back(&inst->max_time, &test_sample_number);
            } else
            {
                /*Ignore this max*/
            }
        } else
        {
            vector_push_back(&inst->max_Amp, &test_value);
            vector_push_back(&inst->max_time, &test_sample_number);
        }

        if (inst->max_time.size > 1)
        {
            size_t _index = &inst->max_time.size - 2;
            double second_from_last = VECTOR_GET_AS(double, &inst->max_time, _index); 
            double last_maxtime = *((double*)vector_back(&inst->max_time));
            double new_heart_rate = (60.0 * kSAMPLE_RATE) / (double) (last_maxtime - second_from_last);

            if (is_heart_rate_valid(new_heart_rate))
            {
                vector_push_back(&inst->heart_rate, &new_heart_rate);
                is_valid = true;
            } else
            {
                /*Into here is a invalid pulse entered*/
                vector_clear(&inst->max_Amp);
                vector_clear(&inst->max_time);
                vector_clear(&inst->min_Amp);
                vector_clear(&inst->min_time);
                vector_clear(&inst->heart_rate);
                vector_push_back(&inst->max_time, &test_sample_number);
                vector_push_back(&inst->max_Amp, &test_value);
                is_valid = false;
            }
        }
    }

    return is_valid;
}

static bool check_maximal(OBPDetection_t *inst)
{
    bool is_valid = false;
    if (inst->oscillation_data.size <= kMIN_DATA_SIZE)
    {
        return false;
    }
    size_t start_index = inst->oscillation_data.size - 2;
    size_t local_domain = inst->oscillation_data.size - 3;
    size_t end_index = inst->oscillation_data.size;
    double data = VECTOR_GET_AS(double, &inst->oscillation_data, start_index);
    if (data < kPROMINENCE)
    {
        return false;
    }

    size_t local_max_idx = start_index;
    double local_max_value;
    for (size_t i = start_index - 1; i < end_index; i++)
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
        size_t first_max = VECTOR_GET_AS(size_t, &inst->max_time, first_time_index);
        size_t last_max = *((size_t*)vector_back(&inst->max_time)); 
        size_t min_pos = first_max;
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
        size_t distance = min_pos - first_max;

        if (inst->min_time.size == (inst->max_time.size - 1))
        {
            vector_pop_back(&inst->min_Amp);
            vector_push_back(&inst->min_Amp, &min_value);
            vector_pop_back(&inst->min_time);
            distance += first_time_index;
            vector_push_back(&inst->min_time, &distance);
        } else
        {
            distance += first_time_index;
            vector_push_back(&inst->min_Amp, &min_value);
            vector_push_back(&inst->min_time, &distance);
        }
    }
}

static bool is_enough_data(OBPDetection_t *inst)
{
    bool is_enough = false;
    if (inst->max_Amp.size < kMIN_NBR_PEAKS)
    {
        return false;
    }

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
    if (max_el > 1.5 && ((last_value < second_from_last && last_value < first_from_last) ||
                            last_value < 2 * kPROMINENCE))
    {
        double cut_off = max_el * (inst->ratio_DBP - kCUTOFF_HYST);

        if ((second_from_last < cut_off) && (first_from_last < cut_off) && (last_value < cut_off))
        {
            is_enough = true;
        }
    }

    return is_enough;
}

static void find_OWME(OBPDetection_t *inst)
{
    Iterator(itime_max1) = vector_begin(&inst->max_time);
    Iterator(iamp_min1) = vector_begin(&inst->min_Amp);
    Iterator(iamp_max1) = vector_begin(&inst->max_Amp);

    for (size_t i = 0; i < inst->min_time.size - 1; i++)
    {
        double amp_min1 = ITERATOR_GET_AS(double, &iamp_min1);
        double amp_max1 = ITERATOR_GET_AS(double, &iamp_max1);
        double time_max1 = ITERATOR_GET_AS(double, &itime_max1);
        double time_max2 = *((double*)iterator_next(&itime_max1));
        double amp_min2 = *((double*)iterator_next(&iamp_min1));
        double amp_max2 = *((double*)iterator_next(&iamp_max1));
        double time_min1 = VECTOR_GET_AS(double, &inst->min_time, i);
        double time_min2 = VECTOR_GET_AS(double, &inst->min_time, i + 1); 

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

/* PUBLIC API */
OBPDetection_t *OBP_get_instance()
{
    ref_counter ++;
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
}

double OBP_get_current_heart_rate(const OBPDetection_t *instance)
{
    return 0;
}

double OBP_get_average_heart_rate(const OBPDetection_t *instance)
{
    return 0;
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

void OBP_reset_instance(OBPDetection_t *instance)
{
    /*Perform reset all measurement variable before next step*/
    instance->result_MAP = 0.0f;
    instance->result_DBP = 0.0f;
    instance->result_SBP = 0.0f;
}
