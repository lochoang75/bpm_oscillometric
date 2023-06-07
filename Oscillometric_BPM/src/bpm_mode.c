#include "bpm_mode.h"
#include "bpm_fir.h"
#include "bpm_adaptive_algorithm.h"

long map(long x, long in_min, long in_max, long out_min, long out_max)
{
    return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

void IPP_Detect()
{
    float pulse_sum = 0;
    float pulse = 0;
    for (int i = 1; i < pulse_index; ++i)
    {
        pulse_sum += (1000.0f / ((index_array[i] - index_array[i - 1]) * (1000.0f / (float)SAMPLING_RATE)) * 60.0f);
    }

    if (pulse_index > 2)
    {
        pulse = pulse_sum / (pulse_index - 1);
    }
    else
    {
        return;
    }

    for (int i = 1; i < pulse_index; ++i)
    {
        float ratio = fabs(1.0f - ((1000.0f / ((index_array[i] - index_array[i - 1]) * (1000.0f / (float)SAMPLING_RATE)) * 60.0f) / pulse));
        if (ratio > IPP_Ratio)
        {
            ++IPP;
        }
    }

    if (((float)IPP / (float)(pulse_index - 1)) > IHB_Ratio)
    {
        ++IHB;
    }
}

void BP_measurement(uint32_t is_usb_mode)
{

    if (calibration_flag)
    {
        if (micros - last_update > (2 * SECOND))
        {
            measure_flag = true;
            //TODO: calibration PID for pump

            calibration_flag = false;

            dc = baseline / init_baseline_count;
        }
        else
        {
        }
    }

    if (ADC3_ready && (measure_flag | calibration_flag) && is_usb_mode == DATA_COLLECTION_MODE)
    {
        ADC3_ready = false;
        unsigned char str[255];

        sprintf(str, "R%d,%d,%d,%d,%d,", micros, ADC3_value[0], ADC3_value[1], ADC3_value[2], ADC3_value[3]);
        if (calibration_flag)
        {
            baseline += FIR_filter(ADC3_value[0], &info2);
            ++init_baseline_count;
        }
        else
        {
            // float value_dc = ADC3_value[0] - dc;
            float value_dc = FIR_filter(ADC3_value[0], &info2) - dc;

            pressure = (a[0] * (value_dc * value_dc)) + (a[1] * value_dc) + a[2]; // 174base

            diff_pressure = pressure - last_pressure;

            Computing(micros);

            if (pressure > PRESSURE_MIN)
            {
                TIM3->CCR2 = map(pwm, PID_PWM_MIN, PID_PWM_MAX, TIMER_PWM_MIN, TIMER_PWM_MAX) / PWM_Freq; // 1~100 = duty cycle 10%~100%
            }
            else
            {
                TIM3->CCR2 = TIMER_PWM_33 / PWM_Freq;
            }

            last_pressure = pressure;
        }
    }
    else if (ADC3_ready && (measure_flag | calibration_flag) && is_usb_mode == BPM_MODE)
    {
        ADC3_ready = false;

        if (calibration_flag)
        {
            baseline += FIR_filter(ADC3_value[0], &info2);
            ++init_baseline_count;
        }
        else if (measure_flag)
        {

            float value = FIR_filter(ADC3_value[1], &info1);
            float value_dc = FIR_filter(ADC3_value[0], &info2) - dc;

            // pressure = (-0.000009f*(value_dc*value_dc))+(0.1361f*value_dc) - 0.661f; //78base
            // pressure = (-0.000001f*(value_dc*value_dc))+(0.1333f*value_dc) - 0.3395f; //121base
            pressure = (a[0] * (value_dc * value_dc)) + (a[1] * value_dc) + a[2]; // 174base

            diff_pressure = pressure - last_pressure;

            Computing(micros);

            if (pressure > PRESSURE_MIN)
            {
                TIM3->CCR2 = map(pwm, PID_PWM_MIN, PID_PWM_MAX, TIMER_PWM_MIN, TIMER_PWM_MAX) / PWM_Freq; // 1~100 = duty cycle 10%~100%
            }
            else
            {
                TIM3->CCR2 = TIMER_PWM_33 / PWM_Freq;
            }

            ++time_to_index;

            static const float CV_LIMIT = 50.0f;
            static const float THRESHOLD_FACTOR = 3.0f;
            float mean = CalculateMean(value);
            float rms = CalculateRootMeanSquare(value);
            float cv = CalculateCoefficientOfVariation(value);
            float threshold;
            if (cv > CV_LIMIT)
            {
                // threshold = rms;
                threshold = dc;
            }
            else
            {
                // threshold = (rms * (cv/100.0f) * THRESHOLD_FACTOR);
                threshold = (dc * (cv / 100.0f) * THRESHOLD_FACTOR);
            }

            bool is_peak;
            SignalPoint_t result;
            result = peak_deatect(value, time_to_index, threshold, &is_peak);
            if (result.index != -1 && (pressure > PRESSURE_MIN && pressure < PRESSURE_MAX))
            {
                if (is_peak)
                {
                    dc_array[pulse_index] = pressure;
                    ac_array[pulse_index] = result.value;
                    index_array[pulse_index] = result.index;
                    ++pulse_index;

                    detect_time = micros;
                }
                else
                {
                }
            }

            last_pressure = pressure;
        }
    }
    else if (ADC3_ready && is_usb_mode == AF_MODE)
    {
        if (!measure_flag && !calibration_flag && (AF_counter > 1) && ((micros - last_AF) > (AF_PERIOD * SECOND)))
        {
            --AF_counter;
            ResetMeasurementParameter();
        }

        ADC3_ready = false;

        if (calibration_flag)
        {
            baseline += FIR_filter(ADC3_value[0], &info2);
            ++init_baseline_count;
        }
        else if (measure_flag)
        {

            float value = FIR_filter(ADC3_value[1], &info1);
            float value_dc = FIR_filter(ADC3_value[0], &info2) - dc;
            // float value = ADC3_value[1];
            // float value_dc = ADC3_value[0] - dc;

            pressure = (a[0] * (value_dc * value_dc)) + (a[1] * value_dc) + a[2];

            diff_pressure = pressure - last_pressure;

            Computing(micros);

            if (pressure > PRESSURE_MIN)
            {
                // TODO: Control pump via pwm
            }
            else
            {
            }

            ++time_to_index;

            static const float CV_LIMIT = 50.0f;
            static const float THRESHOLD_FACTOR = 3.0f;
            float mean = CalculateMean(value);
            float rms = CalculateRootMeanSquare(value);
            float cv = CalculateCoefficientOfVariation(value);
            float threshold;
            if (cv > CV_LIMIT)
            {
                threshold = dc;
            }
            else
            {
                threshold = (dc * (cv / 100.0f) * THRESHOLD_FACTOR);
            }

            bool is_peak;
            SignalPoint_t result;
            result = PeakDetect(value, time_to_index, threshold, &is_peak);
            if (result.index != -1 && (pressure > PRESSURE_MIN && pressure < PRESSURE_MAX))
            {
                if (is_peak)
                {
                    dc_array[pulse_index] = pressure;
                    ac_array[pulse_index] = result.value;
                    index_array[pulse_index] = result.index;
                    ++pulse_index;

                    detect_time = micros;
                }
                else
                {
                }
            }

            last_pressure = pressure;
        }
    }

    if (measure_flag)
    {
        if (micros - last_update > (MEASURMENT_TIME * SECOND) || pressure > STOP_PRESSURE)
        {
            measure_flag = false;
            TIM3->CCR2 = 0;
            GPIO_ResetBits(GPIOB, GPIO_Pin_8);

            if (is_usb_mode == AF_MODE)
            {
                last_AF = micros;
                IPP = 0;
                IPP_Detect();
            }
        }
    }
}

void BP_calculator()
{
    uint32_t MAP;
    uint32_t SBP;
    uint32_t DBP;
    float map_ratio = 0;
    uint32_t map_index = 0;
    for (int i = 0; i < pulse_index; ++i)
    {
        float temp = ac_array[i] / dc;
        if (temp > map_ratio)
        {
            map_ratio = temp;
            map_index = i;

            MAP = dc_array[i];
        }
    }

    float dbp_error = 0;
    for (int i = map_index - 1; i >= 0; --i)
    {
        float temp = ac_array[i] / dc;
        float dbp_ratio = temp / map_ratio;
        if (fabs(dbp_ratio - ad_am_value) < dbp_error || dbp_error == 0)
        {
            dbp_error = fabs(dbp_ratio - ad_am_value);
            DBP = dc_array[i];
        }
    }

    float sbp_error = 0;
    for (int i = map_index + 1; i < pulse_index; ++i)
    {
        float temp = ac_array[i] / dc;
        float sbp_ratio = temp / map_ratio;
        if (fabs(sbp_ratio - as_am_value) < sbp_error || sbp_error == 0)
        {
            sbp_error = fabs(sbp_ratio - as_am_value);
            SBP = dc_array[i];
        }
    }

    float pulse_sum = 0;
    float pulse = 0;
    for (int i = 1; i < pulse_index; ++i)
    {
        pulse_sum += (1000.0f / ((index_array[i] - index_array[i - 1]) * (1000.0f / (float)SAMPLING_RATE)) * 60.0f);
    }
    if (pulse_index > 2)
    {
        pulse = pulse_sum / (pulse_index - 1);
    }
}