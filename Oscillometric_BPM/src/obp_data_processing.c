#include "obp_data_processing.h"
#include "obp_const.h"
#include "obp_detection.h"
#include "filter/filter.h"

/* Local variable */
OBPDetection_t *l_detection = NULL;
BWLowPass *l_Lowpass = NULL;
BWHighPass *l_HightPass = NULL;

/* Private API */
static double convert_mV_to_mmHg(double mmV)
{
    double ambientV = 0.710;
    double mmHg_per_kPa = 7.5006157584566;
    double kPa_per_V = 50;
    double kPa_per_mmHg = 0.133322;
    double corrFact = 2.50;

    double ymmHg = (( mmV - ambientV) * kPa_per_V * corrFact)/ kPa_per_mmHg;
    return ymmHg;
}

/* Public API */
ret_code_t init_obp_processing()
{
    l_detection = OBP_get_instance();
    l_Lowpass = create_bw_low_pass_filter(4, kSAMPLE_RATE, 10);
    l_HightPass = create_bw_high_pass_filter(4, kSAMPLE_RATE, 0.5);
    if (l_detection == NULL || l_Lowpass == NULL || l_HightPass == NULL)
    {
        log_warn("Init failed, detection: %p, lowPass: %p, highpass: %p", l_detection,
                                                                          l_Lowpass,
                                                                          l_HightPass);
        return kRET_FAILED;
    }

    return kRET_OK;
}


ret_code_t obp_process_data(double data)
{
    static bool deflat = false;
    double obp_mmHg = convert_mV_to_mmHg(data);
    double yLp = bw_low_pass(l_Lowpass, obp_mmHg); 
    double yHp = bw_high_pass(l_HightPass, yLp);
    static int counter = 0;
    // printf("%d, %0.4f\n", counter++, yHp);
    if (obp_mmHg >= 180.0f)
    {
        deflat = true;
    }

    if (deflat)
    {
        OBP_process_sample(l_detection, yLp, yHp);
    }
    if (OBP_is_enough_data(l_detection))
    {
        return kRET_OK;
    }

    return kRET_TRY_AGAIN;
}

void obp_get_result(float *oMAP, float *oSBP, float *oDBP)
{
    *oMAP = OBP_get_MAP(l_detection);
    *oSBP = OBP_get_SBP(l_detection);
    *oDBP = OBP_get_DBP(l_detection);
}

void obp_get_current_heart_rate(float *oheart_rate)
{
    *oheart_rate = OBP_get_average_heart_rate(l_detection);
}