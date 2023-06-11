#ifndef OBP_DATA_PROCESSING_H
#define OBP_DATA_PROCESSING_H
#include "log_wrapper.h"
#include "ret_code.h"

ret_code_t init_obp_processing();

ret_code_t obp_process_data(double new_data);

void obp_get_result(float *oMAP, float *oSBP, float *oDBP);

void obp_get_current_heart_rate(float *oheart_rate);

#endif /*OBP_DATA_PROCESSING_H*/