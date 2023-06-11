#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "log_wrapper.h"
#include "obp_data_processing.h"

int main(int argc, char* argv[])
{
    init_obp_processing();
    FILE *data_file = fopen("filtering_result.dat", "r");
    if (data_file == NULL)
    {
        log_error("Unable to open data file");
        return -1;
    }
    
    int nread = 0;
    size_t len = 0;
    char *line = NULL;
    do {
        nread = getline(&line, &len, data_file);
        if (nread < 0)
        {
            break;
        }
        float bpm_value = strtof(line, NULL);
        ret_code_t ret = obp_process_data(bpm_value);
        if (ret == kRET_OK)
        {
            break;
        }
        memset(line, 0, len);
    } while (nread > 0);

    fclose(data_file);
    float MAP, SBP, DBP;
    float heart_rate;
    obp_get_result(&MAP, &SBP, &DBP);
    obp_get_current_heart_rate(&heart_rate);
    log_info("result MAP: %0.4f, SBP: %0.4f, DBP: %0.4f, heart_rate: %0.4f", MAP, SBP, DBP, heart_rate);
    return 0;
}