#include <Iir.h>
#include <stdio.h>

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

int main(int argc, char *argvv[])
{
    Iir::Butterworth::LowPass<4> *iirLP = new Iir::Butterworth::LowPass<4>;
    Iir::Butterworth::HighPass<4> *iirHP = new Iir::Butterworth::HighPass<4>;
    iirLP->setup(1000, 10);
    iirHP->setup(1000, 0.5);
    FILE *data_file = fopen("filtering_result.dat", "r");
    if (data_file == NULL)
    {
        printf("Unable to open data file\n");
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
        double bpm_mmHg = convert_mV_to_mmHg(bpm_value);
        double yLP = iirLP->filter(bpm_mmHg);
        double yHp = iirHP->filter(yLP);
        static int counter = 0;
        printf("%d, %0.4f\n", counter++, yHp);
        memset(line, 0, len);
    } while (nread > 0);

    fclose(data_file);
    return 0;
}