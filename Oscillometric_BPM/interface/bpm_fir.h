#ifndef BPM_FIR_H
#define BPM_FIR_H

typedef struct {
    unsigned int taps;
    float *coeff;
    float *buffer;
    unsigned int offset;
} FIRInfo_t;

float FIR_filter(float, FIRInfo_t *);
void FIR_reset_buffer(FIRInfo_t *); 
#endif /*BPM_FIR_H*/