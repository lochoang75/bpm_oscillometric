#ifndef OBP_DETECTION_H
#define OBP_DETECTION_H 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

/* Private struct */
typedef struct sOBPDetection OBPDetection_t;

/* Warning: only create object use this function, auto allocation may cause segfault*/
OBPDetection_t *OBP_get_instance();

/* Setter */
void OBP_set_SBP_ratio(OBPDetection_t *instance, double val);
void OBP_set_DBP_ratio(OBPDetection_t *instance, double val);
void OBP_set_min_number_of_peaks(OBPDetection_t *instance, int val);

/* Sample processing */
void OBP_process_sample(OBPDetection_t *instance, double pressure, double oscillation);

/* Getter */
double OBP_get_current_heart_rate(const OBPDetection_t *instance);
double OBP_get_average_heart_rate(const OBPDetection_t *instance);
double OBP_get_MAP(const OBPDetection_t *instance);
double OBP_get_SBP(const OBPDetection_t *instance);
double OBP_get_DBP(const OBPDetection_t *instance);

/* Object reset */
void OBP_reset_instance(OBPDetection_t *instance);
#endif /*OBP_DETECTION*/