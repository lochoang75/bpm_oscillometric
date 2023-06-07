set(bpm_oscillometric_SOURCES src/obp_detection.c
                              src/bpm_fir.c
                              src/bpm_mode.c)

add_library(bpm_oscillometric STATIC ${bpm_oscillometric_SOURCES})
target_include_directories(bpm_oscillometric PRIVATE interface)
target_include_directories(bpm_oscillometric PRIVATE include)