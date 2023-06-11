add_subdirectory(3rd_party)
set(bpm_oscillometric_SOURCES src/obp_detection.c
                              src/obp_data_processing.c)

add_library(bpm_oscillometric STATIC ${bpm_oscillometric_SOURCES})
target_link_libraries(bpm_oscillometric PRIVATE vector filter)
include_directories(3rd_party)
target_include_directories(bpm_oscillometric PUBLIC interface)
target_include_directories(bpm_oscillometric PRIVATE include)