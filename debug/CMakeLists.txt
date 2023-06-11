add_library(logging STATIC src/log_wrapper.c)
target_include_directories(logging PRIVATE "include")