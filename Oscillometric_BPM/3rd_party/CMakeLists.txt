add_library(filter STATIC filter/filter.c)
target_link_libraries(filter m)
add_library(vector STATIC vector/c_vector.c)