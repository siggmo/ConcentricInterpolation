./interpolation_demo
./interpolation_demo \
    ../data/directions/dir_D8_N166_from_N83_sym1_s-2.txt \
    ../data/values/values_D8_N166from83_s-2_t01_R5.txt \
    ../data/directions/dir_D8_N10000_EQ.txt \
    ../data/values/values_D8_N10000_EQ_t01_R51.txt \
    outputs_interpolation_demo.txt \
    0.5

# 1st argument is the training directions. see the Section 7.2 of the original paper for explanation of the filenames
# 2nd argument is the original function values along the training directions at the radii specified in interpolation_demo.cxx
# 3rd argument is the validation directions
# 4th argument is the original function values along the validation directions at the validation radii specified in the code
# 5th argument is the name of the output file
# 6th argument is the kernel width parameter gamma. if a gamma search is performed in the code, this value is ignored.
