import sys
sys.path.append('../build/release/python')
from pyconinter import ConcentricData, ConcentricInterpolation

data_root = "../data/"
filename_support_directions = data_root + "directions/D8_N166_from_N83_sym1_s-2.txt"
filename_support_radii = data_root + "radii/radii_4.txt"
filename_support_values = data_root + "values/D8_N166from83_s-2_R4.txt"

support_data = ConcentricData(
    filename_support_directions, filename_support_radii, filename_support_values
)

concentric_interpolant = ConcentricInterpolation(support_data)
#concentric_interpolant.print_info()
print(concentric_interpolant.get_info())

filename_validation_directions = data_root + "directions/D8_N512_from_N256_sym1_s-2.txt"
filename_validation_radii = data_root + "radii/radii_4.txt"
filename_validation_values = data_root + "values/D8_N512from256_s-2_R4.txt"

validation_data = ConcentricData(
    filename_validation_directions,
    filename_validation_radii,
    filename_validation_values,
)

gamma_min = 0.1
gamma_max = 2.0
N_gamma_regular = 10
N_gamma_bisection = 5
concentric_interpolant.optimize_gamma(
    gamma_min,
    gamma_max,
    N_gamma_regular,
    N_gamma_bisection,
    validation_data,
    1,
    0,
    0,
    2,
)

print(f"\nRMS absolute error on whole validation data = {concentric_interpolant.compute_error(validation_data):e}\n\n")