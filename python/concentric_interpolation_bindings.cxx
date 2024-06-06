#include "concentric_interpolation.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <sstream>

// made with ChatGPT, probably needs quite a bit of refinement... but works so far

namespace py = pybind11;

PYBIND11_MODULE(pyconinter, m) {
  py::class_<ConcentricInterpolation>(m, "ConcentricInterpolation")
      .def(py::init<>())
      .def(py::init<ConcentricData &>())
      .def("setup", &ConcentricInterpolation::Setup, py::arg("setup_data"),
           py::arg("gamma") = 0.0, py::arg("sym") = false)
      .def("evaluate",
           [](ConcentricInterpolation &self, py::array_t<double> direction,
              double radius) {
             py::buffer_info buf = direction.request();
             if (buf.ndim != 1)
               throw std::runtime_error("Direction should be a 1D NumPy array");
             auto dir_ptr = static_cast<double *>(buf.ptr);
             std::vector<double> result(buf.shape[0]);
             self.Evaluate(dir_ptr, radius, result.data());
             return py::array(result.size(), result.data());
           })
      .def("compute_error", &ConcentricInterpolation::Error,
           py::arg("validation_data"), py::arg("error_type") = 2,
           py::arg("d_val_start") = -1, py::arg("d_val_end") = -1,
           py::arg("only_radius") = -1)
      .def("optimize_gamma", &ConcentricInterpolation::OptimizeGamma,
           py::arg("gamma_min"), py::arg("gamma_max"),
           py::arg("n_gamma_regular"), py::arg("n_gamma_bisection"),
           py::arg("validation_data"), py::arg("error_type") = 2,
           py::arg("d_val_start") = -1, py::arg("d_val_end") = -1,
           py::arg("only_radius") = -1)
      .def("print_info", [](ConcentricInterpolation &self) {
            self.PrintInfo(stdout);
        }, "Print information about this object")
      .def("get_info", [](ConcentricInterpolation &self) {
            FILE* tmpf = std::tmpfile();
            self.PrintInfo(tmpf);
            std::rewind(tmpf);
            std::string result;
            char buffer[256];
            while (std::fgets(buffer, sizeof(buffer), tmpf) != nullptr) {
                result += buffer;
            }
            std::fclose(tmpf);

            return result;
        }, "Get information about this object, but don't print it. Returns a string.");

  py::class_<ConcentricData>(m, "ConcentricData")
      .def(py::init<char *, char *, char *>());
}
