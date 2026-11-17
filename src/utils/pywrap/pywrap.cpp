// pywrap.cpp
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h> 
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include "basic_simulator.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

// Module/class name: `oilspillsim` for the pip build (setup.py defines it); the CMake target in
// src/utils/pywrap/CMakeLists.txt is called SIMULATOR, so the CMake build gets a module named SIMULATOR
// (the module name must match the file name of the extension for Python to import it).
#ifndef OSS_MODULE_NAME
#define OSS_MODULE_NAME SIMULATOR
#endif

namespace py = pybind11;
constexpr auto byref = py::return_value_policy::reference_internal;

PYBIND11_MODULE(OSS_MODULE_NAME, m) {
    m.doc() = "optional module docstring";
    py::class_<SIMULATOR>(m, MACRO_STRINGIFY(OSS_MODULE_NAME))
    .def(py::init<const std::string &, const double, const double, const double, const double, const double, const int, const int, const int, const int, bool>(),py::arg("_filepath"), py::arg("_dt") = 10., py::arg("_kw") = 0.5, py::arg("_kc") = 1., py::arg("_gamma") = 1., py::arg("_flow") = 50., py::arg("_number_of_sources") = int(3), py::arg("_max_contamination_value") = int(5), py::arg("_source_fuel") = int(1000),  py::arg("_random_seed") = int(-1),  py::arg("_triangular") = false)  
    .def(py::init<Eigen::MatrixXi &, const double, const double, const double, const double, const double, const int, const int, const int, const int, bool>(),py::arg("base_matrix"), py::arg("_dt") = 10., py::arg("_kw") = 0.5, py::arg("_kc") = 1., py::arg("_gamma") = 1., py::arg("_flow") = 50., py::arg("_number_of_sources") = int(3), py::arg("_max_contamination_value") = int(5), py::arg("_source_fuel") = int(1000),  py::arg("_random_seed") = int(-1),  py::arg("_triangular") = false)  
    .def("step", &SIMULATOR::step, py::call_guard<py::gil_scoped_release>())
    .def("reset", &SIMULATOR::reset, py::arg("_seed") = int(-1), py::arg("_source_points_pos") = Eigen::MatrixXi(), py::arg("_v") = Eigen::MatrixXd(), py::arg("_u") = Eigen::MatrixXd(), py::arg("_wind_speed") = Eigen::VectorXd(), py::arg("_contamination_position") = Eigen::MatrixXd(), py::call_guard<py::gil_scoped_release>())
    .def("get_normalized_density", &SIMULATOR::get_normalized_density, py::arg("gaussian") = true, py::call_guard<py::gil_scoped_release>())
    .def_readonly("source_points", &SIMULATOR::source_points, byref)
    .def_readonly("contamination_position", &SIMULATOR::contamination_position, byref)
    .def_readonly("density", &SIMULATOR::density, byref)
    .def_readonly("x", &SIMULATOR::x, byref)
    .def_readonly("y", &SIMULATOR::y, byref)
    .def_readonly("u", &SIMULATOR::u, byref)
    .def_readonly("v", &SIMULATOR::v, byref)
    .def_readonly("wind_speed", &SIMULATOR::wind_speed, byref)
    // ---- v0.1.0: dt-consistent physics by default + optional dynamics (legacy physics: see set_dynamics) ----
    .def("set_dynamics", &SIMULATOR::set_dynamics,
         py::arg("wind_variability") = 0.0, py::arg("wind_relaxation") = 1.0 / 30.0, py::arg("wind_max") = 1.4142135623730951,
         py::arg("source_speed") = 0.0, py::arg("source_turn") = 0.0, py::arg("source_drift") = 0.0,
         py::arg("decay_rate") = 0.0, py::arg("evaporable_fraction") = 1.0, py::arg("diffusion_dt_ref") = 0.3,
         py::arg("emission_dt_ref") = 0.3, py::arg("fix_removal_skip") = true,
         "Configure the dynamics (call after reset; they persist across resets). Calling it resets every option not given "
         "to its default. Legacy physics of the PyPI 0.0.1 wheel: diffusion_dt_ref=-1, emission_dt_ref=-1, fix_removal_skip=False. "
         "wind_variability: OU noise of the wind (0 = constant); wind_relaxation: reversion rate to the prevailing wind; "
         "wind_max: cap on |wind|; source_speed/source_turn: self-propelled moving source (cells/time, rad/sqrt(time)); "
         "source_drift: source drifts with kw*wind + kc*current; decay_rate: evaporation hazard (1/time); "
         "evaporable_fraction: fraction of new particles that can evaporate; diffusion_dt_ref > 0: random walk scaled by sqrt(dt/dt_ref); "
         "emission_dt_ref > 0: emission per unit of time (flow*dt/dt_ref particles per step); "
         "fix_removal_skip: fix the particle that skipped its step after a removal at the border.")
    .def("warmup", &SIMULATOR::warmup, py::arg("n_steps"), py::call_guard<py::gil_scoped_release>(),
         "Run n_steps simulator steps (start the spill before the mission).")
    .def("set_wind", &SIMULATOR::set_wind, py::arg("wind"), py::arg("also_mean") = true)
    .def("set_source_fuel", &SIMULATOR::set_source_fuel, py::arg("fuel"))
    .def("set_density", &SIMULATOR::set_density, py::arg("max_per_cell"), py::arg("log_scale") = false,
         "Crowding cap per cell (a full cell pushes newcomers to a free neighbour) and value that maps to 1 in "
         "get_normalized_density; log_scale: log1p(n)/log1p(max_per_cell) instead of n/max_per_cell. Call after reset, "
         "before stepping (the constructor argument _max_contamination_value sets the initial cap).")
    .def_property_readonly("max_per_cell", &SIMULATOR::get_max_per_cell)
    .def_property_readonly("density_log", &SIMULATOR::get_density_log)
    .def("get_source_position", &SIMULATOR::get_source_position)
    .def_readonly("wind_mean", &SIMULATOR::wind_mean, byref)
    .def_property_readonly("step_count", &SIMULATOR::step_count)
    .def_property_readonly("remaining_fuel", &SIMULATOR::remaining_fuel)
    .def_property_readonly("particle_count", &SIMULATOR::particle_count)
    ;
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}