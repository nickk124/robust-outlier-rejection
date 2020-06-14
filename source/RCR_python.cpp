/*
 Robust Chauvenet Rejection (RCR) Official Codebase
 Active Author: Nick C. Konz
 Former Author: Michael Maples
 See license at https://github.com/nickk124/RCR

This file houses all of the RCR functionality that needs to be exposed to Python
*/
#include "RCR.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;
using namespace RCRLib;

PYBIND11_MODULE(rcr, m) {
    m.doc() = "RCR plugin";

    // enums that need to be exposed to python
    py::enum_<RejectionTechs>(m, "RejectionTechniques", py::arithmetic(), "RCR Rejection Techniques")
        .value("SS_MEDIAN_DL", SS_MEDIAN_DL, "Rejection technique for a symmetric uncontaminated distribution with two-sided contaminants")
        .value("LS_MODE_68", LS_MODE_68, "Rejection technique for a symmetric uncontaminated distribution with one-sided contaminants")
        .value("LS_MODE_DL", LS_MODE_DL, "Rejection technique for a symmetric uncontaminated distribution with a mixture of one-sided and two-sided contaminants")
        .value("ES_MODE_DL", ES_MODE_DL, "Rejection technique for a mildly asymmetric uncontaminated distribution and/or a very low number of data points")
        .export_values();


    // ABOVE: When the special tag py::arithmetic() is specified to the enum_ constructor, 
    // pybind11 creates an enumeration that also supports rudimentary arithmetic 
    // and bit-level operations like comparisons, and, or, xor, negation, etc.

    // results class
    py::class_<RCRResults>(m, "RCRResults", "Results that RCR generates")
        .def_readwrite("mu", &RCRResults::mu, "Recovered mean/median/mode (central value) of uncontaminated data distribution")
        .def_readwrite("stDev", &RCRResults::stDev, "Recovered standard deviation of uncontaminated data distribution")
        .def_readwrite("stDevBelow", &RCRResults::stDevBelow, "Recovered standard deviation below mu (mean/median/mode) of uncontaminated (asymmetric) data distribution")
        .def_readwrite("stDevAbove", &RCRResults::stDevAbove, "Recovered standard deviation above mu (mean/median/mode) of uncontaminated (asymmetric) data distribution")
        .def_readwrite("stDevTotal", &RCRResults::stDevTotal, "Recovered combined standard deviation both above and below mu (mean/median/mode), of uncontaminated (asymmetric) data distribution")
        .def_readwrite("sigma", &RCRResults::sigma, "Recovered robust 68.3-percentile deviation of uncontaminated data distribution")
        .def_readwrite("sigmaBelow", &RCRResults::sigmaBelow, "Recovered robust 68.3-percentile deviation below mu (mean/median/mode) of uncontaminated data distribution")
        .def_readwrite("sigmaAbove", &RCRResults::sigmaAbove, "Recovered robust 68.3-percentile deviation above mu (mean/median/mode) of uncontaminated data distribution")
        .def_readwrite("flags", &RCRResults::flags, "Ordered flags describing outlier status of each corresponding datapoint (true if datapoint is NOT an outlier)")
        .def_readwrite("indices", &RCRResults::indices, "Indices of dataset where are NOT outliers")
        .def_readwrite("cleanW", &RCRResults::cleanW, "Weights of non-outlying datapoints")
        .def_readwrite("cleanY", &RCRResults::cleanY, "Non-outlying datapoints")
        .def_readwrite("rejectedW", &RCRResults::rejectedW, "Weights of outlying datapoints")
        .def_readwrite("rejectedY", &RCRResults::rejectedY, "Outlying datapoints")
        .def_readwrite("originalW", &RCRResults::originalW, "Weights of initial supplied dataset")
        .def_readwrite("originalY", &RCRResults::originalY, "Initial supplied dataset");



    // main class
    py::class_<RCR>(m, "RCR", "Master class used to initialize and run RCR procedures")
        // constructors
        .def(py::init<>())
        .def(py::init<RejectionTechs>())
        
        // main methods
        .def("setRejectionTech", &RCR::setRejectionTech,
            "Set outlier rejection technique", py::arg("RejectionTechnique"))

        .def("performRejection", (void (RCR::*)(std::vector <double> &)) &RCR::performRejection,  // explicitly giving arguments is necessary for overloaded funcs
            "Perform outlier rejection WITHOUT bulk pre-rejection (slower)", py::arg("data"))
        .def("performBulkRejection", (void (RCR::*)(std::vector <double> &)) &RCR::performBulkRejection,
             "Perform outlier rejection WITH bulk pre-rejection (faster)", py::arg("data"))
        .def("performRejection", (void (RCR::*)(std::vector <double> &, std::vector <double> &)) &RCR::performRejection, 
            "Perform outlier rejection WITHOUT bulk pre-rejection (slower)", py::arg("weights"), py::arg("data"))
        .def("performBulkRejection", (void (RCR::*)(std::vector <double> &, std::vector <double> &)) &RCR::performBulkRejection,
             "Perform outlier rejection WITH bulk pre-rejection (faster)", py::arg("weights"), py::arg("data"))

        // .def("setParametricModel", &RCR::setParametricModel)
        // .def("setNonParametricModel", &RCR::setNonParametricModel)

        // .def("setInitialModel", &RCR::setInitialModel)
        // .def("setMuType", &RCR::setMuType);

        // results
        .def_readwrite("result", &RCR::result, "Results from RCR");
}