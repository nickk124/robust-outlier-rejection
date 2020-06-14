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

using namespace pybind11::literals;

// variable argument conversion handlers

// pivot point finding
double pivotFunc_wrapper(std::vector <double> xdata, std::vector <double> weights, std::function<double(double, std::vector<double>)> f, std::vector <double> params, std::function <double(py::kwargs)> pivotFunc_py){
    /*
        Uses python pivot function kwargs to call the function
    */
    py::dict kwargs = py::dict("xdata"_a=xdata, "weights"_a=weights, "f"_a=f, "params"_a=params); // generate list of kwargs to be used by python func

    return pivotFunc_py(kwargs);
}

std::vector <double> pivotFunc_wrapper_ND(std::vector < std::vector <double> > xdata, std::vector <double> weights, std::function<double(std::vector <double>, std::vector<double>)> f, std::vector <double> params, std::function <std::vector <double> (py::kwargs)> pivotFunc_py){
    /*
        Uses python pivot function kwargs to call the function (ND version)
    */
    py::dict kwargs = py::dict("xdata"_a=xdata, "weights"_a=weights, "f"_a=f, "params"_a=params); // generate list of kwargs to be used by python func

    return pivotFunc_py(kwargs);
}

std::function <double(std::vector <double>, std::vector <double>, std::function<double(double, std::vector<double>)>, std::vector <double>)> getPivotFunc_cpp(std::function <double(py::kwargs)> pivotFunc_py){
    /*
        Wraps python pivot point function into cpp version
    */
    std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivot_func_cpp; // use explicit type to not avoid confusion with std::bind
    pivot_func_cpp = std::bind(pivotFunc_wrapper, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, pivotFunc_py);

    return pivot_func_cpp;
}

std::function <std::vector <double> (std::vector < std::vector <double> >, std::vector <double>, std::function<double(std::vector <double>, std::vector<double>)>, std::vector <double>)> getPivotFunc_cpp(std::function <std::vector <double>(py::kwargs)> pivotFunc_py){
    /*
        Wraps python pivot point function into cpp version (ND version)
    */
    std::function < std::vector <double>(std::vector < std::vector <double> >, std::vector <double>, std::function < double(std::vector <double>, std::vector <double>) >, std::vector <double>) > pivot_func_cpp; // use explicit type to not avoid confusion with std::bind
    pivot_func_cpp = std::bind(pivotFunc_wrapper_ND, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, pivotFunc_py);

    return pivot_func_cpp;
}


// factory functions
FunctionalForm getFunctionalFormObject(py::args args, py::kwargs kwargs){
    /*
        Bridges the gap between the single python FunctionalForm constructor with the many overloaded c++ ones

        args: f, xdata, ydata, model_partials, guess
        kwargs: weights, error_y, tol = 1e-6, has_priors = false, pivot_function = None, pivot_guess = None
    */
    FunctionalForm FFobj;

    // extract dimensionality independent types

    // args
    std::vector <double> ydata = py::cast <std::vector <double> >(args[2]);
    std::vector <double> guess = py::cast <std::vector <double> >(args[4]);

    // kwargs
    std::vector <double> error_y;
    std::vector <double> weights(ydata.size(), 1.0); // unweighted by default

    double tol = NAN;
    bool has_priors = false;
    bool find_pivot = false;

    if (kwargs.contains("weights")){
        weights = py::cast <std::vector <double> >(kwargs["weights"]);
    }
    if (kwargs.contains("error_y")){
        error_y = py::cast <std::vector <double> >(kwargs["error_y"]);
    }
    if (kwargs.contains("tol")){
        tol = py::cast <double>(kwargs["tol"]);
    }
    if (kwargs.contains("has_priors")){
        has_priors = py::cast <bool>(kwargs["has_priors"]);
    }
    if (kwargs.contains("pivot_function")){
        find_pivot = true;
    }

    bool isND = false;
    try { // check dimensionality by attempting a type conversion for x data
        py::cast< std::vector <double> >(args[1]);
    } catch (...){
        isND = true;
    }

    bool hasErrs = error_y.size() > 0;
    bool hasTol = !std::isnan(tol);

    const double DEFAULT_TOL = 1e-6;
    if (!hasTol){tol = DEFAULT_TOL;};

    // dimensionality dependent types
    std::function <double(std::vector <double>, std::vector <double> )> f_ND; 
    std::vector < std::vector <double> > xdata_ND; 
    std::vector < std::function <double(std::vector <double>, std::vector <double>)> > model_partials_ND;
    std::function <std::vector <double>(py::kwargs)> pivot_func_py_ND;
    std::function <std::vector <double> (std::vector < std::vector <double> >, std::vector <double>, std::function<double(std::vector <double>, std::vector<double>)>, std::vector <double>)> pivot_func_cpp_ND;
    std::vector <double> pivot_guess_ND;

    std::function <double(double, std::vector <double> )> f; 
    std::vector <double> xdata; 
    std::vector < std::function <double(double, std::vector <double>)> > model_partials;
    std::function <double(py::kwargs)> pivot_func_py;
    std::function <double(std::vector <double>, std::vector <double>, std::function<double(double, std::vector<double>)>, std::vector <double>)> pivot_func_cpp;
    double pivot_guess;

    if (isND){
        f_ND = py::cast< std::function <double(std::vector <double>, std::vector <double> )> >(args[0]);
        xdata_ND = py::cast< std::vector < std::vector <double> > >(args[1]);
        model_partials_ND = py::cast < std::vector < std::function <double(std::vector <double>, std::vector <double>)> > > (args[3]);
        
        if (find_pivot){
            pivot_func_py_ND = py::cast <std::function <std::vector <double>(py::kwargs)>>(kwargs["pivot_function"]);
            pivot_func_cpp_ND = getPivotFunc_cpp(pivot_func_py_ND);
            pivot_guess_ND = py::cast <std::vector <double> >(kwargs["pivot_guess"]);
        }
        
    } else {
        f = py::cast< std::function <double(double, std::vector <double> )> >(args[0]);
        std::vector <double> xdata = py::cast< std::vector <double> >(args[1]);
        model_partials = py::cast < std::vector < std::function <double(double, std::vector <double>)> > > (args[3]);

        if (find_pivot){
            pivot_func_py = py::cast <std::function <double(py::kwargs)>>(kwargs["pivot_function"]);
            pivot_func_cpp = getPivotFunc_cpp(pivot_func_py);
            pivot_guess = py::cast <double>(kwargs["pivot_guess"]);
        }
    }

    // build object
    Priors placeholderPriorsObj;
    if (isND){
        if (find_pivot){
            if (has_priors){
                if (hasErrs){
                    FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess, weights, placeholderPriorsObj, pivot_func_cpp_ND, pivot_guess_ND);
                } else {
                    FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess, weights, placeholderPriorsObj, pivot_func_cpp_ND, pivot_guess_ND);
                }
            } else {
                if (hasErrs){
                    FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess, weights, pivot_func_cpp_ND, pivot_guess_ND);
                } else {
                    FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess, weights, pivot_func_cpp_ND, pivot_guess_ND);
                }
            }
        } else {
            if (has_priors){
                if (hasErrs){
                    FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess, weights, placeholderPriorsObj);
                } else {
                    FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess, weights, placeholderPriorsObj);
                }
            } else {
                if (hasErrs){
                    FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess, weights);
                } else {
                    FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess, weights);
                }
            }
        }
    } else {
        if (find_pivot){
            if (has_priors){
                if (hasErrs){
                    FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess, weights, placeholderPriorsObj, pivot_func_cpp, pivot_guess);
                } else {
                    FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess, weights, placeholderPriorsObj, pivot_func_cpp, pivot_guess);
                }
            } else {
                if (hasErrs){
                    FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess, weights, pivot_func_cpp, pivot_guess);
                } else {
                    FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess, weights, pivot_func_cpp, pivot_guess);
                }
            }
        } else {
            if (has_priors){
                if (hasErrs){
                    FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess, weights, placeholderPriorsObj);
                } else {
                    FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess, weights, placeholderPriorsObj);
                }
            } else {
                if (hasErrs){
                    FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess, weights);
                } else {
                    FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess, weights);
                }
            }
        }
    }
    return FFobj;
}


// python binding functions
PYBIND11_MODULE(rcr, m) { // rcr is module name, m is docstring instance
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



    // main (single value) class
    py::class_<RCR>(m, "RCR", "Master class used to initialize and run RCR procedures")
        // constructors
        .def(py::init<RejectionTechs>(), py::arg("RejectionTechnique"))
        .def(py::init<>())
        
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

        .def("setParametricModel", &RCR::setParametricModel, "Initialize parametric/functional form model with RCR", py::arg("functionalform_model"))
        // .def("setNonParametricModel", &RCR::setNonParametricModel)

        // .def("setInitialModel", &RCR::setInitialModel)
        // .def("setMuType", &RCR::setMuType);

        // results
        .def_readwrite("result", &RCR::result, "Results from RCR");

    
    // functional form/model fitting RCR
    py::class_<FunctionalForm>(m, "FunctionalForm", "Class used to initialize functional form/model-fitting RCR")
    //     // constructors
        .def(py::init(&getFunctionalFormObject), "args: f, xdata, ydata, model_partials, guess\nkwargs: weights = None, error_y = None, tol = 1e-6, has_priors = false, pivot_function = None")
        .def(py::init<>())

        // members
        .def_readwrite("has_priors", &FunctionalForm::has_priors, "Are you applying prior probability distributions to your model parameters?")
        .def_readwrite("parameters", &FunctionalForm::parameters, "Current (or final) estimated model parameters");

    // parameter prior probability distributions
    py::enum_<priorTypes>(m, "priorsTypes", py::arithmetic(), "Parameter prior probability distribution types")
        .value("CUSTOM_PRIORS", CUSTOM_PRIORS, "Custom, function-defined prior probability distribution(s)")
        .value("GAUSSIAN_PRIORS", GAUSSIAN_PRIORS, "Gaussian (normal) prior probability distribution(s)")
        .value("CONSTRAINED_PRIORS", CONSTRAINED_PRIORS, "Bounded/hard-constrained prior probability distribution(s)")
        .value("MIXED_PRIORS", MIXED_PRIORS, "A mixture of gaussian (normal), hard-constrained, and uninformative (uniform/flat) prior probability distributions")
        .export_values();


    py::class_<Priors>(m, "Priors", "Class of prior probability distribution functions that can be applied to model parameters")
        // constructors , std::function <std::vector <double>(std::vector <double>, std::vector <double>)> 
        .def(py::init< priorTypes, std::function <std::vector <double>(std::vector <double>, std::vector <double>)> >()) // custom priors
        .def(py::init< priorTypes, std::vector < std::vector <double> > >()) // only Gaussian or only bounded/hard constraints
        .def(py::init< priorTypes, std::vector < std::vector <double> >, std::vector < std::vector <double> > >()) // mixed priors
        .def(py::init<>())

        // members
        .def_readwrite("priorType", &Priors::priorType, "Type of prior")
        .def_readwrite("p", &Priors::p, "a function that takes in a parameters vector and a weights vector and modifies the weights given the prior probability distirbution")
        .def_readwrite("gaussianParams", &Priors::gaussianParams, "A list that contains a list of mu and sigma for the guassian prior of each param. If no prior, then just use NANs.")
        .def_readwrite("paramBounds", &Priors::paramBounds, "a list that contains lists of the bounds of each param. If not bounded, use NANs, and if there's only one bound, use NAN for the other \"bound\".");
}
