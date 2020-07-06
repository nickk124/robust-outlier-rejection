/*
 Robust Chauvenet Rejection (RCR) Official Codebase
 Active Author: Nick C. Konz
 Former Author: Michael Maples
 See license at https://github.com/nickk124/RCR

This file houses all of the RCR functionality that needs to be exposed to Python.

build cmd: (only tested on mac atm):
    c++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes` RCR_python.cpp RCR.cpp FunctionalForm.cpp MiscFunctions.cpp NonParametric.cpp -o rcr`python3-config --extension-suffix`

*/
#include "RCR.h"
#include <pybind11/pybind11.h> // pybind header files are within ./pybind11/include/pybind11/
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;
using namespace RCRLib;

using namespace pybind11::literals;

// variable argument conversion handlers

// support for variable argument number pivot point functions from python (WIP)
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
void buildFunctionalFormObject( // decision tree for builing a FF obj
    FunctionalForm &FFobj,
    std::function <double(std::vector <double>, std::vector <double> )> &f_ND, 
    std::function <double(double, std::vector <double> )> &f, 
    std::vector < std::vector <double> > &xdata_ND, 
    std::vector <double> &xdata, 
    std::vector <double> &ydata, 
    std::vector <double> &error_y, 
    std::vector < std::function <double(std::vector <double>, std::vector <double>)> > &model_partials_ND, 
    std::vector < std::function <double(double, std::vector <double>)> > &model_partials,
    double tol, 
    std::vector <double> &guess, 
    std::vector <double> &weights, 
    std::function <std::vector <double> (std::vector < std::vector <double> >, std::vector <double>, std::function<double(std::vector <double>, std::vector<double>)>, std::vector <double>)> &pivot_func_cpp_ND, 
    std::function <double(std::vector <double>, std::vector <double>, std::function<double(double, std::vector<double>)>, std::vector <double>)> &pivot_func_cpp, 
    std::vector <double> &pivot_guess_ND,
    double pivot_guess,
    bool has_weights,
    bool isND,
    bool find_pivot,
    bool has_priors,
    bool has_errors
){
    Priors placeholderPriorsObj;
    if (has_weights){
        if (isND){
            if (find_pivot){
                if (has_priors){
                    if (has_errors){
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess, weights, placeholderPriorsObj, pivot_func_cpp_ND, pivot_guess_ND);
                    } else {
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess, weights, placeholderPriorsObj, pivot_func_cpp_ND, pivot_guess_ND);
                    }
                } else {
                    if (has_errors){
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess, weights, pivot_func_cpp_ND, pivot_guess_ND);
                    } else {
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess, weights, pivot_func_cpp_ND, pivot_guess_ND);
                    }
                }
            } else {
                if (has_priors){
                    if (has_errors){
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess, weights, placeholderPriorsObj);
                    } else {
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess, weights, placeholderPriorsObj);
                    }
                } else {
                    if (has_errors){
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess, weights);
                    } else {
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess, weights);
                    }
                }
            }
        } else {
            if (find_pivot){
                if (has_priors){
                    if (has_errors){
                        FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess, weights, placeholderPriorsObj, pivot_func_cpp, pivot_guess);
                    } else {
                        FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess, weights, placeholderPriorsObj, pivot_func_cpp, pivot_guess);
                    }
                } else {
                    if (has_errors){
                        FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess, weights, pivot_func_cpp, pivot_guess);
                    } else {
                        FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess, weights, pivot_func_cpp, pivot_guess);
                    }
                }
            } else {
                if (has_priors){
                    if (has_errors){
                        FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess, weights, placeholderPriorsObj);
                    } else {
                        FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess, weights, placeholderPriorsObj);
                    }
                } else {
                    if (has_errors){
                        FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess, weights);
                    } else {
                        FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess, weights);
                    }
                }
            }
        }
    } else {
        if (isND){
            if (find_pivot){
                if (has_priors){
                    if (has_errors){
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess, placeholderPriorsObj, pivot_func_cpp_ND, pivot_guess_ND);
                    } else {
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess, placeholderPriorsObj, pivot_func_cpp_ND, pivot_guess_ND);
                    }
                } else {
                    if (has_errors){
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess, pivot_func_cpp_ND, pivot_guess_ND);
                    } else {
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess, pivot_func_cpp_ND, pivot_guess_ND);
                    }
                }
            } else {
                if (has_priors){
                    if (has_errors){
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess, placeholderPriorsObj);
                    } else {
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess, placeholderPriorsObj);
                    }
                } else {
                    if (has_errors){
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, error_y, model_partials_ND, tol, guess);
                    } else {
                        FFobj = FunctionalForm(f_ND, xdata_ND, ydata, model_partials_ND, tol, guess);
                    }
                }
            }
        } else {
            if (find_pivot){
                if (has_priors){
                    if (has_errors){
                        FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess, placeholderPriorsObj, pivot_func_cpp, pivot_guess);
                    } else {
                        FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess, placeholderPriorsObj, pivot_func_cpp, pivot_guess);
                    }
                } else {
                    if (has_errors){
                        FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess, pivot_func_cpp, pivot_guess);
                    } else {
                        FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess, pivot_func_cpp, pivot_guess);
                    }
                }
            } else {
                if (has_priors){
                    if (has_errors){
                        FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess, placeholderPriorsObj);
                    } else {
                        FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess, placeholderPriorsObj);
                    }
                } else {
                    if (has_errors){
                        FFobj = FunctionalForm(f, xdata, ydata, error_y, model_partials, tol, guess);
                    } else {
                        FFobj = FunctionalForm(f, xdata, ydata, model_partials, tol, guess);
                    }
                }
            }
        }
    }
}

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
    std::vector <double> error_y, weights;

    const double DEFAULT_TOL = 1e-6;

    // constructor decision booleans
    double tol = DEFAULT_TOL;
    bool has_priors = false;
    bool find_pivot = false;
    bool has_errors = false;
    bool has_weights = false;
    bool isND = false;

    if (kwargs.contains("weights")){
        weights = py::cast <std::vector <double> >(kwargs["weights"]);
        has_weights = true; // must use this boolean instead of, say using equal weights if unweighted, due to the core construction of RCR
    }
    if (kwargs.contains("error_y")){
        error_y = py::cast <std::vector <double> >(kwargs["error_y"]);
        has_errors = true;
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

    try { // check dimensionality by attempting a type conversion for x data
        py::cast< std::vector <double> >(args[1]);
    } catch (...){
        isND = true;
    }

    // dimensionality dependent types
    std::function <double(std::vector <double>, std::vector <double> )> f_ND; // model function
    std::vector < std::vector <double> > xdata_ND; // x data
    std::vector < std::function <double(std::vector <double>, std::vector <double>)> > model_partials_ND; // model function partials
    // std::function <std::vector <double>(py::kwargs)> pivot_func_py_ND;
    std::function <std::vector <double> (std::vector < std::vector <double> >, std::vector <double>, std::function<double(std::vector <double>, std::vector<double>)>, std::vector <double>)> pivot_func_cpp_ND; // pivot point function
    std::vector <double> pivot_guess_ND; // guess for pivot points

    // same as above, but for a single independent ("x") variable
    std::function <double(double, std::vector <double> )> f; 
    std::vector <double> xdata; 
    std::vector < std::function <double(double, std::vector <double>)> > model_partials;
    // std::function <double(py::kwargs)> pivot_func_py;
    std::function <double(std::vector <double>, std::vector <double>, std::function<double(double, std::vector<double>)>, std::vector <double>)> pivot_func_cpp;
    double pivot_guess = NAN;

    // cast from python arguments to cpp types
    if (isND){
        f_ND = py::cast< std::function <double(std::vector <double>, std::vector <double> )> >(args[0]);
        xdata_ND = py::cast< std::vector < std::vector <double> > >(args[1]);
        model_partials_ND = py::cast < std::vector < std::function <double(std::vector <double>, std::vector <double>)> > > (args[3]);
        
        if (find_pivot){
            // pivot_func_py_ND = py::cast <std::function <std::vector <double>(py::kwargs)>>(kwargs["pivot_function"]);
            // pivot_func_cpp_ND = getPivotFunc_cpp(pivot_func_py_ND);
            pivot_func_cpp_ND = py::cast < std::function <std::vector <double> (std::vector < std::vector <double> >, std::vector <double>, std::function<double(std::vector <double>, std::vector<double>)>, std::vector <double>)> >(kwargs["pivot_function"]);
            pivot_guess_ND = py::cast <std::vector <double> >(kwargs["pivot_guess"]);
        }
        
    } else {
        f = py::cast< std::function <double(double, std::vector <double> )> >(args[0]);
        xdata = py::cast< std::vector <double> >(args[1]);
        model_partials = py::cast < std::vector < std::function <double(double, std::vector <double>)> > > (args[3]);

        if (find_pivot){
            // pivot_func_py = py::cast <std::function <double(py::kwargs)>>(kwargs["pivot_function"]);
            // pivot_func_cpp = getPivotFunc_cpp(pivot_func_py);
            pivot_func_cpp = py::cast < std::function <double(std::vector <double>, std::vector <double>, std::function<double(double, std::vector<double>)>, std::vector <double>)> >(kwargs["pivot_function"]);
            pivot_guess = py::cast <double>(kwargs["pivot_guess"]);
        }
    }

    // decision tree for builing a FF obj
    buildFunctionalFormObject(
        FFobj, f_ND, f, xdata_ND, xdata, ydata, error_y, model_partials_ND, model_partials, tol, 
        guess, weights, pivot_func_cpp_ND, pivot_func_cpp, pivot_guess_ND, pivot_guess,
        has_weights, isND, find_pivot, has_priors, has_errors
    );

    return FFobj;
}

// def performRejection_wrapper(py::args args, py::kwargs kwargs){


// }

// python binding functions
PYBIND11_MODULE(rcr, m) { // rcr is module name, m is docstring instance
    m.doc() = "RCR (Robust Chauvenet Outlier Rejection) Package API Details.";

    // MAIN RCR ######################################################################################################################################################

    // enums that need to be exposed to python
    py::enum_<RejectionTechs>(m, "RejectionTechniques", py::arithmetic(), "RCR Standard Rejection Techniques.")
        .value("SS_MEDIAN_DL", SS_MEDIAN_DL, "Rejection technique for a symmetric uncontaminated distribution with two-sided contaminants.")
        .value("LS_MODE_68", LS_MODE_68, "Rejection technique for a symmetric uncontaminated distribution with one-sided contaminants.")
        .value("LS_MODE_DL", LS_MODE_DL, "Rejection technique for a symmetric uncontaminated distribution with a mixture of one-sided and two-sided contaminants.")
        .value("ES_MODE_DL", ES_MODE_DL, "Rejection technique for a mildly asymmetric uncontaminated distribution and/or a very low number of data points.")
        .export_values();


    // ABOVE: When the special tag py::arithmetic() is specified to the enum_ constructor, 
    // pybind11 creates an enumeration that also supports rudimentary arithmetic 
    // and bit-level operations like comparisons, and, or, xor, negation, etc.imp


    // RCR results class
    py::class_<RCRResults>(m, "RCRResults", "Various results from performing outlier rejection with RCR.")
        .def_readwrite("mu", &RCRResults::mu, R"mydelimiter(
            *float*. Mean/median/mode (central value) of uncontaminated data distribution.

            The central value of the uncontaminated part of the provided dataset, recovered from
            performing RCR.
        )mydelimiter")

        .def_readwrite("stDev", &RCRResults::stDev, R"mydelimiter(
            *float*. Standard deviation of uncontaminated data distribution.

            The standard deviation/width :math:`\sigma` of the uncontaminated part of the provided dataset, recovered from
            performing RCR. For the case of a symmetric uncontaminated data distribution.
        )mydelimiter")

        .def_readwrite("stDevBelow", &RCRResults::stDevBelow, R"mydelimiter(
            *float*. Standard deviation below mu (mean/median/mode) of uncontaminated (asymmetric) data distribution.

            The *asymmetric* standard deviation/width :math:`\sigma_-` of the negative side of a mildly asymmetric uncontaminated data distribution, recovered from RCR (for the symmetric case, :math:`\sigma_-=\sigma_+\equiv\sigma`).
        )mydelimiter")

        .def_readwrite("stDevAbove", &RCRResults::stDevAbove, R"mydelimiter(
            *float*. Standard deviation above mu (mean/median/mode) of uncontaminated (asymmetric) data distribution.

            The *asymmetric* standard deviation/width :math:`\sigma_+` of the positive side of a mildly asymmetric uncontaminated data distribution, recovered from RCR (for the symmetric case, :math:`\sigma_+=\sigma_-\equiv\sigma`).
        )mydelimiter")

        .def_readwrite("stDevTotal", &RCRResults::stDevTotal, R"mydelimiter(
            *float*. Combined standard deviation both above and below mu (mean/median/mode) of uncontaminated (asymmetric) data distribution.

            A combination of the *asymmetric* standard deviation/width :math:`\sigma_+` of the positive side of a mildly asymmetric uncontaminated data distribution and the width :math:`\sigma_-` of the negative side 
            of the distribution, recovered from RCR. Can be used to approximate a mildly asymmetric data distribution as symmetric.
        )mydelimiter")

        .def_readwrite("sigma", &RCRResults::sigma, R"mydelimiter(
            *float*. Recovered robust 68.3-percentile deviation of uncontaminated data distribution.

            A more robust (less sensitive to outliers) version of the standard deviation/width :math:`\sigma`
            of the uncontaminated part of the provided dataset (see Section 2.1 of :ref:`papers`), recovered from
            performing RCR. For the case of a symmetric uncontaminated data distribution.
        )mydelimiter")

        .def_readwrite("sigmaBelow", &RCRResults::sigmaBelow, R"mydelimiter(
            *float*. Recovered robust 68.3-percentile deviation below mu (mean/median/mode) of uncontaminated data distribution.

            A more robust (less sensitive to outliers) version of the standard deviation/width :math:`\sigma_-`
            of the negative side of a mildly asymmetric uncontaminated data distribution (see Section 2.1 of :ref:`papers`), recovered from
            performing RCR. (For the symmetric case, :math:`\sigma_-=\sigma_+\equiv\sigma`).
        )mydelimiter")

        .def_readwrite("sigmaAbove", &RCRResults::sigmaAbove, R"mydelimiter(
            *float*. Recovered robust 68.3-percentile deviation abpve mu (mean/median/mode) of uncontaminated data distribution.

            A more robust (less sensitive to outliers) version of the standard deviation/width :math:`\sigma_+`
            of the positive side of a mildly asymmetric uncontaminated data distribution (see Section 2.1 of :ref:`papers`), recovered from
            performing RCR. (For the symmetric case, :math:`\sigma_+=\sigma_-\equiv\sigma`).
        )mydelimiter")

        .def_readwrite("flags", &RCRResults::flags, R"mydelimiter(
            *list of bools*. Ordered flags describing outlier status of each inputted datapoint (True if datapoint is NOT an outlier).

            For example, if a dataset of ``y = [0, 1, -2, 1, 2, 37, 0.5, -100]`` was provided and only the ``37`` and ``-100`` were found to be outliers,
            then ``flags = [True, True, True, True, True, False, True, False]``.
        )mydelimiter")

        .def_readwrite("indices", &RCRResults::indices, R"mydelimiter(
            *list of ints*. A list of indices of datapoints from original inputted dataset that are NOT outliers.

            For example, if a dataset of ``y = [0, 1, -2, 1, 2, 37, 0.5, -100]`` was provided and only the ``37`` and ``-100`` were found to be outliers,
            then ``indices = [0, 1, 2, 3, 4, 6]``.
        )mydelimiter")

        .def_readwrite("cleanW", &RCRResults::cleanW, R"mydelimiter(
            *list of floats*. The user-provided datapoint weights that correspond to NON-outliers in the original dataset.

            For example, if a dataset of ``y = [0, 1, -2, 1, 2, 37, 0.5, -100]`` with weights ``w = [1, 1.1, 0.9, 1.2, 0.8, 0.2, 0.95, 2]``
            was provided, and only the ``37`` and ``-100`` were found to be outliers,
            then ``cleanW = [1, 1.1, 0.9, 1.2, 0.8, 0.95]``.
        )mydelimiter")

        .def_readwrite("cleanY", &RCRResults::cleanY, R"mydelimiter(
            *list of floats*. After performing RCR on some original dataset, these are the datapoints that were NOT found to be outliers.

            For example, if a dataset of ``y = [0, 1, -2, 1, 2, 37, 0.5, -100]`` was provided and only the ``37`` and ``-100`` were found to be outliers,
            then ``cleanY = [0, 1, -2, 1, 2, 0.5]``.
        )mydelimiter")

        .def_readwrite("rejectedW", &RCRResults::rejectedW, R"mydelimiter(
            *list of floats*. The user-provided datapoint weights that correspond to outliers in the original dataset.

            For example, if a dataset of ``y = [0, 1, -2, 1, 2, 37, 0.5, -100]`` with weights ``w = [1, 1.1, 0.9, 1.2, 0.8, 0.2, 0.95, 2]``
            was provided, and only the ``37`` and ``-100`` were found to be outliers,
            then ``rejectedW = [0.2, 2]``.
        )mydelimiter")

        .def_readwrite("rejectedY", &RCRResults::rejectedY, R"mydelimiter(
            *list of floats*. After performing RCR on some original dataset, these are the datapoints that WERE found to be outliers.

            For example, if a dataset of ``y = [0, 1, -2, 1, 2, 37, 0.5, -100]`` was provided and only the ``37`` and ``-100`` were found to be outliers,
            then ``rejectedY = [37, -100]``.
        )mydelimiter")

        .def_readwrite("originalW", &RCRResults::originalW, R"mydelimiter(
            *list of floats*. The user-provided datapoint weights, pre-RCR.

            For example, if a dataset with weights ``w = [1, 1.1, 0.9, 1.2, 0.8, 0.2, 0.95, 2]``
            was provided, then ``originalW = [1, 1.1, 0.9, 1.2, 0.8, 0.2, 0.95, 2]``.
        )mydelimiter")

        .def_readwrite("originalY", &RCRResults::originalY, R"mydelimiter(
            *list of floats*. The user-provided dataset, pre-RCR.

            For example, if a dataset of ``y = [0, 1, -2, 1, 2, 37, 0.5, -100]`` was provided, 
            then ``originalY = [0, 1, -2, 1, 2, 37, 0.5, -100]``.
        )mydelimiter");


    // main (single value) class
    py::class_<RCR>(m, "RCR", "Master class used to initialize and run RCR outlier rejection procedures.")
        // constructors
        .def(py::init<RejectionTechs>(), py::arg("RejectionTechnique"))
        .def(py::init<>())

        // results
        .def_readwrite("result", &RCR::result, R"mydelimiter(
            ``rcr.RCRResults`` *object*. Access various results of RCR with this (see :class:`rcr.RCRResults`).
        )mydelimiter")
        
        // main methods
        .def("setRejectionTech", &RCR::setRejectionTech, R"mydelimiter(
            Modify/set outlier rejection technique to be used with RCR.

            See :ref:`rejectiontechs` for an explanation of each rejection technique, and when to use it.

            Parameters
            ----------
            rejection_technique : :class:`rcr.RejectionTechniques`
                The rejection technique to be used with your instance of :class:`rcr.RCR`.
        )mydelimiter", py::arg("rejection_technique"))

        // explicitly giving arguments is necessary for overloaded funcs
        .def("performRejection", (void (RCR::*)(std::vector <double> &)) &RCR::performRejection, R"mydelimiter(
            Perform outlier rejection WITHOUT the speed-up of bulk pre-rejection (slower; see :ref:`bulk`).

            *Parameters:*
            
            data : list/array_like, 1D
                Dataset to perform outlier rejection (RCR) on. Access results via the ``result`` attribute 
                (:class:`rcr.RCRResults`) of your instance of :class:`rcr.RCR`.
        )mydelimiter", py::arg("data"))

        .def("performBulkRejection", (void (RCR::*)(std::vector <double> &)) &RCR::performBulkRejection, R"mydelimiter(
            Perform outlier rejection WITH the speed-up of bulk pre-rejection (see :ref:`bulk`).

            *Parameters:*
            
            data : list/array_like, 1D
                Dataset to perform outlier rejection (RCR) on. Access results via the ``result`` attribute 
                (:class:`rcr.RCRResults`) of your instance of :class:`rcr.RCR`.
        )mydelimiter", py::arg("data"))

        .def("performRejection", (void (RCR::*)(std::vector <double> &, std::vector <double> &)) &RCR::performRejection, R"mydelimiter(
            Perform outlier rejection WITHOUT the speed-up of bulk pre-rejection (slower; see :ref:`bulk`).

            *Parameters:*
            
            weights : list/array_like, 1D
                Weights for dataset to perform outlier rejection (RCR) on.
            data : list/array_like, 1D
                Dataset to perform outlier rejection (RCR) on. Access results via the ``result`` attribute 
                (:class:`rcr.RCRResults`) of your instance of :class:`rcr.RCR`.
        )mydelimiter", py::arg("weights"), py::arg("data"))

        .def("performBulkRejection", (void (RCR::*)(std::vector <double> &, std::vector <double> &)) &RCR::performBulkRejection, R"mydelimiter(
            Perform outlier rejection WITH the speed-up of bulk pre-rejection (see :ref:`bulk`).

            *Parameters:*
            
            weights : list/array_like, 1D
                Weights for dataset to perform outlier rejection (RCR) on.
            data : list/array_like, 1D
                Dataset to perform outlier rejection (RCR) on. Access results via the ``result`` attribute 
                (:class:`rcr.RCRResults`) of your instance of :class:`rcr.RCR`.
        )mydelimiter", py::arg("weights"), py::arg("data"))

        // functional form/ model-fitting
        .def("setParametricModel", &RCR::setParametricModel, R"mydelimiter(
            Initialize parametric/functional form model to be used with RCR (see :ref:`functional` for a tutorial).

            Parameters
            ----------
            model : :class:`rcr.FunctionalForm`
                :math:`n`-dimensional model to fit data to while performing outlier rejection.
        )mydelimiter", py::arg("model"));


    // FUNCTIONAL FORM/MODEL-FITTING RCR #############################################################################################################################

    // Functional Form RCR results class
    py::class_<FunctionalFormResults>(m, "FunctionalFormResults", "Results from (and unique to) functional form/model-fitting RCR.")
        .def_readwrite("parameters", &FunctionalFormResults::parameters, R"mydelimiter(
            *list of floats*. Best-fit model parameters, post-outlier rejection.

            For example, if you're fitting to some linear model :math:`y(x|b,m)=b+mx`, and you obtain a best fit of :math:`b=1` and :math:`m=2`, then ``parameters = [1, 2]``.
        )mydelimiter")

        .def_readwrite("parameter_uncertainties", &FunctionalFormResults::parameter_uncertainties, R"mydelimiter(
            *list of floats*. Best-fit model parameter uncertainties, post-outlier rejection.

            For example, if you're fitting to some linear model :math:`y(x|b,m)=b+mx`, and you obtain a best fit of :math:`b=1.0\pm0.5` and :math:`m=2\pm 1`, then ``parameter_uncertainties = [0.5, 1]``.
            
            Note that in order for parameter uncertainties to be computed, either/both weights and data error bars/uncertainties must have been provided when constructing the :class:`rcr.FunctionalForm` model.
        )mydelimiter")

        .def_readwrite("pivot", &FunctionalFormResults::pivot, R"mydelimiter(
           *float*. Recovered optimal \"pivot\" point for model that should minimize correlation between the slope and intercept parameters of the linearized model (1D independent variable case).

           See :ref:`pivots`. For example, the pivot point for the model :math:`y(x|b,m) = b + m(x-x_p)` is :math:`x_p`.
        )mydelimiter")
        .def_readwrite("pivot_ND", &FunctionalFormResults::pivot_ND, R"mydelimiter(
           *float* Recovered optimal :math:`n`-dimensional \"pivot\" point for model that should minimize correlation between the slope and intercept parameters of the linearized model (:math:`n`-D independent variable case).

           See :ref:`pivots`. For example, the pivot point for the :math:`n`-dimensional model :math:`y(\vec{x}|\vec{b},\vec{m}) = \vec{b} + \vec{m}^T(\vec{x}-\vec{x}_p)` is :math:`\vec{x}_p`.
        )mydelimiter");


    // main class
    py::class_<FunctionalForm>(m, "FunctionalForm", R"mydelimiter(
            *class*. Class used to initialize functional form/model-fitting RCR (see :ref:`functional`).

            Constructor arguments:

            Parameters
            ----------
            f : function
                Model function :math:`y(\vec{x}|\vec{\theta})` to fit data to while performing outlier rejection, where :math:`\vec{x}` is an :math:`n`-dimensional list/array_like (or float, for 1D models) of independent variables and :math:`\vec{\theta}` is an
                :math:`M`-dimensional list/array_like of model parameters. Arguments for ``f`` must follow this prototype:

                Parameters
                ----------
                x : float or 1D list/array_like
                    Independent variable(s) of model
                params : list/array_like, 1D
                    Parameters of model

                Returns
                -------
                y : float
                    Model evaluated at the corresponding values of ``x`` and ``params``.
            
            xdata : list/array_like, 1D or 2D
                :math:`n`-dimensional independent variable data to fit model to. For 1D models (:math:`n=1`), this will be a 1D list/array_like, while for :math:`n`-D models, this will be a 2D list/array_like where each entry is a list/array_like of length :math:`n`.
            ydata : list/array_like, 1D
                Dependent variable (model function evaluation) data to fit model to. 
            model_partials : list of functions
                A list of functions that return the partial derivatives of the model function ``f`` with respect to each, ordered, model parameter :math:`\vec{\theta}` (See :ref:`functional` for an example). Arguments for each one of these functions must follow this prototype 
                (same as for the model function ``f``):

                Parameters
                ----------
                x : float or 1D list/array_like
                    Independent variable(s) of model
                params : list/array_like, 1D
                    Parameters of model

                Returns
                -------
                y : float
                    Derivative of model (with respect to given model parameter), evaluated at the corresponding values of ``x`` and ``params``.
            
            guess : list/array_like, 1D
                Guess for best fit values of model parameters :math:`\vec{\theta}` (for the fitting algorithm).

            weights : list/array_like, optional, 1D
                Optional weights to be applied to dataset (see :ref:`weighting`).

            error_y : list/array_like, optional, 1D
                Optional error bars/:math:`y`-uncertainties to be applied to dataset (see :ref:`errorbars`).
            
            tol : float, optional
                Default: ``1e-6``. Convergence tolerance of modified Gauss-Newton fitting algorithm.

            has_priors : bool, optional
                Default: ``False``. Set to ``True`` if you're going to apply statistical priors to your model parameters 
                (see :ref:`priors`; you'll also need to create an instance of :class:`rcr.Priors` and set the ``priors`` attribute of this instance of ``FunctionalForm`` equal to it).

            pivot_function : function, optional
                Default: ``None``. Function that returns the pivot point of some linearized model (see :ref:`pivots`). Must be of the form/prototype of:

                Parameters
                ----------
                xdata : list/array_like, 1D or 2D
                    :math:`n`-dimensional independent variable data to fit model to; same as above``xdata``.
                weights : list/array_like, optional, 1D
                    Optional weights to be applied to dataset (see :ref:`weighting`).
                f : function
                    Model function; same as above ``f``.
                params : list/array_like, 1D
                    Parameters of model

                Returns
                -------
                pivot : float or 1D list/array_like
                    Pivot point(s) of the model; (``float`` if you're using a one-dimensional model/independent variable, ``list/array_like`` if :math:`n`-dimensional.)

                However, note that all arguments need to be actually used for the pivot point computation. For example,
                a simple linear model :math:`y(x|b,m) = b + m(x-x_p)` has a pivot point found by :math:`x_p=\sum_iw_ix_i/\sum_iw_i`, where
                :math:`w_i` are the weights of the datapoints.
            
            pivot_guess : float or 1D list/array_like, optional
                Initial guess for the pivot point(s) of the model (``float`` if you're using a one-dimensional model/independent variable, ``list/array_like`` if :math:`n`-dimensional; see :ref:`pivots`).
        )mydelimiter")
        // constructors
        .def(py::init(&getFunctionalFormObject))
        .def(py::init<>())

        // results
        .def_readwrite("result", &FunctionalForm::result, R"mydelimiter(
            ``rcr.FunctionalFormResults`` *object*. Access various results unique to Functional Form RCR with this (see :class:`rcr.FunctionalFormResults`).
        )mydelimiter")

        // members
        .def_readwrite("priors", &FunctionalForm::priors, R"mydelimiter(
            ``rcr.Priors`` *object*. Object describing parameter prior probability distribution(s) applied to :class:`rcr.FunctionalForm` model (see :class:`rcr.Priors`).

            To use priors on model parameters for some :class:`rcr.FunctionalForm` model, this attribute of the model needs to be initialized as some instance of :class:`rcr.Priors`
            (see :ref:`priors`).
        )mydelimiter")
        .def_readwrite("pivot_function", &FunctionalForm::pivotFunc, R"mydelimiter(
            Function used to evaluate pivot point(s) (see ``pivot_function`` optional argument of :class:`rcr.FunctionalForm` model constructor).
        )mydelimiter");

    // parameter prior probability distribution types
    py::enum_<priorTypes>(m, "priorsTypes", py::arithmetic(), "Types of prior probability density functions that can be applied to model parameters.")
        .value("CUSTOM_PRIORS", CUSTOM_PRIORS, "Custom, function-defined prior probability density functions(s).")
        .value("GAUSSIAN_PRIORS", GAUSSIAN_PRIORS, "Gaussian (normal) prior probability density function(s).")
        .value("CONSTRAINED_PRIORS", CONSTRAINED_PRIORS, "Bounded/hard-constrained prior probability density function(s).")
        .value("MIXED_PRIORS", MIXED_PRIORS, "A mixture of gaussian (normal), hard-constrained, and uninformative (uniform/flat) prior probability density functions.")
        .export_values();

    // parameter prior probability distribution class
    py::class_<Priors>(m, "Priors", R"mydelimiter(
            *class*. Class that encapsulates probabalistic priors to be applied to model parameters when using model-fitting/functional form RCR (see :ref:`priors` for an example).

            Constructor arguments:

            Parameters
            ----------
            priorType : :class:`rcr.priorsTypes`
                The type of priors that you're applying to your model (see :class:`rcr.priorsTypes` and :ref:`priorstypes`).

            p : function, optional
                Custom priors function; takes in a vector of model parameters and returns a vector of the prior probability density for each value (see :ref:`priors` for an example).
            
            gaussianParams : 2D list/array_like, optional 2nd argument
                A list that contains lists of mu and sigma for the Gaussian prior of each param. If no prior, then just use NaNs (see :ref:`priors` for an example).
            
            paramBounds : 2D list/array_like, optional 2nd argument (or 3rd, for the case of ``rcr.MIXED_PRIORS``)
                A list that contains lists of the lower and upper hard bounds of each param. If not bounded, use NaNs, and if there's only one bound, use NaN for the other bound (see :ref:`priors` for an example).
        )mydelimiter")
        // constructors
        .def(py::init< priorTypes, std::function <std::vector <double>(std::vector <double>)> >()) // custom priors
        .def(py::init< priorTypes, std::vector < std::vector <double> > >()) // only Gaussian or only bounded/hard constraints
        .def(py::init< priorTypes, std::vector < std::vector <double> >, std::vector < std::vector <double> > >()) // mixed priors
        .def(py::init<>())

        // members
        .def_readwrite("priorType", &Priors::priorType, R"mydelimiter(
            ``rcr.priorsTypes`` *object*. The type of priors that you're applying to your model (see :class:`rcr.priorsTypes` and :ref:`priorstypes`).
        )mydelimiter")
        .def_readwrite("p", &Priors::p, R"mydelimiter(
            *function*. Custom priors function; takes in a vector of model parameters and returns a vector of the prior probability density for each value (see :ref:`priors` for an example).
        )mydelimiter")
        .def_readwrite("gaussianParams", &Priors::gaussianParams, R"mydelimiter(
            *2D list/array_like*. A list that contains lists of mu and sigma for the Gaussian prior of each param. If no prior, then just use NaNs (see :ref:`priors` for an example).
        )mydelimiter")
        .def_readwrite("paramBounds", &Priors::paramBounds, R"mydelimiter(
            A list that contains lists of the lower and upper hard bounds of each param. If not bounded, use NaNs, and if there's only one bound, use NaN for the other bound (see :ref:`priors` for an example).
        )mydelimiter");
}