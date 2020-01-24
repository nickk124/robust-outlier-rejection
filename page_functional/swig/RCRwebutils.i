/* File: RCRwebutils.i */
%module RCRwebutils

%{
#define SWIG_FILE_WITH_INIT
#include "RCRwebutils.h"
%}

%include "std_vector.i"

%include "RCRwebutils.h"
%include "RCR.h"
%include "FunctionalForm.h"
%include "MiscFunctions.h"
%include "NonParametric.h"

/* creates DoubleVector object */
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
}

std::vector <double> requestHandlerUnWeighted(std::vector <double> x, std::vector <double> y, std::vector <double> guess, int fType, int dataSize, int rejTechNo, int priorsCheck, std::vector <double> priorsParams, std::vector <int> hasPriorsVec);
std::vector <double> requestHandlerWeighted(std::vector <double> x, std::vector <double> y, std::vector <double> guess, std::vector <double> w, int fType, int dataSize, int rejTechNo, int priorsCheck, std::vector <double> priorsParams, std::vector <int> hasPriorsVec);