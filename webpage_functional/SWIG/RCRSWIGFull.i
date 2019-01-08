/* File: RCRSWIGFull.i */
%module RCRSWIGFull

%{
#define SWIG_FILE_WITH_INIT
#include "RCRSWIGFull.h"
%}

%include "std_vector.i"

%include "RCRSWIGFull.h"

/* creates DoubleVector object */
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
}

std::vector <double> requestHandlerUnWeighted(std::vector <double> x, std::vector <double> y, std::vector <double> sigma_y, std::vector <double> guess, int fType, int dataSize, int rejTechNo);
std::vector <double> requestHandlerWeighted(std::vector <double> x, std::vector <double> y, std::vector <double> sigma_y, std::vector <double> guess, std::vector <double> w, int fType, int dataSize, int rejTechNo);
