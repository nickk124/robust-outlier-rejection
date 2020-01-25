/* RCRLib.i */
%module RCRLib
%include <std_vector.i>

%{
#define SWIG_FILE_WITH_INIT
#include "RCR.h"
%}

%include "std_vector.i"

%include "RCR.h"
%include "NonParametric.h"
%include "FunctionalForm.h"
%include "MiscFunctions.h"

namespace std {
 %template(DubVec) vector<double>;
 %template(BoolVec) vector<bool>;

}
