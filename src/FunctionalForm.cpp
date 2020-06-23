/*
 Robust Chauvenet Rejection (RCR) Official Codebase
 Active Author: Nick C. Konz
 Former Author: Michael Maples
 See license at https://github.com/nickk124/RCR
 */
#include "FunctionalForm.h"
//std::vector <double>(*p)(std::vector <double>)

double pivot; //see declaration in FunctionalForm.h for explanation as to why this needs to be global.
std::vector <double> pivot_ND;

// Constructors
Priors::Priors(priorTypes priorType, std::function <std::vector <double>(std::vector <double>, std::vector <double>)> p) { //custom priors
	this->priorType = priorType;
	this->p = p;
};

Priors::Priors(priorTypes priorType, std::vector < std::vector <double> > params) { //Gaussian or bounded priors only
	this->priorType = priorType;
	if (priorType == GAUSSIAN_PRIORS) {
		this->gaussianParams = params;
	}
	else if (priorType == CONSTRAINED_PRIORS) {
		this->paramBounds = params;
	}
};

Priors::Priors(priorTypes priorType, std::vector < std::vector <double> > gaussianParams, std::vector < std::vector <double> > paramBounds) { //mixed
	this->priorType = priorType;
	this->paramBounds = paramBounds;
	this->gaussianParams = gaussianParams;
};

//default constructor
Priors::Priors() {

};

FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess)
{
	this->f = f;
	this->x = x;
	this->y = y;
	this->guess = guess;
	this->sigma_y = sigma_y;
	this->w.resize(x.size(), 1.0);
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = true;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector<double> w)
{
	this->f = f;
	this->x = x;
	this->y = y;
	this->w = w;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = true;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess)
{
	this->f_ND = f_ND;
	this->x_ND = x_ND;
	this->y = y;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = true;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector<double> w)
{
	this->f_ND = f_ND;
	this->x_ND = x_ND;
	this->y = y;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->w = w;
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = true;
	this->hasCustomPivot = false;
}

FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess)
{
	this->f = f;
	this->x = x;
	this->y = y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = false;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector<double> w)
{
	this->f = f;
	this->x = x;
	this->y = y;
	this->w = w;
	this->guess = guess;
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = false;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess)
{
	this->f_ND = f_ND;
	this->x_ND = x_ND;
	this->y = y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = false;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector<double> w)
{
	this->f_ND = f_ND;
	this->x_ND = x_ND;
	this->y = y;
	this->guess = guess;
	this->w = w;
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = false;
	this->hasCustomPivot = false;
}
//constructors with priors:
FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, Priors priorsObject)
{
	this->f = f;
	this->priors = priorsObject;
	this->x = x;
	this->y = y;
	this->guess = guess;
	this->sigma_y = sigma_y;
	this->w.resize(x.size(), 1.0);
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = true;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, Priors priorsObject)
{
	this->f = f;
	this->priors = priorsObject;
	this->x = x;
	this->y = y;
	this->w = w;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = true;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, Priors priorsObject)
{
	this->f_ND = f_ND;
	this->priors = priorsObject;
	this->x_ND = x_ND;
	this->y = y;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = true;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, Priors priorsObject)
{
	this->f_ND = f_ND;
	this->priors = priorsObject;
	this->x_ND = x_ND;
	this->y = y;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->w = w;
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = true;
	this->hasCustomPivot = false;
}

FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, Priors priorsObject)
{
	this->f = f;
	this->priors = priorsObject;
	this->x = x;
	this->y = y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = false;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, Priors priorsObject)
{
	this->f = f;
	this->priors = priorsObject;
	this->x = x;
	this->y = y;
	this->w = w;
	this->guess = guess;
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = false;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, Priors priorsObject)
{
	this->f_ND = f_ND;
	this->priors = priorsObject;
	this->x_ND = x_ND;
	this->y = y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = false;
	this->hasCustomPivot = false;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, Priors priorsObject)
{
	this->f_ND = f_ND;
	this->priors = priorsObject;
	this->x_ND = x_ND;
	this->y = y;
	this->guess = guess;
	this->w = w;
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = false;
	this->hasCustomPivot = false;
}

//support for custom pivot/average value:

FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess)
{
	this->f = f;
	this->x = x;
	this->y = y;
	this->guess = guess;
	this->sigma_y = sigma_y;
	this->w.resize(x.size(), 1.0);
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = true;
	this->hasCustomPivot = true;
	this->pivotFunc = pivotFunc;

	pivot = pivot_guess;
}
FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess)
{
	this->f = f;
	this->x = x;
	this->y = y;
	this->w = w;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = true;
	this->hasCustomPivot = true;
	this->pivotFunc = pivotFunc;

	pivot = pivot_guess;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess)
{
	this->f_ND = f_ND;
	this->x_ND = x_ND;
	this->y = y;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = true;
	this->hasCustomPivot = true;
	this->pivotFunc_ND = pivotFunc_ND;

	pivot_ND = pivot_ND_guess;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess)
{
	this->f_ND = f_ND;
	this->x_ND = x_ND;
	this->y = y;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->w = w;
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = true;
	this->hasCustomPivot = true;
	this->pivotFunc_ND = pivotFunc_ND;

	pivot_ND = pivot_ND_guess;
}

FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess)
{
	this->f = f;
	this->x = x;
	this->y = y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = false;
	this->hasCustomPivot = true;
	this->pivotFunc = pivotFunc;

	pivot = pivot_guess;
}
FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess)
{
	this->f = f;
	this->x = x;
	this->y = y;
	this->w = w;
	this->guess = guess;
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = false;
	this->hasCustomPivot = true;
	this->pivotFunc = pivotFunc;

	pivot = pivot_guess;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess)
{
	this->f_ND = f_ND;
	this->x_ND = x_ND;
	this->y = y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = false;
	this->hasCustomPivot = true;
	this->pivotFunc_ND = pivotFunc_ND;

	pivot_ND = pivot_ND_guess;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess)
{
	this->f_ND = f_ND;
	this->x_ND = x_ND;
	this->y = y;
	this->guess = guess;
	this->w = w;
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = false;
	this->hasErrorBars = false;
	this->hasCustomPivot = true;
	this->pivotFunc_ND = pivotFunc_ND;

	pivot_ND = pivot_ND_guess;
}
//constructors with priors:
FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, Priors priorsObject, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess)
{
	this->f = f;
	this->priors = priorsObject;
	this->x = x;
	this->y = y;
	this->guess = guess;
	this->sigma_y = sigma_y;
	this->w.resize(x.size(), 1.0);
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = true;
	this->hasCustomPivot = true;
	this->pivotFunc = pivotFunc;

	pivot = pivot_guess;
}
FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, Priors priorsObject, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess)
{
	this->f = f;
	this->priors = priorsObject;
	this->x = x;
	this->y = y;
	this->w = w;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = true;
	this->hasCustomPivot = true;
	this->pivotFunc = pivotFunc;

	pivot = pivot_guess;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, Priors priorsObject, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess)
{
	this->f_ND = f_ND;
	this->priors = priorsObject;
	this->x_ND = x_ND;
	this->y = y;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = true;
	this->hasCustomPivot = true;
	this->pivotFunc_ND = pivotFunc_ND;

	pivot_ND = pivot_ND_guess;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, Priors priorsObject, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess)
{
	this->f_ND = f_ND;
	this->priors = priorsObject;
	this->x_ND = x_ND;
	this->y = y;
	this->sigma_y = sigma_y;
	this->guess = guess;
	this->w = w;
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = true;
	this->hasCustomPivot = true;
	this->pivotFunc_ND = pivotFunc_ND;

	pivot_ND = pivot_ND_guess;
}

FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, Priors priorsObject, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess)
{
	this->f = f;
	this->priors = priorsObject;
	this->x = x;
	this->y = y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = false;
	this->hasCustomPivot = true;
	this->pivotFunc = pivotFunc;

	pivot = pivot_guess;
}
FunctionalForm::FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector<double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, Priors priorsObject, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess)
{
	this->f = f;
	this->priors = priorsObject;
	this->x = x;
	this->y = y;
	this->w = w;
	this->guess = guess;
	this->partialsvector = partialsvector;
	this->M = (int) partialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = false;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = false;
	this->hasCustomPivot = true;
	this->pivotFunc = pivotFunc;

	pivot = pivot_guess;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, Priors priorsObject, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess)
{
	this->f_ND = f_ND;
	this->priors = priorsObject;
	this->x_ND = x_ND;
	this->y = y;
	this->guess = guess;
	this->w.resize(x.size(), 1.0);
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = false;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = false;
	this->hasCustomPivot = true;
	this->pivotFunc_ND = pivotFunc_ND;

	pivot_ND = pivot_ND_guess;
}
FunctionalForm::FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector <std::vector<double> > x_ND, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector<double> w, Priors priorsObject, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess)
{
	this->f_ND = f_ND;
	this->priors = priorsObject;
	this->x_ND = x_ND;
	this->y = y;
	this->guess = guess;
	this->w = w;
	this->NDpartialsvector = NDpartialsvector;
	this->M = (int) NDpartialsvector.size();
	this->N = (int) y.size();
	this->tolerance = tolerance;
	this->NDcheck = true;
	this->weightedCheck = true;
	this->parameters = guess; //initial 
	this->has_priors = true;
	this->hasErrorBars = false;
	this->hasCustomPivot = true;
	this->pivotFunc_ND = pivotFunc_ND;

	pivot_ND = pivot_ND_guess;
}


//default constructor:
FunctionalForm::FunctionalForm() {

}

/*
void FunctionalForm::initializeBar() {
	if (hasCustomPivot) {
		std::vector <double> w_temp;
		for (int i = 0; i < x.size(); i++) {
			if (!weightedCheck) {
				w_temp.push_back(1.0);
			}
			else {
				w_temp.push_back(w[i]);
			}
		}

		if (NDcheck) { //computes the pivot value for this iteration.
			pivot_ND = pivotFunc_ND(x_ND, w_temp, f_ND, parameters);
		}
		else if (!NDcheck) {

			pivot = pivotFunc(x, w_temp, f, parameters);
		}
	}
};
*/

void FunctionalForm::setTrueVec(std::vector<bool> &flags, std::vector<double> &w, std::vector<double> &y)
{

	int trueCount = 0, currentIndex;
	std::vector<int> indicesVec;
	std::vector<double> trueWVec, trueYVec;
	this->flags = flags;

	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueCount += 1;
		}
	}
	trueWVec.resize(trueCount);
	trueYVec.resize(trueCount);
	indicesVec.resize(trueCount);
	currentIndex = 0;
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueWVec[currentIndex] = (w[i]);
			trueYVec[currentIndex] = (y[i]);
			indicesVec[currentIndex] = i;

			currentIndex += 1;
		}
	}
	trueY = trueYVec;
	trueW = trueWVec;
//	this->trueW = trueW;
	indices = indicesVec;
}
void FunctionalForm::setTrueVec(std::vector<bool> &flags, std::vector<double> &y)
{
	int trueCount = 0, currentIndex;
	std::vector<int> indicesVec;
	std::vector<double> trueYVec;
	this->flags = flags;
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueCount += 1;
		}
	}
	trueYVec.resize(trueCount);
	indicesVec.resize(trueCount);
	currentIndex = 0;
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueYVec[currentIndex] = (y[i]);
			indicesVec[currentIndex] = i;
			currentIndex += 1;
		}
	}
	trueY = trueYVec;
	indices = indicesVec;

}
void FunctionalForm::buildModelSpace()
{
	if (hasCustomPivot) {
		std::vector <double> truex;
		std::vector < std::vector <double> > truexND;
		std::vector <double> truew;

		if (NDcheck) { //computes the pivot value for this iteration.
			for (int i = 0; i < N; i++) {
				if (flags[i]) {
					truexND.push_back(x_ND[i]);
					if (!weightedCheck) {
						truew.push_back(1.0);
					}
					else {
						truew.push_back(w[i]);
					}
				}
			}
			pivot_ND = pivotFunc_ND(truexND, truew, f_ND, parameters);
			result.pivot_ND = pivot_ND;
		}
		else if (!NDcheck) {
			for (int i = 0; i < N; i++) {
				if (flags[i]) {
					truex.push_back(x[i]);
					if (!weightedCheck) {
						truew.push_back(1.0);
					}
					else {
						truew.push_back(w[i]);
					}
				}
			}
			pivot = pivotFunc(truex, truew, f, parameters);
			result.pivot = pivot;
		}
	}

	parameterSpace.clear();
	weightSpace.clear();
	extraParameterSpace.clear();
	extraWeightSpace.clear();
	combos.clear();
	NDcombos.clear();
	combos_indices.clear();
	combosgood_indices.clear();

	//initializes parameterSpace and weightSpace:
	for (int i = 0; i < M; i++) {
		parameterSpace.push_back(innerSpace);
		weightSpace.push_back(innerSpace);
		extraWeightSpace.push_back(innerSpace);
		extraParameterSpace.push_back(innerSpace);
	}

	//std::cout << "Calculating all non-rejected data M-combinations for parameter calculation..." << std::endl;

	if (NDcheck == false) { //creates the vector that has also possible M-combinations of data-points, and a vector of the corresponding indices of these
		getCombos(x, M, 0);
	}
	else if (NDcheck) {
		getCombos(x_ND, M, 0); //ND case
	}
	if (weightedCheck) { //calculates wbar -- average weight of all unrejected data points
		double wsum = 0.0;
		double goodcount = 0.0;
		for (int j = 0; j < N; j++) {
			if (flags[j]) {
				wsum += w[j];
				goodcount += 1.0;
			}
		}
		wbar = wsum / goodcount;
	}
	if (NDcheck == false) {
		for (int i = 0; i < combos.size(); i++) //the following builds the vector of good-flagged combinations, and the vector of their corresponding indices.
		{
			bool check = true;
			for (int j = 0; j < M; j++) {
				if (flags[combos_indices[i][j]] == false) {
					check = false;
				}
			}
			if (check) {
				combosgood_indices.push_back(combos_indices[i]);
			}
		}
	}
	else if (NDcheck) {
		for (int i = 0; i < NDcombos.size(); i++) //the following builds the vector of good-flagged combinations, and the vector of their corresponding indices.
		{
			bool check = true;
			for (int j = 0; j < M; j++) {
				if (flags[combos_indices[i][j]] == false) {
					check = false;
				}
			}
			if (check) {
				combosgood_indices.push_back(combos_indices[i]);
			}
		}
	}
	// the below limits the number of combos used if the criteria from the paper is met (chooses random eigthed draws instead of directly using all of combosgood_indices
	int combolimit = 20000;
	double frac = 0.5;
	int combosgoodcount = (int) combosgood_indices.size();

	std::vector <double> combosweights;
	std::vector <std::vector <int> > chosencomboindices; // (the indices of) combos that will be used, drawn randomly with weights
	double totalcomboweight;

	if ((frac * combosgoodcount) > combolimit)
	{

		
		///std::cout << "Total parameter combination count is very large, switching to weighted random draws..." << std::endl;
		totalcomboweight = 0.0;
		for (int i = 0; i < combosgoodcount; i++) // ith combo; computes weights of each combo
		{
			double comboweight = 0.0;
			for (int j = 0; j < M; j++) {
				comboweight += w[combosgood_indices[i][j]];
			}
			totalcomboweight += comboweight;
			combosweights.push_back(comboweight);
		}
		// drawing the random combos:

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(0, totalcomboweight); // generates psuedo-random numbers between 0 and the sum of weights

																   // doing all random draws:
		for (int j = 0; j < combolimit; j++) //each iteration of this takes ~0.1 seconds 
		{
			double randomnum = dis(gen);
			for (int k = 0; k < combosgoodcount; k++)
			{
				randomnum -= combosweights[k];
				if (randomnum < combosweights[k]) {
					chosencomboindices.push_back(combosgood_indices[k]); // adds the drawn combo (indices) to big vector

					/*
					for (int i = 0; i < combosgood_indices[k].size(); i++) {
					std::cout << combosgood_indices[k][i] << "  ";
					}
					std::cout << std::endl;
					*/
					break;
				}
			}
		}

		//getting rid of duplicates:
		std::set <std::vector <int> > s;
		unsigned size = (int) chosencomboindices.size();
		for (unsigned i = 0; i < size; ++i) s.insert(chosencomboindices[i]);
		chosencomboindices.assign(s.begin(), s.end());

		/*
		testsum = 0;

		for (int i = 0; i < chosencomboindices.size(); i++) {
		for (int j = i + 1; j < chosencomboindices.size(); j++) {
		if ((chosencomboindices[i] == chosencomboindices[j]) && (i != j)) {
		testsum += 1;
		}
		}
		}


		std::cout << std::endl << testsum << " duplicates out of 20,000 random draws" << std::endl;
		*/
		//std::cout << combolimit << " random weighted draws done from " << combosgoodcount << " combinations total..." << std::endl;

		combosgood_indices = chosencomboindices; //now, the random draws will be used.
	}

	//JUST FOR TESTING:
	/*
	int testsum = 0;

	for (int i = 0; i < chosencomboindices.size(); i++) {
	for (int j = i + 1; j < chosencomboindices.size(); j++) {
	if ((chosencomboindices[i] == chosencomboindices[j]) && (i != j)) {
	testsum += 1;
	}
	}
	}

	std::cout << std::endl << testsum << " duplicates out of 20,000 random draws" << std::endl;

	//above took 1,140 seconds
	*/

	/*
	testsum = 0;

	for (int i = 0; i < chosencomboindices.size(); i++) {
	for (int j = i + 1; j < chosencomboindices.size(); j++) {
	if ((chosencomboindices[i] == chosencomboindices[j]) && (i != j)) {
	testsum += 1;
	}
	}
	}


	std::cout << std::endl << testsum << " duplicates out of 20,000 random draws" << std::endl;
	*/


	std::vector <double> comboy, comboparamset, comboparam_uncertainties;

	if (weightedCheck && NDcheck && hasErrorBars) //weighted and >1 dimension of independent variable in model
	{
		std::vector <double> combosigma_y;
		std::vector <std::vector <double> > combox;
		std::vector <double> combow;
		for (int i = 0; i < combosgood_indices.size(); i++) //using each combination
		{
			std::vector <int> combo_indices = combosgood_indices[i]; //the indices of the combo to be used; initializes the combo to be used
			combox.clear();
			comboy.clear();
			combow.clear();
			combosigma_y.clear();
			comboparamset.clear();
			comboparam_uncertainties.clear();

			for (int j = 0; j < M; j++) {
				combox.push_back(x_ND[combo_indices[j]]);
				comboy.push_back(y[combo_indices[j]]);
				combosigma_y.push_back(sigma_y[combo_indices[j]]);
				combow.push_back(w[combo_indices[j]]);
			}


			comboparamset = modifiedGN(f_ND, NDpartialsvector, comboy, combox, parameters, combosigma_y, tolerance, combow); //guess is part of the FunctionalForm Constructor


																												   //next, checks for exceptions

			if (comboparamset.size() == M) //no exceptions triggered
			{
				comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combosigma_y, combow, wbar);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += combow[i] / std::pow(combosigma_y[i], 2.0);
				}

				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						parameterSpace[f].push_back(comboparamset[f]);
					}
					else {
						weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
					}
				}
			}
			else if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
			{
				std::vector <double> semigoodparamvec(M, 0.0);

				for (int j = 0; j < M; j++) {
					semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
				}

				comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combosigma_y, combow, wbar);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += combow[i] / std::pow(combosigma_y[i], 2.0);
				}
				for (int k = 0; k < M; k++) {
					extraParameterSpace[k].push_back(semigoodparamvec[k]);
				}
				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
					}
					else {
						extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
					}
				}
			}
			else if (comboparamset.size() == (M + 2))
			{
				comboparamset = regularGN(f_ND, NDpartialsvector, comboy, combox, parameters, tolerance, combow); //uses regular GN

				if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
				{

					//std::cout << comboparamset[0] << "  " << comboparamset[1] << std::endl;

					std::vector <double> semigoodparamvec(M, 0.0);

					for (int j = 0; j < M; j++) {
						semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
					}

					comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combosigma_y, combow, wbar);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += combow[i] / std::pow(combosigma_y[i], 2.0);
					}
					for (int k = 0; k < M; k++) {
						extraParameterSpace[k].push_back(semigoodparamvec[k]);
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						}
						else {
							extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						}
					}
				}
				else
				{
					comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combosigma_y, combow, wbar);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += combow[i] / std::pow(combosigma_y[i], 2.0);
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
							parameterSpace[f].push_back(comboparamset[f]);
						}
						else {
							weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
							parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
						}
					}
				}
			}
			// otherwise, singlar GN matrix issue due to data; excluded from all calculations


		}
	}
	else if (weightedCheck && !NDcheck && hasErrorBars) //weighted and only 1 dimension of independent variable in model
	{
		std::vector <double> combosigma_y;
		std::vector <double> combox, combow;
		for (int i = 0; i < combosgood_indices.size(); i++) //using each combination
		{
			std::vector <int> combo_indices = combosgood_indices[i]; //the indices of the combo to be used; initializes the combo to be used
			combox.clear();
			comboy.clear();
			combosigma_y.clear();
			comboparamset.clear();
			comboparam_uncertainties.clear();
			combow.clear();


			for (int j = 0; j < M; j++) {
				combox.push_back(x[combo_indices[j]]);
				comboy.push_back(y[combo_indices[j]]);
				combosigma_y.push_back(sigma_y[combo_indices[j]]);
				combow.push_back(w[combo_indices[j]]);
			}

			comboparamset = modifiedGN(f, partialsvector, comboy, combox, parameters, combosigma_y, tolerance, combow); //guess is part of the FunctionalForm Constructor

																											  //next, checks for exceptions

			if (comboparamset.size() == M) //no exceptions triggered
			{
				comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combosigma_y, combow, wbar);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += combow[i] / std::pow(combosigma_y[i], 2.0);
				}

				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						parameterSpace[f].push_back(comboparamset[f]);
					}
					else {
						weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
					}
				}
			}
			else if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
			{
				std::vector <double> semigoodparamvec(M, 0.0);

				for (int j = 0; j < M; j++) {
					semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
				}

				comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combosigma_y, combow, wbar);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += combow[i] / std::pow(combosigma_y[i], 2.0);
				}
				for (int k = 0; k < M; k++) {
					extraParameterSpace[k].push_back(semigoodparamvec[k]);
				}
				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
					}
					else {
						extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
					}
				}
			}
			else if (comboparamset.size() == (M + 2))
			{
				comboparamset = regularGN(f, partialsvector, comboy, combox, parameters, tolerance, combow); //uses regular GN

				if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
				{

					//std::cout << comboparamset[0] << "  " << comboparamset[1] << std::endl;

					std::vector <double> semigoodparamvec(M, 0.0);

					for (int j = 0; j < M; j++) {
						semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
					}

					comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combosigma_y, combow, wbar);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += combow[i] / std::pow(combosigma_y[i], 2.0);
					}
					for (int k = 0; k < M; k++) {
						extraParameterSpace[k].push_back(semigoodparamvec[k]);
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						}
						else {
							extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						}
					}
				}
				else
				{
					comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combosigma_y, combow, wbar);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += combow[i] / std::pow(combosigma_y[i], 2.0);
					}

					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
							parameterSpace[f].push_back(comboparamset[f]);
						}
						else {
							weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
							parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
						}
					}
				}
			}
			// otherwise, singlar GN matrix issue due to data; excluded from all calculations
		}
	}
	else if (!weightedCheck && NDcheck && hasErrorBars) //non-weighted but >1 dimension of independent variable in model
	{
		std::vector <double> combosigma_y;
		std::vector <std::vector <double> > combox;
		for (int i = 0; i < combosgood_indices.size(); i++) //using each combination
		{
			std::vector <int> combo_indices = combosgood_indices[i]; //the indices of the combo to be used; initializes the combo to be used
			combox.clear();
			comboy.clear();
			combosigma_y.clear();
			comboparamset.clear();
			comboparam_uncertainties.clear();

			for (int j = 0; j < M; j++) {
				combox.push_back(x_ND[combo_indices[j]]);
				comboy.push_back(y[combo_indices[j]]);
				combosigma_y.push_back(sigma_y[combo_indices[j]]);
			}


			comboparamset = modifiedGN(f_ND, NDpartialsvector, comboy, combox, parameters, combosigma_y, tolerance); //guess is part of the FunctionalForm Constructor

																												//next, checks for exceptions
			if (comboparamset.size() == M) //no exceptions triggered
			{
				comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combosigma_y);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += 1.0 / std::pow(combosigma_y[i], 2.0);
				}

				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						parameterSpace[f].push_back(comboparamset[f]);
					}
					else {
						weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
					}
				}
			}
			else if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
			{
				std::vector <double> semigoodparamvec(M, 0.0);

				for (int j = 0; j < M; j++) {
					semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
				}

				comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combosigma_y);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += 1.0 / std::pow(combosigma_y[i], 2.0);
				}
				for (int k = 0; k < M; k++) {
					extraParameterSpace[k].push_back(semigoodparamvec[k]);
				}
				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
					}
					else {
						extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
					}
				}
			}
			else if (comboparamset.size() == (M + 2))
			{
				comboparamset = regularGN(f_ND, NDpartialsvector, comboy, combox, parameters, tolerance); //uses regular GN

				if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
				{

					//std::cout << comboparamset[0] << "  " << comboparamset[1] << std::endl;

					std::vector <double> semigoodparamvec(M, 0.0);

					for (int j = 0; j < M; j++) {
						semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
					}

					comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combosigma_y);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += 1.0 / std::pow(combosigma_y[i], 2.0);
					}

					for (int k = 0; k < M; k++) {
						extraParameterSpace[k].push_back(semigoodparamvec[k]);
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						}
						else {
							extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						}
					}
				}
				else
				{
					comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combosigma_y);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += 1.0 / std::pow(combosigma_y[i], 2.0);
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
							parameterSpace[f].push_back(comboparamset[f]);
						}
						else {
							weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
							parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
						}
					}
				}
			}
			// otherwise, singlar GN matrix issue due to data; excluded from all calculations
		}
	}
	else if (!weightedCheck && !NDcheck && hasErrorBars) //non-weighted, 1 dimension of independent variable in model
	{
		std::vector <double> combosigma_y;
		std::vector <double> combox;
		for (int i = 0; i < combosgood_indices.size(); i++) //using each combination
		{
			std::vector <int> combo_indices = combosgood_indices[i]; //the indices of the combo to be used; initializes the combo to be used
			combox.clear();
			comboy.clear();
			combosigma_y.clear();
			comboparamset.clear();
			comboparam_uncertainties.clear();

			for (int j = 0; j < M; j++) {
				combox.push_back(x[combo_indices[j]]);
				comboy.push_back(y[combo_indices[j]]);
				combosigma_y.push_back(sigma_y[combo_indices[j]]);
			}

			comboparamset = modifiedGN(f, partialsvector, comboy, combox, parameters, combosigma_y, tolerance); //parameters is part of the FunctionalForm Constructor			


			//next, checks for exceptions
			if (comboparamset.size() == M) //no exceptions triggered
			{

				comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combosigma_y);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += 1.0 / std::pow(combosigma_y[i], 2.0);
				}

				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						parameterSpace[f].push_back(comboparamset[f]);
					}
					else {
						// FOR TESTING



						weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
					}
				}
			}
			else if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
			{
				//std::cout << comboparamset[0] << "  " << comboparamset[1] << std::endl;

				std::vector <double> semigoodparamvec(M, 0.0);

				for (int j = 0; j < M; j++) {
					semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
				}

				comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combosigma_y);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += 1.0 / std::pow(combosigma_y[i], 2.0);
				}

				for (int k = 0; k < M; k++) {
					extraParameterSpace[k].push_back(semigoodparamvec[k]);
				}
				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
					}
					else {
						extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
					}
				}
			}
			else if (comboparamset.size() == (M + 2)) 
			{
				comboparamset = regularGN(f, partialsvector, comboy, combox, parameters, tolerance); //uses regular GN

				if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
				{

					//std::cout << comboparamset[0] << "  " << comboparamset[1] << std::endl;

					std::vector <double> semigoodparamvec(M, 0.0);

					for (int j = 0; j < M; j++) {
						semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
					}

					comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combosigma_y);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += 1.0 / std::pow(combosigma_y[i], 2.0);
					}
					for (int k = 0; k < M; k++) {
						extraParameterSpace[k].push_back(semigoodparamvec[k]);
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						}
						else {
							extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						}
					}
				}
				else
				{
					comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combosigma_y);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += 1.0 / std::pow(combosigma_y[i], 2.0);
					}

					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
							parameterSpace[f].push_back(comboparamset[f]);
						}
						else {
							weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
							parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
						}
					}
				}
			}
			// otherwise, singlar GN matrix issue due to data; excluded from all calculations
		}
	}

	if (weightedCheck && NDcheck && !hasErrorBars) //weighted and >1 dimension of independent variable in model
	{
		std::vector <std::vector <double> > combox;
		std::vector <double> combow;
		for (int i = 0; i < combosgood_indices.size(); i++) //using each combination
		{
			std::vector <int> combo_indices = combosgood_indices[i]; //the indices of the combo to be used; initializes the combo to be used
			combox.clear();
			comboy.clear();
			combow.clear();
			comboparamset.clear();
			comboparam_uncertainties.clear();

			for (int j = 0; j < M; j++) {
				combox.push_back(x_ND[combo_indices[j]]);
				comboy.push_back(y[combo_indices[j]]);
				combow.push_back(w[combo_indices[j]]);
			}


			comboparamset = modifiedGN(f_ND, NDpartialsvector, comboy, combox, parameters, tolerance, combow); //guess is part of the FunctionalForm Constructor


																												   //next, checks for exceptions

			if (comboparamset.size() == M) //no exceptions triggered
			{
				comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combow, wbar);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += combow[i];
				}

				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						parameterSpace[f].push_back(comboparamset[f]);
					}
					else {
						weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
					}
				}
			}
			else if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
			{
				std::vector <double> semigoodparamvec(M, 0.0);

				for (int j = 0; j < M; j++) {
					semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
				}

				comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combow, wbar);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += combow[i];
				}
				for (int k = 0; k < M; k++) {
					extraParameterSpace[k].push_back(semigoodparamvec[k]);
				}
				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
					}
					else {
						extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
					}
				}
			}
			else if (comboparamset.size() == (M + 2))
			{
				comboparamset = regularGN(f_ND, NDpartialsvector, comboy, combox, parameters, tolerance, combow); //uses regular GN

				if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
				{

					//std::cout << comboparamset[0] << "  " << comboparamset[1] << std::endl;

					std::vector <double> semigoodparamvec(M, 0.0);

					for (int j = 0; j < M; j++) {
						semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
					}

					comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combow, wbar);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += combow[i];
					}
					for (int k = 0; k < M; k++) {
						extraParameterSpace[k].push_back(semigoodparamvec[k]);
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						}
						else {
							extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						}
					}
				}
				else
				{
					comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters, combow, wbar);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += combow[i];
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
							parameterSpace[f].push_back(comboparamset[f]);
						}
						else {
							weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
							parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
						}
					}
				}
			}
			// otherwise, singlar GN matrix issue due to data; excluded from all calculations


		}
	}
	else if (weightedCheck && !NDcheck && !hasErrorBars) //weighted and only 1 dimension of independent variable in model
	{
		std::vector <double> combox, combow;
		for (int i = 0; i < combosgood_indices.size(); i++) //using each combination
		{
			std::vector <int> combo_indices = combosgood_indices[i]; //the indices of the combo to be used; initializes the combo to be used
			combox.clear();
			comboy.clear();
			comboparamset.clear();
			comboparam_uncertainties.clear();
			combow.clear();


			for (int j = 0; j < M; j++) {
				combox.push_back(x[combo_indices[j]]);
				comboy.push_back(y[combo_indices[j]]);
				combow.push_back(w[combo_indices[j]]);
			}

			comboparamset = modifiedGN(f, partialsvector, comboy, combox, parameters, tolerance, combow); //guess is part of the FunctionalForm Constructor

																											  //next, checks for exceptions

			if (comboparamset.size() == M) //no exceptions triggered
			{
				comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combow, wbar);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += combow[i];
				}

				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						parameterSpace[f].push_back(comboparamset[f]);
					}
					else {
						weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
					}
				}
			}
			else if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
			{
				std::vector <double> semigoodparamvec(M, 0.0);

				for (int j = 0; j < M; j++) {
					semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
				}

				comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combow, wbar);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += combow[i];
				}
				for (int k = 0; k < M; k++) {
					extraParameterSpace[k].push_back(semigoodparamvec[k]);
				}
				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
					}
					else {
						extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
					}
				}
			}
			else if (comboparamset.size() == (M + 2))
			{
				comboparamset = regularGN(f, partialsvector, comboy, combox, parameters, tolerance, combow); //uses regular GN

				if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
				{

					//std::cout << comboparamset[0] << "  " << comboparamset[1] << std::endl;

					std::vector <double> semigoodparamvec(M, 0.0);

					for (int j = 0; j < M; j++) {
						semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
					}

					comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combow, wbar);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += combow[i];
					}
					for (int k = 0; k < M; k++) {
						extraParameterSpace[k].push_back(semigoodparamvec[k]);
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						}
						else {
							extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						}
					}
				}
				else
				{
					comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters, combow, wbar);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += combow[i];
					}

					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
							parameterSpace[f].push_back(comboparamset[f]);
						}
						else {
							weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
							parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
						}
					}
				}
			}
			// otherwise, singlar GN matrix issue due to data; excluded from all calculations
		}
	}
	else if (!weightedCheck && NDcheck && !hasErrorBars) //non-weighted but >1 dimension of independent variable in model
	{
		std::vector <std::vector <double> > combox;
		for (int i = 0; i < combosgood_indices.size(); i++) //using each combination
		{
			std::vector <int> combo_indices = combosgood_indices[i]; //the indices of the combo to be used; initializes the combo to be used
			combox.clear();
			comboy.clear();
			comboparamset.clear();
			comboparam_uncertainties.clear();

			for (int j = 0; j < M; j++) {
				combox.push_back(x_ND[combo_indices[j]]);
				comboy.push_back(y[combo_indices[j]]);
			}


			comboparamset = modifiedGN(f_ND, NDpartialsvector, comboy, combox, parameters, tolerance); //guess is part of the FunctionalForm Constructor

																												//next, checks for exceptions
			if (comboparamset.size() == M) //no exceptions triggered
			{
				comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += 1.0;
				}

				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						parameterSpace[f].push_back(comboparamset[f]);
					}
					else {
						weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
					}
				}
			}
			else if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
			{
				std::vector <double> semigoodparamvec(M, 0.0);

				for (int j = 0; j < M; j++) {
					semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
				}

				comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += 1.0;
				}
				for (int k = 0; k < M; k++) {
					extraParameterSpace[k].push_back(semigoodparamvec[k]);
				}
				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
					}
					else {
						extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
					}
				}
			}
			else if (comboparamset.size() == (M + 2))
			{
				comboparamset = regularGN(f_ND, NDpartialsvector, comboy, combox, parameters, tolerance); //uses regular GN

				if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
				{

					//std::cout << comboparamset[0] << "  " << comboparamset[1] << std::endl;

					std::vector <double> semigoodparamvec(M, 0.0);

					for (int j = 0; j < M; j++) {
						semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
					}

					comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += 1.0;
					}

					for (int k = 0; k < M; k++) {
						extraParameterSpace[k].push_back(semigoodparamvec[k]);
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						}
						else {
							extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						}
					}
				}
				else
				{
					comboparam_uncertainties = paramuncertainty(NDpartialsvector, combox, parameters);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += 1.0;
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
							parameterSpace[f].push_back(comboparamset[f]);
						}
						else {
							weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
							parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
						}
					}
				}
			}
			// otherwise, singlar GN matrix issue due to data; excluded from all calculations
		}
	}
	else if (!weightedCheck && !NDcheck && !hasErrorBars) //non-weighted, 1 dimension of independent variable in model
	{
		std::vector <double> combox;
		for (int i = 0; i < combosgood_indices.size(); i++) //using each combination
		{
			std::vector <int> combo_indices = combosgood_indices[i]; //the indices of the combo to be used; initializes the combo to be used
			combox.clear();
			comboy.clear();
			comboparamset.clear();
			comboparam_uncertainties.clear();

			for (int j = 0; j < M; j++) {
				combox.push_back(x[combo_indices[j]]);
				comboy.push_back(y[combo_indices[j]]);
			}

			comboparamset = modifiedGN(f, partialsvector, comboy, combox, parameters, tolerance); //parameters is part of the FunctionalForm Constructor			


			//next, checks for exceptions
			if (comboparamset.size() == M) //no exceptions triggered
			{

				comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += 1.0;
				}

				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						parameterSpace[f].push_back(comboparamset[f]);
					}
					else {
						// FOR TESTING



						weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
					}
				}
			}
			else if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
			{
				//std::cout << comboparamset[0] << "  " << comboparamset[1] << std::endl;

				std::vector <double> semigoodparamvec(M, 0.0);

				for (int j = 0; j < M; j++) {
					semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
				}

				comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters);

				double correctivesum = 0.0;
				for (int i = 0; i < M; i++) {
					correctivesum += 1.0;
				}

				for (int k = 0; k < M; k++) {
					extraParameterSpace[k].push_back(semigoodparamvec[k]);
				}
				for (int f = 0; f < M; f++) {
					double testweight = std::pow(comboparam_uncertainties[f], -2.0);
					if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
						extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
					}
					else {
						extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
					}
				}
			}
			else if (comboparamset.size() == (M + 2))
			{
				comboparamset = regularGN(f, partialsvector, comboy, combox, parameters, tolerance); //uses regular GN

				if (comboparamset.size() == (M + 1)) // "runaway" parameter issue; include for median calc, but not mode
				{

					//std::cout << comboparamset[0] << "  " << comboparamset[1] << std::endl;

					std::vector <double> semigoodparamvec(M, 0.0);

					for (int j = 0; j < M; j++) {
						semigoodparamvec[j] = comboparamset[j]; //takes first M vals of comboparamset
					}

					comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += 1.0;
					}
					for (int k = 0; k < M; k++) {
						extraParameterSpace[k].push_back(semigoodparamvec[k]);
					}
					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							extraWeightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
						}
						else {
							extraWeightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
						}
					}
				}
				else
				{
					comboparam_uncertainties = paramuncertainty(partialsvector, combox, parameters);

					double correctivesum = 0.0;
					for (int i = 0; i < M; i++) {
						correctivesum += 1.0;
					}

					for (int f = 0; f < M; f++) {
						double testweight = std::pow(comboparam_uncertainties[f], -2.0);
						if ((testweight != testweight) || (std::isinf(testweight)) || (testweight == 0.0)) {
							weightSpace[f].push_back(DBL_MIN * correctivesum); //if the weight is NaN , makes it the smallest possible double val
							parameterSpace[f].push_back(comboparamset[f]);
						}
						else {
							weightSpace[f].push_back(testweight * correctivesum); //weighting the calculated parameters
							parameterSpace[f].push_back(comboparamset[f]); // vector of calculated vals for kth parameter
						}
					}
				}
			}
			// otherwise, singlar GN matrix issue due to data; excluded from all calculations
		}
	}


	//Now the parameter space is constructed, apply the priors (if there are any):
	if (has_priors) {

		std::vector <double> innerPostWeightSpace (weightSpace[0].size(), 0.0);
		std::vector <double> innerPostExtraWeightSpace (extraWeightSpace[0].size(), 0.0);
		std::vector <std::vector<double> > postWeightSpace (M, innerPostWeightSpace);
		std::vector <std::vector<double> > postExtraWeightSpace (M, innerPostExtraWeightSpace);
		std::vector <double> holdWeights, holdParams;
		std::vector <double> holdPost;
		std::vector <double> gaussianParamSet;
		std::vector <double> boundsSet;

		switch (priors.priorType) {
		case CUSTOM_PRIORS:
			// weight space
			for (int i = 0; i < weightSpace[0].size(); i++) { // for the ith set of parameters
				holdWeights.clear();
				holdParams.clear();
				for (int j = 0; j < M; j++) { //for each ith set, for the jth parameter of the M params 
					holdWeights.push_back(weightSpace[j][i]);
					holdParams.push_back(parameterSpace[j][i]);
				}
				holdPost = priors.p(holdParams, holdWeights); // takes in a vector of parameter weights, and returns a vector with modified weights
				for (int j = 0; j < M; j++) { //adds new, adjusted weights from priors, to the post weight space
					postWeightSpace[j][i] = holdPost[j];
				}
			}
			// extra weight space
			for (int i = 0; i < extraWeightSpace[0].size(); i++) { // for the ith set of parameters
				holdWeights.clear();
				holdParams.clear();
				for (int j = 0; j < M; j++) { //for each ith set, for the jth parameter of the M params 
					holdWeights.push_back(extraWeightSpace[j][i]);
					holdParams.push_back(extraParameterSpace[j][i]);
				}
				holdPost = priors.p(holdParams, holdWeights);
				for (int j = 0; j < M; j++) { //adds new, adjusted weights from priors, to the post weight space
					postExtraWeightSpace[j][i] = holdPost[j];
				}
			}
			break;
		case GAUSSIAN_PRIORS:
			// weight space
			for (int i = 0; i < weightSpace[0].size(); i++) { // for the ith set of parameters
				holdPost.clear();
				for (int j = 0; j < M; j++) { //for each ith set, for the jth parameter of the M params 
					gaussianParamSet = priors.gaussianParams[j];
					if ((std::isnan(gaussianParamSet[0]) == false) && ( std::isnan(gaussianParamSet[1]) == false) ){ //if neither are nans
						holdPost.push_back(weightSpace[j][i] * gaussian(parameterSpace[j][i], priors.gaussianParams[j][0], priors.gaussianParams[j][1]));
					}
					else { //if they're nans (don't modify priors))                             
						holdPost.push_back(weightSpace[j][i]);
					}
				}
				for (int j = 0; j < M; j++) { //adds new, adjusted weights from priors, to the post weight space
					postWeightSpace[j][i] = holdPost[j];
				}
			}
			// extra weight space
			for (int i = 0; i < extraWeightSpace[0].size(); i++) { // for the ith set of parameters
				holdPost.clear();
				for (int j = 0; j < M; j++) { //for each ith set, for the jth parameter of the M params 
					gaussianParamSet = priors.gaussianParams[j];
					if ((std::isnan(gaussianParamSet[0]) == false) && (std::isnan(gaussianParamSet[1]) == false)) { //if neither are nans
						holdPost.push_back(extraWeightSpace[j][i] * gaussian(parameterSpace[j][i], priors.gaussianParams[j][0], priors.gaussianParams[j][1]));
					}
					else { //if they're nans (don't modify priors))                             
						holdPost.push_back(extraWeightSpace[j][i]);
					}
				}
				for (int j = 0; j < M; j++) { //adds new, adjusted weights from priors, to the post weight space
					postExtraWeightSpace[j][i] = holdPost[j];
				}
			}
			break;
		case CONSTRAINED_PRIORS:
			// weight space
			for (int i = 0; i < weightSpace[0].size(); i++) { // for the ith set of parameters
				holdPost.clear();
				for (int j = 0; j < M; j++) { //for each ith set, for the jth parameter of the M params 
					boundsSet = priors.paramBounds[j];
					if ((std::isnan(boundsSet[0]) == true) && (std::isnan(boundsSet[1]) == true)) { //if both are nans (no bounds on the prior to add)
						holdPost.push_back(weightSpace[j][i]);
					}
					else if ((std::isnan(boundsSet[0]) == false) && (std::isnan(boundsSet[1]) == true)) { // lower bound      
						if (parameterSpace[j][i] < boundsSet[0]) { // is the parameter below the bound?
							holdPost.push_back(0.0); //makes the weight zero if the param is out of bounds
						}
						else {
							holdPost.push_back(weightSpace[j][i]); //continue
						}
					}
					else if ((std::isnan(boundsSet[0]) == true) && (std::isnan(boundsSet[1]) == false)) { // upper bound                          
						if (parameterSpace[j][i] > boundsSet[1]) { // is the parameter above the bound?
							holdPost.push_back(0.0); //makes the weight zero if the param is out of bounds
						}
						else {
							holdPost.push_back(weightSpace[j][i]); //continue
						}
					}
					else if ((std::isnan(boundsSet[0]) == false) && (std::isnan(boundsSet[1]) == false)) { // lower and upper bounds          
						if ( (parameterSpace[j][i] < boundsSet[0] ) || (parameterSpace[j][i] > boundsSet[1]) ) { // is the parameter not within the bounds?
							holdPost.push_back(0.0); //makes the weight zero if the param is out of bounds
						}
						else {
							holdPost.push_back(weightSpace[j][i]); //continue
						}
					}
				}
				for (int j = 0; j < M; j++) { //adds new, adjusted weights from priors, to the post weight space
					postWeightSpace[j][i] = holdPost[j];
				}
			}
			// extra weight space
			for (int i = 0; i < extraWeightSpace[0].size(); i++) { // for the ith set of parameters
				holdPost.clear();
				for (int j = 0; j < M; j++) { //for each ith set, for the jth parameter of the M params 
					boundsSet = priors.paramBounds[j];
					if ((std::isnan(boundsSet[0]) == true) && (std::isnan(boundsSet[1]) == true)) { //if both are nans (no bounds on the prior to add)
						holdPost.push_back(extraWeightSpace[j][i]);
					}
					else if ((std::isnan(boundsSet[0]) == false) && (std::isnan(boundsSet[1]) == true)) { // lower bound      
						if (extraParameterSpace[j][i] < boundsSet[0]) { // is the parameter below the bound?
							holdPost.push_back(0.0); //makes the weight zero if the param is out of bounds
						}
						else {
							holdPost.push_back(extraWeightSpace[j][i]); //continue
						}
					}
					else if ((std::isnan(boundsSet[0]) == true) && (std::isnan(boundsSet[1]) == false)) { // upper bound                          
						if (extraParameterSpace[j][i] > boundsSet[1]) { // is the parameter above the bound?
							holdPost.push_back(0.0); //makes the weight zero if the param is out of bounds
						}
						else {
							holdPost.push_back(extraWeightSpace[j][i]); //continue
						}
					}
					else if ((std::isnan(boundsSet[0]) == false) && (std::isnan(boundsSet[1]) == false)) { // lower and upper bounds          
						if ((extraParameterSpace[j][i] < boundsSet[0]) || (extraParameterSpace[j][i] > boundsSet[1])) { // is the parameter not within the bounds?
							holdPost.push_back(0.0); //makes the weight zero if the param is out of bounds
						}
						else {
							holdPost.push_back(extraWeightSpace[j][i]); //continue
						}
					}
				}
				for (int j = 0; j < M; j++) { //adds new, adjusted weights from priors, to the post weight space
					postExtraWeightSpace[j][i] = holdPost[j];
				}
			}

			break;
		case MIXED_PRIORS:
			// weight space
			for (int i = 0; i < weightSpace[0].size(); i++) { // for the ith set of parameters
				holdPost.clear();
				for (int j = 0; j < M; j++) { //for each ith set, for the jth parameter of the M params 
					gaussianParamSet = priors.gaussianParams[j];
					boundsSet = priors.paramBounds[j];
					//Gaussian:
					if ((std::isnan(gaussianParamSet[0]) == false) && (std::isnan(gaussianParamSet[1]) == false)) { //if neither are nans
						holdPost.push_back(weightSpace[j][i] * gaussian(parameterSpace[j][i], priors.gaussianParams[j][0], priors.gaussianParams[j][1]));
					}
					else { //if they're nans (don't modify priors))                             
						holdPost.push_back(weightSpace[j][i]);
					}
					//Constrained:
					if ((std::isnan(boundsSet[0]) == false) && (std::isnan(boundsSet[1]) == true)) { // lower bound      
						if (parameterSpace[j][i] < boundsSet[0]) { // is the parameter below the bound?
							holdPost[j] = (0.0); //makes the weight zero if the param is out of bounds
						}
					}
					else if ((std::isnan(boundsSet[0]) == true) && (std::isnan(boundsSet[1]) == false)) { // upper bound                          
						if (parameterSpace[j][i] > boundsSet[1]) { // is the parameter above the bound?
							holdPost[j] = (0.0); //makes the weight zero if the param is out of bounds
						}
					}
					else if ((std::isnan(boundsSet[0]) == false) && (std::isnan(boundsSet[1]) == false)) { // lower and upper bounds          
						if ((parameterSpace[j][i] < boundsSet[0]) || (parameterSpace[j][i] > boundsSet[1])) { // is the parameter not within the bounds?
							holdPost[j] = (0.0); //makes the weight zero if the param is out of bounds
						}
					}
				}
				for (int j = 0; j < M; j++) { //adds new, adjusted weights from priors, to the post weight space
					postWeightSpace[j][i] = holdPost[j];
				}
			}
			// extra weight space
			for (int i = 0; i < extraWeightSpace[0].size(); i++) { // for the ith set of parameters
				holdPost.clear();
				for (int j = 0; j < M; j++) { //for each ith set, for the jth parameter of the M params 
					gaussianParamSet = priors.gaussianParams[j];
					boundsSet = priors.paramBounds[j];
					//Gaussian:
					if ((std::isnan(gaussianParamSet[0]) == false) && (std::isnan(gaussianParamSet[1]) == false)) { //if neither are nans
						holdPost.push_back(extraWeightSpace[j][i] * gaussian(extraParameterSpace[j][i], priors.gaussianParams[j][0], priors.gaussianParams[j][1]));
					}
					else { //if they're nans (don't modify priors))                             
						holdPost.push_back(extraWeightSpace[j][i]);
					}
					//Constrained:
					if ((std::isnan(boundsSet[0]) == false) && (std::isnan(boundsSet[1]) == true)) { // lower bound      
						if (extraParameterSpace[j][i] < boundsSet[0]) { // is the parameter below the bound?
							holdPost[j] = (0.0); //makes the weight zero if the param is out of bounds
						}
					}
					else if ((std::isnan(boundsSet[0]) == true) && (std::isnan(boundsSet[1]) == false)) { // upper bound                          
						if (extraParameterSpace[j][i] > boundsSet[1]) { // is the parameter above the bound?
							holdPost[j] = (0.0); //makes the weight zero if the param is out of bounds
						}
					}
					else if ((std::isnan(boundsSet[0]) == false) && (std::isnan(boundsSet[1]) == false)) { // lower and upper bounds          
						if ((extraParameterSpace[j][i] < boundsSet[0]) || (extraParameterSpace[j][i] > boundsSet[1])) { // is the parameter not within the bounds?
							holdPost[j] = (0.0); //makes the weight zero if the param is out of bounds
						}
					}
				}
				for (int j = 0; j < M; j++) { //adds new, adjusted weights from priors, to the post weight space
					postWeightSpace[j][i] = holdPost[j];
				}
			}

			break;
		default:
			break;
		}

		//now that priors have been applied:
		weightSpace = postWeightSpace;
		extraWeightSpace = postExtraWeightSpace;
	};
}

std::vector<double> FunctionalForm::get_bestfit_errorbars(std::vector <double> best_fit_params){
	std::vector <double> goodw, goodsigma_y, goodx, result; // true-flagged x and y data 
	std::vector <std::vector <double> > goodx_ND;

	for (int i = 0; i < N; i++) {
		if (flags[i]) {
			goodw.push_back(w[i]);
			if (hasErrorBars) {
				goodsigma_y.push_back(sigma_y[i]);
			}

			if (NDcheck) {
				goodx_ND.push_back(x_ND[i]);
			}
			if (NDcheck == false) {
				goodx.push_back(x[i]);
			}
		}
	}

	if (weightedCheck) { //calculates wbar -- average weight of all unrejected data points
		double wsum = 0.0;
		double goodcount = 0.0;
		for (int j = 0; j < N; j++) {
			if (flags[j]) {
				wsum += w[j];
				goodcount += 1.0;
			}
		}
		wbar = wsum / goodcount;
	}

	if (weightedCheck && NDcheck) 
	{
		if (hasErrorBars) {
			result = paramuncertainty(NDpartialsvector, goodx_ND, best_fit_params, goodsigma_y, goodw, wbar);
		} else if (!hasErrorBars) {
			result = paramuncertainty(NDpartialsvector, goodx_ND, best_fit_params, goodw, wbar);
		}
	}
	else if (weightedCheck && (NDcheck == false))
	{
		if (hasErrorBars) {
			result = paramuncertainty(partialsvector, goodx, best_fit_params, goodsigma_y, goodw, wbar);
		}
		else if (!hasErrorBars) {
			result = paramuncertainty(partialsvector, goodx, best_fit_params, goodw, wbar);
		}
	}
	else if ((weightedCheck == false) && NDcheck)
	{
		if (hasErrorBars) {
			result = paramuncertainty(NDpartialsvector, goodx_ND, best_fit_params, goodsigma_y);
		}
		else if (!hasErrorBars) {
			// model parameter errors are only defined if you have weights and/or error bars
		}
	}
	else if ((weightedCheck == false) && (NDcheck == false))
	{
		if (hasErrorBars) {
			result = paramuncertainty(partialsvector, goodx, best_fit_params, goodsigma_y);
		}
		else if (!hasErrorBars) {
			// model parameter errors are only defined if you have weights and/or error bars
		}
	}

	return result;
}

std::vector<double> FunctionalForm::regression() //determines best fit params
{
	std::vector <double> goody, goodw, goodsigma_y, goodx, result; // true-flagged x and y data 
	std::vector <std::vector <double> > goodx_ND;

	for (int i = 0; i < N; i++) {
		if (flags[i]) {
			goody.push_back(y[i]);
			goodw.push_back(w[i]);
			if (hasErrorBars) {
				goodsigma_y.push_back(sigma_y[i]);
			}

			if (NDcheck) {
				goodx_ND.push_back(x_ND[i]);
			}
			if (NDcheck == false) {
				goodx.push_back(x[i]);
			}
		}
	}
	// USES REGULAR GN BECAUSE DOESN'T CONVERGE TO CORRECT PARAMS WITH MODIFIED GN; INVESTIGATE LATER
	if (weightedCheck && NDcheck) 
	{
		//result = regularGN(f_ND, NDpartialsvector, goody, goodx_ND, meanstartingpoint, tolerance, w);
		if (hasErrorBars) {
			result = modifiedGN(f_ND, NDpartialsvector, goody, goodx_ND, meanstartingpoint, goodsigma_y, tolerance, goodw);
		}
		else if (!hasErrorBars) {
			result = modifiedGN(f_ND, NDpartialsvector, goody, goodx_ND, meanstartingpoint, tolerance, goodw);
		}
	}
	else if (weightedCheck && (NDcheck == false))
	{
		//result = regularGN(f, partialsvector, goody, goodx, meanstartingpoint, tolerance, w);
		if (hasErrorBars) {
			result = modifiedGN(f, partialsvector, goody, goodx, meanstartingpoint, goodsigma_y, tolerance, goodw);
		}
		else if (!hasErrorBars) {
			result = modifiedGN(f, partialsvector, goody, goodx, meanstartingpoint, tolerance, goodw);
		}
	}
	else if ((weightedCheck == false) && NDcheck)
	{
		//result = regularGN(f_ND, NDpartialsvector, goody, goodx_ND, meanstartingpoint, tolerance);
		if (hasErrorBars) {
			result = modifiedGN(f_ND, NDpartialsvector, goody, goodx_ND, meanstartingpoint, goodsigma_y, tolerance);
		}
		else if (!hasErrorBars) {
			result = modifiedGN(f_ND, NDpartialsvector, goody, goodx_ND, meanstartingpoint, tolerance);
		}
	}
	else if ((weightedCheck == false) && (NDcheck == false))
	{
		//result = regularGN(f, partialsvector, goody, goodx, meanstartingpoint, tolerance);
		if (hasErrorBars) {
			result = modifiedGN(f, partialsvector, goody, goodx, meanstartingpoint, goodsigma_y, tolerance);
		}
		else if (!hasErrorBars) {
			result = modifiedGN(f, partialsvector, goody, goodx, meanstartingpoint, tolerance);
		}
	}

	std::vector <double> new_result = result;
	if (result.size() > M) {
		new_result.clear();
		for (int i = 0; i < M; i++) {
			new_result.push_back(result[i]);
		}
	}

	result = new_result;


	return result; //the final parameter vector (line will be set equal to this)
}
std::vector<double> FunctionalForm::getErrors(std::vector <double> line) //takes in the model function and the vector of parameters
{
	std::vector <double> paramsvec = line;
	/*
	for (int i = 0; i < M; i++) {
		std::cout << paramsvec[i] << "\t";
	}
	std::cout << std::endl;
	*/
	double modeledY;
	std::vector<double> toRet;
	for (int i = 0; i < x.size(); i++)
	{
		if (flags[i])
		{
			modeledY = f(x[i], paramsvec);
			toRet.push_back(y[i] - modeledY);
		}
	}
	return toRet;
}
std::vector<double> FunctionalForm::getErrors_ND(std::vector <double> line) //same, but with the case of >1 "x" type (independent) variables in the function
{
	std::vector <double> paramsvec = line;

	/*
	for (int i = 0; i < M; i++) {
		std::cout << paramsvec[i] << "\t";
	}
	std::cout << std::endl;
	*/
	double modeledY;
	std::vector<double> toRet;
	for (int i = 0; i < x_ND.size(); i++)
	{
		if (flags[i])
		{
			modeledY = f_ND(x_ND[i], paramsvec);
			toRet.push_back(y[i] - modeledY);
		}
	}
	return toRet;
}
void FunctionalForm::printData()
{

}
void FunctionalForm::setModel(std::vector<double> x)
{

}
void FunctionalForm::getCombos(std::vector <double> total, int k, int offset) { //1D case in x

	if (k == M) {

		combos.clear();
		combos_indices.clear();
	}
	if (k == 0) {
		combos.push_back(combination);
		combos_indices.push_back(combination_indices);
		return;
	}
	for (int i = offset; i <= N - k; ++i) {
		combination.push_back(total[i]);
		combination_indices.push_back(i);
		getCombos(total, k - 1, i + 1);
		combination.pop_back();
		combination_indices.pop_back();
	}
}

void FunctionalForm::getCombos(std::vector <std::vector <double> > total, int k, int offset) { //ND case in x

	if (k == M) {

		NDcombos.clear();
		combos_indices.clear();
	}
	if (k == 0) {
		NDcombos.push_back(NDcombination);
		combos_indices.push_back(combination_indices);
		return;
	}
	for (int i = offset; i <= N - k; ++i) {
		NDcombination.push_back(total[i]);
		combination_indices.push_back(i);
		getCombos(total, k - 1, i + 1);
		NDcombination.pop_back();
		combination_indices.pop_back();
	}
}




FunctionalForm::~FunctionalForm()
{
}
