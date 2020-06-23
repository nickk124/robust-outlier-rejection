/*
 Robust Chauvenet Rejection (RCR) Official Codebase
 Active Author: Nick C. Konz
 Former Author: Michael Maples
 See license at https://github.com/nickk124/RCR
 */
#include "MiscFunctions.h"
#include <cfloat>

enum priorTypes {CUSTOM_PRIORS, GAUSSIAN_PRIORS, CONSTRAINED_PRIORS, MIXED_PRIORS};

extern double pivot; //xBar, logxBar; etc. Needs to be defined globally, not just as a member of a FunctionalForm instance, so that it isn't just continually updated, but can be computed within model function and partials without adding it as an extra argument.
extern std::vector <double> pivot_ND; //vector of pivot values, one for each of the multiple indep variables

struct FunctionalFormResults
{
	std::vector <double> parameters;
	std::vector <double> parameter_uncertainties;

	// pivot points / linear parameter correlation minimization
	double pivot;
	std::vector <double> pivot_ND;
};

class Priors
{
public:
	//constructors:
	Priors(priorTypes priorType, std::function <std::vector <double>(std::vector <double>, std::vector <double>)> p); //custom priors
	Priors(priorTypes priorType, std::vector < std::vector <double> > params); //Only Gaussian or only bounded/constrained
	Priors(priorTypes priorType, std::vector < std::vector <double> > gaussianParams, std::vector < std::vector <double> > paramBounds); //mixed
	Priors();

	priorTypes priorType;
	std::function <std::vector <double>(std::vector <double>, std::vector <double>)> p; // a pointer to a function that takes in a parameters vector and a weights vector and modifies it with the priors
	std::vector < std::vector <double> > gaussianParams; // a vector that contains a vector of mu and sigma for the guassian prior of each param. If no prior, then just use NANs.
	std::vector < std::vector <double> > paramBounds; // a vector that contains vectors of the bounds of each param. If not bounded, use NANs, and if there's only one bound, use NAN for the other "bound".
};

class FunctionalForm
{
public:
	/*
	Constructor Organization:
		Non-Priors:
			Non-weighted, with Error Bars:
				1D
				ND
			Weighted, with EB:
				1D
				ND
			**Non-weighted, without Error Bars
				1D
				ND
			**Weighted, without Error Bars
				1D
				ND

		With Priors:
			Non-weighted, with Error Bars:
				1D
				ND
			Weighted, with EB:
				1D
				ND
			**Non-weighted, without Error Bars
				1D
				ND
			**Weighted, without Error Bars
				1D
				ND
	
	
	
	*/



	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w);

	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w);
	//PRIORS support:
	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, Priors priorsObject);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, Priors priorsObject); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, Priors priorsObject);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, Priors priorsObject);

	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, Priors priorsObject);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, Priors priorsObject); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, Priors priorsObject);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, Priors priorsObject);

	//support for xBar/ (logx)Bar/any other custom average value:

	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess);

	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess);
	//PRIORS support:
	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, Priors priorsObject, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, Priors priorsObject, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, Priors priorsObject, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, Priors priorsObject, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess);

	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, Priors priorsObject, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, Priors priorsObject, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(std::function <double(double, std::vector <double> )> f, std::vector <double> x, std::vector<double> y, std::vector <std::function <double(double, std::vector <double>)> > partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, Priors priorsObject, std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc, double pivot_guess);
	FunctionalForm(std::function <double(std::vector <double>, std::vector <double>)> f_ND, std::vector < std::vector <double> > x, std::vector<double> y, std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, Priors priorsObject, std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) >pivotFunc_ND, std::vector <double> pivot_ND_guess);



	//default constructor:
	FunctionalForm();
	~FunctionalForm();

	// results
	FunctionalFormResults result;

	// below need to be public for usage by RCR class
	void buildModelSpace();
	void setTrueVec(std::vector<bool>&, std::vector<double>&, std::vector<double>&);
	void setTrueVec(std::vector<bool>&, std::vector<double>&);
	std::vector<double> regression();
	std::vector<double> get_bestfit_errorbars(std::vector <double> best_fit_params);
	std::vector<double> getErrors(std::vector <double> line);
	std::vector<double> getErrors_ND(std::vector <double> line);
	void setModel(std::vector<double>);

	// priors
	Priors priors;
	bool has_priors;

	// misc
	std::vector <double> trueW;
	std::vector<bool> flags;
	std::vector <double> parameters, meanstartingpoint;
	bool NDcheck;
	std::vector<int> indices;
	std::vector <std::vector<double> > parameterSpace, weightSpace, extraParameterSpace, extraWeightSpace; //generalized x vector of vectors, for >1D cases

	// custom function for computing pivot value (xbar, logxBar, etc)
	std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc;
	// ND case
	std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) > pivotFunc_ND;



private: 
	void printData();
	void getCombos(std::vector <double> total, int k, int offset); //1D case in x
	void getCombos(std::vector <std::vector <double> > total, int k, int offset);//ND case in x


	double wbar; //average of unrejected weights (not constant over time)

	std::vector<double> trueY, x, y, guess, sigma_y, w, modelSpaceW; //params is a vector of the parameters. meanstartingpoint is the initial guess/starting point for generalized mean GN
	std::vector <double> innerSpace;
	std::vector<std::vector<double> > x_ND; //generalized x vector of vectors, for >1D cases
	std::vector <std::function <double(double, std::vector <double>)> > partialsvector;
	std::vector <std::function <double(std::vector <double>, std::vector <double>)> > NDpartialsvector;
	int M, N; //number of params for GN, and the total number of points
	double tolerance = 1e-6; //wanted tolerance for GN
	//1D case in x:
	std::vector < std::vector <double> > combos;
	std::vector < std::vector <int> > combos_indices;
	std::vector <double> combination;
	std::vector <int> combination_indices;
	//ND case in x:
	std::vector < std::vector <std::vector <double > > > NDcombos;
	std::vector < std::vector <double> > NDcombination;
	std::vector < std::vector <int> > combosgood_indices;
	//check to see if this is ND case or not:
	bool weightedCheck, hasCustomPivot;
	//actual modeled function: (these are pointers to the function)
	std::function <double(double, std::vector <double> )> f;
    std::function <double(std::vector <double>, std::vector <double>)> f_ND;
	// //custom function for computing pivot value (xbar, logxBar, etc)
	// std::function < double(std::vector <double>, std::vector <double>, std::function < double(double, std::vector <double>) >, std::vector <double>) > pivotFunc;
	// //ND case
	// std::function <std::vector<double>(std::vector<std::vector<double> >, std::vector <double>, std::function <double(std::vector <double>, std::vector <double>)>, std::vector<double>) > pivotFunc_ND;

	bool hasErrorBars;

	//double function(double x, std::vector <double> params); // 1D case
	//double function_ND(std::vector <double> x, std::vector <double> params); // >1D case


};

