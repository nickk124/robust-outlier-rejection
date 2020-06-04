//#pragma once
#include "MiscFunctions.h"
#include <cfloat>

enum RCRpriorTypes {CUSTOM_PRIORS, GAUSSIAN_PRIORS, CONSTRAINED_PRIORS, MIXED_PRIORS};

extern double bar; //xBar, logxBar; etc. Needs to be defined globally, not just as a member of a FunctionalForm instance, so that it isn't just continually updated, but can be computed within model function and partials without adding it as an extra argument.
extern std::vector <double> bar_ND; //vector of bar values, one for each of the multiple indep variables

class RCRPriors
{
public:
	//constructors:
	RCRPriors(RCRpriorTypes priorType, std::vector <double>(*p)(std::vector <double>, std::vector <double>)); //custom priors
	RCRPriors(RCRpriorTypes priorType, std::vector < std::vector <double> > params); //Only Gaussian or only bounded/constrained
	RCRPriors(RCRpriorTypes priorType, std::vector < std::vector <double> > gaussianParams, std::vector < std::vector <double> > paramBounds); //mixed
	RCRPriors();

	RCRpriorTypes priorType;
	std::vector <double>(*p)(std::vector <double>, std::vector <double>); // a pointer to a function that takes in a parameters vector and a weights vector and modifies it with the priors
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



	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w);

	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w);
	//PRIORS support:
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, RCRPriors priorsObject);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, RCRPriors priorsObject); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, RCRPriors priorsObject);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, RCRPriors priorsObject);

	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, RCRPriors priorsObject);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, RCRPriors priorsObject); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, RCRPriors priorsObject);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, RCRPriors priorsObject);

	//support for xBar/ (logx)Bar/any other custom average value:

	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, double(*gB)(std::vector<double>, std::vector <double>, double(*)(double, std::vector <double>), std::vector<double>), double bar_guess);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector<double>(*gB_ND)(std::vector<std::vector<double> >, std::vector <double>, double(*)(std::vector <double>, std::vector <double>), std::vector<double>), std::vector <double> bar_ND_guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, double(*gB)(std::vector<double>, std::vector <double>, double(*)(double, std::vector <double>), std::vector<double>), double bar_guess);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, std::vector<double>(*gB_ND)(std::vector<std::vector<double> >, std::vector <double>, double(*)(std::vector <double>, std::vector <double>), std::vector<double>), std::vector <double> bar_ND_guess);

	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, double(*gB)(std::vector<double>, std::vector <double>, double(*)(double, std::vector <double>), std::vector<double>), double bar_guess);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector<double>(*gB_ND)(std::vector<std::vector<double> >, std::vector <double>, double(*)(std::vector <double>, std::vector <double>), std::vector<double>), std::vector <double> bar_ND_guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, double(*gB)(std::vector<double>, std::vector <double>, double(*)(double, std::vector <double>), std::vector<double>), double bar_guess);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, std::vector<double>(*gB_ND)(std::vector<std::vector<double> >, std::vector <double>, double(*)(std::vector <double>, std::vector <double>), std::vector<double>), std::vector <double> bar_ND_guess);
	//PRIORS support:
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, RCRPriors priorsObject, double(*gB)(std::vector<double>, std::vector <double>, double(*)(double, std::vector <double>), std::vector<double>), double bar_guess);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, RCRPriors priorsObject, std::vector<double>(*gB_ND)(std::vector<std::vector<double> >, std::vector <double>, double(*)(std::vector <double>, std::vector <double>), std::vector<double>), std::vector <double> bar_ND_guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, RCRPriors priorsObject, double(*gB)(std::vector<double>, std::vector <double>, double(*)(double, std::vector <double>), std::vector<double>), double bar_guess);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, RCRPriors priorsObject, std::vector<double>(*gB_ND)(std::vector<std::vector<double> >, std::vector <double>, double(*)(std::vector <double>, std::vector <double>), std::vector<double>), std::vector <double> bar_ND_guess);

	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, RCRPriors priorsObject, double(*gB)(std::vector<double>, std::vector <double>, double(*)(double, std::vector <double>), std::vector<double>), double bar_guess);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, RCRPriors priorsObject, std::vector<double>(*gB_ND)(std::vector<std::vector<double> >, std::vector <double>, double(*)(std::vector <double>, std::vector <double>), std::vector<double>), std::vector <double> bar_ND_guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, RCRPriors priorsObject, double(*gB)(std::vector<double>, std::vector <double>, double(*)(double, std::vector <double>), std::vector<double>), double bar_guess);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, RCRPriors priorsObject, std::vector<double>(*gB_ND)(std::vector<std::vector<double> >, std::vector <double>, double(*)(std::vector <double>, std::vector <double>), std::vector<double>), std::vector <double> bar_ND_guess);



	//default constructor:
	FunctionalForm();

	~FunctionalForm();


	void buildModelSpace();
	void setTrueVec(std::vector<bool>&, std::vector<double>&, std::vector<double>&);
	void setTrueVec(std::vector<bool>&, std::vector<double>&);
	std::vector<double> regression();
	std::vector<double> getErrors(std::vector <double> line);
	std::vector<double> getErrors_ND(std::vector <double> line);
	void setModel(std::vector<double>);

	//void initializeBar();
	double bar_result; //xBar, logxBar, etc; this is not used by the model function directly (as that needs to be defined globally), so this is just used as the final bar result.
	std::vector <double> bar_ND_result, trueW;
	std::vector<bool> flags;
	std::vector <double> parameters, meanstartingpoint;
	bool NDcheck;
	std::vector<int> indices;

	std::vector<std::vector<double> > parameterSpace, weightSpace, extraParameterSpace, extraWeightSpace; //generalized x vector of vectors, for >1D cases


private: 
	void printData();
	void getCombos(std::vector <double> total, int k, int offset); //1D case in x
	void getCombos(std::vector <std::vector <double> > total, int k, int offset);//ND case in x


	double wbar; //average of unrejected weights (not constant over time)

	std::vector<double> trueY, x, y, guess, sigma_y, w, modelSpaceW; //params is a vector of the parameters. meanstartingpoint is the initial guess/starting point for generalized mean GN
	std::vector <double> innerSpace;
	std::vector<std::vector<double> > x_ND; //generalized x vector of vectors, for >1D cases
	std::vector <double(*)(double, std::vector <double>)> partialsvector;
	std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector;
	int M, N; //number of params for GN, and the total number of points
	double tolerance; //wanted tolerance for GN
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
	bool weightedCheck, customBarCheck;
	//actual modeled function: (these are pointers to the function)
	double(*f)(double, std::vector <double>);
	double(*f_ND)(std::vector <double>, std::vector <double>);
	//custom function for computing bar value (xbar, logxBar, etc)
	double(*gB)(std::vector<double>, std::vector <double>, double(*)(double, std::vector <double>), std::vector<double>);
	//ND case
	std::vector<double>(*gB_ND)(std::vector<std::vector<double> >, std::vector <double>, double(*)(std::vector <double>, std::vector <double>), std::vector<double>);

	RCRPriors priorsObject;
	bool hasPriors;
	bool hasErrorBars;

	//double function(double x, std::vector <double> params); // 1D case
	//double function_ND(std::vector <double> x, std::vector <double> params); // >1D case


};

