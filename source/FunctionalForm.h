#pragma once
#include "MiscFunctions.h"

enum priorTypes {CUSTOM, GAUSSIAN, CONSTRAINED, MIXED};

class Priors
{
public:
	//constructors:
	Priors(priorTypes priorType, std::vector <double>(*p)(std::vector <double>)); //custom priors
	Priors(priorTypes priorType, std::vector < std::vector <double> > params); //Only Gaussian or only bounded/constrained
	Priors(priorTypes priorType, std::vector < std::vector <double> > gaussianParams, std::vector < std::vector <double> > paramBounds); //mixed
	Priors();

	priorTypes priorType;
	std::vector <double>(*p)(std::vector <double>); // a pointer to a function that takes in a weights vector and modifies it with the priors
	std::vector < std::vector <double> > gaussianParams; // a vector that contains a vector of mu and sigma for the guassian prior of each param. If no prior, then just use NANs.
	std::vector < std::vector <double> > paramBounds; // a vector that contains vectors of the bounds of each param. If not bounded, use NANs, and if there's only one bound, use NAN for the other "bound".
};

class FunctionalForm
{
public:
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w);
	//PRIORS support:
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, Priors priorsObject);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, Priors priorsObject); //case of there being <1 indepedent variable (x variable) in the function
	FunctionalForm(double(*f)(double, std::vector <double>), std::vector <double> x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(double, std::vector <double>)> partialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, Priors priorsObject);
	FunctionalForm(double(*f)(std::vector <double>, std::vector <double>), std::vector < std::vector <double> > x, std::vector<double> y, std::vector<double> sigma_y, std::vector <double(*)(std::vector <double>, std::vector <double>)> NDpartialsvector, double tolerance, std::vector <double> guess, std::vector <double> w, Priors priorsObject);

	//default constructor:
	FunctionalForm();

	void setTrueVec(std::vector<bool>&, std::vector<double>&, std::vector<double>&);
	void setTrueVec(std::vector<bool>&, std::vector<double>&);
	void buildModelSpace();
	std::vector<double> regression();
	std::vector<double> getErrors(std::vector <double> line);
	std::vector<double> getErrors_ND(std::vector <double> line);
	void setModel(std::vector<double>);
	void printData();
	void getCombos(std::vector <double> total, int k, int offset); //1D case in x
	void getCombos(std::vector <std::vector <double> > total, int k, int offset);//ND case in x

	~FunctionalForm();

	double wbar; //average of unrejected weights (not constant over time)

	std::vector<bool> flags;
	std::vector<int> indices;
	std::vector<double> trueW, trueY, x, y, guess, sigma_y, w, modelSpaceW, meanstartingpoint; //params is a vector of the parameters. meanstartingpoint is the initial guess/starting point for generalized mean GN
	std::vector <double> innerSpace;
	std::vector <double> parameters;
	std::vector<std::vector<double> > parameterSpace, weightSpace, x_ND, extraParameterSpace, extraWeightSpace; //generalized x vector of vectors, for >1D cases
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
	bool NDcheck, weightedCheck;
	//actual modeled function: (these are pointers to the function)
	double(*f)(double, std::vector <double>);
	double(*f_ND)(std::vector <double>, std::vector <double>);

	Priors priorsObject;
	bool hasPriors;

	//double function(double x, std::vector <double> params); // 1D case
	//double function_ND(std::vector <double> x, std::vector <double> params); // >1D case


};

