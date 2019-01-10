//#include "stdafx.h"
#include "Demo.h"

// EXAMPLE FUNCTIONS (with corresponding partials and vectors of said partials)

/*
Example::Example(double bar) {
	this->bar = bar;
}

Example::Example() {

}

Example exBar;
*/
// LINEAR

double xBar = 4.5;
double lnx_Bar;

double function_linear(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 + a1 * (x - xBar);
}


double partial1_linear(double x, std::vector <double> params) {

	return 1.0;
}

double partial2_linear(double x, std::vector <double> params) {

	return x - xBar;
}

std::vector <double(*)(double, std::vector <double>)> partialsvector_linear = { partial1_linear, partial2_linear };

// QUADRATIC

double function_quadratic(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];
	double a2 = params[2];

	return a0 + a1 * (x - xBar) + a2 * std::pow((x - xBar), 2.0);
}

double partial1_quadratic(double x, std::vector <double> params) {

	return 1.0;
}

double partial2_quadratic(double x, std::vector <double> params) {

	return x-xBar;
}

double partial3_quadratic(double x, std::vector <double> params) {

	return std::pow((x - xBar), 2.0);
}

std::vector <double(*)(double, std::vector <double>)> partialsvector_quadratic = { partial1_quadratic, partial2_quadratic, partial3_quadratic };


// CUBIC

double function_cubic(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];
	double a2 = params[2];
	double a3 = params[3];

	return a0 + a1 * (x - xBar) + a2 * std::pow((x - xBar), 2.0) + a3 * std::pow((x - xBar), 3.0);
}

double partial1_cubic(double x, std::vector <double> params) {

	return 1.0;
}

double partial2_cubic(double x, std::vector <double> params) {

	return x - xBar;
}

double partial3_cubic(double x, std::vector <double> params) {

	return std::pow((x - xBar), 2.0);
}

double partial4_cubic(double x, std::vector <double> params) {

	return std::pow((x - xBar), 3.0);
}

std::vector <double(*)(double, std::vector <double>)> partialsvector_cubic= { partial1_cubic, partial2_cubic, partial3_cubic, partial4_cubic };


// POWER LAW

double function_powerlaw(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * std::pow((x / std::exp(lnx_Bar)), a1);
}

double partial1_powerlaw(double x, std::vector <double> params) {
	double a1 = params[1];

	return std::pow((x / std::exp(lnx_Bar)), a1);
}

double partial2_powerlaw(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * std::pow((x / std::exp(lnx_Bar)), a1) * std::log(x / std::exp(lnx_Bar));
}

std::vector <double(*)(double, std::vector <double>)> partialsvector_powerlaw = { partial1_powerlaw, partial2_powerlaw};

// EXPONENTIAL

double function_exponential(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0*std::exp(a1*(x - xBar));
}

double partial1_exponential(double x, std::vector <double> params) {
	double a1 = params[1];

	return std::exp(a1*(x - xBar));
}

double partial2_exponential(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * (x - xBar) * std::exp(a1*(x - xBar));
}

std::vector <double(*)(double, std::vector <double>)> partialsvector_exponential = { partial1_exponential, partial2_exponential };

// LOGARITHMIC

double function_logarithmic(double x, std::vector <double> params) {
	double a0 = params[0];

	return a0 * std::log(x - xBar);
}

double partial1_logarithmic(double x, std::vector <double> params) {
	double a0 = params[0];

	return std::log(x - xBar);
}

std::vector <double(*)(double, std::vector <double>)> partialsvector_logarithmic = { partial1_logarithmic};