#pragma once

#include <vector>
#include <cmath>

// EXAMPLE FUNCTIONS (with corresponding partials and vectors of said partials)
/*
class Example {
public:
	double bar;

	Example(double bar);
	Example();
};
*/
// LINEAR

double function_linear(double x, std::vector <double> params);


double partial1_linear(double x, std::vector <double> params);

double partial2_linear(double x, std::vector <double> params);



// QUADRATIC

double function_quadratic(double x, std::vector <double> params);


double partial1_quadratic(double x, std::vector <double> params);

double partial2_quadratic(double x, std::vector <double> params);

double partial3_quadratic(double x, std::vector <double> params);




// CUBIC

double function_cubic(double x, std::vector <double> params);


double partial1_cubic(double x, std::vector <double> params);

double partial2_cubic(double x, std::vector <double> params);

double partial3_cubic(double x, std::vector <double> params);

double partial4_cubic(double x, std::vector <double> params);




// POWER LAW

double function_powerlaw(double x, std::vector <double> params);


double partial1_powerlaw(double x, std::vector <double> params);

double partial2_powerlaw(double x, std::vector <double> params);



// EXPONENTIAL

double function_exponential(double x, std::vector <double> params);


double partial1_exponential(double x, std::vector <double> params);

double partial2_exponential(double x, std::vector <double> params);




// LOGARITHMIC

double function_logarithmic(double x, std::vector <double> params);


double partial1_logarithmic(double x, std::vector <double> params);