//#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <set>

//MATRIX MATH
// Function to get cofactor of A[p][q] in temp[][]. n is current
// dimension of A[][]
//void getCofactor(double A[M][M], double temp[M][M], int p, int q, int n)
std::vector<std::vector<double> > getCofactor(std::vector<std::vector <double> > A, int p, int q);
//Recursive function for finding determinant of matrix.
//n is current dimension of A[][].
double determinant(std::vector < std::vector <double> > A);

std::vector<std::vector<double> > adjoint(std::vector < std::vector <double> > A);// Function to get adjoint of A[M][M] in adj[M][M].

std::vector<std::vector <double> > inverse(std::vector < std::vector <double> > A);

std::vector<std::vector<double> > pivotSystem(std::vector <std::vector<double> >A, std::vector<double> b); //makes inverting matrices more numerically accurate

std::vector < std::vector <double> > transpose(std::vector < std::vector <double> > array); //takes the transpose of some input array

std::vector <double> dot(std::vector< std::vector <double> > A, std::vector <double> b);

std::vector< std::vector <double> > dot(std::vector< std::vector <double> > A, std::vector <std::vector <double> > B); // multiplication of two square matrices

std::vector <double> forwardSubstitution(std::vector <std::vector <double> > A, std::vector <double> b); // algorithm used to solve lower-triangular matrix equations																		

std::vector <std::vector <double> > LUInverse(std::vector <std::vector <double> > A);

//std::vector < std::vector <double> > pseudoInverse(std::vector < std::vector <double> > A); // Computes the Moore-Penrose Pseudo-Inverse of some square Matrix A, used for the Gauss-Newton Method

//std::vector < std::vector <double> > weightedPseudoInverse(std::vector < std::vector <double> > A, std::vector < std::vector <double> > W); //

//MODEL FITTING MATH

double chiSquared(double (*f)(double, std::vector <double>), std::vector <double> y, std::vector <double> x, std::vector <double> params, std::vector <double> sigma_y); // computes the chi-squared value. takes the function y=... as an argument (via pointer) 

double chiSquared(double (*f)(std::vector <double>, std::vector <double>), std::vector <double> y, std::vector <std::vector<double> > x, std::vector <double> params, std::vector <double> sigma_y); // multi-dimensional independent variables case

double chiSquared(double(*f)(double, std::vector <double>), std::vector <double> y, std::vector <double> x, std::vector <double> params, std::vector <double> w, int K); // computes the chi-squared value. takes the function y=... as an argument (via pointer) 

double chiSquared(double(*f)(std::vector <double>, std::vector <double>), std::vector <double> y, std::vector <std::vector<double> > x, std::vector <double> params, std::vector <double> w, int K); // multi-dimensional independent variables case


std::vector <double> residuals(double (*f)(double, std::vector <double>), std::vector <double> y, std::vector <double> x, std::vector <double> params); // computes the residuals vector needed in the GN algorithm

std::vector <double> residuals(double (*f)(std::vector <double>, std::vector <double>), std::vector <double> y, std::vector <std::vector<double> > x, std::vector <double> params); // multi-dimensional independent variables case

std::vector < std::vector <double> > jacobian(std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> x, std::vector <double> params); //Jacobian matrix, CASE of 1 independent var, takes argument of vector of parameters. USER INPUTTED.

std::vector < std::vector <double> > jacobian(std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector< std::vector <double> > x, std::vector <double> params); //Jacobian matrix, CASE of >1 independent var, takes argument of vector of parameters. USER INPUTTED.
//weighted:
std::vector <double> paramuncertainty(std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> x, std::vector <double> params, std::vector <double> sigma_y, std::vector <double> w, double wbar); //1D case

std::vector <double> paramuncertainty(std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector< std::vector <double> > x, std::vector <double> params, std::vector <double> sigma_y, std::vector <double> w, double wbar); // >1D case

//non weighted:
std::vector <double> paramuncertainty(std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> x, std::vector <double> params, std::vector <double> sigma_y); //1D case

std::vector <double> paramuncertainty(std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector< std::vector <double> > x, std::vector <double> params, std::vector <double> sigma_y); //>1D case

//w, no sy, ND
std::vector <double> paramuncertainty(std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> x, std::vector <double> params, std::vector <double> w, double wbar); //1D case
//w, no sy, ND
std::vector <double> paramuncertainty(std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector< std::vector <double> > x, std::vector <double> params, std::vector <double> w, double wbar); // >1D case

																																																							   // no w or sy
std::vector <double> paramuncertainty(std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> x, std::vector <double> params); //1D case
// no w or sy, ND
std::vector <double> paramuncertainty(std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector< std::vector <double> > x, std::vector <double> params); //>1D case
//GAUSS-NEWTON ALGORITHMS

//MODIFIED
//NON-WEIGHTED
std::vector <double> modifiedGN(double (*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, std::vector <double> sigma_y, double tolerance); //the case of 1 independent (x) variables in the function

std::vector <double> modifiedGN(double (*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, std::vector <double> sigma_y, double tolerance); //the case of >1 independent (x) variables in the function
//WEIGHTED
std::vector <double> modifiedGN(double (*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, std::vector <double> sigma_y, double tolerance, std::vector <double> w); //the case of 1 independent (x) variables in the function

std::vector <double> modifiedGN(double (*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, std::vector <double> sigma_y, double tolerance, std::vector <double> w); //the case of >1 independent (x) variables in the function

//NON-WEIGHTED, without error bars:
std::vector <double> modifiedGN(double(*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, double tolerance); //the case of 1 independent (x) variables in the function

std::vector <double> modifiedGN(double(*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, double tolerance); //the case of >1 independent (x) variables in the function

//WEIGHTED, without error bars:
std::vector <double> modifiedGN(double(*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, double tolerance, std::vector <double> w); //the case of 1 independent (x) variables in the function

std::vector <double> modifiedGN(double(*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, double tolerance, std::vector <double> w); //the case of >1 independent (x) variables in the function


//REGULAR GAUSS NEWTON 
//NON-WEIGHTED
std::vector <double> regularGN(double(*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, double tolerance); //the case of 1 independent (x) variables in the function

std::vector <double> regularGN(double(*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, double tolerance); //the case of >1 independent (x) variables in the function
//WEIGHTED:
std::vector <double> regularGN(double(*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, double tolerance, std::vector <double> w); //the case of 1 independent (x) variables in the function

std::vector <double> regularGN(double(*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, double tolerance, std::vector <double> w); //the case of >1 independent (x) variables in the function

// misc
double factorial(double n);

double fRand(double fMin, double fMax);

double gaussian(double x, double mu, double sig);

double getAvg(std::vector<double> x, std::vector <double> w, double(*f)(double, std::vector <double>), std::vector<double> params);

double getAvg_Exp(std::vector<double> x, std::vector <double> w, double(*f)(double, std::vector <double>), std::vector<double> params);

double getLogXBar_PowerLaw(std::vector<double> x, std::vector <double> w, double(*f)(double, std::vector <double>), std::vector<double> params);

double getLogXBar_Log(std::vector<double> x, std::vector <double> w, double(*f)(double, std::vector <double>), std::vector<double> params);
