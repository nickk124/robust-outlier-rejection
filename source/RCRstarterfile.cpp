//#include "stdafx.h"
#include "RCR.h"
//#include "RCRwebutils.h"
#include <iostream>
#include <iterator>
#include <fstream>

// SEE STARTER GUIDE ATTACHED IN .ZIP FILE FOR BASIC HOW-TO INSTRUCTIONS

// EXAMPLE CODE FOR FUNCTIONAL / MODEL-FITTING RCR COMMENTED OUT BELOW:

/*
// ENTER YOUR FUNCTION HERE (either with one independent variable, or multiple, in a vector.
//also takes vector of parameters in as argument
//function arguments need to be defined in that order:
//example function given: y = b*exp(m(x-0.5))

double func(double x, std::vector <double> params) {
    double b = params[0];
    double m = params[1];

    return b * std::exp(m *(x - 0.5));
}

// ENTER THE PARTIAL DERIVATIVES HERE, IN ORDER OF PARAMETERS
// EXAMPLE PARTIALS FOR y=b*exp(mx) shown

double partial1(double x, std::vector <double> params) //example of case of 1 indep var
{
    //here you may declare the needed parameters. MAKE SURE TO USE THE CORRECT ORDERING OF PARAMS AS THAT OF YOUR PARAM VECTOR
    double b = params[0];
    double m = params[1];

    //now return the partial with everything plugged in.
    return std::exp(m *(x - 0.5));

};

double partial2(double x, std::vector <double> params) //example of case of 1 indep var
{
    //here you may declare the needed parameters. MAKE SURE TO USE THE CORRECT ORDERING OF PARAMS AS THAT OF YOUR PARAM VECTOR
    double b = params[0];
    double m = params[1];

    //now return the partial with everything plugged in.
    return (x-0.5) * b * std::exp(m *(x - 0.5));

};
*/
// EXAMPLE FOR >1 INDEP VARIABLES SHOWN for y = a0 + a1*x1 + a2*x2
// (in this case, these independent "x" variables are in a vector. In the below example, x = {x1, x2})

//example function for >1 indep. variable, y = a_0 + a_1*x_1 + a_2*x_2 (plane)
double func_ND(std::vector <double> x, std::vector <double> params) {
    double a0 = params[0];
    double a1 = params[1];
    double a2 = params[2];

    double x1 = x[0];
    double x2 = x[1];

    return a0 + a1 * x1 + a2 * x2;

}

double NDpartial1(std::vector <double> x, std::vector <double> params) {

    return 1.0;
}

double NDpartial2(std::vector <double> x, std::vector <double> params) {

    double x1 = x[0];

    return x1;
}

double NDpartial3(std::vector <double> x, std::vector <double> params) {

    double x2 = x[1];

    return x2;
}

//Next create your vector of partial derivatives for use in the jacobian (use the lower definition if you have >1 independent variable):

// CASE OF ONLY ONE INDEP VARIABLE:
//std::vector <double(*)(double, std::vector <double>)> partialsvector = { partial1, partial2 };

// CASE OF >1 INDEP VARIABLES:
std::vector <double(*)(std::vector <double>, std::vector <double>)> partialsvector_ND = { NDpartial1, NDpartial2, NDpartial3 };

//example 1D linear function

//double xBAR;

double func_linear(double x, std::vector <double> params) {
    double a0 = params[0];
    double a1 = params[1];

    return a0 + a1 * (x - bar);
}


double partial1_lin(double x, std::vector <double> params) {

    return 1.0;
}

double partial2_lin(double x, std::vector <double> params) {

    return x - bar;
}

std::vector <double(*)(double, std::vector <double>)>  partialsvector_lin = { partial1_lin, partial2_lin };

/*
double func(double x, std::vector <double> params) {
    double sigma = params[0];
    double a = params[1];

    return std::exp((-1.0*std::pow((x - a), 2.0)) / (2.0*std::pow(sigma, 2.0))) / (std::sqrt(2.0*PI) * sigma);
}

double partial1(double x, std::vector <double> params) //example of case of 1 indep var
{
    //here you may declare the needed parameters. MAKE SURE TO USE THE CORRECT ORDERING OF PARAMS AS THAT OF YOUR PARAM VECTOR
    double sigma = params[0];
    double a = params[1];

    //now return the partial with everything plugged in.
    return std::exp((-1.0*std::pow((x - a), 2.0)) / (2.0*std::pow(sigma, 2.0))) * (std::pow((x - a), 2.0) / (std::sqrt(2 * PI) * std::pow(sigma, 4.0)));

};

double partial2(double x, std::vector <double> params) //example of case of 1 indep var
{
    //here you may declare the needed parameters. MAKE SURE TO USE THE CORRECT ORDERING OF PARAMS AS THAT OF YOUR PARAM VECTOR
    double sigma = params[0];
    double a = params[1];

    //now return the partial with everything plugged in.
    return std::exp((-1.0*std::pow((x - a), 2.0)) / (2.0*std::pow(sigma, 2.0))) * ((x - a) / (std::sqrt(2.0*PI) * std::pow(sigma, 3.0)));

};

// CASE OF ONLY ONE INDEP VARIABLE:
std::vector <double(*)(double, std::vector <double>)> partialsvector = { partial1, partial2 };
*/


// Example of custom priors function
/*
std::vector <double> custom_priors(std::vector <double> weights)
{
    std::vector <double> hold;
    for (int i = 0; i < weights.size(); i++) {
        hold.push_back(weights[i]*2.0);
    };

    return hold;
};
*/

int main()
{
    /*
    // SEE STARTER GUIDE ATTACHED IN .ZIP FILE FOR BASIC HOW-TO INSTRUCTIONS

    //  ----------EXAMPLE CODE FOR EXECUTION OF RCR COMMENTED OUT BELOW: ----------


    // EXAMPLES SHOWN BELOW USE SS_MEDIAN_DL AND BULK REJECTION TECHNIQUES


    // ----------BASIC RCR ----------

    // UNWEIGHTED DATA
    RCR rcr = RCR(SS_MEDIAN_DL);
    rcr.performBulkRejection(y);

    // WEIGHTED DATA
    RCR rcr = RCR(SS_MEDIAN_DL);
    rcr.performBulkRejection(weights, y);

    //  ----------NON-PARAMETRIC RCR ----------

    // UNWEIGHTED DATA
    NonParametric model = NonParametric();
    RCR rcr = RCR(SS_MEDIAN_DL);
    rcr.setNonParametricModel(model);
    rcr.performBulkRejection(weights, y);

    // WEIGHTED DATA
    NonParametric model = NonParametric();
    RCR rcr = RCR(SS_MEDIAN_DL);
    rcr.setNonParametricModel(model);
    rcr.performBulkRejection(weights, y);

    //  ----------PARAMETRIC/FUNCTIONAL FORM RCR ----------
    // (SEE ABOVE int.main() FOR EXMAPLES OF FUNCTION, PARTIALS, ETC. DEFINITIONS


    // UNWEIGHTED DATA, SINGLE INDEPENDENT VARIABLE IN MODEL FUNCTION
    FunctionalForm model = FunctionalForm(func, x, y, sigma_y, partialsvector, tolerance, guess);
    RCR rcr = RCR(SS_MEDIAN_DL);
    rcr.setParametricModel(model);
    rcr.performBulkRejection(weights, y);

    // UNWEIGHTED DATA, MULTIPLE INDEPENDENT VARIABLES IN MODEL FUNCTION
    FunctionalForm model = FunctionalForm(func, x_ND, y, sigma_y, partialsvector_ND, tolerance, guess);
    RCR rcr = RCR(SS_MEDIAN_DL);
    rcr.setParametricModel(model);
    rcr.performBulkRejection(weights, y);

    // WEIGHTED DATA, SINGLE INDEPENDENT VARIABLE IN MODEL FUNCTION
    FunctionalForm model = FunctionalForm(func, x, y, sigma_y, partialsvector, tolerance, guess, weights);
    RCR rcr = RCR(SS_MEDIAN_DL);
    rcr.setParametricModel(model);
    rcr.performBulkRejection(weights, y);


    // WEIGHTED DATA, MULTIPLE INDEPENDENT VARIABLES IN MODEL FUNCTION
    FunctionalForm model = FunctionalForm(func, x_ND, y, sigma_y, partialsvector_ND, tolerance, guess, weights);
    RCR rcr = RCR(SS_MEDIAN_DL);
    rcr.setParametricModel(model);
    rcr.performBulkRejection(weights, y);
    */

    //EXAMPLE with weighted linear data and Priors:

    /*

    std::vector <double(*)(double, std::vector <double>)> partialsvector_linear = { partial1_linear, partial2_linear };


    double tolerance = 0.01;
    std::vector <double> x, y, w, guess, sigma_y;
    guess = { 2.1, 0.9 };
    y = { 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 20.0, 11.0 };
    sigma_y = { 0.1, 0.2, 0.4, 0.2, 0.3, 0.5, 0.10, 0.11, 1.0, 0.2 };
    x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
    w = { 0.5, 1.0, 1.2, 3.0, 1.0, 0.6, 0.2, 1.0, 15.0, 2.0 };

    std::vector < std::vector <double> > gaussianParams = { {NAN, NAN}, {1.0, 2.0} }; //params specifiying the Gaussian priors on the model function parameters
    std::vector < std::vector <double> > boundedParams = { {0.0, NAN}, {NAN, NAN} };  //params specifiying the bounds (priors) on the model function parameters

    Priors testPriors = Priors(MIXED, gaussianParams, boundedParams); //construction of priors object

    FunctionalForm model = FunctionalForm(function_linear, x, y, sigma_y, partialsvector_linear, tolerance, guess, w, testPriors); //setting up for functional form
    RCR rcr = RCR(SS_MEDIAN_DL); //setting up for RCR with this rejection technique
    rcr.setParametricModel(model);
    rcr.performBulkRejection(w, y); //running Bulk rejection RCR

    std::vector <double> final_parameters = model.parameters; // final calculated model function parameters

    std::cout << "\n\nparameters: \n"; //prints out calculated model function parameters
    for (int i = 0; i < final_parameters.size(); i++) {
        std::cout << final_parameters[i] << std::endl;
    }
    std::cout << "\n\n\n";

    std::vector <double> rejected_data = rcr.result.rejectedY;
    std::vector <double> nonrejected_data = rcr.result.cleanY;
    std::vector <bool> flags = rcr.result.flags; // booleans specifying whether each point, in order, was rejected; true for nonrejected, false for rejected.




    //example of usage with >1 independent ("x") variable (weighted data, with priors)

    double tolerance_ND = 0.01;
    std::vector <double> y_ND, w_ND, guess_ND, real_params, sigma_y_ND;
    real_params = { 1.0, 2.0, 3.0 };
    guess_ND = { 1.1, 1.9, 3.3 };

    std::vector <double> innerx_ND;
    std::vector <std::vector<double> > x_ND(25, innerx_ND);

    int count = 0;

    //building x data
    for (int i = -2; i < 3; i++) {
        for (int j = -2; j < 3; j++) {
            x_ND[count].push_back(i);
            x_ND[count].push_back(j);

            y_ND.push_back(func_ND(x_ND[count], real_params));
            sigma_y_ND.push_back(0.2);
            w_ND.push_back(1.0);

            count += 1;
        }
    };

    y_ND[6] = 2200.0; //intentionally bad data point
    sigma_y_ND[6] = 5.0;
    w_ND[6] = 0.1;

    std::vector < std::vector <double> > gaussianParams_ND = { {NAN, NAN}, {1.0, 2.0}, {NAN, NAN} }; //params specifiying the Gaussian priors on the model function parameters
    std::vector < std::vector <double> > boundedParams_ND = { {0.0, NAN}, {NAN, NAN}, {NAN, NAN} };  //params specifiying the bounds (priors) on the model function parameters

    Priors testPriors_ND = Priors(MIXED, gaussianParams_ND, boundedParams_ND); //construction of priors object

    FunctionalForm model_ND = FunctionalForm(func_ND, x_ND, y_ND, sigma_y_ND, partialsvector_ND, tolerance_ND, guess_ND, w_ND, testPriors_ND);
    RCR rcr_ND = RCR(SS_MEDIAN_DL); //setting up for RCR with this rejection technique
    rcr_ND.setParametricModel(model_ND);
    rcr_ND.performBulkRejection(w_ND, y_ND); //running Bulk rejection RCR

    std::vector <double> final_parameters_ND = model_ND.parameters; // final calculated model function parameters

    std::cout << "\n\nND parameters: \n"; //prints out calculated model function parameters
    for (int i = 0; i < final_parameters_ND.size(); i++) {
        std::cout << final_parameters_ND[i] << std::endl;
    }
    std::cout << "\n\n\n";

    std::vector <double> rejected_data_ND = rcr_ND.result.rejectedY;
    std::vector <double> nonrejected_data_ND = rcr_ND.result.cleanY;
    std::vector <bool> flags_ND = rcr_ND.result.flags; // booleans specifying whether each point, in order, was rejected; true for nonrejected, false for rejected.

    */

    //example of usage with >1 independent ("x") variable (weighted data, with priors, but no error bars)

    
    

//    std::vector <double> result;
//    std::vector <double> xt = { 1,2,3,4,5,6,7,8,9 };
//    std::vector <double> yt = { 1.0, 2.95, 5.55, 8.69, 12.31, 16.37, 20.81, 25.63, 51.0 };
//    std::vector <double> wt = { 1,1.1,0.9,1,1.1,0.9,1.2,1,0.1 };
//    double bt = 1.1;
//    std::vector <double> gt = { 1.0,2.0};
//    std::vector <double> pP;
//    std::vector <int> hP;
//    int funcInt = 4;
//    int rejTech = 2;
//    result = requestHandlerWeighted(xt, yt, gt, wt, funcInt, xt.size(), rejTech, 0, pP, hP, bt);




    double tolerance_ND = 0.01;
    std::vector <double> y_ND, w_ND, guess_ND, real_params;
    real_params = { 1.0, 2.0, 3.0 };
    guess_ND = { 1.1, 1.9, 3.3 };

    std::vector <double> innerx_ND;
    std::vector <std::vector<double> > x_ND(25, innerx_ND);

    int count = 0;

    //building x data
    for (int i = -2; i < 3; i++) {
        for (int j = -2; j < 3; j++) {
            x_ND[count].push_back(i);
            x_ND[count].push_back(j);

            y_ND.push_back(func_ND(x_ND[count], real_params));
            w_ND.push_back(1.0);

            count += 1;
        }
    };

    y_ND[6] = 2200.0; //intentionally bad data point
    w_ND[6] = 0.1;

    std::vector < std::vector <double> > gaussianParams_ND = { {NAN, NAN}, {1.0, 2.0}, {NAN, NAN} }; //params specifiying the Gaussian priors on the model function parameters
    std::vector < std::vector <double> > boundedParams_ND = { {0.0, NAN}, {NAN, NAN}, {NAN, NAN} };  //params specifiying the bounds (priors) on the model function parameters

    Priors testPriors_ND = Priors(MIXED, gaussianParams_ND, boundedParams_ND); //construction of priors object

    FunctionalForm model_ND = FunctionalForm(func_ND, x_ND, y_ND, partialsvector_ND, tolerance_ND, guess_ND, w_ND, testPriors_ND);
    RCR rcr_ND = RCR(SS_MEDIAN_DL); //setting up for RCR with this rejection technique
    rcr_ND.setParametricModel(model_ND);
    rcr_ND.performBulkRejection(w_ND, y_ND); //running Bulk rejection RCR

    std::vector <double> final_parameters_ND = model_ND.parameters; // final calculated model function parameters

    std::cout << "\n\nND parameters: \n"; //prints out calculated model function parameters
    for (int i = 0; i < final_parameters_ND.size(); i++) {
        std::cout << final_parameters_ND[i] << std::endl;
    }
    std::cout << "\n\n\n";

    std::vector <double> rejected_data_ND = rcr_ND.result.rejectedY;
    std::vector <double> nonrejected_data_ND = rcr_ND.result.cleanY;
    std::vector <bool> flags_ND = rcr_ND.result.flags; // booleans specifying whether each point, in order, was rejected; true for nonrejected, false for rejected.

    

    //EXAMPLE with priors, but weighted and no error bars (1D case)

    double tolerance = 0.01;
    std::vector <double> x, y, guess, w;
    guess = { 2.1, 0.9 };
    y = { 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 20.0, 11.0 };
    x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
    w = { 0.5, 1.0, 1.2, 3.0, 1.0, 0.6, 0.2, 1.0, 0.1, 2.0 };
    //std::vector <double> w(10, 1.0);

    double xb_guess = 4.5;

    std::vector < std::vector <double> > gaussianParams = { {NAN, NAN}, {1.0, 2.0} }; //params specifiying the Gaussian priors on the model function parameters
    std::vector < std::vector <double> > boundedParams = { {0.0, NAN}, {NAN, NAN} };  //params specifiying the bounds (priors) on the model function parameters

    Priors testPriors = Priors(MIXED, gaussianParams, boundedParams); //construction of priors object

    FunctionalForm model = FunctionalForm(func_linear, x, y, partialsvector_lin, tolerance, guess, w, testPriors, getAvg, xb_guess);// , getAvg, xb_guess); //setting up for functional form
    RCR rcr = RCR(SS_MEDIAN_DL); //setting up for RCR with this rejection technique
    rcr.setParametricModel(model);
    rcr.performBulkRejection(w, y); //running Bulk rejection RCR

    std::vector <double> final_parameters = model.parameters; // final calculated model function parameters

    std::cout << "\n\nparameters: \n"; //prints out calculated model function parameters
    for (int i = 0; i < final_parameters.size(); i++) {
        std::cout << final_parameters[i] << std::endl;
    }
    std::cout << "\n\n\n";

    std::vector <double> rejected_data = rcr.result.rejectedY;
    std::vector <double> nonrejected_data = rcr.result.cleanY;
    std::vector <bool> flags = rcr.result.flags; // booleans specifying whether each point, in order, was rejected; true for nonrejected, false for rejected.
};
