//#include "stdafx.h"
#include "RCRwebutils.h"

std::vector <double(*)(double, std::vector <double>)> partialsvector_linear = { partial1_linear, partial2_linear };
std::vector <double(*)(double, std::vector <double>)> partialsvector_quadratic = { partial1_quadratic, partial2_quadratic, partial3_quadratic };
std::vector <double(*)(double, std::vector <double>)> partialsvector_cubic = { partial1_cubic, partial2_cubic, partial3_cubic, partial4_cubic };
std::vector <double(*)(double, std::vector <double>)> partialsvector_powerlaw = { partial1_powerlaw, partial2_powerlaw };
std::vector <double(*)(double, std::vector <double>)> partialsvector_exponential = { partial1_exponential, partial2_exponential };
std::vector <double(*)(double, std::vector <double>)> partialsvector_logarithmic = { partial1_logarithmic };

double xBar;
double lnx_Bar;

// EXAMPLE FUNCTIONS (with corresponding partials and vectors of said partials)

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

//std::vector <double(*)(double, std::vector <double>)> partialsvector_linear = { partial1_linear, partial2_linear };

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

	return x - xBar;
}

double partial3_quadratic(double x, std::vector <double> params) {

	return std::pow((x - xBar), 2.0);
}

//std::vector <double(*)(double, std::vector <double>)> partialsvector_quadratic = { partial1_quadratic, partial2_quadratic, partial3_quadratic };


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

//std::vector <double(*)(double, std::vector <double>)> partialsvector_cubic= { partial1_cubic, partial2_cubic, partial3_cubic, partial4_cubic };


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

//std::vector <double(*)(double, std::vector <double>)> partialsvector_powerlaw = { partial1_powerlaw, partial2_powerlaw};

// EXPONENTIAL

double function_exponential(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * std::exp(a1*(x - xBar));
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

//std::vector <double(*)(double, std::vector <double>)> partialsvector_exponential = { partial1_exponential, partial2_exponential };

// LOGARITHMIC

double function_logarithmic(double x, std::vector <double> params) {
	double a0 = params[0];

	return a0 * std::log(x - xBar);
}

double partial1_logarithmic(double x, std::vector <double> params) {
	double a0 = params[0];

	return std::log(x - xBar);
}

Priors getPriors(int priorsCheck, std::vector <double> priorsParams, std::vector <int> hasPriorsVec, std::vector <double> guess) {
	int M = guess.size();
	std::vector < std::vector <double> > gaussianParams;
	std::vector < std::vector <double> > paramBounds;
	std::vector <double> temp_params;
	Priors priorsResult;

	switch (priorsCheck) {
	case 1: //gaussian
		for (int i = 0; i < M; i++) {
			temp_params.clear();

			for (int j = 0; j < 2; j++) {
				if (hasPriorsVec[i * 2 + j] == 1) {//has a prior
					temp_params.push_back(priorsParams[i * 2 + j]);
				}
				else if (hasPriorsVec[i * 2 + j] == 0) { //no prior
					temp_params.push_back(NAN);
				}
			}
			gaussianParams.push_back(temp_params);
		}
		priorsResult = Priors(GAUSSIAN, gaussianParams);

		break;
	case 2: //constrained
		for (int i = 0; i < M; i++) {
			temp_params.clear();

			for (int j = 0; j < 2; j++) {
				if (hasPriorsVec[i * 2 + j] == 1) {//has a prior
					temp_params.push_back(priorsParams[i * 2 + j]);
				}
				else if (hasPriorsVec[i * 2 + j] == 0) { //no prior
					temp_params.push_back(NAN);
				}
			}
			paramBounds.push_back(temp_params);
		}
		priorsResult = Priors(CONSTRAINED, paramBounds);

		break;
	case 3: //mixed
		for (int i = 0; i < M; i++) {
			temp_params.clear();

			for (int j = 0; j < 2; j++) {
				if (hasPriorsVec[i * 2 + j] == 1) {//has a prior
					temp_params.push_back(priorsParams[i * 2 + j]);
				}
				else if (hasPriorsVec[i * 2 + j] == 0) { //no prior
					temp_params.push_back(NAN);
				}
			}

			gaussianParams.push_back(temp_params);
		}

		for (int i = 0; i < M; i++) {
			temp_params.clear();

			for (int j = 0; j < 2; j++) {
				if (hasPriorsVec[2*M + i * 2 + j] == 1) {//has a prior
					temp_params.push_back(priorsParams[i * 2 + j]);
				}
				else if (hasPriorsVec[2*M + i * 2 + j] == 0) { //no prior
					temp_params.push_back(NAN);
				}
			}
			paramBounds.push_back(temp_params);
		}
		priorsResult = Priors(MIXED, gaussianParams, paramBounds);

		break;

	default:
		break;
	}


	return priorsResult;
}

std::vector <double> requestHandlerUnWeighted(std::vector <double> x, std::vector <double> y, std::vector <double> guess, int fType, int dataSize, int rejTechNo, int priorsCheck, std::vector <double> priorsParams, std::vector <int> hasPriorsVec)
{
	//partialsvector_linear = { partial1_linear, partial2_linear };
	//partialsvector_quadratic = { partial1_quadratic, partial2_quadratic, partial3_quadratic };
	//partialsvector_cubic = { partial1_cubic, partial2_cubic, partial3_cubic, partial4_cubic };
	//partialsvector_powerlaw = { partial1_powerlaw, partial2_powerlaw };
	//partialsvector_exponential = { partial1_exponential, partial2_exponential };
	//partialsvector_logarithmic = { partial1_logarithmic };


	//tolerance:
	double tolerance = 0.01;

	//sets rejection technique:
	RejectionTechs rejTech;

	switch (rejTechNo) {
	case 1:
		rejTech = SS_MEDIAN_DL;
		break;
	case 2:
		rejTech = LS_MODE_68;
		break;
	case 3:
		rejTech = LS_MODE_DL;
		break;
	case 4:
		rejTech = ES_MODE_DL;
		break;
	}

	/*
	  // converts the data arrays to vectors:
	  std::vector <double> x(xarray, xarray + dataSize);
	  std::vector <double> y(yarray, yarray + dataSize);
	  std::vector <double> sigma_y(sigma_yarray, sigma_yarray + dataSize);
	*/

	// gets the xBar needed (or lnx_bar):
	if (fType == 4) { // power law case
		lnx_Bar = getLnX_Bar(x);
	}
	else { // otherwise
		xBar = getAvg(x);
	}

	if (fType == 1) { // linear

					  //converts guess array to vector
		/*
		std::vector <double> guess(guessarray, guessarray + 2);
	*/

	// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(y);



		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 2; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;

	}
	if (fType == 2) { // quadratic

					  //converts guess array to vector
		//std::vector <double> guess(guessarray, guessarray + 3);

		// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(y);

		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 3; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;
	}
	if (fType == 3) { // cubic

					  //converts guess array to vector
		//std::vector <double> guess(guessarray, guessarray + 4);

		// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(y);

		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 4; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;
	}
	if (fType == 4) { // power law

					  //converts guess array to vector
		//std::vector <double> guess(guessarray, guessarray + 2);

		// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(y);

		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 2; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;
	}
	if (fType == 5) { // exponential

					  //converts guess array to vector
		//std::vector <double> guess(guessarray, guessarray + 2);

		// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(y);

		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 2; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;
	}
	if (fType == 6) { // logarithmic

					  //converts guess array to vector
		//std::vector <double> guess(guessarray, guessarray + 1);

		// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(y);

		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 1; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;
	}
}

std::vector <double> requestHandlerWeighted(std::vector <double> x, std::vector <double> y, std::vector <double> guess, std::vector <double> w, int fType, int dataSize, int rejTechNo, int priorsCheck, std::vector <double> priorsParams, std::vector <int> hasPriorsVec)// handles the request for the functionalForm. fType is an character corresponding to a certain function type. takes in (C?) arrays that are then converted to usable cpp vectors.
{
	partialsvector_linear = { partial1_linear, partial2_linear };
	partialsvector_quadratic = { partial1_quadratic, partial2_quadratic, partial3_quadratic };
	partialsvector_cubic = { partial1_cubic, partial2_cubic, partial3_cubic, partial4_cubic };
	partialsvector_powerlaw = { partial1_powerlaw, partial2_powerlaw };
	partialsvector_exponential = { partial1_exponential, partial2_exponential };
	partialsvector_logarithmic = { partial1_logarithmic };


	//tolerance:
	double tolerance = 0.01;

	//sets rejection technique:
	RejectionTechs rejTech;

	switch (rejTechNo) {
	case 1:
		rejTech = SS_MEDIAN_DL;
		break;
	case 2:
		rejTech = LS_MODE_68;
		break;
	case 3:
		rejTech = LS_MODE_DL;
		break;
	case 4:
		rejTech = ES_MODE_DL;
		break;
	}


	// converts the data arrays to vectors:
	/*
	std::vector <double> x(xarray, xarray + dataSize);
	std::vector <double> y(yarray, yarray + dataSize);
	std::vector <double> sigma_y(sigma_yarray, sigma_yarray + dataSize);
	std::vector <double> w(warray, warray + dataSize);
  */

  //gets xBar / lnX_Bar
	if (fType == 4) { // power law case
		lnx_Bar = getLnX_Bar(x, w);
	}
	else { // otherwise
		xBar = getAvg(x, w);
	}

	if (fType == 1) { // linear

					  //converts guess array to vector
		//std::vector <double> guess(guessarray, guessarray + 2);

		// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(w, y);

		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 2; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;

	}
	if (fType == 2) { // quadratic

					  //converts guess array to vector
		//std::vector <double> guess(guessarray, guessarray + 3);

		// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(w, y);

		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 3; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;
	}
	if (fType == 3) { // cubic

					  //converts guess array to vector
		//std::vector <double> guess(guessarray, guessarray + 4);

		// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(w, y);

		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 4; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;
	}
	if (fType == 4) { // power law

					  //converts guess array to vector
		//std::vector <double> guess(guessarray, guessarray + 2);

		// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(w, y);

		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 2; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;
	}
	if (fType == 5) { // exponential

					  //converts guess array to vector
		//std::vector <double> guess(guessarray, guessarray + 2);

		// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(w, y);

		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 2; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;
	}
	if (fType == 6) { // logarithmic

					  //converts guess array to vector
		//std::vector <double> guess(guessarray, guessarray + 1);

		// performs functionalRCR:
		FunctionalForm model;
		if (priorsCheck == 0) { //no priors

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w);
		}
		else { //are priors
			Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, guess);

			model = FunctionalForm(function_linear, x, y, partialsvector_linear, tolerance, guess, w, priorsObj);
		}
		RCR rcr = RCR(rejTech);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(w, y);

		std::vector <double> result;
		//first, lists the final calculated params
		for (int j = 0; j < 1; j++) {
			result.push_back(model.parameters[j]);
		}
		result.push_back(0.0);
		for (int k = 0; k < dataSize; k++) { //adding booleans to the file, showing whether the point was rejected or not
			result.push_back(rcr.result.flags[k]);
		}

		return result;
	}
}
