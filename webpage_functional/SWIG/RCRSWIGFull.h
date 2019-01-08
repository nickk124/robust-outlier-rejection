#include <fstream>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <float.h>
#include <iterator>
#include <vector>
#include <iostream>
#include <random>
#include <set>

std::vector <double> requestHandlerUnWeighted(std::vector <double> x, std::vector <double> y, std::vector <double> sigma_y, std::vector <double> guess, int fType, int dataSize, int rejTechNo);
std::vector <double> requestHandlerWeighted(std::vector <double> x, std::vector <double> y, std::vector <double> sigma_y, std::vector <double> guess, std::vector <double> w, int fType, int dataSize, int rejTechNo);
