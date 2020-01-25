//#include "stdafx.h"
#include "MiscFunctions.h"

// Function to get cofactor of A[p][q] in temp[][]. n is current
// dimension of A[][]
//void getCofactor(double A[M][M], double temp[M][M], int p, int q, int n)
std::vector<std::vector<double> > getCofactor(std::vector<std::vector <double> > A, int p, int q)
{
	int i = 0, j = 0;
	int M = (int) A[0].size();
	std::vector <double> innertemp(M - 1, 0.0);
	std::vector<std::vector <double> > temp(M - 1, innertemp);


	// Looping for each element of the matrix
	for (int row = 0; row < M; row++)
	{
		for (int col = 0; col < M; col++)
		{
			//  Copying into temporary matrix only those element
			//  which are not in given row and column
			if (row != p && col != q)
			{
				temp[i][j++] = A[row][col];

				// Row is filled, so increase row index and
				// reset col index
				if (j == M - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
	return temp;
}
//Recursive function for finding determinant of matrix.
//n is current dimension of A[][].
double determinant(std::vector < std::vector <double> > A)
{
	int M = (int) A[0].size();
	double D = 0; // Initialize result

				  //  Base case : if matrix contains single element
				  //std::cout << "n/M = " << n << std::endl;
	if (M == 1)
		return A[0][0];


	//std::vector<double> innertemp(M, 0.0);
	//std::vector <std::vector <double> > temp(M, innertemp); // To store cofactors

	int sign = 1;  // To store sign multiplier

				   // Iterate for each element of first row
	for (int f = 0; f < M; f++)
	{
		// Getting Cofactor of A[0][f]
		std::vector <std::vector <double> > temp = getCofactor(A, 0, f);
		D += sign * A[0][f] * determinant(temp);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}
// Function to get adjoint of A[M][M] in adj[M][M].
std::vector<std::vector<double> > adjoint(std::vector < std::vector <double> > A)
{
	int M = (int) A[0].size();
	if (M == 1)
	{
		std::vector<std::vector<double> > result = { { 1.0 } };
		return result;
	}

	std::vector <double> innertemp(M - 1, 0.0);
	std::vector<std::vector <double> > temp(M - 1, innertemp);

	std::vector <double> inneradj(M, 0.0);
	std::vector<std::vector <double> > adj(M, inneradj);

	// temp is used to store cofactors of A[][]
	int sign = 1;

	//std::vector <std::vector <double> > temp;

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < M; j++)
		{
			// Get cofactor of A[i][j]
			std::vector <std::vector <double> > temp = getCofactor(A, i, j);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j][i] = (sign)*(determinant(temp));
		}
	}
	return adj;
}
// Function to calculate and store inverse, returns false if
// matrix is singular
std::vector<std::vector <double> > inverse(std::vector < std::vector <double> > A)
{
	int M = (int) A[0].size();
	// Find determinant of A[][]
	double det = determinant(A);
	//std::cout << "Determinant = " << det << std::endl;
	if (det == 0.0)
	{
		std::cout << "Singular matrix, can't find its inverse";
		return std::vector<std::vector <double>> { {0} };
	}

	// Find adjoint

	std::vector <std::vector <double> > adj = adjoint(A);


	std::vector <double> innerinverse(M, 0.0);
	std::vector<std::vector <double> > inverse(M, innerinverse);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i = 0; i < M; i++)
		for (int j = 0; j < M; j++)
			inverse[i][j] = adj[i][j] / float(det);

	return inverse;
}

std::vector<std::vector<double> > pivotSystem(std::vector <std::vector<double> >Amatrix, std::vector<double> b) // columnCount is number of columns in the matrix A
{
	int columnCount = (int) Amatrix[0].size();

	std::vector <double> A; //vector form of the matrix A

	for (int i = 0; i < columnCount; i++) {
		for (int j = 0; j < columnCount; j++) {
			A.push_back(Amatrix[i][j]);
		}
	}

	int pivotIndex;
	double pivot, swap, multiplier;
	std::vector<std::vector<double> > tensorHold;

	for (int i = columnCount - 1; i > 0; i--)
	{
		pivot = std::abs(A[columnCount * i + i]);
		pivotIndex = i;
		for (int j = i; j >= 0; j--)
		{
			if (std::abs(A[columnCount * j + i]) > pivot)
			{
				pivot = std::abs(A[columnCount * j + i]);
				pivotIndex = j;
			}
		}
		if (pivotIndex == i)
		{
			for (int j = 0; j < columnCount; j++)
			{
				A[columnCount * i + j] = A[columnCount * i + j] / pivot;
			}
			b[i] = b[i] / pivot;
		}
		else
		{
			for (int j = 0; j < columnCount; j++)
			{
				swap = A[columnCount * i + j];
				A[columnCount * i + j] = A[columnCount * pivotIndex + j] / pivot;
				A[columnCount * pivotIndex + j] = swap;
			}
			swap = b[i];
			b[i] = b[pivotIndex] / pivot;
			b[pivotIndex] = swap;
		}
		for (int j = i - 1; j >= 0; j--)
		{
			multiplier = A[j*columnCount + i] / A[i*columnCount + i];
			for (int k = 0; k < columnCount; k++)
			{
				A[j*columnCount + k] -= multiplier*A[i*columnCount + k];
			}
			b[j] -= multiplier*b[i];
		}
	}

	tensorHold.resize(2);
	tensorHold[0].resize(A.size());
	tensorHold[1].resize(b.size());

	tensorHold[0] = A;
	tensorHold[1] = b;

	return tensorHold;
}

std::vector < std::vector <double> > transpose(std::vector < std::vector <double> > array) { //takes the transpose of some input array

	int n = (int) array.size();
	int m = (int) array[0].size();

	std::vector<double> innertransposedArray(n, 0.0);
	std::vector <std::vector <double> > transposedArray(m, innertransposedArray);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			transposedArray[j][i] = array[i][j];
	return transposedArray;
};

std::vector <double> dot(std::vector< std::vector <double> > A, std::vector <double> b) // simple square matrix-vector multiplication function
{
	int n = (int) A.size(); //dimensionality of square matrix/vector
	int m = (int) b.size();
	std::vector <double> c; //vector result

	for (int i = 0; i < n; i++) {
		double sum = 0.0;
		for (int j = 0; j < m; j++) {
			sum += A[i][j] * b[j];
		}
		c.push_back(sum);
	}
	return c;
}

std::vector< std::vector <double> > dot(std::vector< std::vector <double> > A, std::vector <std::vector <double> > B) // multiplication of two matrices (can be non-square)
{
	int n = (int) A.size();
	int m = (int) A[0].size();
	int p = (int) B[0].size();
	std::vector <double> innerC(p, 0.0);
	std::vector< std::vector <double> > C(n, innerC); //Matrix result C
	double sum;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < p; j++) {
			sum = 0.0;
			for (int k = 0; k < m; k++) {
				double added = (A[i][k]) * (B[k][j]);
				sum += added;
			}
			C[i][j] = sum;
		}
	}
	return C;
}

std::vector <double> forwardSubstitution(std::vector <std::vector <double> > A, std::vector <double> b)
{
	std::vector <double> x;
	int N = (int) b.size();

	x.push_back(b[0] / A[0][0]);
	for (int m = 1; m < N; m++) {
		double sum = 0.0;
		for (int n = 0; n < m; n++) {
			sum += (A[m][n] * x[n]);
		}
		x.push_back((b[m] - sum) / A[m][m]);
	}

	return x;
}

std::vector <std::vector <double> > LUInverse(std::vector <std::vector <double> > A)
{
	int m = (int) A.size();

	std::vector <double> innerresult_transposed;
	std::vector< std::vector <double > > result_transposed(m, innerresult_transposed);

	for (int i = 0; i < m; i++) {

		std::vector <double> b(m, 0.0);

		b[i] = 1.0;

		result_transposed[i] = forwardSubstitution(A, b);

	}

	std::vector< std::vector <double > > result = transpose(result_transposed);

	return result;


}

/*
std::vector < std::vector <double> > pseudoInverse(std::vector < std::vector <double> > A) // Computes the Moore-Penrose Pseudo-Inverse of some square Matrix A, used for the Gauss-Newton Method
{


std::vector < std::vector <double> > Atrans = transpose(A);
std::vector < std::vector <double> > dotted = dot(Atrans, A);
std::vector < std::vector <double> > inv = inverse(dotted);
std::vector < std::vector <double> > result = dot(inv, Atrans);


std::vector < std::vector <double> > Auw = dot(transpose(A), A);
buw

return result;
};

std::vector < std::vector <double> > weightedPseudoInverse(std::vector < std::vector <double> > A, std::vector < std::vector <double> > W) {

std::vector < std::vector <double> > Atrans = transpose(A);
std::vector < std::vector <double> > inver = inverse(dot(Atrans, dot(W, A)));
std::vector < std::vector <double> > result = dot(inver, Atrans);

return result;
}

*/
double chiSquared(double(*f)(double, std::vector <double>), std::vector <double> y, std::vector <double> x, std::vector <double> params, std::vector <double> sigma_y) // computes the chi-squared value. takes the function y=... as an argument (via pointer) 
{ //FUNCTION MUST BE USER-PROVIDED
	int N = (int) y.size();
	double sum = 0.0;

	for (int i = 0; i < N; i++) {
		sum += std::pow(((y[i] - (*f)(x[i], params) ) / sigma_y[i]), 2.0);
	}
	return sum;
}

double chiSquared(double(*f)(std::vector <double>, std::vector <double>), std::vector <double> y, std::vector <std::vector<double> > x, std::vector <double> params, std::vector <double> sigma_y) // multi-dimensional independent variables case
{ //FUNCTION MUST BE USER-PROVIDED
	int N = (int) y.size();
	double sum = 0.0;

	for (int i = 0; i < N; i++) {
		sum += std::pow(((y[i] - (*f)(x[i], params) )/ sigma_y[i]), 2.0);
	}
	return sum;

}

double chiSquared(double(*f)(double, std::vector <double>), std::vector <double> y, std::vector <double> x, std::vector <double> params, std::vector <double> w, int K) // computes the chi-squared value. takes the function y=... as an argument (via pointer) 
{ //FUNCTION MUST BE USER-PROVIDED
	//K = 0 if nonweighted, 1 if weighted
	int N = (int) y.size();

	if (K == 0) {
		w.clear();
		for (int i = 0; i < N; i++) {
			w.push_back(1.0);
		}
	}
	double sum = 0.0;

	for (int i = 0; i < N; i++) {
		sum += w[i] * std::pow(y[i] - (*f)(x[i], params), 2.0);
	}
	return sum;
}

double chiSquared(double(*f)(std::vector <double>, std::vector <double>), std::vector <double> y, std::vector <std::vector<double> > x, std::vector <double> params, std::vector <double> w, int K) // multi-dimensional independent variables case
{ //FUNCTION MUST BE USER-PROVIDED
	//K = 0 if nonweighted, 1 if weighted
	int N = (int) y.size();

	if (K == 0) {
		w.clear();
		for (int i = 0; i < N; i++) {
			w.push_back(1.0);
		}
	}
	double sum = 0.0;

	for (int i = 0; i < N; i++) {
		sum += w[i] * std::pow(y[i] - (*f)(x[i], params), 2.0);
	}
	return sum;
}

std::vector <double> residuals(double(*f)(double, std::vector <double>), std::vector <double> y, std::vector <double> x, std::vector <double> params) // computes the residuals vector needed in the GN algorithm
{ //FUNCTION MUST BE USER-PROVIDED
	int N = (int) y.size();
	std::vector <double> result(N, 0.0); // placeholder

	for (int i = 0; i < N; i++) {
		result[i] = ((*f)(x[i], params) - y[i]);
	}
	return result;
}

std::vector <double> residuals(double(*f)(std::vector <double>, std::vector <double>), std::vector <double> y, std::vector <std::vector<double> > x, std::vector <double> params) // multi-dimensional independent variables case
{ //FUNCTION MUST BE USER-PROVIDED
	int N = (int) y.size();
	std::vector <double> result(N, 0.0); // placeholder

	for (int i = 0; i < N; i++) {
		result[i] = ((*f)(x[i], params) - y[i]);
	}
	return result;
}

std::vector < std::vector <double> > jacobian(std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> x, std::vector <double> params) //Jacobian matrix, CASE of 1 independent var, takes argument of vector of parameters. USER INPUTTED.
{
	int M = (int) params.size(); //dimensionality
	int N = (int) x.size();

	std::vector <double> jacobrow(M, 0.0); //initializing the Jacobian
	std::vector < std::vector <double> > jacob(N, jacobrow);


	for (int i = 0; i < N; i++) { //looping through rows
		for (int j = 0; j < M; j++) { //looping through columns
			jacob[i][j] = parsvector[j](x[i], params);
		}
	}

	return jacob;

};

std::vector < std::vector <double> > jacobian(std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector< std::vector <double> > x, std::vector <double> params) //Jacobian matrix, CASE of >1 independent var, takes argument of vector of parameters. USER INPUTTED.
{
	int M = (int) params.size(); //dimensionality
	int N = (int) x.size();

	std::vector <double> jacobrow(M, 0.0); //initializing the Jacobian
	std::vector < std::vector <double> > jacob(N, jacobrow);


	for (int i = 0; i < N; i++) { //looping through rows
		for (int j = 0; j < M; j++) { //looping through columns
			jacob[i][j] = parsvector[j](x[i], params);
		}
	}

	return jacob;

};
//sy, w
std::vector <double> paramuncertainty(std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> x, std::vector <double> params, std::vector <double> sigma_y, std::vector <double> w, double wbar) //1D case
{
	int M = (int) params.size();

	std::vector <double> J_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > J_pivoted(M, J_pivoted_inner);

	std::vector <double> sigma_params(M, 0.0);
	std::vector <double> sig_y(M, 0.0);

	for (int j = 0; j < M; j++) {
		sig_y[j] = std::sqrt(wbar / w[j]) * sigma_y[j];
	}

	std::vector < std::vector <double> > J = jacobian(parsvector, x, params);

	std::vector < std::vector <double> > result = pivotSystem(J, sig_y);

	std::vector <double> J_pivoted_vec = result[0];
	std::vector <double> sig_y_pivoted = result[1];


	int d = 0; //turns the outputted pivoted A vector back into an A matrix
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			J_pivoted[i][j] = J_pivoted_vec[d];
			d += 1;
		}
	}

	std::vector <std::vector <double> > invJ_pivoted = LUInverse(J_pivoted);

	double sm;
	for (int k = 0; k < M; k++) {
		sm = 0.0;
		for (int j = 0; j < M; j++) {
			sm += std::pow((invJ_pivoted[k][j] * sig_y_pivoted[j]), 2.0); //sums in quadrature
		}
		sigma_params[k] = std::sqrt(sm);
	}
	return sigma_params;
}
//sy, w ND
std::vector <double> paramuncertainty(std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector< std::vector <double> > x, std::vector <double> params, std::vector <double> sigma_y, std::vector <double> w, double wbar) // >1D case
{
	int M = (int) params.size();

	std::vector <double> J_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > J_pivoted(M, J_pivoted_inner);

	std::vector <double> sigma_params(M, 0.0);
	std::vector <double> sig_y(M, 0.0);

	for (int j = 0; j < M; j++) {
		sig_y[j] = std::sqrt(wbar / w[j]) * sigma_y[j];
	}

	std::vector < std::vector <double> > J = jacobian(parsvector, x, params);

	std::vector < std::vector <double> > result = pivotSystem(J, sig_y);

	std::vector <double> J_pivoted_vec = result[0];
	std::vector <double> sig_y_pivoted = result[1];


	int d = 0; //turns the outputted pivoted A vector back into an A matrix
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			J_pivoted[i][j] = J_pivoted_vec[d];
			d += 1;
		}
	}

	std::vector <std::vector <double> > invJ_pivoted = LUInverse(J_pivoted);

	double sm;
	for (int k = 0; k < M; k++) {
		sm = 0.0;
		for (int j = 0; j < M; j++) {
			sm += std::pow((invJ_pivoted[k][j] * sig_y_pivoted[j]), 2.0); //sums in quadrature
		}
		sigma_params[k] = std::sqrt(sm);
	}
	return sigma_params;
}
//sy
std::vector <double> paramuncertainty(std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> x, std::vector <double> params, std::vector <double> sigma_y) //1D case
{
	int M = (int) params.size();

	std::vector <double> J_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > J_pivoted(M, J_pivoted_inner);

	std::vector <double> sigma_params(M, 0.0);

	std::vector <double> sig_y = sigma_y;

	std::vector < std::vector <double> > J = jacobian(parsvector, x, params);

	std::vector < std::vector <double> > result = pivotSystem(J, sig_y);

	std::vector <double> J_pivoted_vec = result[0];
	std::vector <double> sig_y_pivoted = result[1];


	int d = 0; //turns the outputted pivoted A vector back into an A matrix
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			J_pivoted[i][j] = J_pivoted_vec[d];
			d += 1;
		}
	}

	std::vector <std::vector <double> > invJ_pivoted = LUInverse(J_pivoted);


	double sm;
	for (int k = 0; k < M; k++) {
		sm = 0.0;
		for (int j = 0; j < M; j++) {
			sm += std::pow((invJ_pivoted[k][j] * sig_y_pivoted[j]), 2.0); //sums in quadrature
		}
		sigma_params[k] = std::sqrt(sm);
	}
	return sigma_params;
}
// sy, ND
std::vector <double> paramuncertainty(std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector< std::vector <double> > x, std::vector <double> params, std::vector <double> sigma_y) //>1D case
{
	int M = (int) params.size();

	std::vector <double> J_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > J_pivoted(M, J_pivoted_inner);

	std::vector <double> sigma_params(M, 0.0);

	std::vector <double> sig_y = sigma_y;

	std::vector < std::vector <double> > J = jacobian(parsvector, x, params);

	std::vector < std::vector <double> > result = pivotSystem(J, sig_y);

	std::vector <double> J_pivoted_vec = result[0];
	std::vector <double> sig_y_pivoted = result[1];


	int d = 0; //turns the outputted pivoted A vector back into an A matrix
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			J_pivoted[i][j] = J_pivoted_vec[d];
			d += 1;
		}
	}

	std::vector <std::vector <double> > invJ_pivoted = LUInverse(J_pivoted);

	double sm;
	for (int k = 0; k < M; k++) {
		sm = 0.0;
		for (int j = 0; j < M; j++) {
			sm += std::pow((invJ_pivoted[k][j] * sig_y_pivoted[j]), 2.0); //sums in quadrature
		}
		sigma_params[k] = std::sqrt(sm);
	}
	return sigma_params;
}

//w, no sy, ND
std::vector <double> paramuncertainty(std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> x, std::vector <double> params, std::vector <double> w, double wbar) //1D case
{
	int M = (int) params.size();

	std::vector <double> J_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > J_pivoted(M, J_pivoted_inner);

	std::vector <double> sigma_params(M, 0.0);
	std::vector <double> sig_y(M, 0.0);

	for (int j = 0; j < M; j++) {
		sig_y[j] = std::sqrt(wbar / w[j]);
	}

	std::vector < std::vector <double> > J = jacobian(parsvector, x, params);

	std::vector < std::vector <double> > result = pivotSystem(J, sig_y);

	std::vector <double> J_pivoted_vec = result[0];
	std::vector <double> sig_y_pivoted = result[1];


	int d = 0; //turns the outputted pivoted A vector back into an A matrix
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			J_pivoted[i][j] = J_pivoted_vec[d];
			d += 1;
		}
	}

	std::vector <std::vector <double> > invJ_pivoted = LUInverse(J_pivoted);

	double sm;
	for (int k = 0; k < M; k++) {
		sm = 0.0;
		for (int j = 0; j < M; j++) {
			sm += std::pow((invJ_pivoted[k][j] * sig_y_pivoted[j]), 2.0); //sums in quadrature
		}
		sigma_params[k] = std::sqrt(sm);
	}
	return sigma_params;
}
//w, no sy, ND
std::vector <double> paramuncertainty(std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector< std::vector <double> > x, std::vector <double> params, std::vector <double> w, double wbar) // >1D case
{
	int M = (int) params.size();

	std::vector <double> J_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > J_pivoted(M, J_pivoted_inner);

	std::vector <double> sigma_params(M, 0.0);
	std::vector <double> sig_y(M, 0.0);

	for (int j = 0; j < M; j++) {
		sig_y[j] = std::sqrt(wbar / w[j]);
	}

	std::vector < std::vector <double> > J = jacobian(parsvector, x, params);

	std::vector < std::vector <double> > result = pivotSystem(J, sig_y);

	std::vector <double> J_pivoted_vec = result[0];
	std::vector <double> sig_y_pivoted = result[1];


	int d = 0; //turns the outputted pivoted A vector back into an A matrix
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			J_pivoted[i][j] = J_pivoted_vec[d];
			d += 1;
		}
	}

	std::vector <std::vector <double> > invJ_pivoted = LUInverse(J_pivoted);

	double sm;
	for (int k = 0; k < M; k++) {
		sm = 0.0;
		for (int j = 0; j < M; j++) {
			sm += std::pow((invJ_pivoted[k][j] * sig_y_pivoted[j]), 2.0); //sums in quadrature
		}
		sigma_params[k] = std::sqrt(sm);
	}
	return sigma_params;
}
// no w or sy
std::vector <double> paramuncertainty(std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> x, std::vector <double> params) //1D case
{
	int M = (int) params.size();

	std::vector <double> J_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > J_pivoted(M, J_pivoted_inner);

	std::vector <double> sigma_params(M, 0.0);

	std::vector <double> sig_y(M, 1.0);

	std::vector < std::vector <double> > J = jacobian(parsvector, x, params);

	std::vector < std::vector <double> > result = pivotSystem(J, sig_y);

	std::vector <double> J_pivoted_vec = result[0];
	std::vector <double> sig_y_pivoted = result[1];


	int d = 0; //turns the outputted pivoted A vector back into an A matrix
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			J_pivoted[i][j] = J_pivoted_vec[d];
			d += 1;
		}
	}

	std::vector <std::vector <double> > invJ_pivoted = LUInverse(J_pivoted);


	double sm;
	for (int k = 0; k < M; k++) {
		sm = 0.0;
		for (int j = 0; j < M; j++) {
			sm += std::pow((invJ_pivoted[k][j] * sig_y_pivoted[j]), 2.0); //sums in quadrature
		}
		sigma_params[k] = std::sqrt(sm);
	}
	return sigma_params;
}
// no w or sy, ND
std::vector <double> paramuncertainty(std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector< std::vector <double> > x, std::vector <double> params) //>1D case
{
	int M = (int) params.size();

	std::vector <double> J_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > J_pivoted(M, J_pivoted_inner);

	std::vector <double> sigma_params(M, 0.0);

	std::vector <double> sig_y(M, 1.0);

	std::vector < std::vector <double> > J = jacobian(parsvector, x, params);

	std::vector < std::vector <double> > result = pivotSystem(J, sig_y);

	std::vector <double> J_pivoted_vec = result[0];
	std::vector <double> sig_y_pivoted = result[1];


	int d = 0; //turns the outputted pivoted A vector back into an A matrix
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			J_pivoted[i][j] = J_pivoted_vec[d];
			d += 1;
		}
	}

	std::vector <std::vector <double> > invJ_pivoted = LUInverse(J_pivoted);

	double sm;
	for (int k = 0; k < M; k++) {
		sm = 0.0;
		for (int j = 0; j < M; j++) {
			sm += std::pow((invJ_pivoted[k][j] * sig_y_pivoted[j]), 2.0); //sums in quadrature
		}
		sigma_params[k] = std::sqrt(sm);
	}
	return sigma_params;
}
//NON-WEIGHTED, with error bars:
std::vector <double> modifiedGN(double(*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, std::vector <double> sigma_y, double tolerance) //the case of 1 independent (x) variables in the function
{
	int M = (int) parsvector.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_uw, result, A_uw_pivoted_inv;
	std::vector <double> res, b_uw, b_uw_pivoted, A_uw_pivoted_vec, incrementvect;
	std::vector <double> A_uw_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_uw_pivoted(M, A_uw_pivoted_inner);
	double newChiSq, oldChiSq;
	double multiplier = 1.0;
	int checkcount;
//	int itercount = 0;

	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_uw = dot(transpose(J), J);
		b_uw = dot(transpose(J), res);
		result = pivotSystem(A_uw, b_uw);
		A_uw_pivoted_vec = result[0];
		b_uw_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_uw_pivoted[i][j] = A_uw_pivoted_vec[d];
				d += 1;
			}
		}

		//A_uw_pivoted_inv = inverse(A_uw_pivoted);
		//incrementvect = dot(A_uw_pivoted_inv, b_uw_pivoted);
		incrementvect = forwardSubstitution(A_uw_pivoted, b_uw_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - multiplier * incrementvect[i];
		}

		oldChiSq = chiSquared((*f), y, x, params_old, sigma_y);
		newChiSq = chiSquared((*f), y, x, params_new, sigma_y);

		int samecheck = 0;
		for (int j = 0; j < M; j++) {
			if (params_old[j] != guess[j]) {
				samecheck += 1;
			}
		}

		if ((newChiSq != newChiSq) && (samecheck == 0)) //if the iteration messes up due to singular matrices (this is only true if newChiSq is NAN
		{
			std::vector <double> badvec;

			return badvec; //returns an empty vector if this worst-case is true
		}

		else if ((newChiSq != newChiSq) && (samecheck >= 0)) // due to the data, one or more of the parameters runs off to "infinity"
		{
			params_old.push_back(0.0); //add an extra element to params_old to signify that this happened

			return params_old;
		}

		if (std::abs(newChiSq - oldChiSq) <= tolerance) {

			int sametcheck = 0;
			for (int j = 0; j < M; j++) {
				if (std::abs(params_new[j] - guess[j]) <= tolerance) {
					sametcheck += 1;
				}
			}
			if (sametcheck == M) { //another bad case
				params_old.push_back(0.0);
				params_old.push_back(0.0);
				return params_old;
			}
			return params_new;
		}

		if (newChiSq <= oldChiSq) //good case
		{
			params_old = params_new; // applies increment vector
		}
		else if (newChiSq > oldChiSq && (std::abs(newChiSq - oldChiSq) > tolerance)) //bad case
		{
			multiplier *= 0.5;
		}

		//itercount += 1;
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};

std::vector <double> modifiedGN(double(*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, std::vector <double> sigma_y, double tolerance) //the case of >1 independent (x) variables in the function
{
	int M = (int) parsvector.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_uw, result, A_uw_pivoted_inv;
	std::vector <double> res, b_uw, b_uw_pivoted, A_uw_pivoted_vec, incrementvect;
	std::vector <double> A_uw_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_uw_pivoted(M, A_uw_pivoted_inner);
	double newChiSq, oldChiSq;
	double multiplier = 1.0;
	int checkcount;
//	int itercount = 0;

	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_uw = dot(transpose(J), J);
		b_uw = dot(transpose(J), res);
		result = pivotSystem(A_uw, b_uw);
		A_uw_pivoted_vec = result[0];
		b_uw_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_uw_pivoted[i][j] = A_uw_pivoted_vec[d];
				d += 1;
			}
		}

		//A_uw_pivoted_inv = inverse(A_uw_pivoted);
		//incrementvect = dot(A_uw_pivoted_inv, b_uw_pivoted);
		incrementvect = forwardSubstitution(A_uw_pivoted, b_uw_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - multiplier * incrementvect[i];
		}

		oldChiSq = chiSquared((*f), y, x, params_old, sigma_y);
		newChiSq = chiSquared((*f), y, x, params_new, sigma_y);

		int samecheck = 0;
		for (int j = 0; j < M; j++) {
			if (params_old[j] != guess[j]) {
				samecheck += 1;
			}
		}

		if ((newChiSq != newChiSq) && (samecheck == 0)) //if the iteration messes up due to singular matrices (this is only true if newChiSq is NAN
		{
			std::vector <double> badvec;

			return badvec; //returns an empty vector if this worst-case is true
		}

		else if ((newChiSq != newChiSq) && (samecheck >= 0)) // due to the data, one or more of the parameters runs off to "infinity"
		{
			params_old.push_back(0.0); //add an extra element to params_old to signify that this happened

			return params_old; //returns an empty vector if this worst-case is true
		}

		if (std::abs(newChiSq - oldChiSq) <= tolerance) {

			int sametcheck = 0;
			for (int j = 0; j < M; j++) {
				if (std::abs(params_new[j] - guess[j]) <= tolerance) {
					sametcheck += 1;
				}
			}
			if (sametcheck == M) { //another bad case
				params_old.push_back(0.0);
				params_old.push_back(0.0);
				return params_old;
			}
			return params_new;
		}

		if (newChiSq <= oldChiSq) //good case
		{
			params_old = params_new; // applies increment vector
		}
		else if (newChiSq > oldChiSq && (std::abs(newChiSq - oldChiSq) > tolerance)) //bad case
		{
			multiplier *= 0.5;
		}
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};
//WEIGHTED, with error bars:
std::vector <double> modifiedGN(double(*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, std::vector <double> sigma_y, double tolerance, std::vector <double> w) //the case of 1 independent (x) variables in the function
{
	int Nr = (int) y.size();
	int M = (int) parsvector.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_w, result, A_w_pivoted_inv;
	std::vector <double> res, b_w, b_w_pivoted, A_w_pivoted_vec, incrementvect;
	std::vector <double> A_w_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_w_pivoted(M, A_w_pivoted_inner);
	double newChiSq, oldChiSq;
	double multiplier = 1.0;
	int checkcount;
//	int itercount = 0;

	// creating weight matrix W:
	std::vector <double> innerW(Nr, 0.0);
	std::vector <std::vector <double> > W(Nr, innerW);
	for (int k = 0; k < Nr; k++) {
		W[k][k] = w[k];
	}
	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_w = dot(transpose(J), dot(W, J));
		b_w = dot(transpose(J), dot(W, res));
		result = pivotSystem(A_w, b_w);
		A_w_pivoted_vec = result[0];
		b_w_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_w_pivoted[i][j] = A_w_pivoted_vec[d];
				d += 1;
			}
		}

		//A_w_pivoted_inv = inverse(A_w_pivoted);
		//incrementvect = dot(A_w_pivoted_inv, b_w_pivoted);
		incrementvect = forwardSubstitution(A_w_pivoted, b_w_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - multiplier * incrementvect[i];
		}

		oldChiSq = chiSquared((*f), y, x, params_old, sigma_y);
		newChiSq = chiSquared((*f), y, x, params_new, sigma_y);

		int samecheck = 0;
		for (int j = 0; j < M; j++) {
			if (params_old[j] != guess[j]) {
				samecheck += 1;
			}
		}

		if ((newChiSq != newChiSq) && (samecheck == 0)) //if the iteration messes up due to singular matrices (this is only true if newChiSq is NAN
		{
			std::vector <double> badvec;

			return badvec; //returns an empty vector if this worst-case is true
		}

		else if ((newChiSq != newChiSq) && (samecheck >= 0)) // due to the data, one or more of the parameters runs off to "infinity"
		{
			params_old.push_back(0.0); //add an extra element to params_old to signify that this happened

			return params_old; //returns an empty vector if this worst-case is true
		}

		if (std::abs(newChiSq - oldChiSq) <= tolerance) {

			int sametcheck = 0;
			for (int j = 0; j < M; j++) {
				if (std::abs(params_new[j] - guess[j]) <= tolerance) {
					sametcheck += 1;
				}
			}
			if (sametcheck == M) { //another bad case
				params_old.push_back(0.0);
				params_old.push_back(0.0);
				return params_old;
			}
			return params_new;
		}

		if (newChiSq <= oldChiSq) //good case
		{
			params_old = params_new; // applies increment vector
		}
		else if (newChiSq > oldChiSq && (std::abs(newChiSq - oldChiSq) > tolerance)) //bad case
		{
			multiplier *= 0.5;
		}
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};

std::vector <double> modifiedGN(double(*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, std::vector <double> sigma_y, double tolerance, std::vector <double> w) //the case of >1 independent (x) variables in the function
{
	int Nr = (int) y.size();
	int M = (int) parsvector.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_w, result, A_w_pivoted_inv;
	std::vector <double> res, b_w, b_w_pivoted, A_w_pivoted_vec, incrementvect;
	std::vector <double> A_w_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_w_pivoted(M, A_w_pivoted_inner);
	double newChiSq, oldChiSq;
	double multiplier = 1.0;
	int checkcount;
//	int itercount = 0;

	// creating weight matrix W:
	std::vector <double> innerW(Nr, 0.0);
	std::vector <std::vector <double> > W(Nr, innerW);
	for (int k = 0; k < Nr; k++) {
		W[k][k] = w[k];
	}
	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_w = dot(transpose(J), dot(W, J));
		b_w = dot(transpose(J), dot(W, res));
		result = pivotSystem(A_w, b_w);
		A_w_pivoted_vec = result[0];
		b_w_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_w_pivoted[i][j] = A_w_pivoted_vec[d];
				d += 1;
			}
		}

		//A_w_pivoted_inv = inverse(A_w_pivoted);
		//incrementvect = dot(A_w_pivoted_inv, b_w_pivoted);
		incrementvect = forwardSubstitution(A_w_pivoted, b_w_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - multiplier * incrementvect[i];
		}

		oldChiSq = chiSquared((*f), y, x, params_old, sigma_y);
		newChiSq = chiSquared((*f), y, x, params_new, sigma_y);

		int samecheck = 0;
		for (int j = 0; j < M; j++) {
			if (params_old[j] != guess[j]) {
				samecheck += 1;
			}
		}

		if ((newChiSq != newChiSq) && (samecheck == 0)) //if the iteration messes up due to singular matrices (this is only true if newChiSq is NAN
		{
			std::vector <double> badvec;

			return badvec; //returns an empty vector if this worst-case is true
		}

		else if ((newChiSq != newChiSq) && (samecheck >= 0)) // due to the data, one or more of the parameters runs off to "infinity"
		{
			params_old.push_back(0.0); //add an extra element to params_old to signify that this happened

			return params_old; //returns an empty vector if this worst-case is true
		}

		if (std::abs(newChiSq - oldChiSq) <= tolerance) {

			int sametcheck = 0;
			for (int j = 0; j < M; j++) {
				if (std::abs(params_new[j] - guess[j]) <= tolerance) {
					sametcheck += 1;
				}
			}
			if (sametcheck == M) { //another bad case
				params_old.push_back(0.0);
				params_old.push_back(0.0);
				return params_old;
			}
			return params_new;
		}

		if (newChiSq <= oldChiSq) //good case
		{
			params_old = params_new; // applies increment vector
		}
		else if (newChiSq > oldChiSq && (std::abs(newChiSq - oldChiSq) > tolerance)) //bad case
		{
			multiplier *= 0.5;
		}

		//itercount += 1;
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};
//NON-WEIGHTED, without error bars:
std::vector <double> modifiedGN(double(*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, double tolerance) //the case of 1 independent (x) variables in the function
{
	int M = (int) parsvector.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_uw, result, A_uw_pivoted_inv;
	std::vector <double> res, b_uw, b_uw_pivoted, A_uw_pivoted_vec, incrementvect;
	std::vector <double> A_uw_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_uw_pivoted(M, A_uw_pivoted_inner);
	double newChiSq, oldChiSq;
	double multiplier = 1.0;
	int checkcount;
//	int itercount = 0;

	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_uw = dot(transpose(J), J);
		b_uw = dot(transpose(J), res);
		result = pivotSystem(A_uw, b_uw);
		A_uw_pivoted_vec = result[0];
		b_uw_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_uw_pivoted[i][j] = A_uw_pivoted_vec[d];
				d += 1;
			}
		}

		//A_uw_pivoted_inv = inverse(A_uw_pivoted);
		//incrementvect = dot(A_uw_pivoted_inv, b_uw_pivoted);
		incrementvect = forwardSubstitution(A_uw_pivoted, b_uw_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - multiplier * incrementvect[i];
		}

		std::vector <double> w_temp;

		oldChiSq = chiSquared((*f), y, x, params_old, w_temp, 0);
		newChiSq = chiSquared((*f), y, x, params_new, w_temp, 0);

		int samecheck = 0;
		for (int j = 0; j < M; j++) {
			if (params_old[j] != guess[j]) {
				samecheck += 1;
			}
		}

		if ((newChiSq != newChiSq) && (samecheck == 0)) //if the iteration messes up due to singular matrices (this is only true if newChiSq is NAN
		{
			std::vector <double> badvec;

			return badvec; //returns an empty vector if this worst-case is true
		}

		else if ((newChiSq != newChiSq) && (samecheck >= 0)) // due to the data, one or more of the parameters runs off to "infinity"
		{
			params_old.push_back(0.0); //add an extra element to params_old to signify that this happened

			return params_old;
		}

		if (std::abs(newChiSq - oldChiSq) <= tolerance) {

			int sametcheck = 0;
			for (int j = 0; j < M; j++) {
				if (std::abs(params_new[j] - guess[j]) <= tolerance) {
					sametcheck += 1;
				}
			}
			if (sametcheck == M) { //another bad case
				params_old.push_back(0.0);
				params_old.push_back(0.0);
				return params_old;
			}
			return params_new;
		}

		if (newChiSq <= oldChiSq) //good case
		{
			params_old = params_new; // applies increment vector
		}
		else if (newChiSq > oldChiSq && (std::abs(newChiSq - oldChiSq) > tolerance)) //bad case
		{
			multiplier *= 0.5;
		}

		//itercount += 1;
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};

std::vector <double> modifiedGN(double(*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, double tolerance) //the case of >1 independent (x) variables in the function
{
	int M = (int) parsvector.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_uw, result, A_uw_pivoted_inv;
	std::vector <double> res, b_uw, b_uw_pivoted, A_uw_pivoted_vec, incrementvect;
	std::vector <double> A_uw_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_uw_pivoted(M, A_uw_pivoted_inner);
	double newChiSq, oldChiSq;
	double multiplier = 1.0;
	int checkcount;
//	int itercount = 0;

	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_uw = dot(transpose(J), J);
		b_uw = dot(transpose(J), res);
		result = pivotSystem(A_uw, b_uw);
		A_uw_pivoted_vec = result[0];
		b_uw_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_uw_pivoted[i][j] = A_uw_pivoted_vec[d];
				d += 1;
			}
		}

		//A_uw_pivoted_inv = inverse(A_uw_pivoted);
		//incrementvect = dot(A_uw_pivoted_inv, b_uw_pivoted);
		incrementvect = forwardSubstitution(A_uw_pivoted, b_uw_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - multiplier * incrementvect[i];
		}

		std::vector <double> w_temp;

		oldChiSq = chiSquared((*f), y, x, params_old, w_temp, 0);
		newChiSq = chiSquared((*f), y, x, params_new, w_temp, 0);



		int samecheck = 0;
		for (int j = 0; j < M; j++) {
			if (params_old[j] != guess[j]) {
				samecheck += 1;
			}
		}

		if ((newChiSq != newChiSq) && (samecheck == 0)) //if the iteration messes up due to singular matrices (this is only true if newChiSq is NAN
		{
			std::vector <double> badvec;

			return badvec; //returns an empty vector if this worst-case is true
		}

		else if ((newChiSq != newChiSq) && (samecheck >= 0)) // due to the data, one or more of the parameters runs off to "infinity"
		{
			params_old.push_back(0.0); //add an extra element to params_old to signify that this happened

			return params_old; //returns an empty vector if this worst-case is true
		}

		if (std::abs(newChiSq - oldChiSq) <= tolerance) {

			int sametcheck = 0;
			for (int j = 0; j < M; j++) {
				if (std::abs(params_new[j] - guess[j]) <= tolerance) {
					sametcheck += 1;
				}
			}
			if (sametcheck == M) { //another bad case
				params_old.push_back(0.0);
				params_old.push_back(0.0);
				return params_old;
			}
			return params_new;
		}

		if (newChiSq <= oldChiSq) //good case
		{
			params_old = params_new; // applies increment vector
		}
		else if (newChiSq > oldChiSq && (std::abs(newChiSq - oldChiSq) > tolerance)) //bad case
		{
			multiplier *= 0.5;
		}
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};
//WEIGHTED, without error bars:
std::vector <double> modifiedGN(double(*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, double tolerance, std::vector <double> w) //the case of 1 independent (x) variables in the function
{
	int Nr = (int) y.size();
	int M = (int) parsvector.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_w, result, A_w_pivoted_inv;
	std::vector <double> res, b_w, b_w_pivoted, A_w_pivoted_vec, incrementvect;
	std::vector <double> A_w_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_w_pivoted(M, A_w_pivoted_inner);
	double newChiSq, oldChiSq;
	double multiplier = 1.0;
	int checkcount;
//	int itercount = 0;

	// creating weight matrix W:
	std::vector <double> innerW(Nr, 0.0);
	std::vector <std::vector <double> > W(Nr, innerW);
	for (int k = 0; k < Nr; k++) {
		W[k][k] = w[k];
	}
	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_w = dot(transpose(J), dot(W, J));
		b_w = dot(transpose(J), dot(W, res));
		result = pivotSystem(A_w, b_w);
		A_w_pivoted_vec = result[0];
		b_w_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_w_pivoted[i][j] = A_w_pivoted_vec[d];
				d += 1;
			}
		}

		//A_w_pivoted_inv = inverse(A_w_pivoted);
		//incrementvect = dot(A_w_pivoted_inv, b_w_pivoted);
		incrementvect = forwardSubstitution(A_w_pivoted, b_w_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - multiplier * incrementvect[i];
		}

		oldChiSq = chiSquared((*f), y, x, params_old, w, 1);
		newChiSq = chiSquared((*f), y, x, params_new, w, 1);

		int samecheck = 0;
		for (int j = 0; j < M; j++) {
			if (params_old[j] != guess[j]) {
				samecheck += 1;
			}
		}

		if ((newChiSq != newChiSq) && (samecheck == 0)) //if the iteration messes up due to singular matrices (this is only true if newChiSq is NAN
		{
			std::vector <double> badvec;

			return badvec; //returns an empty vector if this worst-case is true
		}

		else if ((newChiSq != newChiSq) && (samecheck >= 0)) // due to the data, one or more of the parameters runs off to "infinity"
		{
			params_old.push_back(0.0); //add an extra element to params_old to signify that this happened

			return params_old; //returns an empty vector if this worst-case is true
		}

		if (std::abs(newChiSq - oldChiSq) <= tolerance) {

			int sametcheck = 0;
			for (int j = 0; j < M; j++) {
				if (std::abs(params_new[j] - guess[j]) <= tolerance) {
					sametcheck += 1;
				}
			}
			if (sametcheck == M) { //another bad case
				params_old.push_back(0.0);
				params_old.push_back(0.0);
				return params_old;
			}
			return params_new;
		}

		if (newChiSq <= oldChiSq) //good case
		{
			params_old = params_new; // applies increment vector
		}
		else if (newChiSq > oldChiSq && (std::abs(newChiSq - oldChiSq) > tolerance)) //bad case
		{
			multiplier *= 0.5;
		}
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};

std::vector <double> modifiedGN(double(*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, double tolerance, std::vector <double> w) //the case of >1 independent (x) variables in the function
{
	int Nr = (int) y.size();
	int M = (int) parsvector.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_w, result, A_w_pivoted_inv;
	std::vector <double> res, b_w, b_w_pivoted, A_w_pivoted_vec, incrementvect;
	std::vector <double> A_w_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_w_pivoted(M, A_w_pivoted_inner);
	double newChiSq, oldChiSq;
	double multiplier = 1.0;
	int checkcount;

	// creating weight matrix W:
	std::vector <double> innerW(Nr, 0.0);
	std::vector <std::vector <double> > W(Nr, innerW);
	for (int k = 0; k < Nr; k++) {
		W[k][k] = w[k];
	}
	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_w = dot(transpose(J), dot(W, J));
		b_w = dot(transpose(J), dot(W, res));
		result = pivotSystem(A_w, b_w);
		A_w_pivoted_vec = result[0];
		b_w_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_w_pivoted[i][j] = A_w_pivoted_vec[d];
				d += 1;
			}
		}

		//A_w_pivoted_inv = inverse(A_w_pivoted);
		//incrementvect = dot(A_w_pivoted_inv, b_w_pivoted);
		incrementvect = forwardSubstitution(A_w_pivoted, b_w_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - multiplier * incrementvect[i];
		}

		oldChiSq = chiSquared((*f), y, x, params_old, w, 1);
		newChiSq = chiSquared((*f), y, x, params_new, w, 1);

		int samecheck = 0;
		for (int j = 0; j < M; j++) {
			if (params_old[j] != guess[j]) {
				samecheck += 1;
			}
		}

		if ((newChiSq != newChiSq) && (samecheck == 0)) //if the iteration messes up due to singular matrices (this is only true if newChiSq is NAN
		{
			std::vector <double> badvec;

			return badvec; //returns an empty vector if this worst-case is true
		}

		else if ((newChiSq != newChiSq) && (samecheck >= 0)) // due to the data, one or more of the parameters runs off to "infinity"
		{
			params_old.push_back(0.0); //add an extra element to params_old to signify that this happened

			return params_old; //returns an empty vector if this worst-case is true
		}

		if (std::abs(newChiSq - oldChiSq) <= tolerance) {

			int sametcheck = 0;
			for (int j = 0; j < M; j++) {
				if (std::abs(params_new[j] - guess[j]) <= tolerance) {
					sametcheck += 1;
				}
			}
			if (sametcheck == M) { //another bad case
				params_old.push_back(0.0);
				params_old.push_back(0.0);
				return params_old;
			}
			return params_new;
		}

		if (newChiSq <= oldChiSq) //good case
		{
			params_old = params_new; // applies increment vector
		}
		else if (newChiSq > oldChiSq && (std::abs(newChiSq - oldChiSq) > tolerance)) //bad case
		{
			multiplier *= 0.5;
		}

		//itercount += 1;
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};




//for generalized mean part of code, modifications removed:
//NON-WEIGHTED:
std::vector <double> regularGN(double(*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, double tolerance) //the case of 1 independent (x) variables in the function
{
	int M = (int) parsvector.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_uw, result, A_uw_pivoted_inv;
	std::vector <double> res, b_uw, b_uw_pivoted, A_uw_pivoted_vec, incrementvect;
	std::vector <double> A_uw_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_uw_pivoted(M, A_uw_pivoted_inner);
	int checkcount;
	int itercount = 0;

	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_uw = dot(transpose(J), J);
		b_uw = dot(transpose(J), res);

		result = pivotSystem(A_uw, b_uw);
		A_uw_pivoted_vec = result[0];
		b_uw_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_uw_pivoted[i][j] = A_uw_pivoted_vec[d];
				d += 1;
			}
		}

		//A_uw_pivoted_inv = inverse(A_uw_pivoted);
		//incrementvect = dot(A_uw_pivoted_inv, b_uw_pivoted);
		incrementvect = forwardSubstitution(A_uw_pivoted, b_uw_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - incrementvect[i];
		}

		for (int i = 0; i < M; i++) // checks for NaNs/ blowing up params
		{
			if ((params_new[i] != params_new[i]) || std::isinf(params_new[i])) {
				params_old.push_back(0.0);
				return params_old;
			}
		}

		for (int j = 0; j < M; j++) // stops the iteration if the difference between each element in the latest param vector with it's previous counterpart are all less than the specificed tolerance
		{
			if (std::abs(params_new[j] - params_old[j]) <= tolerance) {
				checkcount += 1;
			}

		}

		if (checkcount == M) {
			tolcheck = false;
		}
		params_old = params_new; // applies increment vector

		itercount += 1;
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};

std::vector <double> regularGN(double(*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, double tolerance) //the case of >1 independent (x) variables in the function
{
	int M = (int) parsvector.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_uw, result, A_uw_pivoted_inv;
	std::vector <double> res, b_uw, b_uw_pivoted, A_uw_pivoted_vec, incrementvect;
	std::vector <double> A_uw_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_uw_pivoted(M, A_uw_pivoted_inner);
	int checkcount;
	int itercount = 0;

	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_uw = dot(transpose(J), J);
		b_uw = dot(transpose(J), res);
		result = pivotSystem(A_uw, b_uw);
		A_uw_pivoted_vec = result[0];
		b_uw_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_uw_pivoted[i][j] = A_uw_pivoted_vec[d];
				d += 1;
			}
		}

		//A_uw_pivoted_inv = inverse(A_uw_pivoted);
		//incrementvect = dot(A_uw_pivoted_inv, b_uw_pivoted);
		incrementvect = forwardSubstitution(A_uw_pivoted, b_uw_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - incrementvect[i];
		}

		for (int i = 0; i < M; i++) // checks for NaNs/ blowing up params
		{
			if ((params_new[i] != params_new[i]) || std::isinf(params_new[i])) {
				params_old.push_back(0.0);
				return params_old;
			}
		}

		for (int j = 0; j < M; j++) // stops the iteration if the difference between each element in the latest param vector with it's previous counterpart are all less than the specificed tolerance
		{
			if (std::abs(params_new[j] - params_old[j]) <= tolerance) {
				checkcount += 1;
			}

		}

		if (checkcount == M) {
			tolcheck = false;
		}
		params_old = params_new; // applies increment vector

		itercount += 1;
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};
//WEIGHTED:
std::vector <double> regularGN(double(*f)(double, std::vector <double>), std::vector <double(*)(double, std::vector <double>)> parsvector, std::vector <double> y, std::vector <double> x, std::vector <double> guess, double tolerance, std::vector <double> w) //the case of 1 independent (x) variables in the function
{
	int M = (int) parsvector.size();
	int Nr = (int) y.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_w, result, A_w_pivoted_inv;
	std::vector <double> res, b_w, b_w_pivoted, A_w_pivoted_vec, incrementvect;
	std::vector <double> A_w_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_w_pivoted(M, A_w_pivoted_inner);
	int checkcount;
	int itercount = 0;

	// creating weight matrix W:
	std::vector <double> innerW(Nr, 0.0);
	std::vector <std::vector <double> > W(Nr, innerW);
	for (int k = 0; k < Nr; k++) {
		W[k][k] = w[k];
	}
	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_w = dot(transpose(J), dot(W, J));
		b_w = dot(transpose(J), dot(W, res));
		result = pivotSystem(A_w, b_w);
		A_w_pivoted_vec = result[0];
		b_w_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_w_pivoted[i][j] = A_w_pivoted_vec[d];
				d += 1;
			}
		}

		//A_w_pivoted_inv = inverse(A_w_pivoted);
		//incrementvect = dot(A_w_pivoted_inv, b_w_pivoted);
		incrementvect = forwardSubstitution(A_w_pivoted, b_w_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - incrementvect[i];
		}

		for (int i = 0; i < M; i++) // checks for NaNs/ blowing up params
		{
			if ((params_new[i] != params_new[i]) || std::isinf(params_new[i])) {
				params_old.push_back(0.0);
				return params_old;
			}
		}

		for (int j = 0; j < M; j++) // stops the iteration if the difference between each element in the latest param vector with it's previous counterpart are all less than the specificed tolerance
		{
			if (std::abs(params_new[j] - params_old[j]) <= tolerance) {
				checkcount += 1;
			}

		}

		if (checkcount == M) {
			tolcheck = false;
		}
		params_old = params_new; // applies increment vector

		itercount += 1;
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};

std::vector <double> regularGN(double(*f)(std::vector <double>, std::vector <double>), std::vector <double(*)(std::vector <double>, std::vector <double>)> parsvector, std::vector <double> y, std::vector< std::vector <double> > x, std::vector <double> guess, double tolerance, std::vector <double> w) //the case of >1 independent (x) variables in the function
{
	int M = (int) parsvector.size();
	int Nr = (int) y.size();
	bool tolcheck = true;
	std::vector <double> params_old = guess;
	std::vector <double> params_new(M, 0.0); //placeholder
	std::vector < std::vector <double> > J, A_w, result, A_w_pivoted_inv;
	std::vector <double> res, b_w, b_w_pivoted, A_w_pivoted_vec, incrementvect;
	std::vector <double> A_w_pivoted_inner(M, 0.0);
	std::vector < std::vector <double> > A_w_pivoted(M, A_w_pivoted_inner);
	int checkcount;
	int itercount = 0;

	// creating weight matrix W:
	std::vector <double> innerW(Nr, 0.0);
	std::vector <std::vector <double> > W(Nr, innerW);
	for (int k = 0; k < Nr; k++) {
		W[k][k] = w[k];
	}
	while (tolcheck) {
		checkcount = 0;
		J = jacobian(parsvector, x, params_old);
		res = residuals((*f), y, x, params_old);

		A_w = dot(transpose(J), dot(W, J));
		b_w = dot(transpose(J), dot(W, res));
		result = pivotSystem(A_w, b_w);
		A_w_pivoted_vec = result[0];
		b_w_pivoted = result[1];

		int d = 0; //turns the outputted pivoted A vector back into an A matrix
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_w_pivoted[i][j] = A_w_pivoted_vec[d];
				d += 1;
			}
		}

		//A_w_pivoted_inv = inverse(A_w_pivoted);
		//incrementvect = dot(A_w_pivoted_inv, b_w_pivoted);
		incrementvect = forwardSubstitution(A_w_pivoted, b_w_pivoted);

		for (int i = 0; i < M; i++) // iterates over each row of the GN matrix equation
		{
			params_new[i] = params_old[i] - incrementvect[i];
		}

		for (int i = 0; i < M; i++) // checks for NaNs/ blowing up params
		{
			if ((params_new[i] != params_new[i]) || std::isinf(params_new[i])) {
				params_old.push_back(0.0);
				return params_old;
			}
		}

		for (int j = 0; j < M; j++) // stops the iteration if the difference between each element in the latest param vector with it's previous counterpart are all less than the specificed tolerance
		{
			if (std::abs(params_new[j] - params_old[j]) <= tolerance) {
				checkcount += 1;
			}

		}

		if (checkcount == M) {
			tolcheck = false;
		}
		params_old = params_new; // applies increment vector

		itercount += 1;
		//std::cout << params_old[0] << "  " << params_old[1] << std::endl;
	}
	//std::cout << "Modified Gauss-Newton completed in " << itercount << " iterations.\n";

	return params_new;

};

//MISC
double factorial(double n)
{
	if (n == 0)
		return 1;
	return n * factorial(n - 1);
}

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

double pi = 3.1415926535897932384626434;

double gaussian(double x, double mu, double sig) {
	return (std::exp((-0.5) * (std::pow((x - mu), 2.0) / (2.0 * std::pow(sig, 2.0)))) / std::sqrt(2.0*pi*std::pow(sig, 2.0)));
}

double getAvg(std::vector<double> x, std::vector <double> w, double(*f)(double, std::vector <double>), std::vector<double> params) {
	double uppersum = 0.0;
	double lowersum = 0.0;

	for (int i = 0; i < x.size(); i++) {
		uppersum += (x[i] * w[i]);
		lowersum += w[i];
	}

	return (uppersum / lowersum);
}

double getAvg_Exp(std::vector<double> x, std::vector <double> w, double(*f)(double, std::vector <double>), std::vector<double> params) {
	double uppersum = 0.0;
	double lowersum = 0.0;

	for (int i = 0; i < x.size(); i++) {
		double y_i = f(x[i], params);
		uppersum += w[i]*x[i]*std::pow(y_i, 2.0);
		lowersum += w[i]*std::pow(10.0, -2.0*y_i);
	}

	return (uppersum / lowersum);
}

double getLogXBar_PowerLaw(std::vector<double> x, std::vector <double> w, double(*f)(double, std::vector <double>), std::vector<double> params) {
	double uppersum = 0.0;
	double lowersum = 0.0;
	for (int i = 0; i < x.size(); i++) {
		double y_i = f(x[i], params);
		uppersum += w[i] * std::log10(x[i]) * std::pow(y_i, 2.0);
		lowersum += w[i] * std::pow(10.0, -2.0*y_i);
	}

	return (uppersum / lowersum);
}

double getLogXBar_Log(std::vector<double> x, std::vector <double> w, double(*f)(double, std::vector <double>), std::vector<double> params) {
	double uppersum = 0.0;
	double lowersum = 0.0;

	for (int i = 0; i < x.size(); i++) {
		double y_i = f(x[i], params);
		uppersum += w[i] * std::log10(x[i]) * std::pow(10.0, -2.0*y_i);
		lowersum += w[i] * std::pow(10.0, -2.0*y_i);
	}

	return (uppersum / lowersum);
}
