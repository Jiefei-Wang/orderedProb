#include <Rcpp.h>
#include <vector>
#include <chrono>
using namespace Rcpp;

std::vector<double> myFactorial;

double getPoisson(double k, double rate) {
	double result = k * std::log(rate) - rate - myFactorial[k];
	return std::exp(result);
}

void computeFactorialUpTo(R_xlen_t n){
	R_xlen_t oldN = myFactorial.size();
	myFactorial.reserve(n+1);
	for (R_xlen_t i = oldN; i < n+1; ++i) {
		if (i == 0 || i == 1) {
			myFactorial[i] = 0;
		}
		else {
			myFactorial[i]= myFactorial[i-1]+std::log(i);
		}
	}
}


#define SETVALUE(x,n,value)\
for(R_xlen_t I = 0;I<n;++I){\
(x)[I]=value;\
}

template<class T>
SEXP allocWithInit(int type, R_xlen_t len, T value) {
	SEXP R_res = Rf_allocVector(type, len);
	T* res = (double*)DATAPTR(R_res);
	for (R_xlen_t i = 0; i < len; ++i) {
		res[i] = value;
	}
	return R_res;
}


// [[Rcpp::export]]
double compute_prob(R_xlen_t m, SEXP R_g_value,SEXP R_h_value,
	R_xlen_t n_t, SEXP R_diff_t) {
	computeFactorialUpTo(m);
	double* Q = new double[m + 1];
	SETVALUE(Q, m + 1, 0);
	Q[0] = 1;
	double* newQ = new double[m + 1];
	double* g_value = (double*)DATAPTR(R_g_value);
	double* h_value = (double*)DATAPTR(R_h_value);
	double* diff_t = (double*)DATAPTR(R_diff_t);
	for (R_xlen_t i = 0; i < n_t - 1; ++i) {
		double gt_i = g_value[i];
		double gt_i_plus = g_value[i + 1];
		double ht_i_plus = h_value[i + 1];
		R_xlen_t maxRange = ht_i_plus - gt_i_plus - 1;
		SETVALUE(newQ, maxRange, 0);
		for (R_xlen_t j = 0; j < maxRange; ++j) {
			R_xlen_t curM = j + gt_i_plus + 1;
			for (R_xlen_t l = gt_i + 1; l <= curM; ++l) {
				double curQ = Q[l];
				double curPi = getPoisson(curM - l, m * diff_t[i]);
				newQ[j] += curQ * curPi;
			}
		}
		for (R_xlen_t j = gt_i_plus + 1; j <= ht_i_plus - 1; ++j) {
			Q[j] = newQ[j - (R_xlen_t)gt_i_plus -1];
		}

	}
	double result = Q[m] / getPoisson(m, m);
	delete[] Q;
	delete[] newQ;
	return result;
}

R_xlen_t findMaxMemSize(R_xlen_t n, double* g_value, double* h_value) {
	R_xlen_t maxSize = 0;
	for (R_xlen_t i = 0; i < n-1; ++i) {
		R_xlen_t curSize = h_value[i + 1] - g_value[i] - 1;
		maxSize = maxSize > curSize ? maxSize : curSize;
	}
	maxSize *= 2;
	double logSize = std::log2((double)maxSize);
	logSize = std::ceil(logSize);
	maxSize = std::pow(2, logSize);
	return maxSize;
}
R_xlen_t findCurMemSize(double g, double h) {
	R_xlen_t curSize = h - g - 1;
	curSize *= 2;
	double logSize = std::log2((double)curSize);
	logSize = std::ceil(logSize);
	curSize = std::pow(2, logSize);
	return curSize;
}
void convolution(R_xlen_t N,
	double* x_real, double* x_img,
	double* y_real, double* y_img);
// [[Rcpp::export]]
double compute_prob_fft(R_xlen_t m, SEXP R_g_value, SEXP R_h_value,
	R_xlen_t n_t, SEXP R_diff_t) {
	computeFactorialUpTo(m);
	double* Q = new double[m + 1];
	SETVALUE(Q, m + 1, 0);
	Q[0] = 1;
	double* g_value = (double*)DATAPTR(R_g_value);
	double* h_value = (double*)DATAPTR(R_h_value);
	double* diff_t = (double*)DATAPTR(R_diff_t);

	R_xlen_t max_size = findMaxMemSize(n_t, g_value, h_value);
	double* x_real = new double[max_size];
	double* y_real = new double[max_size];
	double* x_img = new double[max_size];
	double* y_img = new double[max_size];
	auto start = std::chrono::high_resolution_clock::now();
	for (R_xlen_t i = 0; i < n_t - 1; ++i) {
		double gt_i = g_value[i];
		double ht_i_plus = h_value[i + 1];
		R_xlen_t maxRange = ht_i_plus - gt_i - 1;
		R_xlen_t curSize = findCurMemSize(gt_i, ht_i_plus);
		//Rprintf("Curlen:%d,MaxMem:%d\n", maxRange, curSize);
		SETVALUE(x_real + maxRange, curSize - maxRange, 0);
		SETVALUE(y_real + maxRange, curSize - maxRange, 0);
		SETVALUE(x_img, curSize, 0);
		SETVALUE(y_img, curSize, 0);
		
		for (R_xlen_t j = 0; j < maxRange; j++) {
			x_real[j] = Q[j + (R_xlen_t)gt_i + 1];
			y_real[j] = getPoisson(j, m * diff_t[i]);
		}
		convolution(curSize, x_real, x_img, y_real, y_img);
		for (R_xlen_t j = 0; j < maxRange; j++) {
			Q[j + (R_xlen_t)gt_i + 1] = x_real[j];
			//Rprintf("%f,", x_real[j]);
		}
		//Rprintf("\n");
	}
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	Rprintf("%f\n", elapsed.count());
	double result = Q[m] / getPoisson(m, m);
	delete[] Q;
	delete[] x_real;
	delete[] x_img;
	delete[] y_real;
	delete[] y_img;
	return result;
}
