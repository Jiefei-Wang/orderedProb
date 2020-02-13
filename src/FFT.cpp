//The code is from :https://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B with minor changes

#include <complex>
#include <iostream>
#include <valarray>
#include <Rcpp.h>
#include "kiss_fft.h"
#include <vector>
using std::vector;

// [[Rcpp::export]]
SEXP test() {
	const int nfft = 256;
	kiss_fft_cfg fwd = kiss_fft_alloc(nfft, 0, NULL, NULL);
	kiss_fft_cfg inv = kiss_fft_alloc(nfft, 1, NULL, NULL);

	vector<std::complex<float>> x(nfft, 0.0);
	vector<std::complex<float>> fx(nfft, 0.0);

	x[0] = 1;
	x[1] = std::complex<float>(0, 3);

	kiss_fft(fwd, (kiss_fft_cpx*)& x[0], (kiss_fft_cpx*)& fx[0]);
	for (int k = 0; k < nfft; ++k) {
		fx[k] = fx[k] * conj(fx[k]);
		fx[k] *= 1. / nfft;
	}
	kiss_fft(inv, (kiss_fft_cpx*)& fx[0], (kiss_fft_cpx*)& x[0]);
	std::cout
		<< x[0] << ","
		<< x[1] << ","
		<< x[2] << ","
		<< x[3] << " ... " << std::endl;
	kiss_fft_free(fwd);
	kiss_fft_free(inv);
	return 0;
}


#define complex_mul_real(x_real,x_img,y_real,y_img) (x_real*y_real - x_img*y_img)
#define complex_mul_img(x_real,x_img,y_real,y_img) (x_real*y_img + x_img*y_real)

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;


void fft(R_xlen_t N,double* x_real, double* x_img)
{
	// DFT
	R_xlen_t k = N, n;
	double thetaT = PI / N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				//Complex t = x[a] - x[b];
				Complex t = Complex(x_real[a] - x_real[b], x_img[a] - x_img[b]);
				//x[a] += x[b];
				x_real[a] += x_real[b];
				x_img[a] += x_img[b];
				//x[b] = t * T;
				Complex tmp = t * T;
				x_real[b] = tmp.real();
				x_img[b] = tmp.imag();
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			//Complex t = x[a];
			//x[a] = x[b];
			//x[b] = t;
			std::swap(x_real[a], x_real[b]);
			std::swap(x_img[a], x_img[b]);
		}
	}
}

void fft_negative_img(R_xlen_t N,double* x_real, double* x_img)
{
	// DFT
	R_xlen_t k = N, n;
	double thetaT = PI / N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				//Complex t = x[a] - x[b];
				Complex t = Complex(x_real[a] - x_real[b], -x_img[a] + x_img[b]);
				//x[a] += x[b];
				x_real[a] += x_real[b];
				x_img[a] += x_img[b];
				//x[b] = t * T;
				Complex tmp = t * T;
				x_real[b] = tmp.real();
				x_img[b] = -tmp.imag();
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			//Complex t = x[a];
			//x[a] = x[b];
			//x[b] = t;
			auto tmp = x_real[a];
			x_real[a] = x_real[b];
			x_real[b] = tmp;

			tmp = x_img[a];
			x_img[a] = x_img[b];
			x_img[b] = tmp;

			//std::swap(x_real[a], x_real[b]);
			//std::swap(x_img[a], x_img[b]);
		}
	}
}

// inverse fft (in-place)
void ifft(R_xlen_t N, double* x_real, double* x_img)
{
	for (R_xlen_t i = 0; i < N; ++i) {
		x_img[i] = -x_img[i];
	}
    // conjugate the complex numbers
    //x = x.apply(std::conj);
    

    // forward fft
	fft(N,x_real, x_img);
    
    // conjugate the complex numbers again
    //x = x.apply(std::conj);
	for (R_xlen_t i = 0; i < N; ++i) {
		x_img[i] = -x_img[i]/N;
		x_real[i] /= N;
	}
    // scale the numbers
    //x /= x.size();
}

// inverse fft (in-place)
void ifft_no_img(R_xlen_t N, double* x_real, double* x_img)
{
	// forward fft
	fft_negative_img(N,x_real, x_img);

	// conjugate the complex numbers again
	//x = x.apply(std::conj);
	for (R_xlen_t i = 0; i < N; ++i) {
		x_real[i] /= N;
	}
	// scale the numbers
	//x /= x.size();
}

//Results will be in x
void convolution(R_xlen_t N, 
	double* x_real, double* x_img,
	double* y_real, double* y_img) {
	fft(N, x_real, x_img);
	fft(N, y_real, y_img);

	for (R_xlen_t i = 0; i < N; ++i) {
		double tmp = complex_mul_real(x_real[i], x_img[i], y_real[i], y_img[i]);
		x_img[i] = complex_mul_img(x_real[i], x_img[i], y_real[i], y_img[i]);
		x_real[i] = tmp;
	}

	ifft_no_img(N, x_real, x_img);
}

// [[Rcpp::export]]
void test(SEXP x_real, SEXP x_img, SEXP y_real, SEXP y_img) {
	R_xlen_t N = Rf_xlength(x_real);
	convolution(N, (double*)DATAPTR(x_real), (double*)DATAPTR(x_img),
		(double*)DATAPTR(y_real), (double*)DATAPTR(y_img));
}



// [[Rcpp::export]]
SEXP testfft(SEXP x_real) {
	R_xlen_t N = Rf_xlength(x_real);
	SEXP res = Rf_allocVector(REALSXP, N);
	
	fft(N,(double*)DATAPTR(x_real), (double*)DATAPTR(res));

	return res;
}



// [[Rcpp::export]]
int main1(int N)
{
    double test_real[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0, 0.0 };
	double test_img[] = { 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0 };
    
    // forward fft
    fft(N,test_real, test_img);
    
    std::cout << "fft" << std::endl;
    for (int i = 0; i < N; ++i)
    {
        std::cout << test_real[i]<<","<< test_img [i]<< std::endl;
    }
    
    // inverse fft
	ifft(N,test_real, test_img);
    
    std::cout << std::endl << "ifft" << std::endl;
    for (int i = 0; i < N; ++i)
    {
        std::cout << test_real[i] << "," << test_img[i] << std::endl;
    }
    return 0;
}