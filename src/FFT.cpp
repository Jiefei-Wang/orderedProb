//The code is from :https://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B with minor changes

#include <Rcpp.h>
#include "kiss_fftr.h"
#define complex_mul_real(x_real,x_img,y_real,y_img) (x_real*y_real - x_img*y_img)
#define complex_mul_img(x_real,x_img,y_real,y_img) (x_real*y_img + x_img*y_real)
// [[Rcpp::export]]
void testfft(Rcpp::NumericVector x) {
	int nfft = XLENGTH(x);
	size_t required_size;
	kiss_fftr_alloc(nfft, 0, NULL, &required_size);
	Rprintf("%llu\n", required_size);
	void* buffer = malloc(required_size);
	kiss_fftr_cfg fwd = kiss_fftr_alloc(nfft, 0, buffer, &required_size);

	kiss_fft_cpx* res = new kiss_fft_cpx[nfft/2+1];
	kiss_fftr(fwd, (kiss_fft_scalar*)DATAPTR(x), res);
	for (int i = 0; i < nfft / 2 + 1; i++) {
		Rprintf("%f,%f\n", res[i].r, res[i].i);
	}
	free(buffer);
	delete[] res;
}



// [[Rcpp::export]]
void test() {
	const int nfft = 8;
	kiss_fft_cfg fwd = kiss_fft_alloc(nfft, 0, NULL, NULL);
	kiss_fft_cfg inv = kiss_fft_alloc(nfft, 1, NULL, NULL);

	kiss_fft_cpx* x=new kiss_fft_cpx[nfft];
	kiss_fft_cpx* fx = new kiss_fft_cpx[nfft];
	for (int i = 0; i < nfft; i++) {
		if (i < 4) {
			x[i].r = 1;
		}
		else {
			x[i].r = 0;
		}
		x[i].i = 0;
	}
	kiss_fft(fwd, (kiss_fft_cpx*)& x[0], (kiss_fft_cpx*)& fx[0]);
	Rprintf("before:\n");
	for (int i = 0; i < nfft; i++) {
		Rprintf("%f,%f\n", x[i].r, x[i].i);
	}
	Rprintf("After:\n");
	for (int i = 0; i < nfft; i++) {
		Rprintf("%f,%f\n", fx[i].r, fx[i].i);
	}
	kiss_fft_free(fwd);
	kiss_fft_free(inv);
	delete[] x;
	delete[] fx;
}
void convolution(void* buffer, R_xlen_t N,
	kiss_fft_scalar* x_in, kiss_fft_cpx* x_out,
	kiss_fft_scalar* y_in, kiss_fft_cpx* y_out) {

	size_t required_size;
	kiss_fftr_alloc(N, 0, NULL, &required_size);
	kiss_fftr_cfg fwd = kiss_fftr_alloc(N, 0, buffer, &required_size);
	kiss_fftr(fwd, x_in, x_out);
	kiss_fftr(fwd, y_in, y_out);
	for (R_xlen_t i = 0; i < N / 2 + 1; ++i) {
		double tmp = complex_mul_real(x_out[i].r, x_out[i].i, y_out[i].r,  y_out[i].i);
		x_out[i].i = complex_mul_img(x_out[i].r,  x_out[i].i, y_out[i].r,  y_out[i].i);
		x_out[i].r = tmp;
	}
	kiss_fftr_cfg inv = kiss_fftr_alloc(N, 1, buffer, &required_size);
	kiss_fftri(inv, x_out, x_in);
}
// [[Rcpp::export]]
void testConv(Rcpp::NumericVector x, Rcpp::NumericVector y) {
	int nfft = XLENGTH(x);
	size_t required_size;
	kiss_fftr_alloc(nfft, 0, NULL, &required_size);
	void* buffer = malloc(required_size);
	size_t len = nfft / 2 + 1;
	kiss_fft_cpx* x_out = new kiss_fft_cpx[len];
	kiss_fft_cpx* y_out = new kiss_fft_cpx[len];

	convolution(buffer,nfft,
		(kiss_fft_scalar*)DATAPTR(x), x_out, 
		(kiss_fft_scalar*)DATAPTR(y), y_out);
	free(buffer);
	delete[] x_out;
	delete[] y_out;
}