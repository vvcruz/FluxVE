#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>


int run_all_fft(double *bcoefre,double *bcoefim, double *X, double *Y, double *xknot, double *yknot, double *tknot,int kx,double ti,double tf,double steptspl,int nxr,int nyr,int nf,int NTG, int *fpr, double width, double *WPERe,double *WPEIm, double tmax,double stren,double t0);

int gendata_complex(fftw_complex *work, double *bcoefre, double *bcoefim, double *xknot, double *yknot, double *tknot, int kx, int nxr, int nyr, int nf, double x, double y, double ti,double tf, double steptspl, int NTG, double width, double tmax,double stren,double t0);

int run_fftw(fftw_complex *in, fftw_complex *out, int N);

int center_fft(fftw_complex *out,int N);

int gendata_function(fftw_complex *work,double tf,double steptspl,int NTG,double width);

int gendata(double *work, double *bcoef, double *xknot, double *yknot, double *tknot, int nxr, int nyr, int nf, double x, double y, double tf, double steptspl, int NTG, double width, int type);

double gausswp(double x, double x0, double hwhm, double norm);

int print_bcoef(double *bcoefre, double *bcoefim, int nx,int ny,int nz,char *filenam, char *type);

double gausswp_k(double k, double k0, double hwhm);

