#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>

#include "splinesurf.h"
#define one 1.0E+0

/*
  gendata:
  this routine generates a full data set from its [0,tf] spline representation, then mirrors it back to [-tf,0].
  it is also applied a gaussian window to the data set, according to the width input value.
  the main purpose is to prepare the data for fourier transform.
  the type variable should be 1 for even functions f(-x)=f(x) and 2 for odd functions f(-x)=-f(x).
 */

int gendata(double *work, double *bcoef, double *xknot, double *yknot, double *tknot, int nxr, int nyr, int nf, double x, double y, double tf, double steptspl, int NTG, double width, int type){
  int i,j,nt,korder;
  double t,val,window,taux;
  korder = 3.0;
  
  taux = pow(tf,2)/(log(one/width)/log(M_E));

  for(i=0;i<NTG;i++){
    t = -tf + i*steptspl;
    window = exp(-pow(t,2)/taux);

    if(type==1 || (type==2 && t>0.00e+0) ){
      t = fabs(t);
      work[i] = window*dbs3vl_ (&y,&x,&t,&korder,&korder,&korder,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoef);
    }else if(type==2 && t<0.00e+0){
      t = fabs(t);
      work[i] = -window*dbs3vl_ (&y,&x,&t,&korder,&korder,&korder,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoef);
    }

  }

}


/*
  gendata_complex:

  this routine generates a full complex data set from its [0,tf] spline representation, then mirrors it back to [-tf,0].
  it is also applied a supergaussian window to the data set, according to the width input value.
  the main purpose is to prepare the data for fourier transform.
 */

int gendata_complex(fftw_complex *work, double *bcoefre, double *bcoefim, double *xknot, double *yknot, double *tknot, int kx, int nxr, int nyr, int nf, double x, double y, double ti,double tf, double steptspl, int NTG, double width, double tmax,double stren,double t0){
  int i,j,nt,m;
  double t,val,window,taux;
  //FILE *wind=fopen("window.dat","w");

  taux = pow(tf,2)/(log(one/width)/log(M_E));
  m=6;
  taux = 12000;
  tmax = 1500;
  

  for(i=0;i<NTG;i++){
    //t = -tf + i*steptspl;
    t = ti + i*steptspl;
    //window = exp(-pow((t - tmax),2)/taux);
    //window = exp(-pow((t - tmax)/taux,2));
    //window = 1.0/(1.0 + exp(-2.0*0.01*(t - tmax - 100)));
    //window = exp(-pow((t - tmax)/taux,2*m));
    window = 1.000 - 1.000/(1 + exp(-2.0*stren*(t - t0)));
    //fprintf(wind,"%E %E\n",t,window);
    //window = 1.0000e+0;
    work[i][0] = window*dbs3vl_ (&y,&x,&t,&kx,&kx,&kx,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoefre);
    work[i][1] = window*dbs3vl_ (&y,&x,&t,&kx,&kx,&kx,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoefim);
  }
}

/*
 shifts zero frequency to the center of the array
*/

int center_fft(fftw_complex *out,int N){
  int i;
  double work;

  //centering the fourier transform
  for(i=0;i<N/2;i++){
    work= out[N/2 +i][0];
    out[N/2 +i][0] = out[i][0];
    out[i][0] = work;

    work= out[N/2 +i][1];
    out[N/2 +i][1] = out[i][1];
    out[i][1] = work;  
  }
}



/*
  run_all_fft

  This routine converts a wavepacket propagation from |f(x,y,t)> to |f(x,y,E)> , from a previously calculated cubic spline representation, using the fftw library routines

  main input: bcoefre, bcoefim. real and complex part spline matrices of |f(x,y,t)>
  main output: WPERe, WPEIm. real and complex part vectors of |f(x,y,E)>

 */
int run_all_fft(double *bcoefre,double *bcoefim, double *X, double *Y,double *xknot, double *yknot, double *tknot, int kx,double ti,double tf,double steptspl,int nxr,int nyr,int nf,int NTG, int *fpr, double width,double *WPERe,double *WPEIm, double tmax,double stren,double t0){
  int i,j,l,ll,k,nE;
  double x,y,E,ke,norm,stepxspl;
  fftw_complex *workin,*workout;
  fftw_plan p;
  double Ei,stepE;
  FILE *filin=fopen("fft_check_in.dat","w");
  FILE *filout=fopen("fft_check_out.dat","w");

  stepE = 2*M_PI/(NTG*steptspl);
  Ei = -NTG*stepE/(2.0E+0);
  nE =  NTG; //NTG/2
  printf("\nEnergy step: %E \n\n",stepE);

  //--- allocating work vectors
  workin = fftw_malloc(sizeof(fftw_complex) * NTG);
  workout = fftw_malloc(sizeof(fftw_complex) * NTG);


  for(i=0;i<nxr;i++){
    for(j=0;j<nyr;j++){
      //generates data set from spline coeff.
      gendata_complex(workin,bcoefre,bcoefim,xknot,yknot,tknot,kx,nxr,nyr,nf,X[i],Y[j],ti,tf,steptspl,NTG,width,tmax,stren,t0);

      //centers t=0 at the zero frequency position of the fft input vector(first) due to periodicity requirement of fft
      //center_fft(workin,NTG);

      //--- do forward fft
      // FFTW_BACKWARD -> sign in the exponent = 1.
      p = fftw_plan_dft_1d(NTG,workin,workout,FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(p);  

      //shifts the fft result back to the center of the array
      center_fft(workout,NTG);

      //debug ---
      if(i==5 && j==26){
	//center_fft(workin,NTG);
	fprintf(filin,"# %E %E\n",X[i],Y[j]);
	//for(k=0;k<NTG;k++)fprintf(filin,"%E %E %E\n",-tf + k*steptspl,workin[k][0],workin[k][1]);
	for(k=0;k<NTG;k++)fprintf(filin,"%E %E %E\n",k*steptspl,workin[k][0],workin[k][1]);
	fprintf(filout,"# %E %E\n",X[i],Y[j]);
	//for(k=0;k<NTG;k++)fprintf(filout,"%E %E %E\n",Ei + k*stepE,steptspl*workout[k][0],steptspl*workout[k][1]);
	for(k=0;k<NTG;k++)fprintf(filout,"%E %E %E\n",Ei + k*stepE,(1.0/sqrt(2*M_PI))*steptspl*workout[k][0],(1.0/sqrt(2*M_PI))*steptspl*workout[k][1]);
      }
      //----------
      
      //copies fft results to permanent array
      for(l=0;l<NTG;l++){ 
	//result of discrete FT must be multiplied by the time step (see num. recipes)
	ll = j + i*nyr + l*nxr*nyr;
	//printf("%d %d %d -> %d \n",j,i,k,ll);
	//fprintf(fil," %E %E %E\n",Ei + l*stepE,steptspl*workout[l][0],steptspl*workout[l][1]);
	WPERe[ll] = (1.0/sqrt(2*M_PI))*steptspl*workout[l][0];
	WPEIm[ll] = (1.0/sqrt(2*M_PI))*steptspl*workout[l][1];
      }
      //fprintf(fil,"\n");
    }
  }
  //----------


  fftw_destroy_plan(p);
  fftw_free(workin); 
  fftw_free(workout);
  return 0;
}


/*
 gausswp

 this routine generates the initial state spatial distribution

 */

double gausswp(double x, double x0, double hwhm, double norm){
  double gauss,X,ln_2;
  X = x - x0;
  ln_2 = log(2);
  gauss = norm*exp(-ln_2*pow(X/hwhm,2));
  return gauss;
}

/*
 gausswp_k

 this routine generates the initial state momentum distribution

 */

double gausswp_k(double k, double k0, double hwhm){
  double gauss,norm,K,ln2;
  K = k - k0;
  // M_LN2 = 6.931471805599453E-01
  ln2 = 0.69314718e+0; //value used in eSPec
  //norm = sqrt(sqrt( hwhm*hwhm/(2.0*M_PI*M_LN2)));
  //gauss = norm*exp(-pow(hwhm*K,2)/(4.0*M_LN2));
  norm = sqrt(sqrt( hwhm*hwhm/(2.0*M_PI*ln2)));
  gauss = norm*exp(-pow(hwhm*K,2)/(4.0*ln2));
  return gauss;
}


/*
  print_bcoef

  this routine prints the spline coefficients into a file for a later calculation.

 */

int print_bcoef(double *bcoefre, double *bcoefim, int nx,int ny,int nz,char *filenam ,char *type){
  int i,nt;
  FILE *out=fopen(filenam,"w");
  nt = nx*ny*nz;

  fprintf(out,"%s\n",type);
  fprintf(out,"%d %d %d\n",ny,nx,nz);

  for(i=0;i<nt;i++){
    fprintf(out,"%E.15 %E.15 ",bcoefre[i],bcoefim[i]);
  }

  fclose(out);

  return 0;

}



/*
  gendata_function
  routine to test the fourier transform
*/
int gendata_function(fftw_complex *work,double tf,double steptspl,int NTG,double width){
  int i;
  double t,t0,window,taux,w0;
  //FILE *wind=fopen("window.dat","w");
  t0 = 0.0*tf/2;//0.0e+0;//tf/2;
  w0 = 5.0e+0;
  width = 1.0e-2;
  taux = pow(tf/(log(one/width)/log(M_E)),2);
  for(i=0;i<NTG;i++){
    t = -tf + i*steptspl;
    window = exp(-pow(t,2)/taux);
    //fprintf(wind,"%E %E\n",t,window);
    //window = 1.0e+0;
    //work[i][0]=window*(cos(2*t)+cos(t)+5.0*cos(5*t));
    work[i][0]=exp(-M_PI*pow(t-t0,2))*cos(w0*t);
    work[i][1]=exp(-M_PI*pow(t-t0,2))*sin(w0*t); //0.0;//window*sin(t);//0.0E+0;
  }
  return 0;
}


/*

OLD BACKUP


/*
  gendata_complex:

  this routine generates a full complex data set from its [0,tf] spline representation, then mirrors it back to [-tf,0].
  it is also applied a supergaussian window to the data set, according to the width input value.
  the main purpose is to prepare the data for fourier transform.

int gendata_complex(fftw_complex *work, double *bcoefre, double *bcoefim, double *xknot, double *yknot, double *tknot, int nxr, int nyr, int nf, double x, double y, double tf, double steptspl, int NTG, double width){
  int i,j,nt,korder;
  double t,val,window,taux;
  korder = 3.0;
  
  taux = pow(tf,2)/(log(one/width)/log(M_E));

  for(i=0;i<NTG;i++){
    t = -tf + i*steptspl;
    window = exp(-pow(t,2)/taux);
    window = 1.0000;
    if(t>=0.00e+0){
      t = fabs(t);
      work[i][0] = window*dbs3vl_ (&y,&x,&t,&korder,&korder,&korder,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoefre);
      work[i][1] = window*dbs3vl_ (&y,&x,&t,&korder,&korder,&korder,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoefim);
    }else{
      t = fabs(t);
      work[i][0] =  window*dbs3vl_ (&y,&x,&t,&korder,&korder,&korder,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoefre);
      work[i][1] = -window*dbs3vl_ (&y,&x,&t,&korder,&korder,&korder,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoefim);
    }
  }
}






 */
