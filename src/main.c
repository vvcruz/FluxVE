#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#include<time.h>
#include<string.h>

#include "splinesurf.h"
#include "jacobi.h"
#include "Jflux.h"
#include "readf.h"
#include "morse_vib.h"
#include "fourier.h"

/*
Program: FluxVE

Author: Vinicius Vaz da Cruz - viniciusvcruz@gmail.com

History: Based on FluxV written by Vinicius

Purpose: Compute energy resolved reaction probabilities from 2D collision wavepacket simulations

Method: Perform FFT from |f(R,r,t)> to |f(R,r,E)>, then apply flux operator.
        see ref: Neuhauser et al. comp. phys. comm. 63(1991)460-481
	         Meijer et al. Chem. Phys. Lett. 293(1998) 270-276

Goiania, 22nd of september of 2014
 */

#define NAV   6.0221415E+23  //number of entities in 1 mol
#define Me    9.10938188E-31 //mass of electron in Kg
#define FSAU(x) (x)*41.3413745758e+0 ///0.024189E+0 //fs to au units of time conversion                        mtrxdiag_ (char *dim, int *il, int *iu, int *info, int *mxdct, int *n, double *abstol, int *iwork, int *np, double *eigvl, double *shm, double *vpot, double *work, double *wk, double *eigvc);

int rdinput(FILE *arq,int *nx, int *ny, double *xi, double *xf, double *yi, double *yf,char *file,int *stfil,int *endfil, double *ti, double *tf, double *m1, double *m2, double *fluxline, char *axis, double *nv, int *nf, double *stepxpls, double *stepyspl, int *twopow, double *a, double *x0, double *k0, double *Evib, double *potshift, int *prtRMSwhole, int *type, int *noreg, double *ma, double *mb, double *mc, int *wholegrid,double *stren,double *t0);

void readallspl_(char *file,int *iflth,double *X,double *Y,double *T,double *bcoefRe,double *bcoefIm,double *xknot,double *yknot,double *tknot,int *nx, int *ny, int *nf, int *kx, int *nxr, int *nyr, int *fpr, double *norm, int *prtRMSwhole, int *noreg, int *wholegrid, int *stfil);

int main(){
  int i,j,k,l,m,nx,ny,nf,nxr,nyr,nfr,xs,ys,ks,ldf,mdf,nw,nz,nwr,nzr;
  int fp[10],fpr[10],stfil,endfil,iflth,nrmse,type,noreg;
  int wholegrid,twopow,prtwpE,prtRMSwhole;
  double fluxline[10],nv,mem,xrange[2],yrange[2],flux,rflux[2],Prob,norm;
  double ti,tf,xi,xf,yi,yf,stept,stepx,stepy,width,potshift,maxF[2],maxaE[2],stepw,stepz;
  double *X,*Y,*T,rmse,xwork,*W,*Z;
  double m1,m2,x,y,pk,pki,pkf,steppk,ansk,val[5];
  char axis,file[30],wp_Enam[20],num[20];
  FILE *arq=fopen("fluxe.inp","r");
  FILE *deb=fopen("debug.dat","w");
  //--- spline variables
  int NXG,NYG,NTG,fpspl[4],ftinfo,kx;
  double *bcoefre,*bcoefim,*xknot,*yknot,*tknot;
  double stepxspl,stepyspl,steptspl,Xspl[3],Yspl[3],Respl[3],Imspl[3];
  //-- initial wavepacket parameters
  double a,x0,k0;
  //--- fourier variables
  int nE;
  double *WPERe,*WPEIm,*E,Ei,aE,stepE,*eknot,Evib;
  double kmin,kmax,Emin,Emax,stren,t0;
  fftw_plan p;
  FILE *probf = fopen("probability.dat","w");
  FILE *tdprobf = fopen("td-flux.dat","w");
  FILE *wp_E;
  //------ Jacobi variables
  double stepR,stepr,*wzbcoefre,*wzbcoefim,wrange[2],zrange[2];
  double Ri,Rf,ri,rf,m1beta,m2beta,ma,mb,mc,R,r,wR,zr;
  //------------------
  clock_t begint,endt;
  double timediff;


  //--- Default values
  type=0;
  noreg=1;
  axis='x';
  nv=1.0E+0;
  wholegrid=1;
  stepxspl=1.0E-3;
  stepyspl=1.0E-3;
  NYG = 1024;
  t0=0.0e+0;
  stren = 1.5e-3;
  twopow=9; //11;
  width = 1.0e-9;
  kx=8.0e+0;
  prtwpE=1; // =0 to print transformed wp into files -> NEED TO ADD TO rdinput
  prtRMSwhole=1; // = 0 to compute RMSE in relation to whole data -> NEED TO ADD TO rdinput

  //-------------------------------------------------//
  printf("\n ============================================\n");
  printf("||     -> -> -> -> -> -> -> -> -> -> -> ->   ||\n");
  printf("||     -> -> -> -> -> -> -> -> -> -> -> ->   ||\n");
  printf("||                                           ||\n");
  printf("||                  FluxVE                   ||\n");
  printf("||                                           ||\n");
  printf("||    -> -> -> -> -> -> -> -> -> -> -> ->    ||\n");
  printf("||    -> -> -> -> -> -> -> -> -> -> -> ->    ||\n");
  printf(" ============================================\n");
  printf("Vinicius V. Cruz and Freddy F. Guimaraes\n\n");
  printf("Reading input parameters...\n\n");

  //--- Read Input File
  rdinput(arq,&nx,&ny,&xi,&xf,&yi,&yf,file,&stfil,&endfil,&ti, &tf,&m1, &m2,fluxline,&axis,&nv,&nf,&stepxspl,&stepyspl,&twopow,&a,&x0,&k0,&Evib,&potshift,&prtRMSwhole,&type,&noreg,&ma,&mb,&mc,&wholegrid,&stren,&t0);
  fclose(arq);
  //-----


  //--- Units Conversion
  //m1=(m1*1.0E-3)/(NAV * Me);
  //m2=(m2*1.0E-3)/(NAV * Me);
  //conversion below used in eSPec

  //heaviside window center
  printf("LOL 2 %lf %lf \n",stren,t0);
  if(t0 == 0.0e+0) t0 = tf - 0.15 * (tf - ti); 
  t0 = FSAU(t0);

  stepx=(xf -xi)/(nx -1.0);
  stepy=(yf -yi)/(ny -1.0);
  stept=(tf - ti)/(nf -1.0);
  m1 = 1.836152667539e+3*m1;
  m2 = 1.836152667539e+3*m2;
  if(type == 0){
    fp[0]= round((fluxline[0] - xi)/stepx) + 2;
    fp[1]= round((fluxline[1] - xi)/stepx) + 2; //-2? if flux was along y we might neef to change this  
    fp[2]= round((fluxline[2] - yi)/stepy) + 2;
    fp[3]= round((fluxline[3] - yi)/stepy) - 2; 
  }else if(type == 1){
    ma = 1.836152667539e+3*ma;
    mb = 1.836152667539e+3*mb;
    mc = 1.836152667539e+3*mc;
    m1beta = ma*(mb + mc)/(ma + mb + mc);
    m2beta = mb*mc/(mb + mc);
    Ri = fluxline[0];
    Ri = fluxline[1];
    ri = fluxline[2];
    rf = fluxline[3];
    stepr = (rf - ri)/(NYG - 1.000);
    stepR = 1.0e-7;
  }

  //NTG = pow(2,twopow+1);
  NTG = pow(2,twopow);
  //steptspl=(2.0*tf)/(NTG); //-1.0e+0);
  steptspl=(tf - ti)/(NTG);
  //--- spline parameters determination 

  //--- checks if fluxlines are within the x and y range
  if(type==0){
    
    if( xi+fp[0]*stepx > xf ||  xi+fp[0]*stepx < xi){
      printf("error! fluxline x %lf must be within the x range \n", xi+fp[0]*stepx);
      return 66;
    } else if(xi+fp[1]*stepx > xf ||  xi+fp[1]*stepx < xi){
      printf("error! fluxline x %lf must be within the x range \n", xi+fp[1]*stepx);
      return 66;
    } else if(yi+fp[2]*stepy > yf ||  yi+fp[2]*stepy < yi){
      printf("error! fluxline x %lf must be within the x range \n", yi+fp[2]*stepy);
      return 66;
    } else if(yi+fp[3]*stepy > yf ||  yi+fp[3]*stepy < yi){
      printf("error! fluxline x %lf must be within the x range \n", yi+fp[3]*stepy);
      return 66;
    }
    
  }else if(type==1){
    //--- checks if fluxline is within the provided range
    
    Jacobitrans(&wR,&zr,ma,mb,mc,&Ri,&ri,-1);
    if(wR > xf || wR < xi || zr > yf || zr < yi){
      printf("Error! transformed fluxline not within the provided grid! \n");
      printf("R = %lf ri = %lf ->  x = %lf, y = %lf \n",Ri,ri,wR,zr);
      return 666;
    }
    Jacobitrans(&wR,&zr,ma,mb,mc,&Ri,&rf,-1);
    if(wR > xf || wR < xi || zr > yf || zr < yi){
      printf("Error! transformed fluxline not within the provided grid! \n");
      printf("R = %lf rf = %lf ->  x = %lf, y = %lf \n",Ri,rf,wR,zr);
      return 666;
    }

  }

  //--- Input Parameters -----------------------------------------------------------------------//
  printf("\n<< Starting input section >>\n\n");

  printf("> Run type:\n");
  if(type==0){
    printf("This is a non-reactive calculation\n");
  }else if(type==1){
    printf("This is a reactive flux calculation\n");
    printf("a Jacobi coordinate transformation will be performed \n\n");
     wholegrid = 0;
  }

  printf(">Grid parameters:\n");
  printf("npoints x %d \n", nx);
  printf("npoints y %d \n", ny);
  //--- these limits are set in the splinesurf.f file
  if(nx > 200 || ny > 420 || nf > 2050){
    printf("error! at least one of the number of points is larger than the limit\n");
    printf(" nx_max = 200; nx = %d \n",nx);
    printf(" ny_max = 420; nx = %d \n",ny);
    printf(" nf_max = 2050; nf = %d \n",nf);
    printf("\n you may adjust these limits in the file splinesurf.f, and recompile the program \n");
    return 666;
  }

  printf("xrange: %lf %lf\n", xi, xf);
  printf("yrange %lf %lf\n", yi, yf);
  printf("step X %.7lf Bohr\n",stepx);
  printf("step Y %.7lf Bohr\n", stepy);
  
  printf("\n>Wavepacket files:\n");
  printf("File name: %s\n",file);
  printf("ti %lf tf %lf step %lf fs nfile: %d \n",ti,tf,stept,nf);
  printf("masses m1 = %lf a.u.\n", m1);
  printf("       m2 = %lf a.u.\n", m2);
  
  //printf("\nMemory required for 3D spline: %E GB\n",mem);
  if(type==1){
    printf("\n>Jacobi Transformation\n");
    printf("masses ma = %lf a.u., mb = %lf a.u., mc = %lf a.u. \n", ma,mb,mc);
    printf("       m1beta = %lf a.u.\n", m1beta);
    printf("       m2beta = %lf a.u.\n", m2beta);
  }
  
  printf("\n>Numerical methods:\n");
  printf("Derivative method: centered difference O(h^4)\n");
  printf("Integration method: Simpson's 3/8 rule\n");
  printf("Fourier Tranform: FFTW\n");
  printf("Cubic spline interpolation will be used.\n");


  printf("\n>Flux parameters:\n");
  if(type == 0){
  printf("flux will be computed along the %c axis \n", axis);
  printf(" fluxline: x: %lf %lf \n", xi+fp[0]*stepx, xi+fp[1]*stepx);
  printf("           y: %lf %lf \n", yi+fp[2]*stepy, yi+fp[3]*stepy);
  } else if(type == 1){
    printf("\n Reactive fluxline (in transformed Jacobi coordinates) \n");
    Jacobitrans(&wR,&zr,ma,mb,mc,&Ri,&ri,-1);
    printf("         R = %lf ri = %lf ->  x = %lf, y = %lf \n",Ri,ri,wR,zr);
    Jacobitrans(&wR,&zr,ma,mb,mc,&Ri,&rf,-1);
    printf("         R = %lf rf = %lf ->  x = %lf, y = %lf \n\n",Ri,rf,wR,zr);    
  }
  // CHECK THIS ---------------------------- >>
  fluxline[0] = xi+fp[0]*stepx;
  fluxline[1] = xi+fp[1]*stepx;
  fluxline[2] = yi+fp[2]*stepy;
  fluxline[3] = yi+fp[3]*stepy;


  

  if(wholegrid == 1){
    printf("\n>Reduced spline grid:\n");

    if(fp[0]>=4 && xi + (fp[0] - 4)*stepx >= xi ) fpr[0] = fp[0] - 4; 
    else if(fp[0]>=3) fpr[0] = fp[0] - 3;
    else if(fp[0]>=2) fpr[0] = fp[0] - 2;
    else if(fp[0]>=1) fpr[0] = fp[0] - 1;
    else fpr[0] = fp[0];

    if(fp[1]+4<=nx  && xi + (fp[0] + 4)*stepx <= xf) fpr[1] = fp[1] + 4;
    else if(fp[1]+3<=nx) fpr[1] = fp[1] + 3;
    else if(fp[1]+2<=nx) fpr[1] = fp[1] + 2;
    else if(fp[1]+1<=nx) fpr[1] = fp[1] + 1;
    else fpr[1] = fp[1];

    if(fp[2]>=4) fpr[2] = fp[2] - 4; 
    else if(fp[2]>=3) fpr[2] = fp[2] - 3;
    else if(fp[2]>=2) fpr[2] = fp[2] - 2;
    else if(fp[2]>=1) fpr[2] = fp[2] - 1;
    else fpr[2] = fp[2];

    if(fp[3]+4<=ny)  fpr[3] = fp[3] + 4;
    else if(fp[3]+3<=ny)fpr[3] = fp[3] + 3;
    else if(fp[3]+2<=ny)fpr[3] = fp[3] + 2;
    else if(fp[3]+1<=ny)fpr[3] = fp[3] + 1;
    else fpr[3] = fp[3];

    xrange[0]=xi + fpr[0]*stepx;
    xrange[1]=xi + fpr[1]*stepx;
    yrange[0]=yi + fpr[2]*stepy;
    yrange[1]=yi + fpr[3]*stepy;

  
    nxr = fpr[1]-fpr[0]+1;
    nyr = fpr[3]-fpr[2]+1;
  }else if(wholegrid==0){
    fpr[0] = 0;
    fpr[1] = nx;
    fpr[2] = 0;
    fpr[3] = ny;
    xrange[0]=xi + fpr[0]*stepx;
    xrange[1]=xi + fpr[1]*stepx;
    yrange[0]=yi + fpr[2]*stepy;
    yrange[1]=yi + fpr[3]*stepy;

  
    nxr = nx;
    nyr = ny;
    printf("the whole grid provided will be used for spline \n");
    printf("%d -> %d , %d -> %d \n\n",nx,nxr,ny,nyr);
  }
  //------------------------------------------->>

  //NXG = (xrange[1]-xrange[0])/stepxspl;
  //NYG = 1 + (yrange[1]-yrange[0])/stepyspl;
  stepyspl = (yrange[1]-yrange[0])/(NYG - 1.000);

  printf("Grid upon which the spline will be applied:\n");
  printf("x: %lf %lf, %d points\n",xrange[0],xrange[1],nxr);
  printf("y: %lf %lf, %d -> %d points\n",yrange[0],yrange[1],nyr,NYG);
  if(kx > nxr){
    printf("spline order changed from %d to %d \n",kx,nxr);
    kx = nxr;
  }

  printf("\n>Initial collision wavepacket parameters:\n");
  printf("a = %lf a.u., x0 =  %lf a.u., k0 = %lf a.u.\n",a,x0,k0);
  printf("Evib = %E a.u., potential shift =  %E a.u.\n",Evib,potshift);
  //momentum and Energy interval
  kmin = k0 - 3.0*(sqrt(2*M_LN2)/a);
  kmax = k0 + 3.0*(sqrt(2*M_LN2)/a);
  if(kmin >= 0.0e+0){
    Emin = kmin*kmin/(2.0*m1) + Evib + potshift;
  }else if(kmin < 0){
    Emin = 0.0e+0 + Evib + potshift;
  }
  Emax = kmax*kmax/(2.0*m1) + Evib + potshift;
  printf("\nE0 = k0^2/2m: %.11E,\n",k0*k0/(2.0*m1));
  printf("Momentum Interval: k_min = %E, k_max = %E \n",kmin,kmax);
  printf("Energy Interval: E_min = %E, E_max = %E \n\n",Emin,Emax);

  //--------------- 

  printf("\n>Fourier transformparameters:\n");
  printf("grid is 2^%d, time step = %E \n",twopow,steptspl);
  stepE = 2*M_PI/(NTG*FSAU(steptspl));
  Ei = -NTG*stepE/(2.0E+0);
  nE = NTG;
  printf("Energy step = %E a.u., Ef = %E a.u. \n\n",stepE,-Ei);
  printf("FFT window: heaviside function \n");
  printf("k = %lf, t0 = %lf \n\n",stren,t0);


  mem = (4.0*(nxr*nyr*nf*sizeof(double))/(1024*1024*1024) + 2*(nxr+nyr+nf+NTG)*sizeof(double) + 2.0*(nxr*nyr*NTG*sizeof(double)))/(1024*1024*1024);
  printf("Memory requirement estimation: %E GB\n",mem);
  
  
  printf("\nFinished input section!\n");
  printf("------------------------------------------------------\n\n\n");
    
  // if(proj==0){
  //printf("%d Flux Projections will be computed\n", nst);
  //if(vib==0)  printf("Morse vibrational parameters will be computed \n");
  //}
  //--------------------------------------------------------------------------------------------//

  //--- Allocate arrays
  bcoefre = malloc(nyr*nxr*nf*sizeof(double));
  bcoefim = malloc(nyr*nxr*nf*sizeof(double));
  X = malloc(nxr*sizeof(double));
  Y = malloc(nyr*sizeof(double));
  T = malloc(nf*sizeof(double));
  xknot = malloc((nxr+kx)*sizeof(double));
  yknot = malloc((nyr+kx)*sizeof(double));
  tknot = malloc((nf+kx)*sizeof(double));
  WPERe = malloc(nxr*nyr*(NTG)*sizeof(double)); 
  WPEIm = malloc(nxr*nyr*(NTG)*sizeof(double)); 
  //aE = fftw_malloc(sizeof(fftw_complex) * NTG);

  printf("\n<< Starting reading and spline section >>\n");
  
  // eSPec output files have unitary euclidean norm, in order to renormalize it to the integral norm:
  norm = 1.0e+0/sqrt(stepx*stepy);
  //norm = 1.0;

  // convert time variables from fs to au
  ti = FSAU(ti);
  tf = FSAU(tf);
  stept = FSAU(stept);
  steptspl = FSAU(steptspl);

  //--- Read all wavepackets and generate spline coefficient files
  iflth = strlen(file);
  readallspl_(file,&iflth,X,Y,T,bcoefre,bcoefim,xknot,yknot,tknot,&nx,&ny,&nf,&kx,&nxr,&nyr,fpr,&norm,&prtRMSwhole,&noreg,&wholegrid,&stfil);

  //-------TD flux
  Prob = 0.0e+0;
  printf("\nTime dependent flux\n");
  maxF[1]=0.0e+0;

  for(k=0;k<nf;k++){
    if(type == 0) flux = nv*spl_integflux_x(bcoefre,bcoefim,xknot,yknot,tknot,kx,fluxline[0],Y,T[k],nxr,nyr,nyr,nf,1.0e-7,stepy,fp,m1);
    else if(type == 1) flux = nv*spl_jac_integflux_x(bcoefre,bcoefim,xknot,yknot,tknot,kx,Ri,ri,T[k],nxr,nyr,NYG,nf,1.0e-7,stepr,ma,mb,mc,m1beta);
    // find maximum--
    if(flux >= maxF[1]){
      maxF[1] = flux;
      maxF[2] = T[k];
    }//---
    Prob = Prob + flux*stept;
    printf("%E %E %E \n",T[k],flux,Prob);
    fprintf(tdprobf,"%E %E %E \n",T[k],flux,Prob);
  }
  fclose(tdprobf);
  printf("\n P_infty = %E",Prob);
  printf("\nFlux maximum at T = %E\n",maxF[2]);

  printf("\n Finished TD Flux section\n");
  printf("------------------------------------------------------\n\n\n");

  if(type==2){
  printf("\n# End of Calculation!\n");
  printf("# FluxVE terminated successfully!\n");
    return 100;
  }

  //if(type==1)return 5;
  //-----------* /

  /*
  printf("vect check!\n\n");
  printf("X\n");
  for(i=0;i<nxr;i++)printf("%E %E \n",xrange[0] + i*stepx,X[i]);
  printf("Y\n");
  for(i=0;i<nyr;i++)printf("%E %E \n",yrange[0] + i*stepy,Y[i]);
  printf("T\n");
  for(i=0;i<nf;i++)printf("%E %E \n",ti + i*stept,T[i]);
  printf("double T\n");
  for(i=0;i<NTG;i++)printf("%E \n",-tf + i*steptspl);
  */
  printf("\n<< Starting Fourier section >>\n");

  /* it does not improve
  //---- trying to improve R derivative for high energies using spline in the time domain.
  nxr = nxr*2.0;
  stepx = (X[nxr] - X[0])/(nxr - 1);
  xwork = X[0];
  free(X);
  X = malloc(nxr*sizeof(double));
  for(i=0;i<nxr;i++) X[i] = xwork + i*stepx; 
  
  WPERe = malloc(nxr*nyr*(NTG)*sizeof(double)); 
  WPEIm = malloc(nxr*nyr*(NTG)*sizeof(double)); 
  */
  //----------

  begint = clock();
  ftinfo=run_all_fft(bcoefre,bcoefim,X,Y,xknot,yknot,tknot,kx,ti,tf,steptspl,nxr,nyr,nf,NTG,fpr,width,WPERe,WPEIm,maxF[0],stren,t0);
  endt = clock();
  timediff = (double)(endt - begint)/CLOCKS_PER_SEC;

  stepE = 2*M_PI/(NTG*steptspl);
  Ei = -NTG*stepE/(2.0E+0);
  nE = NTG;
  
  if(ftinfo!=0){
    printf("\n An error occured while computing the wavepacket's Fourier transform: The program will stop!\n");
    return 666;
  }
  printf("\n Wavepacket has been trasformed to the Energy domain!\n");
  printf("\n it took %lf seconds \n",timediff);


   //checkmatrix_ (WPERe,WPEIm,&nxr,&nyr,&nE,X,Y,E);

  printf("\nGenerating 3D spline representation of the fourier transformed wavepacket\n");
  free(bcoefre);
  free(bcoefim);
  //free(T); can I free T ???
  free(xknot);
  free(yknot);
  free(tknot);
  xknot = malloc((nxr+kx)*sizeof(double));
  yknot = malloc((nyr+kx)*sizeof(double));
  eknot = malloc((nE +kx)*sizeof(double));
  E = malloc(nE*sizeof(double));
  bcoefre = malloc(nyr*nxr*nE*sizeof(double));
  bcoefim = malloc(nyr*nxr*nE*sizeof(double));

  for(k=0;k<nE;k++){
    E[k] = Ei + k*stepE;
  }

  //---
  begint = clock();
  dbsnak_ (&nyr, Y, &kx, yknot);
  dbsnak_ (&nxr, X, &kx, xknot);
  dbsnak_ (&nE, E, &kx, eknot);
  //----------
  dbs3in_ (&nyr,Y,&nxr,X,&nE,E,WPERe,&nyr,&nxr,&kx,&kx,&kx,yknot,xknot,eknot,bcoefre);
  dbs3in_ (&nyr,Y,&nxr,X,&nE,E,WPEIm,&nyr,&nxr,&kx,&kx,&kx,yknot,xknot,eknot,bcoefim);
  endt = clock();
  timediff = (double)(endt - begint)/CLOCKS_PER_SEC;
  //----------

  //free(WPERe);
  //free(WPEIm);

  printf("3D spline matrices calculated!\n");
  printf("it took %lf seconds\n",timediff);

  prtwpE=1;
  //-----print energy wavepackets ----------------------
  if(prtwpE==0){
    l=0;
    for(k=nE/2;k<nE;k=k + 3){
      strcpy(wp_Enam,"wpE_");
      sprintf(num,"%d.dat",l+1);
      strcat(wp_Enam,num);
      wp_E = fopen(wp_Enam,"w");
      val[2] = Ei + k*stepE;
      fprintf(wp_E,"# Energy = %E \n",val[2]);
      pk = sqrt(2.0*m1*(val[2] - Evib - potshift));
      fprintf(wp_E,"# Momentum = %E \n",pk);

      for(i=0;i<nxr;i++){
	for(j=0;j<nyr;j++){
	  val[0] = dbs3vl_ (&Y[j],&X[i],&val[2],&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);
	  val[1] = dbs3vl_ (&Y[j],&X[i],&val[2],&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
	  fprintf(wp_E,"%E %E %E %E\n",X[i],Y[j],val[0],val[1]);
	}
	fprintf(wp_E,"\n");
      }
      fclose(wp_E);
      l=l+1;
    }
  }
  //---------------------------------------------------
  //---DEBUG--
  //--- print only cuts of the energy transformed wavepackets along r = rmax
  /*
  for(k=nE/2 + 1;k<nE;k++){
    val[2] = Ei + k*stepE;
    pk = sqrt(2.0*m1*(val[2] - Evib - potshift));
    aE = (m1/pk)*gausswp_k(pk,k0,a)*gausswp_k(pk,k0,a);
    
    strcpy(wp_Enam,"cutwpE_");
    sprintf(num,"%.4lf.dat",pk);
    strcat(wp_Enam,num);
    wp_E = fopen(wp_Enam,"w");
    
    fprintf(wp_E,"# Energy = %E \n",val[2]);
    fprintf(wp_E,"# Momentum = %E \n",pk);
    fprintf(wp_E,"# r = %E \n",Y[26]);
    fprintf(wp_E,"# a(E) = %E \n",aE);

    for(i=0;i<nxr;i++){
      val[0] = dbs3vl_ (&Y[26],&X[i],&val[2],&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);
      val[1] = dbs3vl_ (&Y[26],&X[i],&val[2],&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
      fprintf(wp_E,"%E %E %E\n",X[i],val[0],val[1]);
    }
    fclose(wp_E); 
  }*/
    

  //---------------------------------------------------


  printf("\n Finished Fourier section\n");
  printf("------------------------------------------------------\n\n\n");
  

  printf("\nStarting Flux and Probability section\n\n");

  maxF[1]=0.0e+0;
  maxaE[1]=0.0e+0;

  stepx = 1.0e-7;
  //stepy = stepyspl;
  rmse = 0.0e+0;
  nrmse = 0.0e+0;
  //NYG = nyr;
  printf("\nchanging stepx to %E\n",stepx);
  printf("\nchanging stepy to %E\n",stepy);

  fprintf(probf,"#    Collision Energy  ,   Probability  \n");

  for(k=0;k<nE;k++){

    if(E[k] >= Emin && E[k] <= Emax){
      //flux = spl_integflux_x(bcoefre,bcoefim,xknot,yknot,eknot,kx,fluxline[0],Y,E[k],nxr,nyr,NYG,nE,stepx,stepy,fp,m1);
      if(type == 0) flux = nv*spl_integflux_x2(WPERe,WPEIm,X,fluxline[0],Y,E[k],k,nxr,nyr,NYG,nE,stepx,stepy,fp,m1);
      else if(type == 1) flux = nv*spl_jac_integflux_x(bcoefre,bcoefim,xknot,yknot,eknot,kx,Ri,ri,E[k],nxr,nyr,NYG,nE,stepR,stepr,ma,mb,mc,m1beta);

      pk = sqrt(2.0*m1*(E[k] - Evib - potshift));
      aE = gausswp_k(pk,k0,a);
      Prob = flux/((m1/pk)*aE*aE);
      nrmse = nrmse + 1;
      rmse = rmse + (Prob - 1.0e+0)*(Prob - 1.0e+0);
      // finding maxima ---
      if((m1/pk)*aE*aE >= maxaE[1]){
	maxaE[1] = (m1/pk)*aE*aE;
	maxaE[0] = E[k];
      }
      if(flux >= maxF[1]){
	maxF[1] = flux;
	maxF[0] = E[k];
      }
      //-------probf
      fprintf(probf,"%E %E \n", 27.2114*(E[k] - Evib - potshift),Prob);
      fprintf(deb,"%E %E %E %E \n",pk,aE*aE,(pk/m1)*flux,Prob);
      printf("%E %E %E %E \n",E[k],(m1/pk)*aE*aE,flux,Prob);
    }
  }
  fclose(probf);
  fclose(deb);
    
  rmse = sqrt(rmse/nrmse);
  printf("\n Root mean square error in relation to unity probability: %.15E \n\n",rmse);
  pk = sqrt(2.0*m1*(maxF[0] - Evib - potshift));
  printf("maximum Flux: E = %.15E, k = %E, F = %.15E \n",maxF[0],pk,maxF[1]);
  pk = sqrt(2.0*m1*(maxaE[0] - Evib - potshift));
  printf("maximum coeff: E = %.15E, k = %E, F = %.15E \n",maxaE[0],pk,maxaE[1]);
  pk = 0.5e+0*( k0 + sqrt(k0*k0 - 4*M_LN2/(a*a)) );
  maxaE[0] = pk*pk/(2.0*m1) + Evib + potshift;
  maxaE[1] = (m1/pk)*gausswp_k(pk,k0,a)*gausswp_k(pk,k0,a);
  printf("analytical maximum for coeff: E = %.15E, k = %E, F = %.15E \n",maxaE[0],pk,maxaE[1]);

  
  printf("\n Finished Flux and Probability section\n");
  printf("------------------------------------------------------\n\n\n");
  
  printf("\n# End of Calculation!\n");
  printf("# FluxVE terminated successfully!\n");
  return 0;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
