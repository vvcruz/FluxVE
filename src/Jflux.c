#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "jacobi.h"
#include "Jflux.h"
#include "splinesurf.h"
#include "conversion.h"

/*
Flux Vector Routines

Author: Vinícius V. Cruz

Purpose: Compute the J vector, and the flux along a line, either on the original wavefuntion or applying a spline procedure. Also computes the J vector performing a jacobi coordinate transformation.

Date: 7/05/2012

References: Cohen Tannoudji, Quantum Mechanics; Griffiths, Introduction to Quantum Mechanics;

 */


//DERIVATION METHODS
double deriv(double y1, double y2, double step)
{//numerical central difference of order O(h^2)
  double dydx;
  dydx= (y2 - y1)/(2.0*step);
  return dydx;
}

double derivO4(double y1, double y2, double y3, double y4, double step)
{//numerical central difference of order O(h^4)
  double dydx;
  dydx= (-y4 + 8.0*y3 - 8.0*y2 + y1)/(12.0*step);
  return dydx;
}



//FLUX ROUTINES
double prob_current1D(double *J, double *Re, double *Im, int nx, double stepx, int dyn, double m){
  // J = hbar/mass(Re(psi* d/dx psi))
  int i,j;
  double DRe,DIm,work;

    for(i=1;i<nx-1;i++){
      DRe=deriv(Re[i-1],Re[i+1],stepx);
      DIm=deriv(Im[i-1],Im[i+1],stepx);
      J[i] = (1/m)*(Re[i]*(DIm) - Im[i]*DRe);
      //      printf("%E %E \n",DRe,DIm);
      //      printf("%lf \n", (1/m)*(Re[i][0]*(DIm) - Im[i][0]*DRe));
    }
}

double prob_current2D(double **Jx, double **Jy, double **Re, double **Im, int nx, int ny, double stepx, double stepy, double m1, double m2, int order){
  // J = hbar/mass(Re(psi* (gradient) psi))
  int i,j;
  double DRex,DImx,work,DRey,DImy;
  FILE *debug;

  if(order==2){

    for(i=1;i<nx-1;i++){
      for(j=1;j<ny-1;j++){
	DRex=deriv(Re[i-1][j],Re[i+1][j],stepx);
	DImx=deriv(Im[i-1][j],Im[i+1][j],stepx);
	Jx[i][j] = (1/m1)*(Re[i][j]*(DImx) - Im[i][j]*DRex);
      }
    }
    for(i=1;i<nx-1;i++){
      for(j=1;j<ny-1;j++){
	DRey=deriv(Re[i][j-1],Re[i][j+1],stepy);
	DImy=deriv(Im[i][j-1],Im[i][j+1],stepy);
	Jy[i][j] = (1/m2)*(Re[i][j]*(DImy) - Im[i][j]*DRey);
      }
    }

  }else if(order==4){

    for(i=2;i<nx-2;i++){
      for(j=2;j<ny-2;j++){
	DRex=derivO4(Re[i-2][j],Re[i-1][j],Re[i+1][j],Re[i+2][j],stepx);
	DImx=derivO4(Im[i-2][j],Im[i-1][j],Im[i+1][j],Im[i+2][j],stepx);
	Jx[i][j] = (1/m1)*(Re[i][j]*(DImx) - Im[i][j]*DRex);
      }
    } 
    for(i=2;i<nx-2;i++){
      for(j=2;j<ny-2;j++){
	DRey=derivO4(Re[i][j-2],Re[i][j-1],Re[i][j+1],Re[i][j+2],stepy);
	DImy=derivO4(Im[i][j-2],Im[i][j-1],Im[i][j+1],Im[i][j+2],stepy);
	Jy[i][j] = (1/m2)*(Re[i][j]*(DImy) - Im[i][j]*DRey);
      }
    }
  }else{
    printf("Invalid derivation Method \n");
    return 9;
  }
}

double integflux_x(double **Re,double **Im,int fp0,int fp1,int fp2,double stepx, double stepy,double mx,double yi,int order){
  //calculates flux across a line for: X=x0, yi <= Y <=yf
  // does not compute entire J vector, only the points regarding the path of integration
  int i,j;
  double DRex,DImx,flux;
  //FILE *debug=fopen("debugderiv.dat","w");
  flux=0.0E+0;
  for(j=fp1;j<fp2;j++){
    if(order==2){
      DRex=deriv(Re[fp0-1][j],Re[fp0+1][j],stepx);
      DImx=deriv(Im[fp0-1][j],Im[fp0+1][j],stepx);
    }else if(order==4){
      DRex=derivO4(Re[fp0-2][j],Re[fp0-1][j],Re[fp0+1][j],Re[fp0+2][j],stepx);
      DImx=derivO4(Im[fp0-2][j],Im[fp0-1][j],Im[fp0+1][j],Im[fp0+2][j],stepx);
    }else{
      printf("invalid derivation method \n");
      return 1;
    }
    flux=flux+ (1/mx)*(Re[fp0][j]*DImx - Im[fp0][j]*DRex)*stepy;
  }
  return flux;
}

double integflux_y(double **Re,double **Im,int fp0,int fp1,int fp2,double stepx, double stepy,double my,double xi,int order){
  //calculates flux across a line for: Y=x0, xi <= X <=xf
  // does not compute entire J vector, only the points regarding the path of integration
  int i,j;
  double DRey,DImy,flux;
  //FILE *debug=fopen("debugderiv.dat","w");
  flux=0.0E+0;
  for(j=fp1;j<fp2;j++){
    if(order==2){
      DRey=deriv(Re[j][fp0-1],Re[j][fp0+1],stepy);
      DImy=deriv(Im[j][fp0-1],Im[j][fp0+1],stepy);
    }else if(order==4){
      DRey=derivO4(Re[j][fp0-2],Re[j][fp0-1],Re[j][fp0+1],Re[j][fp0+2],stepy);
      DImy=derivO4(Im[j][fp0-2],Im[j][fp0-1],Im[j][fp0+1],Im[j][fp0+2],stepy);
    }else{
      printf("invalid derivation method \n");
      return 1;
    }
    flux=flux+ (1/my)*(Re[j][fp0]*DImy - Im[j][fp0]*DRey)*stepx;
  }
  return flux;
}

double integflux_x_jac(double R0, double ri,int fp1, int fp2,double *XKNOT, double *YKNOT,int nx,int ny,double *BCOEFRe, double *BCOEFIm, double stepR, double stepr,double mx,int order,double ma,double mb,double mc){
  //calculates flux across a line for: X=x0, yi <= Y <=yf  performing Jacobi transformation
  // does not compute entire J vector, only the points regarding the path of integration
  int i,j,KXORD,KYORD;
  double DRex,DImx,flux,X[5],Y[5],Re[5],Im[5],R[5],r;
  //FILE *debug=fopen("debugderivspl.dat","w");
  flux=0.0E+0;
  KXORD=3.0;
  KYORD=3.0;  
  R[0]=R0 - stepR;
  R[1]=R0;
  R[2]=R0 + stepR;
  for(j=fp1;j<fp2;j++){
    r=ri + j*stepr;
    Jacobitrans(&X[0], &Y[0],ma,mb,mc,&R[0],&r,-1);
    Jacobitrans(&X[1], &Y[1],ma,mb,mc,&R[1],&r,-1);
    Jacobitrans(&X[2], &Y[2],ma,mb,mc,&R[2],&r,-1);
    Re[0] = dbs2vl_ (&Y[0],&X[0],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
    Re[1] = dbs2vl_ (&Y[1],&X[1],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
    Re[2] = dbs2vl_ (&Y[2],&X[2],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
    Im[0] = dbs2vl_ (&Y[0],&X[0],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
    Im[1] = dbs2vl_ (&Y[1],&X[1],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
    Im[2] = dbs2vl_ (&Y[2],&X[2],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
    if(order==2){
      DRex=deriv(Re[0],Re[2],stepR);
      DImx=deriv(Im[0],Im[2],stepR);
    }else if(order==4){
      R[3]=R0 + 2*stepR;
      R[4]=R0 - 2*stepR;      
      Jacobitrans(&X[3], &Y[3],ma,mb,mc,&R[3],&r,-1);
      Jacobitrans(&X[4], &Y[4],ma,mb,mc,&R[4],&r,-1);     

      Re[3] = dbs2vl_ (&Y[3],&X[3],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
      Re[4] = dbs2vl_ (&Y[4],&X[4],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
      Im[3] = dbs2vl_ (&Y[3],&X[3],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
      Im[4] = dbs2vl_ (&Y[4],&X[4],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
      DRex=derivO4(Re[4],Re[0],Re[2],Re[3],stepR);
      DImx=derivO4(Im[4],Im[0],Im[2],Im[3],stepR);
    }else{
      printf("invalid derivation method \n");
      return 1;
    }

    flux=flux+ (1/mx)*(Re[1]*DImx - Im[1]*DRex)*stepr;
    //if(f==15) fprintf(debug,"%E %E %E %E %E \n", ri+j*stepr,DRex,DImx,Re[1],Im[1]);
  }
  return flux;
}

double integflux_x_spl(double x0, double yi,int fp1, int fp2,double *XKNOT, double *YKNOT,int nx,int ny,double *BCOEFRe, double *BCOEFIm, double stepx, double stepy,double mx,int order){
  //calculates flux across a line for: X=x0, yi <= Y <=yf
  // does not compute entire J vector, only the points regarding the path of integration
  int i,j,KXORD,KYORD;
  double DRex,DImx,flux,X[5],Y,Re[5],Im[5];
  //FILE *debug=fopen("debugderivspl.dat","w");
  flux=0.0E+0;
  KXORD=3.0;
  KYORD=3.0;  
  X[0]= x0 - stepx;
  X[1]= x0;
  X[2]= x0 + stepx;
  for(j=fp1;j<fp2;j++){
    Y= yi + j*stepy;
    Re[0] = dbs2vl_ (&Y,&X[0],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
    Re[1] = dbs2vl_ (&Y,&X[1],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
    Re[2] = dbs2vl_ (&Y,&X[2],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
    Im[0] = dbs2vl_ (&Y,&X[0],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
    Im[1] = dbs2vl_ (&Y,&X[1],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
    Im[2] = dbs2vl_ (&Y,&X[2],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
    if(order==2){
      DRex=deriv(Re[0],Re[2],stepx);
      DImx=deriv(Im[0],Im[2],stepx);
    }else if(order==4){
      X[3]= x0 + 2*stepx;
      X[4]= x0 - 2*stepx;
      Re[3] = dbs2vl_ (&Y,&X[3],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
      Re[4] = dbs2vl_ (&Y,&X[4],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
      Im[3] = dbs2vl_ (&Y,&X[3],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
      Im[4] = dbs2vl_ (&Y,&X[4],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
      DRex=derivO4(Re[4],Re[0],Re[2],Re[3],stepx);
      DImx=derivO4(Im[4],Im[0],Im[2],Im[3],stepx);
    }else{
      printf("invalid derivation method \n");
      return 1;
    }
    flux=flux+ (1/mx)*(Re[1]*DImx - Im[1]*DRex)*stepy;
    
    //if(f==15) fprintf(debug,"%E %E %E %E %E \n", yi+j*stepy,DRex,DImx,Re[1],Im[1]);
  }
  return flux;
}


double integflux_y_spl(double y0, double xi,int fp1, int fp2,double *XKNOT, double *YKNOT,int nx,int ny,double *BCOEFRe, double *BCOEFIm, double stepx, double stepy,double my,int order){
  //calculates flux across a line for: Y=y0, xi <= X <=xf
  // does not compute entire J vector, only the points regarding the path of integration
  int i,j,KXORD,KYORD;
  double DRey,DImy,flux,X,Y[5],Re[5],Im[5];
  //FILE *debug=fopen("debugderivspl.dat","w");
  flux=0.0E+0;
  KXORD=3.0;
  KYORD=3.0;  
  Y[0]= y0 - stepy;
  Y[1]= y0;
  Y[2]= y0 + stepy;
  for(j=fp1;j<fp2;j++){
    X= xi + j*stepx;
    Re[0] = dbs2vl_ (&Y[0],&X,&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
    Re[1] = dbs2vl_ (&Y[1],&X,&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
    Re[2] = dbs2vl_ (&Y[2],&X,&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
    Im[0] = dbs2vl_ (&Y[0],&X,&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
    Im[1] = dbs2vl_ (&Y[1],&X,&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
    Im[2] = dbs2vl_ (&Y[2],&X,&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
    if(order==2){
      DRey=deriv(Re[0],Re[2],stepy);
      DImy=deriv(Im[0],Im[2],stepy);
    }else if(order==4){
      Y[3]= y0 + 2*stepy;
      Y[4]= y0 - 2*stepy;
      Re[3] = dbs2vl_ (&Y[3],&X,&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
      Re[4] = dbs2vl_ (&Y[4],&X,&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
      Im[3] = dbs2vl_ (&Y[3],&X,&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
      Im[4] = dbs2vl_ (&Y[4],&X,&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
      DRey=derivO4(Re[4],Re[0],Re[2],Re[3],stepy);
      DImy=derivO4(Im[4],Im[0],Im[2],Im[3],stepy);
    }else{
      printf("invalid derivation method \n");
      return 1;
    }
    flux=flux+ (1/my)*(Re[1]*DImy - Im[1]*DRey)*stepx;
    
    //if(f==15) fprintf(debug,"%E %E %E %E %E \n", yi+j*stepy,DRex,DImx,Re[1],Im[1]);
  }
  return flux;
}



/*
 spl_integflux_x
 computes flux across a line for: X=x0, yi <= Y <=yf, given a wavepacket's spline representation
 */
double spl_integflux_x(double *bcoefre,double *bcoefim, double *xknot, double *yknot, double *eknot, int kx,double x0, double *Y, double E,int nxr,int nyr, int NYG, int nE,double stepx,double stepy, int *fp,double mx){
  int i,j,l,k;
  double DRex,DImx,Jflux[NYG],flux,Re[5],Im[5],X[5],val[5];
  double YSPL;
  //-----
  //FILE *debug=fopen("debugderiv.dat","w");
  //double ansDRex,ansDImx;
  //fprintf(debug,"#X Y DRe ansDRe DIm ansDIm \n");
  //printf("NYG = %d Y[fp[2]] = %E \n",NYG,Y[fp[2]]);
  flux=0.0E+0;
  //kx=3.0;
  X[0]= x0;
  X[1]= x0 - 2*stepx;
  X[2]= x0 - stepx;
  X[3]= x0 + stepx;
  X[4]= x0 + 2*stepx;
  //  for(j=fp[2];j<fp[3];j++){
  for(j=0;j<NYG;j++){
    YSPL = Y[fp[2]] + j*stepy;
    //fprintf(debug,"%E \n",YSPL);
    if(YSPL > Y[fp[3]]) break;

    Re[0] = dbs3vl_ (&YSPL,&X[0],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);
    Re[1] = dbs3vl_ (&YSPL,&X[1],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);
    Re[2] = dbs3vl_ (&YSPL,&X[2],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);
    Re[3] = dbs3vl_ (&YSPL,&X[3],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);
    Re[4] = dbs3vl_ (&YSPL,&X[4],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);

    Im[0] = dbs3vl_ (&YSPL,&X[0],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
    Im[1] = dbs3vl_ (&YSPL,&X[1],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
    Im[2] = dbs3vl_ (&YSPL,&X[2],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
    Im[3] = dbs3vl_ (&YSPL,&X[3],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
    Im[4] = dbs3vl_ (&YSPL,&X[4],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
    
    DRex=derivO4(Re[1],Re[2],Re[3],Re[4],stepx);
    DImx=derivO4(Im[1],Im[2],Im[3],Im[4],stepx);
    Jflux[j]=(1.0e+0/mx)*(Re[0]*DImx - Im[0]*DRex);
    flux = flux + Jflux[j]*stepy;
    //printf(" %E  %E %E %E %E \n",YSPL,DRex,DImx,Jflux[j],flux);
    //fprintf(debug,"%E %E %E %E %E \n", Y[j],DRex,DImx,Re[1],Im[1]);
  }
  //flux = comp_simpson38(Jflux,stepy,NYG);
  //printf("%E \n",flux);
  return flux;
}


//-----------------------------

/*
 spl_integflux_x
 computes reactve flux across a line for: R=R0, ri <= r <=rf, given a wavepacket's spline representation
 and the respective Jacobi coordinate transformation
 */
double spl_jac_integflux_x(double *bcoefre,double *bcoefim, double *xknot, double *yknot, double *eknot, int kx,double R0, double ri, double E,int nxr,int nyr, int NYG, int nE,double stepR,double stepr,double ma,double mb, double mc, double mx){
  int i,j,l,k,m;
  double DRex,DImx,Jflux[NYG],flux,Re[5],Im[5],X[5],val[5],R[5],rspl;
  double Y[5];
  //-----
  //FILE *debug=fopen("debugderiv.dat","w");
  //double ansDRex,ansDImx;
  //fprintf(debug,"#X Y DRe ansDRe DIm ansDIm \n");
  //printf("NYG = %d Y[fp[2]] = %E \n",NYG,Y[fp[2]]);
  flux=0.0E+0;
  //kx=3.0;
  R[0]= R0;
  R[1]= R0 - 2*stepR;
  R[2]= R0 - stepR;
  R[3]= R0 + stepR;
  R[4]= R0 + 2*stepR;
  //  for(j=fp[2];j<fp[3];j++){
  for(j=0;j<NYG;j++){
    rspl = ri + j*stepr;
    //fprintf(debug,"%E \n",YSPL);

    for(m=0;m<5;m++) Jacobitrans(&X[m],&Y[m],ma,mb,mc,&R[m],&rspl,-1);

    Re[0] = dbs3vl_ (&Y[0],&X[0],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);
    Re[1] = dbs3vl_ (&Y[1],&X[1],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);
    Re[2] = dbs3vl_ (&Y[2],&X[2],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);
    Re[3] = dbs3vl_ (&Y[3],&X[3],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);
    Re[4] = dbs3vl_ (&Y[4],&X[4],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefre);

    Im[0] = dbs3vl_ (&Y[0],&X[0],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
    Im[1] = dbs3vl_ (&Y[1],&X[1],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
    Im[2] = dbs3vl_ (&Y[2],&X[2],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
    Im[3] = dbs3vl_ (&Y[3],&X[3],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
    Im[4] = dbs3vl_ (&Y[4],&X[4],&E,&kx,&kx,&kx,yknot,xknot,eknot,&nyr,&nxr,&nE,bcoefim);
    
    DRex=derivO4(Re[1],Re[2],Re[3],Re[4],stepR);
    DImx=derivO4(Im[1],Im[2],Im[3],Im[4],stepR);
    Jflux[j]=(1.0e+0/mx)*(Re[0]*DImx - Im[0]*DRex);
    flux = flux + Jflux[j]*stepr;
    //printf(" %E  %E %E %E %E \n",YSPL,DRex,DImx,Jflux[j],flux);
    //fprintf(debug,"%E %E %E %E %E \n", Y[j],DRex,DImx,Re[1],Im[1]);
  }
  //flux = comp_simpson38(Jflux,stepy,NYG);
  //printf("%E \n",flux);
  return flux;
}

//---------------------------



/*
  perform integration using the composite simpson's 3/8 rule

  if the number of intervals is not divisible by 3, then the first few 
  values are estimated using the rectangles method, and the 3*m points 
  are computed using simpson's 3/8 formula

*/

double comp_simpson38(double *fx,double stepx,int n){
  int i,j,m,remd,ini;
  double val[4],integ;

  integ = 0.0e+0;
  ini = 0;
  remd = (n - 1)%3;
  
  ini = remd; 
  m = (n -1 - remd)/3;
  for(i=0;i<remd;i++) integ = integ + fx[i]*stepx;

  for(j=1;j<m+1;j++){
    val[0] = fx[3*j - 3 + ini];
    val[1] = fx[3*j - 2 + ini];
    val[2] = fx[3*j - 1 + ini];
    val[3] = fx[3*j + ini];
    integ = integ + simpson38(val,stepx);
  }

  return integ;
}


/*
  perform integration using the simpson's 3/8 rule
*/
double simpson38(double *y,double stepx){
  double integ;

  integ = (3.0e+0*stepx/8.0e+0)*(y[0] + 3*y[1] + 3*y[2] + y[3]);
  return integ;
}



/*
  perform integration using the simpson's 1/3 rule
*/
double simpson13(double *y,double stepx){
  double integ;

  integ = (stepx/3.0e+0)*(y[0] + 4*y[1] + y[2]);
  return integ;
}


//----------------------

/*
 spl_integflux_x2
 computes flux across a line for: X=x0, yi <= Y <=yf,
 using splines to increase the accuracy of the derivative
 */
double spl_integflux_x2(double *WPERe,double *WPEIm,double *XX,double x0, double *Y, double E, int kE,int nxr,int nyr, int NYG, int nE,double stepx,double stepy, int *fp,double mx){
  int i,j,l,k,ll,kx;
  double DRex,DImx,Jflux[NYG],flux,Re[5],Im[5],X[5],val[5],stepxx;
  kx = 5.0;
  double YSPL,xknot[nxr+kx],valRe[nxr],valIm[nxr];
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *bcoefre,*bcoefim;
  //-----
  //FILE *debug=fopen("debugderiv.dat","w");
  //double ansDRex,ansDImx;
  //fprintf(debug,"#X Y DRe ansDRe DIm ansDIm \n");
  //printf("NYG = %d Y[fp[2]] = %E \n",NYG,Y[fp[2]]);
  flux=0.0E+0;
  //kx=5.0;
  X[0]= x0;
  X[1]= x0 - 2*stepx;
  X[2]= x0 - stepx;
  X[3]= x0 + stepx;
  X[4]= x0 + 2*stepx;

  dbsnak_ (&nxr,XX,&kx,xknot);

  //for(i=0;i<5;i++)printf("%E ",X[i]);
  //printf("\n");
  //  for(j=fp[2];j<fp[3];j++){
  //for(j=0;j<NYG;j++){
  for(j=0;j<nyr;j++){
    YSPL = Y[j];
    //fprintf(debug,"%E \n",YSPL);
    if(YSPL > Y[fp[3]]) break;

    for(i=0;i<nxr;i++){
      ll = j + i*nyr + kE*nxr*nyr;
      //printf("%E \n",XX[i]);
      valRe[i] = WPERe[ll];
      valIm[i] = WPEIm[ll];
      //if(kE==600)printf("%E %E \n",valRe[i],valIm[i]);
    }
    //bcoefre = malloc(nxr*sizeof(double));
    //bcoefim = malloc(nxr*sizeof(double));
    bcoefre = gsl_spline_alloc (gsl_interp_cspline, nxr);
    bcoefim = gsl_spline_alloc (gsl_interp_cspline, nxr);

    gsl_spline_init (bcoefre, XX, valRe, nxr);
    gsl_spline_init (bcoefim, XX, valIm, nxr);
    //dbsint_ (&nxr,XX,valRe,&kx,xknot,bcoefre);
    //dbsint_ (&nxr,XX,valIm,&kx,xknot,bcoefim);

    Re[0] = gsl_spline_eval (bcoefre, X[0], acc);
    //Re[0] = dbsval_ (&X[0],&kx,xknot,&nxr,bcoefre);
    //Re[1] = dbsval_ (&X[1],&kx,xknot,&nxr,bcoefre);
    //Re[2] = dbsval_ (&X[2],&kx,xknot,&nxr,bcoefre);
    //Re[3] = dbsval_ (&X[3],&kx,xknot,&nxr,bcoefre);
    //Re[4] = dbsval_ (&X[4],&kx,xknot,&nxr,bcoefre);

    Im[0] = gsl_spline_eval (bcoefim, X[0], acc);
    //Im[0] = dbsval_ (&X[0],&kx,xknot,&nxr,bcoefim);
    //Im[1] = dbsval_ (&X[1],&kx,xknot,&nxr,bcoefim);
    //Im[2] = dbsval_ (&X[2],&kx,xknot,&nxr,bcoefim);
    //Im[3] = dbsval_ (&X[3],&kx,xknot,&nxr,bcoefim);
    //Im[4] = dbsval_ (&X[4],&kx,xknot,&nxr,bcoefim);

    //if(kE==600)printf("> %E %E %E %E %E %E %E %E\n",Re[0],Re[1],Re[2],Re[3],Im[0],Im[1],Im[2],Im[3]);
    
    //DRex=derivO4(Re[1],Re[2],Re[3],Re[4],stepx);
    //DImx=derivO4(Im[1],Im[2],Im[3],Im[4],stepx);
    DRex = gsl_spline_eval_deriv (bcoefre,X[0],acc);
    DImx = gsl_spline_eval_deriv (bcoefim,X[0],acc);

    Jflux[j]=(1.0e+0/mx)*(Re[0]*DImx - Im[0]*DRex);
    flux = flux + Jflux[j]*stepy;
    //printf(" %E  %E %E %E %E \n",YSPL,DRex,DImx,Jflux[j],flux);
    //fprintf(debug,"%E %E %E %E %E \n", Y[j],DRex,DImx,Re[1],Im[1]);
    free(bcoefre);
    free(bcoefim);
  }
  //flux = comp_simpson38(Jflux,stepy,NYG);
  //printf("%E \n",flux);
  return flux;
}
