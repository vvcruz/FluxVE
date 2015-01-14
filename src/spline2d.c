#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "nrutil.h"

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a);
void splin2(float x1a[], float x2a[], float **ya, float **y2a, int m, int n, float x1, float x2, float *y);

//-----------------------

void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a)
/*Given an m by n tabulated function ya[1..m][1..n] , and tabulated independent variables
x2a[1..n], this routine constructs one-dimensional natural cubic splines of the rows of ya
and returns the second-derivatives in the array y2a[1..m][1..n] . (The array x1a[1..m] is
ncluded in the argument list merely for consistency with routine splin2.)*/
{
  printf("inside.. toasty \n");
    int j;
    for (j=1;j<=m;j++){
         spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);  
    }   

    printf("END of splie2");
}  

void splin2(float x1a[], float x2a[], float **ya, float **y2a, int m, int n, float x1, float x2, float *y)
/*Given x1a, x2a, ya, m, n as described in splie2 and y2a as produced by that routine; and
given a desired interpolating point x1,x2; this routine returns an interpolated function value y
by bicubic spline interpolation.*/
{
    int j;
    float *ytmp,*yytmp;
    ytmp=vector(1,m);
    yytmp=vector(1,m);                 //Perform m evaluations of the row splines constructed by
    for (j=1;j<=m;j++){                     //splie2, using the one-dimensional spline evaluator
      splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);               //splint.
      spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);                 //Construct the one-dimensional col-
      splint(x1a,yytmp,ytmp,m,x1,y);                              //umn spline and evaluate it.
      free_vector(yytmp,1,m);
      free_vector(ytmp,1,m);
    }
}

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[])
/*Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f (xi ), with
x1 < x2 < . . . < xN , and given values yp1 and ypn for the first derivative of the interpolating
function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
the second derivatives of the interpolating function at the tabulated points xi . If yp1 and/or
ypn are equal to 1 × 1030 or larger, the routine is signaled to set the corresponding boundary
condition for a natural spline, with zero second derivative on that boundary.*/
{
     int i,k;
     float p,qn,sig,un,*u;
     u=vector(1,n-1);
     if (yp1 > 0.99e30){                    //The lower boundary condition is set either to be "nat-
       y2[1]=u[1]=0.0;  
     }                     // ural"
     else {                                //or else to have a specified first derivative.
         y2[1] = -0.5;
         u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
     }

     for (i=2;i<=n-1;i++) {          //This is the decomposition loop of the tridiagonal al-
       sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);           //gorithm. y2 and u are used for tem-
       p=sig*y2[i-1]+2.0;                           //porary storage of the decomposed
       y2[i]=(sig-1.0)/p;                           //factors.
       u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
       u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
     }
     if (ypn > 0.99e30) qn=un=0.0; //The upper boundary condition is set either to be "natural"
     else {                          //or else to have a specified first derivative.
       qn=0.5;
       un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
     }
     y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
     for (k=n-1;k>=1;k--){            //This is the backsubstitution loop of the tridiagonal
       y2[k]=y2[k]*y2[k+1]+u[k];       //algorithm.
     }     
     free_vector(u,1,n-1);
}


void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)
/*Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai 's in order),
and given the array y2a[1..n], which is the output from spline above, and given a value of
x, this routine returns a cubic-spline interpolated value y.*/
{
     int klo,khi,k;
     float h,b,a;
     klo=1;                            //We will find the right place in the table by means of
     khi=n;                                 //bisection. This is optimal if sequential calls to this
     while (khi-klo > 1) {                  //routine are at random values of x. If sequential calls
       k=(khi+klo) >> 1;                 //are in order, and closely spaced, one would do better
       if (xa[k] > x) khi=k;             //to store previous values of klo and khi and test if
       else klo=k;                      // they remain appropriate on the next call.
     }                                 //klo and khi now bracket the input value of x.
     h=xa[khi]-xa[klo];
     if (h == 0.0) nrerror("Bad xa input to routine splint"); //The xa's must be dis-
     a=(xa[khi]-x)/h;                                                             //tinct.
     b=(x-xa[klo])/h;                  //Cubic spline polynomial is now evaluated.
     *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
