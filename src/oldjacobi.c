#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/*

Jacobi transformation routine

Author: Vinícius V. Cruz

Purpose: Transform between Jacobi coordinate sets as following:

    Alpha arrangement                       Beta arrangement
                           
    A     B       C       -(s=-1)->        A     B       C
     __r__                <-(s=1)-               ___r___
       ____R______                          ____R______


    (  Rbeta )  =  ( -Xi  lambda*Xi -1 )( Ralpha )
    (  rbeta )  =  (  1     -lambda    )( ralpha )

    ( Ralpha )  =  ( -lambda  1 -lambda*Xi )( Rbeta )
    ( ralpha )  =  (   -1           -Xi    )( rbeta )

    lambda = ma/(ma+mb)
    Xi = mc/(mb+mc)

Date: 20/04/2012

References: 

*/


void Jacobitrans(double *Ralpha, double *ralpha,double ma, double mb, double mc, double *Rbeta, double *rbeta, int s){

  double Xi,lambda;
  Xi= mc/(mb + mc);
  lambda= ma/(ma + mb);

  //printf("Xi: %lf lambda: %lf \n",Xi,lambda);

  if(s==1){
    *Rbeta= -Xi * (*Ralpha) + (lambda * Xi -1)* (*ralpha);
    *rbeta= (*Ralpha) -lambda * (*ralpha);
  } else if(s==-1){
    *Ralpha= -lambda * (*Rbeta) + (1 -lambda * Xi)* (*rbeta);
    *ralpha= -(*Rbeta) -Xi * (*rbeta);
  }else{
    printf("invalid option of transformation! \n");
  }

}
