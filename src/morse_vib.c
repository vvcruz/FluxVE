#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

/*
  Vibrational Parameters
  
  Purpose: Compute vibrational parameters of diatomic molecules
           applying a Linear Morse Parametrization of the eigenvalues:
	   (E_n+1 - E_n)/hc = -2*WeXe*n + (We -2*WeXe)

 Author: Vinicius V. Cruz

 Date: fist version: february of 2011
       this version: 25/05/2012 (Towel day)

References: Mario F. Triola, Elementary Statistics, 10th ed. 2008
            

 */

int vibrationalp(int n, double E[n])
{
  // n is the number of eigenvalues
  // a partir de seus autovalores.  
  //printf("entre com numero de auto estados \n");
  int m,i,j,k,lixo,l;
  double worka,workb,workc,workd,worke;
  double sumx, sumy, sumy2,sumxy,sumx2,p,sum2x,sum2y,R,Xm,Ym,a,b;
  double F,B,deltE[n-1],We,WeXe;
  F=4.3597482E-18;
  //h=6.626E-34 C=2.998E8, B = h*C
  B=1.9864748E-25;
  m=n-1;
  sumx=0;
  sumy=0;
  sumx2=0;
  sumy2=0;
  sumxy=0;
  for(i=0;i<n;i++){
    //fscanf(arq, "%d %lf", &lixo, &E[i]);
    E[i]=(E[i]*F)/B;
  }
  for(k=0; k<n-1;k++){
    deltE[k]= E[k+1] - E[k];
    //conversão para cm-1
    deltE[k]=deltE[k]/100;
    // parametros da regressão linear
    sumx=sumx + k;
    sumx2=sumx2+(k * k);
    sumxy= sumxy + (deltE[k] * k);
    sumy=sumy+deltE[k];
    sumy2=sumy2+(deltE[k] * deltE[k]);
  } 
  Xm=sumx/m;
  Ym=sumy/m;
  sum2x=sumx * sumx;
  sum2y=sumy * sumy;
  
  //coeficiente de correlação
  
  worka=sqrt((m * sumx2) - sum2x);
  workb=sqrt((m * sumy2) - sum2y);
  workc=( (m * sumxy)- ((sumx)*(sumy)) ) ;
  R= workc/((worka) * (workb));
  // reta de ajuste: Y = ax + b
  workd=(m * sumxy) - (sumx * sumy) ;
  worke=(m * sumx2) - sum2x;
  a=workd/worke;
  b= Ym - (a * Xm); 

  // deltE(n) = We(1-2Xe) - 2WeXe*n  
  WeXe=((-1)*(a))/2;
  We= (2*WeXe)+b;

  //debug of statistical parameters
  //  printf(">>>>> %lf %E %E %E %E %E %E \n", sumx, sumy, sumxy, sumx2, sumy2, sum2x, sum2y);
  //  printf("X médio = %E Y médio = %E\n", Xm, Ym);
  
  printf("\n------------------------------ \n");
   printf("    Vibrational Analysis   \n   Morse Potential Parametrization \n");
  printf("------------------------------\n\n");
  printf("eigenstates included: %d \n\n", n);
  printf("Linear Correlation Coefficient: %lf \n", R);
  printf("Regression Line: Y = %lf*X + %lf \n",a,b);
  printf("Vibrational Parameters: \n We= %lf \n WeXe= %lf \n\n",We, WeXe);
  // printf("n   deltE /cm-1 \n");  
  //for(l=0;l<n-1;l++)printf("%d %E\n", l, deltE[l]);

  return 0;
}
