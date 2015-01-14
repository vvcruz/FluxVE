#include<stdio.h>
#include<stdlib.h>
#include<math.h>


//double projector(double *u, double **psi, int nx, int ny, double stepx, double stepy){
main(){
  int i,j,k,l,m,nx,ny;
  FILE *arq=fopen("matrix","r");
  nx=3;
  ny=4;
  double u[20],v[nx][ny];
  double projmatrx[ny][ny],Pi,work;

  Pi=0;
  
  for(i=0;i<ny;i++){
    for(j=0;j<ny;j++){
      // projmatrx[i][j]=u[i]*u[j];
    }
  }
  for(i=0;i<nx;i++){
    k=0;
    for(l=0;l<ny;l++){
      k=k+1;
      for(j=0;j<nx;j++){
	//v[i][k]=v[i][k]+(projmatrix[i][j]*psi[j][l]);
	//printf(" v[%d][%d] = pmatrix[%d][%d]  psi[%d][%d]",i,k,i,j,j,l);
	printf(" pmatrix[%d][%d]  psi[%d][%d]",i,j,j,l);
      }
      printf("\n");
    }
  }
  
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      //Pi=Pi+(psi[i][j]*v[i][j])*stepx*stepy;
    }
  }

  return Pi; 
}
