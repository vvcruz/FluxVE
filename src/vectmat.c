#include <stdio.h>
#include <stdlib.h>

main(){
  int i,j,l,nx,ny;
  nx=3;
  ny=3;
  double B[nx*ny];

  for(i=0;i<(nx*ny);i++){
    scanf("%lf",&B[i]);
  }

  for(i=0;i<nx;i++){
    l=i*nx -1;
    for(j=0;j<ny;j++){
      l=l+1;
      A[i][j]=B[l];
      printf("%lf \n", A[i][j]);
    }
  }

}
