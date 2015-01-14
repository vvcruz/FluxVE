#include<stdio.h>


main(){
  int nx,ny,neig,i,j,k;
  double x,y,eig[10];
  FILE *arq;
  char nam[20];
  FILE *out=fopen("plot_eig.out","w");

  printf("enter nx ny\n");
  scanf("%d %d",&nx,&ny);
  printf("enter file");
  scanf("%s",nam);
  arq=fopen(nam,"r");

  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      fscanf(arq,"%lf %lf %lf %lf",&x,&y,&eig[0],&eig[1]);
      fprintf(out,"%E %E %E %E\n",x,y,eig[0],eig[1]);
    }
    fprintf(out,"\n");
  }

}
