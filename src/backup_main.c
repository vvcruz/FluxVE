#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include "spline2d.h"

double prob_current1D(double *J, double *Re, double *Im, int nx, double stepx, int dyn, double m);
double prob_current2D(double **Jx, double **Jy, double **Re, double **Im, int nx, int ny, double stepx, double stepy, double m1, double m2);
double deriv(double y1, double y2, double step);
void read_ (double *X, double *Y, char *fil, int *nx, int *ny);
void readspl_ (double *X, double *Y, char *fil, int *nx, int *ny, int *ldz);
void readpot_ (double *X, char *fil, int *nx, int *ny, int *ldz);
double real_normalize(double *X, int n, double step);
double project(double *u, double *Re, double *Im,int n, double step);
void mtrxdiag_ (char *dim, int *il, int *iu, int *info, int *mxdct, int *n, double *abstol, int *iwork, int *np, double *eigvl, double *shm, double *vpot, double *work, double *wk, double *eigvc);
double projector(double *u, double **psi, int nx, int ny, double stepx, double stepy);

void dbsnak_ (int *nx,double *xvec, int *kxord, double *xknot);
void dbs2in_ (int *nx,double *xvec,int *ny,double *yvec,double *xydata,int *ldf,int *kx,int *ky, double *xknot,double *yknot,double *bcoef);
void coefcal_ (int *nx,double *xvec,int *ny,double *yvec,double *xydata,int *ldf,int *kx,int *ky, double *xknot,double *yknot,double *bcoef);
void splvalue_ (double *x,double *y,int *kx,int *ky,double *xknot,double *yknot, int *nx,int *ny,double *bcoef,double *Val);
double dbs2vl_ (double *x,double *y,int *kx,int *ky,double *xknot,double *yknot, int *nx,int *ny,double *bcoef);


void main ()

{
  int i,j,k,l,d,m,dim,narq,nx,ny,nf,dyn,fp[4],prtj,f,proj,nst,npp,stfil,endfil,itt,prtpot,info,*iwork,iil,nptu,mxdct,np,prteig,spl,jactrans,LDZ;
  double *Jx, **Jx2, **Jy2,*Re, *Im,**Re2, **Im2, **eigenst/*, **Deigenst,**Ji*/,probs,**projmat, *pot1d; // J[nx][ny], Re[nx][ny] ...
  double xi,xf,yi,yf,stepx,stepy,ti,tf,stept,m1,m2,stepp,projet[20],shm,abstol,eigvl[20],*wk1,*wk2,*gst,dotp,wor;
  FILE *arq, *out, *wp,*fout,*jout,*poten,*eig,*eigg,*projfil,*derivatives,*check;
  char workk[20],file[10],filen[13],filep[13],num[5],outt[20],numout[5],axis,potf[20],dimmm[3];
  double  NAV, Me, flux, prob,fluxline[4],norm,norm2,projRe[20],projIm[20],Pii[20],ReD,ImD,projDRe[20],projDIm[20],projsum,fluxproj;
  //spline variables
  LDZ=1200;
  int KXORD,KYORD,NXG,NYG,fpspl[4];
  double *BCOEFRe,*BCOEFIm,*BCOEFP,*X,*Y,*XKNOT,*YKNOT,*Z;
  double stepxspl,stepyspl,Xspl[3],Yspl[3],Respl[3],Imspl[3];
  double *R, *r, Ri, Rf,ri, rf,ma, mb, mc;
  KXORD=3;
  KYORD=3;
  //---------------
  NAV=6.0221415E+23;
  Me=9.10938188E-31;

  arq=fopen("flux.inp","r");
  out=fopen("output","w");
  fout=fopen("fluxx.out","w");
  projfil=fopen("projections.dat","w");
  derivatives=fopen("derivs.dat","w");
  //default values-----------------------------------//
  dyn=1;
  stept=45.0;
  ny=1;
  m2=1;
  prob=0;
  k=0;
  prtj=1;
  prtpot=1;
  prteig=1;
  proj=1;
  axis='x';
  flux=0;
  norm=0;
  itt=0;
  info=0;
  abstol=1E-12;
  iil=1;
  jactrans=1;
  spl=1;
  strcpy(filep,"veff_0001.dat");
  //-------------------------------------------------//

  printf("      FLUXV       \n\n");
  printf("Vinicius V. Cruz\n\n");
  printf("Reading input parameters...\n\n");
  //rdinput
   while(fscanf(arq,"%s", &workk)!=EOF){
     if(strcmp(workk,"Dimension:")==0) fscanf(arq,"%d", &dim);
     if(strcmp(workk,"xnpoints:")==0) fscanf(arq,"%d", &nx);
     if(strcmp(workk,"xrange:")==0) fscanf(arq,"%lf %lf", &xi, &xf);
     if(dim==1){
       d=0;
     }
     else{
       if(strcmp(workk,"ynpoints:")==0) fscanf(arq,"%d", &ny);
       if(strcmp(workk,"yrange:")==0) fscanf(arq,"%lf %lf", &yi, &yf);  
     }
     if(strcmp(workk,"filename:")==0) fscanf(arq, "%s", &file);
     if(strcmp(workk,"Dynamics")==0) dyn=0;
     if(strcmp(workk,"Single")==0) dyn=1;
     if(strcmp(workk,"nfiles:")==0 && dyn==0) fscanf(arq, "%d %d %d", &nf,&stfil, &endfil);
     if(strcmp(workk,"timeinterval:")==0 && dyn==0) fscanf(arq, "%lf %lf", &ti, &tf);
     if(strcmp(workk,"Mass:")==0 && dim==1) fscanf(arq, "%lf", &m1);
     if(strcmp(workk,"Mass:")==0 && dim==2) fscanf(arq, "%lf %lf", &m1, &m2);
     if(strcmp(workk,"Fluxline:")==0){
       fscanf(arq, "%lf %lf", &fluxline[0], &fluxline[1]);
       if(dim==2){
	 fscanf(arq,"%lf %lf %c", &fluxline[2], &fluxline[3],&axis);
       }
    }
     if(strcmp(workk,"PrintJ")==0) prtj=0;
     if(strcmp(workk,"PrintPot")==0) prtpot=0;
     if(strcmp(workk,"PrintEig")==0) prteig=0;
     if(strcmp(workk,"Spline")==0){
       spl=0;
       fscanf(arq,"%s %d %s %d",&workk,&NXG,&workk,&NYG);
     }

     if(strcmp(workk,"Jacobi")==0){
       jactrans=0;
       fscanf(arq,"%s", &workk);
       if(strcmp(workk,"Masses:")==0) fscanf(arq,"%lf %lf %lf",&ma,&mb,&mc);
       else if(strcmp(workk,"Range:")==0) fscanf(arq,"%lf %lf %lf %lf",&Ri,&Rf,&ri,&rf);
       else printf("wrong input for Jacobi section!");
       fscanf(arq,"%s", &workk);
       if(strcmp(workk,"Range:")==0) fscanf(arq,"%lf %lf %lf %lf",&Ri,&Rf,&ri,&rf);
       else if(strcmp(workk,"Masses:")==0) fscanf(arq,"%lf %lf %lf a.u.",&ma,&mb,&mc);
       else printf("wrong input for Jacobi section!");
     }

     if(strcmp(workk,"Projections")==0){
       fscanf(arq,"%d %s", &nst, &potf);
       poten=fopen(potf,"r");
       eigenst=malloc(nst*sizeof(double));
       if(axis=='x' && dim==2){
	 pot1d=malloc(ny*sizeof(double));
	 for(i=0;i<ny;i++){
	   eigenst[i]=malloc(ny*sizeof(double));
	 }
	 npp=ny;
       }else if(axis=='y' && dim==2){
	 for(i=0;i<nx;i++){
	   eigenst[i]=malloc(nx*sizeof(double));
	   pot1d=malloc(nx*sizeof(double));
	 }
	 npp=nx;
       }
       proj=0;
     }  
   }

   //Mass conversion, steps calculation and flux points determination.

   if( jactrans==0 && spl==1 ){
     //it is strictly necessary to perform the spline procedure in order to transform jacobi coordinates transformation. 
     //If the user does not enter spline parameters, the program will use the same grid
     spl=0;
     NXG=nx;
     NYG=ny;  
   }

   m1=(m1*1.0E-3)/(NAV * Me);
   stepx=(xf -xi)/(nx -1);
   if(dyn==0) stept=(tf - ti)/(nf-1);
   
   fp[0]=((fluxline[0] - xi)/stepx) + 1;
   fp[1]=((fluxline[1] - xi)/stepx) + 1;                          
   if(dim==2){
     stepy=(yf -yi)/(ny -1);
     fp[2]=((fluxline[2] - yi)/stepy) + 1;
     fp[3]=((fluxline[3] - yi)/stepy) + 1; 
    m2=(m2*1.0E-3)/(NAV * Me);
   }

   //spline flux lines and steps
   if(spl==0 && jactrans==1){
     stepxspl=(xf -xi)/(NXG -1);
     fpspl[0]=((fluxline[0] - xi)/stepxspl) + 1;
     fpspl[1]=((fluxline[1] - xi)/stepxspl) + 1;                          
     if(dim==2){
       stepyspl=(yf -yi)/(NYG -1);
       fpspl[2]=((fluxline[2] - yi)/stepyspl) + 1;
       fpspl[3]=((fluxline[3] - yi)/stepyspl) + 1; 
     }
   }

   if(jactrans==0){
     ma=(ma*1.0E-3)/(NAV * Me);
     mb=(mb*1.0E-3)/(NAV * Me);
     mc=(mc*1.0E-3)/(NAV * Me);
      stepxspl=(Rf -Ri)/(NXG -1);
     fpspl[0]=((fluxline[0] - Ri)/stepxspl) + 1;
     fpspl[1]=((fluxline[1] - Ri)/stepxspl) + 1;                          
     stepyspl=(rf -ri)/(NYG -1);
     fpspl[2]=((fluxline[2] - ri)/stepyspl) + 1;
     fpspl[3]=((fluxline[3] - ri)/stepyspl) + 1; 
   }
   //----spline.-------

   //-------------------------------------
   

     printf("INPUT SETTINGS: \n");
     printf("Dimension %d npoints x %d, npoints y %d \n", dim, nx, ny);
     printf(" %lf %lf %lf %lf \n", xi, xf, yi, yf);
     printf("step X %lf , step Y %lf Bohr\n",stepx, stepy);
     printf("dyn %d , %s\n", dyn,file);
     if(dyn==0) printf("ti %lf tf %lf step %lf fs nfile: %d \n",ti,tf,stept,nf);
     printf("step default test %lf \n",stept);
     printf("masses m1 = %lf a.u., m2=%lf a.u.\n", m1, m2);
     printf(" fluxline: x: %lf %lf \n", xi+fp[0]*stepx, xi+fp[1]*stepx);
     if(dim==2) printf("           y: %lf %lf \n", yi+fp[2]*stepy, yi+fp[3]*stepy);
     if(proj==0) printf("%d Flux Projections will be computed\n", nst);

     if(proj==0 && jactrans==0){
       printf("Jacobi transformation: \nWill be performed a Jacobi coordinates transformation from the Alpha arrangement (AB + C) to the Beta arrangement (A + BC).\n");
       printf("ma: %lf mb: %lf mc: %lf \n",ma,mb,mc);
       printf("Transformed range: R %lf %lf, r %lf %lf \n",Ri,Rf,ri,rf);
     }

     if(proj==0 && spl==0){
       printf("Flux Projections will be computed applying a cubic spline interpolation. \n");
       printf("spline parameters: \n");
       printf("npoints x %d, npoints y %d \n", NXG, NYG);
       printf("step X %lf , step Y %lf Bohr\n",stepxspl, stepyspl);
       printf(" fluxline: x: %lf %lf \n", xi+fpspl[0]*stepxspl, xi+fpspl[1]*stepxspl);
       if(dim==2) printf("           y: %lf %lf \n", yi+fpspl[2]*stepyspl, yi+fpspl[3]*stepyspl);
     }
     printf("Printing options: \n");
     if(prtj==0)    printf("Print Current Density Vector J.\n");
     if(prtpot==0)  printf("Print Potential. \n");
     if(prteig==0)  printf("Print Eigenvectors. \n");
     if(prtj==1 && prtpot==1 && prteig==1)  printf("None of the printing options are active.\n");
     
//------------------------------------------------------- 
  //allocating vectors
     Jx=  malloc(nx*ny*sizeof(double));
     Jx2= malloc(nx*sizeof(double));
     Jy2= malloc(nx*sizeof(double));
     Re=  malloc(LDZ*LDZ*sizeof(double));
     Im=  malloc(LDZ*LDZ*sizeof(double));
     Re2= malloc(nx*sizeof(double));
     Im2= malloc(nx*sizeof(double));

     if(jactrans==0){
       R=malloc(LDZ*sizeof(double));
       r=malloc(LDZ*sizeof(double));
     }

     if(spl==0){
       X=malloc(LDZ*sizeof(double));
       Y=malloc(LDZ*sizeof(double));
       XKNOT=malloc(LDZ*sizeof(double));
       YKNOT=malloc(LDZ*sizeof(double));
       BCOEFRe= malloc(LDZ*LDZ*sizeof(double));
       BCOEFIm= malloc(LDZ*LDZ*sizeof(double));
       BCOEFP = malloc(LDZ*LDZ*sizeof(double));
       Z=malloc(LDZ*LDZ*sizeof(double));
     }

     if(dim==2){
       for(i=0;i<nx;i++){
	 Jx2[i]= malloc(ny*sizeof(double));
	 Jy2[i]= malloc(ny*sizeof(double));
	 Re2[i]= malloc(ny*sizeof(double));
	 Im2[i]= malloc(ny*sizeof(double));
	 //projmat[i]= malloc(ny*sizeof(double));
       }
     }
     //--------------------------------------------------
      
   if(proj==0 && dim==2){
     //reading potential----------------------------------

     if(spl==1){
       fscanf(poten,"%s %s %s", &workk, &workk, &workk);
       for(i=0;i<nx;i++){
	 for(j=0;j<ny;j++){
	   fscanf(poten,"%s %s %lf", &workk, &workk, &Jx2[i][j]);
	   //fprintf(out,"%E %E %E \n", xi + i*stepx, yi + j*stepy,Jx2[i][j]);
	 }
	 //fprintf(out," \n");
       }
       //printf(">>>> toatsy 2 \n\n");
       
       if(axis=='x'){
	 j=0;
	 for(i=fp[2];i<fp[3];i++){
	   pot1d[j]=Jx2[fp[0]][i];
	   //fprintf(out,"%lf %lf \n",yi+i*stepy,Jx2[fp[0]][i]);
	   j=j+1;
	 }
       }else if(axis=='y'){
	 //not yet
       }
     //----< calculating eigenvectors
     shm=( 1.0/(24*m2*stepy*stepy) );
     mxdct=ny*nx;
     nptu=fp[3]-fp[2]; //nptu=fp[3]-fp[2]+1;
     wk1=malloc(mxdct*sizeof(double));
     wk2=malloc(mxdct*sizeof(double));
     gst=malloc(20*mxdct*sizeof(double));
     iwork=malloc(mxdct*sizeof(int));
     strcpy(dimmm,".1D");
     //printf(">>>>%s \n", dimmm);
     np=ny;
     iil=1;
     mtrxdiag_(dimmm,&iil,&nst,&info,&mxdct,&nptu,&abstol,iwork,&np,eigvl,&shm,pot1d,wk1,wk2,gst);

     } else if(spl==0){
 
       for(i=0;i<nx;i++){
	 X[i]=xi+i*stepx;
       }
       for(j=0;j<ny;j++){
	 Y[j]=yi+j*stepy;
       }
       dbsnak_ (&nx, X, &KXORD, XKNOT);
       dbsnak_ (&ny, Y, &KYORD, YKNOT);
       //reading potential-------
       readpot_(Z,filep,&nx,&ny,&LDZ);

       /*
       check=fopen("checkpotential.dat","w");
       k=0;
       for(i=0;i<nx;i++){
	 for(j=0;j<ny;j++){
	   Jx2[i][j]=Z[k];
	   fprintf(check,"%E %E %E \n",X[i],Y[j],Jx2[i][j]);
	   k=k+1;
	 }
	 fprintf(check,"\n");
	 }
       */

       //splines coeficients
       dbs2in_ (&ny, Y, &nx, X, Z, &LDZ , &KYORD, &KXORD, YKNOT, XKNOT, BCOEFP);
      
      

       if(axis=='x'){
	 j=0;	
	 for(i=fpspl[2];i<fpspl[3];i++){
	   if(jactrans==0){
	     R[0]=Ri + fpspl[0]*stepxspl;
	     r[i]=ri+i*stepyspl;
	     Jacobitrans(&Xspl[0], &Yspl[0],ma,mb,mc,&R[0],&r[i],-1);
	     printf("%d, X %lf, Y %lf \n",i,Xspl[0],Yspl[0]);
	   }else if(jactrans==1){
	     Xspl[0]=xi+fpspl[0]*stepxspl;
	     Yspl[0]=yi+i*stepyspl;
	   }

	   //pot1d[j]=dbs2vl_ (&Xspl,&Yspl,&KXORD,&KYORD,XKNOT,YKNOT,&nx,&ny,BCOEFP);
	   pot1d[j]=dbs2vl_ (&Yspl[0],&Xspl[0],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFP);
	   //splvalue_ (&Y[i],&X[fpspl[0]],&KXORD,&KYORD,YKNOT,XKNOT,&ny,&nx,BCOEFP,&pot1d[j]);
	   //fprintf(out,"%lf %lf \n",yi+i*stepy,Jx2[fp[0]][i]);
	   j=j+1;
	 }
       }else if(axis=='y'){
	 //not yet
       }

       printf("TOASTY!!! 3\n"); 
      //----< calculating eigenvectors
     shm=( 1.0/(24*m2*stepyspl*stepyspl) );
     mxdct=NYG*NXG;
     nptu=fpspl[3]-fpspl[2]; //nptu=fp[3]-fp[2]+1;
     wk1=malloc(mxdct*sizeof(double));
     wk2=malloc(mxdct*sizeof(double));
     gst=malloc(20*mxdct*sizeof(double));
     iwork=malloc(mxdct*sizeof(int));
     strcpy(dimmm,".1D");
     //printf(">>>>%s \n", dimmm);
     np=NYG;
     iil=1;
     mtrxdiag_(dimmm,&iil,&nst,&info,&mxdct,&nptu,&abstol,iwork,&np,eigvl,&shm,pot1d,wk1,wk2,gst);



     }
     //-----END SPLINE IF--------------------------
     
     printf("\n> Eigenvalues/Hartree \n");
     for(i=0;i<nst;i++){
       printf("%d %E \n", i, eigvl[i]);
     }



     // for(i=0;i<20*mxdct;i++){
       //fprintf(out,"%E \n", gst[i] );
       //}
     if(prtpot==0 && prteig==0){
       eig=fopen("eigv_pot.dat","w");
     } else if(prtpot==0 && prteig==1){
       eig=fopen("potential.dat","w");
     }else if(prtpot==1 && prteig==0){
       eig=fopen("eigenvectors.dat","w");
     }
     //--------testtt
       for(i=0;i<nptu;i++){
	 if(spl==1)fprintf(eig,"%E ", yi + fp[2]*stepy + i*stepy);
	 if(spl==0)fprintf(eig,"%E ", yi + fpspl[2]*stepyspl + i*stepyspl);
	 for(l=0;l<nst;l++){
	   eigenst[l][i]=gst[l*mxdct + i];
	    if(prteig==0) fprintf(eig,"%E ", eigenst[l][i]);
	 }
	 if(prtpot==0){
	   fprintf(eig,"%E \n", pot1d[i]);
	 } else if(prteig==0) fprintf(eig," \n");
       }

       //normalization of eigenvectors------------------------
       if(spl==1){
	 if(proj==0 && axis=='x') {
	   for(i=0;i<nst;i++) real_normalize(eigenst[i],nptu,stepy);
	   stepp=stepy;
	 }else if(proj==0 && axis=='y'){
	   for(i=0;i<nst;i++) real_normalize(eigenst[i],nptu,stepx);
	   stepp=stepx; 
	 }
       }else if(spl==0){
	 if(proj==0 && axis=='x') {
	   for(i=0;i<nst;i++) real_normalize(eigenst[i],nptu,stepyspl);
	   stepp=stepyspl;
	 }else if(proj==0 && axis=='y'){
	   for(i=0;i<nst;i++) real_normalize(eigenst[i],nptu,stepxspl);
	   stepp=stepxspl; 
	 }
	 printf("Eigenvectors are normalized \n");
       }
       //--------------------------------------------------
       
       
       //close(eig);
       //printf(">>>> toatsy \n\n");
     free(gst);
     free(wk1);
     free(wk2);
     free(iwork);
 
     //----END 2D 
   }
   else if(proj==0 && dim==1){
     pot1d=malloc(nx*sizeof(double));
     fscanf(poten,"%s %s",&workk,&workk);
     for(i=0;i<nst;i++){
       eigenst[i]=malloc(nx*sizeof(double));
	 }
     //pot1d=malloc(nx*sizeof(double));
     for(i=0;i<nx;i++){
       fscanf(poten,"%s %lf",&workk,&pot1d[i]);
       //printf("%s %lf \n",workk,pot1d[i]);
     }
     
     
     //calculating eigenvectors
     shm=( 1.0/(24*m1*stepx*stepx) );
     mxdct=nx*nx;
     nptu=nx;
     wk1=malloc(mxdct*sizeof(double));
     wk2=malloc(mxdct*sizeof(double));
     gst=malloc(20*mxdct*sizeof(double));
     iwork=malloc(mxdct*sizeof(int));
     strcpy(dimmm,".1D");
     //printf(">>>>%s \n", dimmm);
     np=nx;
     iil=1;
     mtrxdiag_(dimmm,&iil,&nst,&info,&mxdct,&nptu,&abstol,iwork,&np,eigvl,&shm,pot1d,wk1,wk2,gst);
   
     
     printf("\n> Eigenvalues/Hartree \n");
     for(i=0;i<nst;i++){
       printf("%d %E \n", i, eigvl[i]);
     }

     
      if(prtpot==0 && prteig==0){
       eig=fopen("eigv_pot.dat","w");
     } else if(prtpot==0 && prteig==1){
       eig=fopen("potential.dat","w");
     }else if(prtpot==1 && prteig==0){
       eig=fopen("eigenvectors.dat","w");
     }

      //
      for(i=0;i<nx;i++){
	fprintf(eig,"%E ", xi + i*stepx);
	for(l=0;l<nst;l++){
	  eigenst[l][i]=gst[l*mxdct + i];
	  if(prteig==0) fprintf(eig,"%E ", eigenst[l][i]);
	}
	if(prtpot==0){
	  fprintf(eig,"%E \n", pot1d[i]);
	} else if(prteig==0) fprintf(eig," \n");
      }
     
       //printf(">>>> toatsy \n\n");
     free(gst);
     free(wk1);
     free(wk2);
     free(iwork); 
     
     for(i=0;i<nst;i++) real_normalize(eigenst[i],nx,stepx);
     //---
   }
   /*
   eigg=fopen("HCleigen0","r");
   fprintf(out,"\n");
   printf("toasty");
   for(l=0;l<nst;l++){
     for(i=0;i<ny;i++){
       fscanf(eig,"%s %lf",&workk, &eigenst[l][i]);
       fprintf(out,"%lf %lf \n", yi+i*stepy,eigenst[l][i]);
     }
   }*/
   //normalization of eigenstates
   //------------------------------------

   if(proj==0){
     for(l=0;l<nst;l++){
       Pii[l]=0.0E+0;
     }
   }


  if(dyn==0){    
//----------generating wavepacket files's names------
    for(f=stfil;f<endfil+1;f++){
      strcpy(filen,file);
      strcpy(outt,"Jvect");
      if(f<10){
	sprintf(num,"000%d.dat",f);
	strcat(filen,num);
      }
      else if( f>=10 && f<100){
	sprintf(num,"00%d.dat",f);
	strcat(filen,num);
      }
      else if(f>=100 && f<1000){
	sprintf(num,"0%d.dat",f);
	strcat(filen,num);
      }
      else if(f>=1000){
	sprintf(num,"%d.dat",f);
	strcat(filen,num);
      }
      sprintf(numout,"%d",f);
      strcat(outt,numout);
      // printf("%s %s \n", file, filen);
      
      wp=fopen(filen,"r");

      if(wp==NULL){
	printf("## Error! could not open file %d %s \n", f, filen);
	return;
      }
//-------------------------------------------------

      fscanf(wp,"%s %s %s %s", &workk, &workk, &workk, &workk);

      if(dim==1){ //-----1D------------------------------------------------------
      for(i=0;i<nx;i++){
	fscanf(wp,"%lf %lf %lf", &wor, &Re[i], &Im[i]);
      //fprintf(out,"%E %E %E \n", wor, Re[i][0], Im[i][0]);
	norm=norm + ( (Re[i]*Re[i]) + (Im[i]*Im[i]) )*stepx;
      }
      //printf("\n >>>> norm: %E \n", norm);
       //normalization-----------------
      for(i=0;i<nx;i++){
	Re[i]=Re[i]/sqrt(norm);
	Im[i]=Im[i]/sqrt(norm);
      }
      norm=0;
      //--------------------------------

      prob_current1D(Jx,Re,Im,nx,stepx,dyn,m1);
      //print flux vector J
      if(prtj==0){
	jout=fopen(outt,"w");
	for(i=1;i<nx-1;i++){
	  fprintf(jout,"%E %E \n", xi+i*stepx, Jx[i]);
	}
      }

//-----------------------------PROJECTIONS 1D------------

      if(proj==0){
	fprintf(projfil,"%lf ",ti+itt*stept);
	for(l=0;l<nst;l++){
	  dotp=project(eigenst[l],Re,Im,nx,stepx);
	  fprintf(projfil,"%lf ",dotp);
	}
	fprintf(projfil,"\n");
      }
//-----------------------------------------------------

//-------------------Flux through point--------------------------------
      flux=Jx[fp[0]];
	//atomic unit of time  1 ta = 2.41889E-2 fs
      prob= prob + flux*stept/0.024189;
      fprintf(fout,"%E %E %E \n", ti+itt*stept, flux, prob);
      //printf("%E %E %E \n",ti+itt*stept, flux, prob);
//----------------------------------------------------------------
      itt=itt+1;
      } //-----END--1D---------------------------------------------

      if(dim==2){//---2D---------------------------------------------
	printf("\n %s \n", filen);
	read_(Re,Im,filen,&nx,&ny);

	//transcripting----------------------------------------------
	k=0;
	for(i=0;i<nx;i++){
	  for(j=0;j<ny;j++){
	    //printf("%E %E \n",Re[j],Im[j]);
	    Re2[i][j]=Re[k];
	    Im2[i][j]=Im[k];
	    k=k+1;
	  }
	}
	
	//--------------------------------------------
	//printf("toatsy>>> \n");
	if(f==stfil){
	  //----------------------
	  for(i=0;i<nx;i++){
	    for(j=0;j<ny;j++){
	      //fprintf(out,"%E %E %E %E\n",(xi + i*stepx),(yi+j*stepy),Re2[i][j], Im2[i][j]);
	      //fprintf(out,"%E %E %E \n",(xi + i*stepx),(yi+j*stepy),(Re2[i][j] * Re2[i][j]) + (Im2[i][j] * Im2[i][j]));
	      norm=norm+ ( (Re2[i][j] * Re2[i][j]) + (Im2[i][j] * Im2[i][j]) ) * stepx * stepy;
	    }
	    //fprintf(out,"\n");
	  }
	  //--------------
	}

	printf(" norm: %E \n", norm);
	//normalization-----------------
	for(i=0;i<nx;i++){
	  for(j=0;j<ny;j++){
	    Re2[i][j]=Re2[i][j]/sqrt(norm);
	    Im2[i][j]=Im2[i][j]/sqrt(norm);
	    norm2=norm2+ ( (Re2[i][j] * Re2[i][j]) + (Im2[i][j] * Im2[i][j]) ) * stepx * stepy;
	  }
	}
	printf(" normalized: %E \n", norm2);
	//norm=0;
	norm2=0;

	//-------------------Flux Vector Calc.--------------------------------
	prob_current2D(Jx2,Jy2,Re2,Im2,nx,ny,stepx,stepy,m1,m2);
	//print flux vector J -------------------
	if(prtj==0){
	  jout=fopen(outt,"w");
	  for(i=1;i<nx-1;i++){
	    for(j=1;j<ny-1;j++){
	      fprintf(jout,"%E %E %E %E \n", (xi + i*stepx),(yi+j*stepy),Jx2[i][j],Jy2[i][j]);
	    }
	  }
	}
	//--------------------------------------
	//-------------------Flux through curve--------------------------------//

	if(axis=='x'){
	  for(i=fp[2];i<fp[3];i++){
	    flux = flux + Jx2[fp[0]][i]*stepy;
	    //fprintf(out,"%E %E %E \n", yi + i *stepy ,Jx2[fp[0]][i],flux);
	  }
	  prob=prob + flux*stept/0.024189;
	  fprintf(fout,"%E %E %E \n", ti + itt*stept, flux, prob);
	  printf("\n f: %d time: %lf %E %E \n",f,ti+itt*stept, flux, prob);
	  //flux=0;
	}
	/*else if(axis=='y'){
	  for(i=fp[0];i<fp[1];i++){
	  flux = flux + Jy2[i][fp[2]]*stepx;
	  }
	  prob=prob + flux*stept/0.024189;
	  fprintf(fout,"%E %E %E \n", (f-1)*stept, flux, prob);
	  printf(" \n f: %d time: %lf",f,(f-1)*stept);
	  }*/

	//---------------------------------------------------------------------
	//------------Projections----------------------------------------------
	if(proj==0){
	  //check=fopen("checkeigen.dat","w");
	  if(spl==0){
	    readspl_ (Re,Im,filen,&nx,&ny,&LDZ);

	    for(i=0;i<LDZ*LDZ;i++){
	      Re[i]=Re[i]/sqrt(norm);
	      Im[i]=Im[i]/sqrt(norm);
	    }
	    
	    dbs2in_ (&ny, Y, &nx, X, Re, &LDZ , &KYORD, &KXORD, YKNOT, XKNOT, BCOEFRe);
	    dbs2in_ (&ny, Y, &nx, X, Im, &LDZ , &KYORD, &KXORD, YKNOT, XKNOT, BCOEFIm);	    	    
	    
	    /*check=fopen("checkeigen.dat","w");
	    for(i=0;i<nx;i++){
	      Xspl[0]=xi+i*stepx;
	      for(j=0;j<ny;j++){
		Yspl[0]=yi+j*stepy; 
		Respl[0]=dbs2vl_ (&Yspl[0],&Xspl[0],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
		Imspl[0]=dbs2vl_ (&Yspl[0],&Xspl[0],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
		//fprintf(check,"%E %E %E \n",Xspl[0],Yspl[0],Respl[0]*Respl[0] + Imspl[0]*Imspl[0] );
	      }
	      fprintf(check,"\n");
	      }*/
	    fluxproj=0;
	    for(l=0;l<nst;l++){
	      //starting vectors
	      projRe[l]  =  0.0E+0;
	      projIm[l]  =  0.0E+0;
	      projDRe[l] =  0.0E+0;
	      projDIm[l] =  0.0E+0;
	      projet[l]  =  0.0E+0;
	      //--------------
	      i=0; 
	      
	      for(j=fpspl[2];j<fpspl[3];j++){
		if(axis=='x'){

		  if(jactrans==0){
		    printf("not implemented yet \n");
		  }else if(jactrans==1){
		    Xspl[0]= xi + (fpspl[0]-1)*stepxspl;
		    Xspl[1]= xi + (fpspl[0]  )*stepxspl;
		    Xspl[2]= xi + (fpspl[0]+1)*stepxspl;  
		    Yspl[0]= yi + j*stepyspl; 
		  }
		  //if(l==0 && j==fpspl[2]) printf("x0 %E x1 %E x2 %E y %E \n",Xspl[0],Xspl[1],Xspl[2],Yspl[0]);
		 

		  Respl[0]=dbs2vl_ (&Yspl[0],&Xspl[0],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
		  Respl[1]=dbs2vl_ (&Yspl[0],&Xspl[1],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);
		  Respl[2]=dbs2vl_ (&Yspl[0],&Xspl[2],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFRe);

		  Imspl[0]=dbs2vl_ (&Yspl[0],&Xspl[0],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
		  Imspl[1]=dbs2vl_ (&Yspl[0],&Xspl[1],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
		  Imspl[2]=dbs2vl_ (&Yspl[0],&Xspl[2],&KYORD,&KXORD,YKNOT,XKNOT,&ny,&nx,BCOEFIm);
		  
		  ReD=deriv(Respl[0],Respl[2],stepxspl);
		  ImD=deriv(Imspl[0],Imspl[2],stepxspl);
		  //fprintf(check,"%E %E %E \n",Yspl[0],ReD,ImD); 

		  projRe[l]  =  projRe[l]   +  ( eigenst[l][i] * Respl[1] ) * stepyspl;
		  projIm[l]  =  projIm[l]   +  ( eigenst[l][i] * Imspl[1] ) * stepyspl;
		  projDRe[l] =  projDRe[l]  +  ( eigenst[l][i] * ReD      ) * stepyspl;
		  projDIm[l] =  projDIm[l]  +  ( eigenst[l][i] * ImD      ) * stepyspl;
		  //if(l==0 && f==20) fprintf(check,"%E %E \n",Yspl[0],eigenst[l][i]);
		  i=i+1;
		}/* else if(axis=='y'){
		    projRe= projRe + ( eigenst[l][i] * Re2[i][j] ) * stepx;
		    projIm= projIm + ( eigenst[l][i] * Im2[i][j] ) * stepx; 
		    }*/
		//-------------------
	      }
	      projet[l]  =  (1.0/m1)*( projRe[l]*projDIm[l] - projIm[l]*projDRe[l] );
	      // printf("flux%d %E ", l,(1.0/m1)*( projRe[l]*projDIm[l] - projIm[l]*projDRe[l] ));
	      fluxproj=fluxproj+projet[l];
	    }
	  } 
	  else if(spl==1){
	    
	    fluxproj=0;
	    for(l=0;l<nst;l++){
	      //starting vectors
	      projRe[l]  =  0.0E+0;
	      projIm[l]  =  0.0E+0;
	      projDRe[l] =  0.0E+0;
	      projDIm[l] =  0.0E+0;
	      projet[l]  =  0.0E+0;
	      //--------------
	      i=0; 
	      for(j=fp[2];j<fp[3];j++){
		if(axis=='x'){
		  ReD=deriv(Re2[fp[0]-1][j],Re2[fp[0]+1][j],stepx);
		  ImD=deriv(Im2[fp[0]-1][j],Im2[fp[0]+1][j],stepx);
		  projRe[l]  =  projRe[l]   +  ( eigenst[l][i] * Re2[fp[0]][j] ) * stepy;
		  projIm[l]  =  projIm[l]   +  ( eigenst[l][i] * Im2[fp[0]][j] ) * stepy;
		  projDRe[l] =  projDRe[l]  +  ( eigenst[l][i] * ReD           ) * stepy;
		  projDIm[l] =  projDIm[l]  +  ( eigenst[l][i] * ImD           ) * stepy;
		  i=i+1;
		  //fprintf(check,"%E %E %E \n",yi+j*stepy,ReD,ImD);
		}/* else if(axis=='y'){
		    projRe= projRe + ( eigenst[l][i] * Re2[i][j] ) * stepx;
		    projIm= projIm + ( eigenst[l][i] * Im2[i][j] ) * stepx; 
		    }*/
		//-------------------
	      }
	      projet[l]  =  (1.0/m1)*( projRe[l]*projDIm[l] - projIm[l]*projDRe[l] );
	      // printf("flux%d %E ", l,(1.0/m1)*( projRe[l]*projDIm[l] - projIm[l]*projDRe[l] ));
	      fluxproj=fluxproj+projet[l];
	    }
	  }
	  
	  
	  //printf("%E %E Pi: %E \n",projRe,projIm,projet[0]);
	  fprintf(out,"%E ", ti+itt*stept);
	  fprintf(out,"%E ", fluxproj);
	  printf("time: %E ",ti+itt*stept);
	  projsum=0.0E+0;
	  for(l=0;l<nst;l++){
	    Pii[l]=Pii[l]+projet[l]*stept/0.024189;
	    fprintf(out,"%E ", Pii[l]);
	    printf("Pi%d: %E ",l, Pii[l]);
	    projsum=projsum + Pii[l];
	  }
	  printf("Total: %E ",projsum);
	  fprintf(out,"\n");
	  printf("\n");
	  //                                                 
	  //            |\                             / /   
	  //            | \                           / /    
	  //      bow   |  )   >>--arrow-->   knee   (_ )    
	  //            | /                           \ \    
	  //            |/                            _\_\   
	  //                                         (____|  
	  //                                                 
	}
	//---------------------------------------------------------------------
	flux=0;
	itt=itt+1;	
      }//-----END--2D-----------------------------------------------------------//   
     
    }
  }
  //---------end checking files
  
  /* else if(dyn==1){
     wp=fopen(file,"r");
     for(i=0;i<nx;i++){
     fscanf(wp,"%lf %lf %lf", &wor, &Re[i], &Im[i]);
     fprintf(out,"%E %E %E \n", wor, Re[i][0], Im[i][0]);
     }
     prob_current1D(Jx,Re,Im,nx,stepx,dyn,m1);
     printf("%d %lf \n",xi+i*stepx,J[i][1]);
     }*/ 
//-----------------------------------------------------------
  free(Jx);
  free(Jx2);
  free(Jy2);
  free(Re);
  free(Im);
  free(Re2);
  free(Im2);
 
}

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

double prob_current2D(double **Jx, double **Jy, double **Re, double **Im, int nx, int ny, double stepx, double stepy, double m1, double m2){
  // J = hbar/mass(Re(psi* (gradient) psi))
  int i,j;
  double DRex,DImx,work,DRey,DImy;
  FILE *debug;

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

}

double project(double *u, double *Re, double *Im,int n, double step){
  int i,j,k;
  double dotprod,dotr,doti;
  dotprod=0.0;
  dotr=0.0;
  doti=0.0;
  for(i=0;i<n;i++){
    dotr=dotr + u[i]*Re[i]*step;
    doti=doti + u[i]*Im[i]*step;
  }
  dotprod=dotr*dotr + doti*doti;
  return dotprod;
}

double deriv(double y1, double y2, double step)
{
  double dydx;
  dydx= (y2 - y1)/(2*step);
  return dydx;
}

double real_normalize(double *X, int n, double step){
  int i;
  double norm,norm2;
  //printf("inside the matrix \n\n");
  norm=0;
  for(i=0;i<n;i++){
    norm=norm+X[i]*X[i]*step;
  }
  //printf("eigenstate norm %E \n",norm);
  for(i=0;i<n;i++){
    X[i]=X[i]/sqrt(norm);
    norm2=norm2+X[i]*X[i]*step;
  }
  //printf("normalized eigenstate %E \n",norm2);
  norm2=0;
}
//-----------------------------------------------------------------------------------------
