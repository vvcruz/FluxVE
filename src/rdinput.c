#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

int rdinput(FILE *arq,int *nx, int *ny, double *xi, double *xf, double *yi, double *yf,char *file,int *stfil,int *endfil, double *ti, double *tf, double *m1, double *m2, double *fluxline, char *axis, double *nv, int *nf, double *stepxspl, double *stepyspl, int *twopow,double *a, double *x0, double *k0, double *Evib, double *potshift, int *prtRMSwhole, int *type, int *noreg, double *ma, double *mb, double *mc, int *wholegrid, double *stren, double t0){
  int i,j,k,spl;
  int NXG,NYG,jcheck1,jcheck2;
  char workk[50],coment[150],typenam[20];

  //rdinput
  while(fscanf(arq,"%s", workk)!=EOF){
    //MAIN INPUT GROUP
    // printf("3> %s \n", workk);
  check:
    if(strncasecmp(workk,"#Main",5)==0){   
      //printf("inside #Main \n");
      for(i=0;;i++){
	fscanf(arq,"%s", workk);
	while(workk[0]=='!'){
	  fgets(coment,90,arq);
	  fscanf(arq,"%s", workk);
	  if(feof(arq)) break;
	}
	if(feof(arq) || workk[0]=='#') goto check; //break;
	if(strncasecmp(workk,"runtype",7)==0){
	  fscanf(arq, "%s", typenam);
	  if(strncasecmp(typenam,"non-reac",8)==0)*type=0;
	  else if(strncasecmp(typenam,"reac",4)==0)*type=1;
	  else if(strncasecmp(typenam,"td-only",6)==0)*type=2;
	  else{
	    printf("invalid run type entered: %s",typenam);
	    return 666;
	  }
	}
	if(strncasecmp(workk,"wholegrid",9)==0) *wholegrid = 0;
	if(strncasecmp(workk,"noregion",8)==0) *noreg = 0;
	if(strncasecmp(workk,"xrange",6)==0) fscanf(arq,"%lf %lf", xi, xf);
	if(strncasecmp(workk,"xnpoints",8)==0) fscanf(arq,"%d", nx);
	if(strncasecmp(workk,"ynpoints",8)==0) fscanf(arq,"%d", ny);
	if(strncasecmp(workk,"yrange",6)==0) fscanf(arq,"%lf %lf", yi, yf);  	
	if(strncasecmp(workk,"filename",8)==0) fscanf(arq, "%s", file);
      }
    }
    //FLUX GROUP INPUT 
    if(strncasecmp(workk,"#Flux",5)==0){
      // printf("inside #Flux \n");
      for(i=0;;i++){
	fscanf(arq,"%s", workk);
	//coment structure
	while(workk[0]=='!'){
	  fgets(coment,90,arq);
	  fscanf(arq,"%s", workk);
	  if(feof(arq)) break;
	}
	if(feof(arq) || workk[0]=='#') goto check;// break;
	if(strncasecmp(workk,"nfiles",6)==0){
	  fscanf(arq, "%d %d",stfil, endfil);
	  *nf=(*endfil-*stfil) + 1 ;
	}
	if(strncasecmp(workk,"timeinterval",12)==0) fscanf(arq, "%lf %lf", ti, tf);
	if(strncasecmp(workk,"Mass",4)==0) fscanf(arq, "%lf %lf", m1, m2);
	if(strncasecmp(workk,"Fluxline",8)==0){
	  fscanf(arq, "%lf %lf", &fluxline[0], &fluxline[1]);
	  fscanf(arq,"%lf %lf %c", &fluxline[2], &fluxline[3],axis);
	}
	if(strncasecmp(workk,"NormalVector",12)==0)  fscanf(arq, "%lf",nv);
	if(strncasecmp(workk,"Spline",6)==0){
	  spl=0;
	  fscanf(arq,"%s %lf %s %lf",workk,stepxspl,workk,stepxspl);
	}
	if(strncasecmp(workk,"Fourier",7)==0){
	  spl=0;
	  fscanf(arq,"%d",twopow);
	}
	if(strncasecmp(workk,"window",6)==0)fscanf(arq,"%lf %lf",stren,t0);
	if(strncasecmp(workk,"InitWP",6)==0){
	  fscanf(arq,"%lf %lf %lf",a,x0,k0);
	}
	if(strncasecmp(workk,"Evib",4)==0){
	  fscanf(arq,"%lf",Evib);
	}
	if(strncasecmp(workk,"PotShift",8)==0){
	  fscanf(arq,"%lf",potshift);
	}

      }
    }
    /*//PRINT INPUT GROUP
    if(strncasecmp(workk,"#Print",6)==0){
      // printf("inside #Print \n");
      for(i=0;;i++){
	fscanf(arq,"%s", &workk);
	while(workk[0]=='!'){
	  fgets(coment,90,arq);
	  fscanf(arq,"%s", &workk);
	  if(feof(arq)) break;
	}
	if(feof(arq) || workk[0]=='#') goto check;//break;
	if(strncasecmp(workk,"PrintPot",8)==0) prtpot=0;
	if(strncasecmp(workk,"PrintEig",8)==0) prteig=0;
	if(strncasecmp(workk,"PrintPotjac",11)==0)  prtjacpot=0;
	if(strncasecmp(workk,"PrintFile",9)==0)  prtfil=0;
	if(strncasecmp(workk,"Debug",5)==0)  debug=0;
	if(strncasecmp(workk,"DebugF",5)==0)  debugf=0;
      }
    }*/
    //JACOBI INPUT GROUP
    if(strncasecmp(workk,"#Jacobi",7)==0){
      //printf("inside #Jacobi \n");
      for(i=0;;i++){
	fscanf(arq,"%s", &workk);
	//Coment structure
	while(workk[0]=='!'){
	  fgets(coment,90,arq);
	  fscanf(arq,"%s", &workk);
	  if(feof(arq)) break;
	}	
	if(feof(arq) || workk[0]=='#') goto check;// break;
	//jactrans=0;
	/*if(spl==1){ //it is strictly necessary to perform the spline procedure in order to transform jacobi coordinates. 
	  //If the user does not enter spline parameters, the program will use the parameters of the original grid
	  NXG=nx;
	  NYG=ny;
	  }*/
	if(strncasecmp(workk,"Masses",6)==0){
	  fscanf(arq,"%lf %lf %lf",ma,mb,mc);
	  jcheck1=0;
	  jcheck2=0;
	}
	if(strncasecmp(workk,"Fluxline",5)==0){
	  fscanf(arq,"%lf %lf %lf %lf",&fluxline[0],&fluxline[1],&fluxline[2],&fluxline[3]);
	  jcheck2=0;
	}
	if(jcheck1==1 && jcheck2==1){
	  printf("wrong input for Jacobi section! \n");
	  return 65;
	}
	if(strncasecmp(workk,"Spline",6)==0){
	  spl=0;
	  fscanf(arq,"%s %d %s %d",&workk,&NXG,&workk,&NYG);
	}
	if(strncasecmp(workk,"Potfile",7)==0){
	  //fscanf(arq,"%s",&potf);
	  //poten=fopen(potf,"r");
	}
      }
    }
    /*//PROJECTIONS INPUT GROUP
    if(strncasecmp(workk,"#Projections",12)==0){
      proj=0;
      //printf("inside #Projections \n");
      for(i=0;;i++){
	fscanf(arq,"%s", &workk);
	//Coment structure
	while(workk[0]=='!'){
	  fgets(coment,90,arq);
	  fscanf(arq,"%s", &workk);
	  if(feof(arq)) break;
	}	
	if(feof(arq) || workk[0]=='#') goto check;// break;
	if(strncasecmp(workk,"Potfile",7)==0){
	  fscanf(arq,"%s",&potf);
	  poten=fopen(potf,"r");
	}
	if(strncasecmp(workk,"Nstates",7)==0) fscanf(arq,"%d",&nst);
	if(strncasecmp(workk,"Maxstates",9)==0){
	  fscanf(arq,"%lf",&pottool);
	  maxst=0;
	}
	if(strncasecmp(workk,"Printstates",11)==0){
	  fscanf(arq,"%d",&prtst);
	  prcheck=0;
	}
	if(strncasecmp(workk,"Vibration",9)==0)	vib=0;
      }
      }*/
  }
  //-----END OF INPUT
 }
