//---------TRASH

  //---- generating coefficients (fourier transform of initial spatial distribution)-------------------
  //stepxspl = (2*xf)/NTG;
  //printf(">> %E %E %E %E %E\n",xf,stepxspl,a,x0,k0);
  //norm = 0.0e+0;
  //x0=0;
  //k0=0;
  /*
  for(i=0;i<NTG;i++){
    x = -xf + i*stepxspl;
    norm = norm + pow(gausswp(x,x0,a,1.0e+0),2)*stepxspl;
  }
  norm = 1.0/sqrt(norm);
  printf("norm1 = %E\n",norm);
  */

  /*
  pki = 0.0e+0;
  steppk = 30.0/(NTG - 1);
  for(i=0;i<NTG;i++){
    pk = pki + i*steppk;
    //divide by twopi to keep the transformed wavefunction normalized
    aE[i][0] = gausswp_k(pk,k0,a);
    aE[i][1] = 0.0e+0;
    fprintf(g_k,"%E %E %E %E %E\n",pk, aE[i][0],aE[i][1],pow(aE[i][0],2) + pow(aE[i][1],2),pk*pk/(2*m1));
  }
  fclose(g_r);fclose(g_k);
  */

  /* this part is gonna be replaced by the analytical expression for the momentum distribution

  norm = sqrt(sqrt(2*log(2)/(a*a*M_PI))); // gaussian normalization factor
  //printf("norm1 = %E\n",norm);
  for(i=0;i<NTG;i++){
    x = -xf + i*stepxspl;
    workk[i][0] =gausswp(x,x0,a,norm)*cos(k0*x);//gausswp(x,x0,a,norm)*cos(k0*x); exp(-0.02*x*x)*cos(k0*x);
    workk[i][1] = gausswp(x,x0,a,norm)*sin(k0*x);//gausswp(x,x0,a,norm)*sin(k0*x);//0.0e+0;-exp(-0.02*x*x)*sin(k0*x);
    fprintf(g_r,"%E %E %E %E\n",x,workk[i][0],workk[i][1],pow(workk[i][0],2) + pow(workk[i][1],2));
  }
  center_fft(workk,NTG);
  p = fftw_plan_dft_1d(NTG,workk,aE,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(p);
  free(workk);
  center_fft(aE,NTG);
  steppk = 2*M_PI/(NTG*stepxspl);
  pki = -NTG*steppk/(2.0E+0);
  printf(">> %E %E \n",pki,steppk);
  for(i=0;i<NTG;i++){
    pk = pki + i*steppk;
    //divide by twopi to keep the transformed wavefunction normalized
    aE[i][0] = aE[i][0]*stepxspl/(sqrt(2.0*M_PI));
    aE[i][1] = aE[i][1]*stepxspl/(sqrt(2.0*M_PI));
    //ansk = sqrt(sqrt( 2*M_PI*a*a/M_LN2))*exp( -pow(a*(pk-k0),2)/(4*M_LN2) );//sqrt(2.0*M_PI); //a*sqrt(M_PI/log(2)) //(1.0/a)*sqrt(sqrt(2)*log(2)/M_PI)
    fprintf(g_k,"%E %E %E %E %E\n",pk,aE[i][0],aE[i][1],pow(aE[i][0],2) + pow(aE[i][1],2),pk*pk/(2*m1));
  }
  fclose(g_r);fclose(g_k);
  //---------------------------------------------------------------------------------
  */

 /*/TEST FFT ----------------------------------------------
  NTG=1024;
  tf = 20.1*M_PI;
  steptspl = (2*tf)/(NTG);
  stepE = 2*M_PI/(NTG*steptspl);
  //printf(">> Energy step: %E \n",stepE);
  Ei = -NTG*stepE/(2.0E+0);

//------- generate test function
  gendata_function(workin,tf,steptspl,NTG,width); 

  for(i=0;i<NTG;i++) fprintf(func,"%E %E\n",-tf + i*steptspl,workin[i][0]);
  fprintf(func,"\n");
  for(i=0;i<NTG;i++) fprintf(func,"%E %E\n",-tf + i*steptspl,workin[i][1]);
  
  //calling this routine here centers t=0 at the zero frequency position of the fft input vector
  center_fft(workin,NTG);

  //------- do forward fft
  p = fftw_plan_dft_1d(NTG,workin,workout,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(p);
  
  //calling this routine here shifts the fft result back to the center of the array
  center_fft(workout,NTG);

  for(i=0;i<NTG;i++) fprintf(trans,"%E %E\n",Ei + i*stepE,steptspl*workout[i][0]);
  fprintf(trans,"\n");
  for(i=0;i<NTG;i++) fprintf(trans,"%E %E\n",Ei + i*stepE,steptspl*workout[i][1]);
  fprintf(trans,"\n");
  /*
  for(i=0;i<NTG;i++) fprintf(trans,"%E %E\n",Ei + i*stepE,cos(tf/2*(Ei + i*stepE))*exp(-pow((Ei - 5.0 + i*stepE),2)/(4*M_PI)));
  fprintf(trans,"\n");
  for(i=0;i<NTG;i++) fprintf(trans,"%E %E\n",Ei + i*stepE,-sin(tf/2*(Ei + i*stepE))*exp(-pow((Ei -5.0 + i*stepE),2)/(4*M_PI)));
  * /

  //------- do backward fft
  center_fft(workout,NTG);
  p = fftw_plan_dft_1d(NTG,workout,workk,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(p);
  center_fft(workk,NTG);

  for(i=0;i<NTG;i++) fprintf(itrans,"%E %E\n",-tf + i*steptspl,workk[i][0]/(NTG*1.0e+0));
  fprintf(itrans,"\n");
  for(i=0;i<NTG;i++) fprintf(itrans,"%E %E\n",-tf + i*steptspl,workk[i][1]/(NTG*1.0e+0));
  //END OF TEST FFT ------------------------------------------------- */


//return 1;

/*/------------ UTRA DEBUG

  for(l=0;l<30;l++){

    //k0 = 18.000 + l*(1.50000/29.000);
    potshift = 0.0000 + l*(0.008/29.000);
    fprintf(deb,"# %.11E\n",potshift);
    for(k=0;k<nE;k++){
      if(E[k] >= Emin && E[k] <= Emax){
	pk = sqrt(2.0*m1*(E[k] - Evib - potshift));
	aE[k][0] = gausswp_k(pk,k0,a);
	// finding maxima ---
	if((m1/pk)*aE[k][0]*aE[k][0] >= maxaE[1]){
	  maxaE[1] = (m1/pk)*aE[k][0]*aE[k][0];
	  maxaE[0] = E[k];
	}
	//-------
	fprintf(deb,"%E %E \n",E[k],(m1/pk)*aE[k][0]*aE[k][0]);
	//-----
      }
    }
    fprintf(deb,"\n");

    pk = sqrt(2.0*m1*(maxF[0] - Evib - potshift));
    printf("maximum Flux: E = %.15E, k = %E, F = %.15E \n",maxF[0],pk,maxF[1]);
    pk = sqrt(2.0*m1*(maxaE[0] - Evib - potshift));
    printf("maximum coeff: E = %.15E, k = %E, F = %.15E \n",maxaE[0],pk,maxaE[1]);
    pk = 0.5e+0*( k0 + sqrt(k0*k0 - 4*M_LN2/(a*a)) );
    maxaE[0] = pk*pk/(2.0*m1) + Evib + potshift;
    maxaE[1] = (m1/pk)*gausswp_k(pk,k0,a)*gausswp_k(pk,k0,a);
    printf("analytical maximum for coeff: E = %.15E, k = %E, F = %.15E \n",maxaE[0],pk,maxaE[1]);
  }
  //---------------------*/
