/* Interpret the input data.*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alloc.h"
#include "defmol.h"

#define ZERO          0.0E+0

extern char    *toc;
extern char    *title;
extern char    *daltonpath;
extern char    *gamesspath;
extern char    *potential;
extern char    *print_pop;
extern char    *rmoutfile;
extern char    *rmdatfile;
extern char    *optstruct;
extern char    *loc;
extern char    *psdpot;
extern char    *bs;
extern char    eeecp;
extern int     natoms;
extern int     nindiv;
extern int     noas;
extern int     mltpl;
extern int     icharg;
extern int     GA;
extern int     RS;
extern int     absenrgy;
extern int     nstep;
extern int     percmut;
extern int     mutmin;
extern int     mutmax;
extern int     typmat;
extern int     ngen;
extern int     *conv;
extern int     *genafconv;
extern int     nrotgp;
extern int     typga;
extern double  bohr2angstron;
extern double  *energy;
extern long    nseed;
extern defmol  *PMolC;
extern defmol  *MPop;
extern Txyz    *XYZ;
extern mzmat   *cnct;
extern mzmat   *zmat;
extern Trtt    *rtt;
extern Tatom   *rttatm;

int RdInput(const char *Infile, const int nchar)
{
  FILE     *fl0, *fl1;

  int      ng = 0, nsg = 0, nmat = 0, eof = 1, lock = 0, nn = 0, nos=1;
  int      i, j, k, t, length, kai, kaf, status = 0, nga = 0;
  int      act, ang, die, natmrot, center, plane, atm, count;
  double   distfactor = 1.0E+0;

  char     fname[11], rdaux[nchar];

  defmol   *pPMolC;
  Txyz     *pxyz;

  /* Check and open the input file */
  if ( ! (fl0=fopen(Infile,"r")) )
    {
      printf("|\n\nError [RdInput]: The file \"%s\" is not available!\n",Infile);
      exit(EXIT_FAILURE);
    }

  fscanf(fl0, "%s", &rdaux);
  printf("#");

  while ( !feof(fl0) && eof == 1 )
    {
      /* Uncomment bellow to debug */
      /*printf ("\n\n <<<>>> rdaux: %s <<<>>> \n\n ", rdaux);*/   
      if( !strncasecmp(rdaux,"**Main", 6) ) 
	ng = 1; 
      else if( !strncasecmp(rdaux,"**Molecule", 10) ) 
	ng = 2; 
      else if( !strncasecmp(rdaux,"**GA", 4) ) 
	ng = 3; 
      else if( !strncasecmp(rdaux,"**Force_Field", 13) ) 
	ng = 4;
      else if( !strncasecmp(rdaux,"**Dalton_Input", 14) ) 
	ng = 5;
      else if( !strncasecmp(rdaux,"**Gamess_Input", 14) ) 
	ng = 6;
      /*else if( !strncasecmp(rdaux,"**RS", 4) ) 
	ng = 7;*/
      else if( !strncasecmp(rdaux,"**End", 5) ) 
	ng = 1000;
      else
	{
	  length = strlen(rdaux);
	  printf("#|\n\nError: Group \"");
	  for (i = 0; i <= length - 1; i++) printf("%c",rdaux[i]);
	  printf("\" is not available!\n");
	  exit(EXIT_FAILURE);
	}
 
      switch (ng)
	{
	case 1: /* "**Main" */
	  for ( i = 0; ; i++ )
	    {
	      fscanf(fl0, "%s", &rdaux);
	      /* Uncomment bellow to debug */
	      /*printf ("\n\n\n <<<>>> rdaux: %s <<<>>>  \n\n\n", rdaux);*/
	      printf("#");

	      nsg = 0;
	      if( feof(fl0) || rdaux[1] == '*' ) break;
	      /* else if ( rdaux == "/$" ) break; */
	      else if ( !strncasecmp(rdaux,"*Type_of_calculation", 20) ) 
		nsg = 11;
	      else if ( !strncasecmp(rdaux,"*Title", 6) )
		nsg = 12;
	      else if ( !strncasecmp(rdaux,"*Dalton_path", 11) ) 
		nsg = 13;
	      else if ( !strncasecmp(rdaux,"*Gamess_path", 11) ) 
		nsg = 14;
	      else if ( !strncasecmp(rdaux,"*Angstron", 9) ) 
		nsg = 15;
	      else if ( !strncasecmp(rdaux,"*Print_Population", 17) ) 	      
		nsg = 16;
	      else if ( !strncasecmp(rdaux,"*Remove_DAT_Files", 24) ) 
		nsg = 17;
	      else if ( !strncasecmp(rdaux,"*ABS_Energy_Comparison", 22) ) 
		nsg = 18;
	      else if ( !strncasecmp(rdaux,"*Remove_OUT_Files", 17) )
		nsg = 19;
	      else if ( !strncasecmp(rdaux,"*Optimize_Structures", 20) )
		nsg = 20;
	      else if ( !strncasecmp(rdaux,"*NSTEP", 6) )
		nsg = 21;
	      
	      switch (nsg) 
		{
		case 11: /* "*Type of Calculation" */
		  /* fgets(toc, nchar, fl0);
		     fscanf does not allow to include comments in the input file*/
		  fscanf(fl0, "%s", &rdaux);
		  toc = Char_alloc(rdaux);
		  break;

		case 12: /* "*Title" */
		  fscanf(fl0, "\n", &rdaux); 
		  fgets(rdaux, nchar, fl0);
		  title = Char_alloc(rdaux);
		  /*printf("Title=%s\n", title);*/
		  break; 

		case 13: /* "*Dalton_path" */
		  fscanf(fl0, "%s", &rdaux);
		  daltonpath = Char_alloc(rdaux);
		  break;
		  
		case 14: /* "*Gamess_path" */
		  fscanf(fl0, "%s", &rdaux);
		  gamesspath = Char_alloc(rdaux);
		  break;

		case 15: /* "*Angstron" */
		  distfactor = bohr2angstron;
		  break;

		case 16: /* "Print_Population" */
		  fscanf(fl0, "%s", &rdaux);
		  print_pop = Char_alloc(rdaux);
		  break;
		  
		case 17: /* "Remove_DAT_Files" */
		  fscanf(fl0, "%s", &rdaux);
		  rmdatfile = Char_alloc(rdaux);
		  break;
		  
		case 18: /* "*ABS_Energy_Comparison" */
		  fscanf(fl0, "%lf", &absenrgy);
		  break;
		
		case 19: /* "Remove_OUT_Files" */
		  fscanf(fl0, "%s", &rdaux);
		  rmoutfile = Char_alloc(rdaux);
		  break;

		case 20: /* "Optimize_Structures" */
		  fscanf(fl0, "%s", &rdaux);
		  optstruct = Char_alloc(rdaux);
		  break;

		case 21: /* "NSTEP" */
		  fscanf(fl0, "%d", &nstep);
		  break;

		default:
		  exitg("**Main",rdaux);
		}
	    } 
	  break;

	case 2: /* "**Molecule" */
	  for ( i = 0; ; i++ )
	    {
	      fscanf(fl0, "%s", &rdaux);
	      /* Uncomment bellow to debug */
	      /*printf ("\n\n\n <<<>>> rdaux: %s <<<>>>  \n\n\n", rdaux);*/ 
	      if( feof(fl0) || rdaux[1] == '*' ) break;	       
	      /* else if ( rdaux[0] == '\0' ) break; */
	      else if ( !strncasecmp(rdaux,"*Number_of_atoms", 16) ) nsg = 21;
	      else if ( !strncasecmp(rdaux,"*ZMAT", 5) ) nsg = 22;

	      switch (nsg)
		{
		case 21: /* "*Number of atoms" */
		  fscanf(fl0, "%d", &natoms);
		  if( nindiv == ZERO ){
		    printf("The \"**GA\" parameters have to be declared before the \"**Molecule\"");
		    exit(1);
		  }
		  if (&GA != NULL ) nn = 2 * nindiv * natoms; 
		  else nn = nindiv * natoms; 
		  /*printf("\nnn= %d",nn);*/
		  MPop = (defmol *) calloc( nn, sizeof(defmol) );
		  XYZ  = (Txyz *) calloc( nn, sizeof(Txyz) );
		  energy = (double *) calloc( nn, sizeof(double) );
		  
		  if (&GA != NULL ){ 
		    cnct = (mzmat *) calloc( 3*natoms, sizeof(mzmat) );
		    zmat = (mzmat *) calloc( 9*natoms, sizeof(mzmat) );
		    conv = (int *) calloc ( nn, sizeof(int) );
		  }
		  break;

		case 22: /* "*ZMAT" */
		  /*FMolC = (defmol *) calloc(natoms, sizeof(defmol));*/ 
		 
		  for ( j = 1; j >= natoms; j++ ) {
		    fscanf(fl0, "%s\n", &rdaux);
		    scanf(rdaux);
		  }
		  break;
		  
		default:
		  exitg("**Molecule",rdaux);
		}
	    }
	  break;



	case 3: /* "**GA" */
	  for ( i = 0; ; i++ )
	    {
	      fscanf(fl0, "%s", &rdaux);
	      /* Uncomment bellow to debug */
	      /*printf ("\n\n\n <<<>>> rdaux: %s <<<>>>  \n\n\n", rdaux); /**/ 
	      if( feof(fl0) || rdaux[1] == '*' ) break;	       
	      /* else if ( rdaux[0] == '\0' ) break; */
	      else if ( !strncasecmp(rdaux,"*Type_of_GA", 11) ) nsg = 30;
	      else if ( !strncasecmp(rdaux,"*Population_Size", 16) ) nsg = 31;
	      else if ( !strncasecmp(rdaux,"*Generations", 12) ) nsg = 32;
	      else if ( !strncasecmp(rdaux,"*Generation_After_Convergence", 29) ) nsg = 33;
	      else if ( !strncasecmp(rdaux,"*Seed", 5) ) nsg = 34;
	      else if ( !strncasecmp(rdaux,"*Number_Rotation_Groups", 23) ) nsg = 35;
	      else if ( !strncasecmp(rdaux,"*Rotation_Mutation_Range", 24) ) nsg = 36;
	      else if ( !strncasecmp(rdaux,"*Type_of_Mating", 15) ) nsg = 37;
	      else if ( !strncasecmp(rdaux,"*Percentage_of_Mutation", 23) ) nsg = 38;
	      else if ( !strncasecmp(rdaux,"*Groups_to_be_Rotate", 20) ) nsg = 39;
	      
	      switch (nsg)
		{
		case 30: /* Type_of_GA */
		  fscanf(fl0, "%s", &rdaux);
		  /* Uncomment bellow to debug */
		  /*printf ("\n\n\n <<<>>> rdaux: %s <<<>>>  \n\n\n", rdaux); /**/ 
		  if( feof(fl0) || rdaux[1] == '*' ) break;
		  else if ( !strncasecmp(rdaux,"random-search", 13) ) nga = 301;
		  else if ( !strncasecmp(rdaux,"elitist", 7) ) nga = 302;
		  else if ( !strncasecmp(rdaux,"non-elitist", 11) ) nga = 303;
		  
		  switch (nga){
		  case 301: /* "random-seach" */
		    RS = 1;
		    break;
		    
		  case 302: /* "elitist" */
		    GA = 1;
		    break;
		    
		  case 303: /* "non-elitist" */
		    GA = 1;
		    typga = 1;
		    break;
		    
		  default:
		    exitg("*Type_of_GA",rdaux);
		  }
		  break;


		case 31: /* "*Population_Size" */
		  fscanf(fl0, "%d", &nindiv);
		  //printf("\nPop=%d", nindiv);
		  break;

		case 32: /* "*Generations" */
		  fscanf(fl0, "%d", &ngen);
		  break;

		case 33: /* "*Generation_After_Convergence" */
		  fscanf(fl0, "%d", &genafconv);
		  //printf("gen = %d", genafconv);
		  break;

		case 34: /* "*Seed" */
		  fscanf(fl0, "%li", &nseed);
		  break;

		case 35: /* "*Number_Rotation_Groups" */
		  fscanf(fl0, "%d", &nrotgp);
		  //printf("nrotgp = %d\n", nrotgp);
		  break;

		case 36: /* "*Rotation_Mutation_Range" */
		  fscanf(fl0, "%d %d", &mutmin, &mutmax);
		  break;

		case 37: /* "*Type_of_Mating" */
		    fscanf(fl0, "%s", &rdaux);
		    /* Uncomment bellow to debug */
		    /*printf ("\n\n\n <<<>>> rdaux: %s <<<>>>  \n\n\n", rdaux); /**/ 
		    if( feof(fl0) || rdaux[1] == '*' ) break;
		    //else if ( rdaux[0] == '\0' ) break; 
		    else if ( !strncasecmp(rdaux,"single", 6) ) nmat = 371;
		    else if ( !strncasecmp(rdaux,"double", 6) ) nmat = 372;
		    else if ( !strncasecmp(rdaux,"mixed", 5) ) nmat = 373;
		    
		    switch (nmat){
		    case 371: /* "single" */
		      typmat = 1;
		      break;
		      
		    case 372: /* "double" */
		      typmat = 2;
		      break;
		      
		    case 373: /* "mixed" */
		      typmat = 3;
		      break;

		    default:
		      exitg("*Type_of_Mating",rdaux);
		    }
		  break;
		  
		case 38: /* "*Percentage_of_Mutation" */
		  fscanf(fl0, "%d", &percmut);
		  break;

		case 39: /* "*Groups_to_be_Rotate" */
		  rtt  = (Trtt *) calloc( 3*nrotgp, sizeof(Trtt) );
		  rttatm  = (Tatom *) calloc( 500*nrotgp, sizeof(Tatom) );

		  count=0;
		  for( i = 0; i < nrotgp; i++ ){
		    fscanf(fl0, "%d", &natmrot);
		    rtt[i].natmrot = natmrot;
		    
		    fscanf(fl0, "%d", &center);
		    rtt[i].center = center;
		    
		    fscanf(fl0, "%d", &plane);
		    rtt[i].plane = plane;
		    
		    for( j = 0; j < natmrot; j++){
		      fscanf(fl0, "%d", &atm);
		      rttatm[count].atom = atm;
		      count=count+1;
		    }
		  }
		  break;


		default:
		  exitg("**GA",rdaux);
		}
	    }
	  break;


	case 4: /* "**Force Field" */
	  for ( i = 0; ; i++ )
	    {
	      fscanf(fl0, "%s", &rdaux);
	      /* Uncomment bellow to debug */
	      /*printf ("\n\n\n <<<>>> rdaux: %s <<<>>>  \n\n\n", rdaux);*/
	      if( feof(fl0) || rdaux[1] == '*' ) break;
	      else if ( !strncasecmp(rdaux,"*Number_of_partitions", 21) )
		nsg = 41;
	      else if ( !strncasecmp(rdaux,"*Levels_of_calculation", 22) )
		nsg = 42;
	      else if ( !strncasecmp(rdaux,"*Basis_sets", 11) )
		nsg = 43;
	      else if ( !strncasecmp(rdaux,"*Multiplicity", 13) )
		nsg = 44;
	      else if ( !strncasecmp(rdaux,"*External_EEC_programs", 22) )
		nsg = 45;
	      else if ( !strncasecmp(rdaux,"*Molecular_charge", 17) )
		nsg = 46;
	      else if ( !strncasecmp(rdaux,"*Potential", 10) )
		nsg = 47;
	      else
		nsg = 999;

	      switch (nsg)
		{
		case 41: /* "*Number of Sections" */
		  fscanf(fl0, "%d", &nos);
		  /* Uncomment bellow to debug */
		  /*printf("nos_rdinput2 = %d",*nos); /**/

		  noas = *Int_alloca(nos);

		  kai = 0;
		  for ( i = 0; i < nos; i++ )
		    {
		      fscanf(fl0, "%s", &rdaux);
		      if ( strncasecmp(rdaux,".Partition", 8) )
			{
			  printf("|\n\n ** Error: Syntax declaration inside section **\n");
			  exit (EXIT_FAILURE);
			}

		      fscanf(fl0, "%d %s", (noas + i), &rdaux);

		      kaf = kai + (noas + i);
		      /*  if ( pPMolC == NULL )
			PMolC = (defmol *) calloc(kaf, sizeof(defmol));
		      else
			PMolC = realloc (PMolC, kaf*sizeof(defmol));

		      pPMolC = PMolC;

		      getpmol(fl0, kai, kaf, &PMolC, &FMolC, distfactor); */

		      kai = kaf;
		    }
		  /* PMolC[11].xyz.x = 9.9989; */
		  break;

		case 42: /* "*Level_of_calculation" */
		  fscanf(fl0, "%s", &rdaux);

		  lock = 0;
		  if ( !strncasecmp(rdaux,"All", 3) ) lock = 1;

		  loc = (char *) malloc (nos * sizeof(char));
		  psdpot = (char *) malloc (nos * sizeof(char));

		  fscanf(fl0, "%s", &rdaux);
		  for (i = 0; i < nos; i++)
		    {
		      if ( lock == 0 && i != 0 )
			fscanf(fl0, "%s %s", rdaux, rdaux);
		      loc = Char_alloc(rdaux);
		      
		      if ( !strncasecmp(&loc[i], "DFT", 3) && lock == 0 ){
			fscanf(fl0, "%s", &rdaux);
			psdpot = Char_alloc(rdaux);
		      }
		      else if ( !strncasecmp(&loc[i], "DFT", 3) && lock == 1
				&& i == nos - 1 ){
			fscanf(fl0, "%s", &rdaux);
			for ( j = 0;  j < nos; j++ ) 
			  psdpot[j] = *Char_alloc(rdaux);
		      }
		      else psdpot = NULL; 
		    }
		  lock = 0;
		  break;

		case 43: /* "*Basis sets" */
		  fscanf(fl0, "%s", &rdaux);
		  lock = 0;
		  if ( !strncasecmp(rdaux,"All", 3) ) lock = 1;

		  bs = (char *) malloc (nos * sizeof(char));/**/
		  fscanf(fl0, "%s", &rdaux);
		  //fgets(rdaux, 27, fl0);
		  for (i = 0; i < nos; i++)
		    {
		      if ( lock == 0 && i != 0 ) {
			fscanf(fl0, "%s %s", &rdaux, &rdaux);
			fscanf(fl0, "%s", &rdaux);
			fgets(rdaux, 27, fl0);
		      }
		      bs = Char_alloc(rdaux);
		    }
		  lock = 0;
		  break;
 
		case 44: /* "*Multiplicity" */
		  fscanf(fl0, "%s", &rdaux);
		  
		  lock = 0;
		  if ( !strncasecmp(rdaux,"All", 3) ) lock = 1;
		  
		  mltpl = *Int_alloca(nos);
		  
		  fscanf(fl0, "%d", &mltpl);
		  /*
		  for (i = 1; i < nos; i++)
		    {
		      if ( lock == 0 )
			fscanf(fl0, "%s %d", &rdaux, (mltpl + i));
		      else
			(mltpl + i) = mltpl[0];
		    }
		  */
		  lock = 0;
		  

		  break;

		case 45: /* "*External_EEC_Programs" */
		  fscanf(fl0, "%s", &rdaux);
		  lock = 0;
		  if ( !strncasecmp(rdaux,"All", 3) ) lock = 1;

		  /*eeecp = (char **) malloc (nos * sizeof(char)) */;/**/
		  fscanf(fl0, "%s", &rdaux);
		  for (i = 0; i < nos; i++)
		    {
		      if ( lock == 0 && i != 0 )
			fscanf(fl0, "%s %s", rdaux, rdaux);
		      /*  eeecp[i] = Char_alloc(rdaux); */
		    }
		  lock = 0;
		  break;

		case 46:
		  fscanf(fl0, "%s", &rdaux);
		  if ( !strncasecmp(rdaux,"All", 3) ) lock = 1;

		  icharg = *Int_alloca(nos);
		  
		  fscanf(fl0, "%d", &icharg);
		  /*
		  for (i = 1; i < nos; i++)
		    {
		      if ( lock == 0 )
			fscanf(fl0, "%s %d", &rdaux, (icharg + i));
		      else
			(icharg + i) = icharg[0];
		    }
		  */
		  lock = 0;
		  
		  break;

		case 47:
		  fscanf(fl0, "%s", &rdaux);
		  potential = Char_alloc(rdaux);
		  /*printf("potential inside rdinput -> %s", potential);*/
		  break;

		default:
		  exitg("**Force_Field",rdaux);
		}
	    }
	  break;

	case 5:
	  fl1 = fopen("dalton.dal", "w");
	  while( !strncasecmp(rdaux,"*End_Dalton", 11) )
	    {
	      fscanf(fl0, "\n", &rdaux);
	      fgets(rdaux, nchar, fl0);
	      fprintf(fl1, "%s", rdaux);
	    }
	  break;

	case 6:
	  strcpy(fname, "gamess0.inp");
	  
	  fl1 = fopen(fname, "w");
	  
	  fgets(rdaux, nchar, fl0);
	  fprintf(fl1, "%s", rdaux);
	  
	  do
	    {
	      fscanf(fl0, "\n", &rdaux);
	      fgets(rdaux, nchar, fl0);
	      fprintf(fl1, "\n%s", rdaux);
	    }while( !strncasecmp(rdaux,"$DATA", 5) );
	  
	  if ( loc != NULL || bs != NULL )fprintf(fl1,"\r", *loc);
	  close(fl1);
	  
	  /*	  for(i = 0; i < nos; i++)
		  {
		  fname[7] = i;
		  fl1 = fopen(fname, "w");
		  
		  fgets(rdaux, nchar, fl0);
		  fprintf(fl1, "%s", rdaux);
		  
		  if ( loc != NULL )
		  {
		  fprintf(fl1,"\r$CONTRL SCFTYP=%s $END", *loc);
		  fprintf(fl1,"\n$DATA");
		  }
		  
		  if ( bs != NULL )
		  {
		  fprintf(fl1,"\r$BASIS GBASIS=%s $END", *bs);
		  fprintf(fl1,"\n$DATA");
		  }*/
	  
	  break;

	  
	case 1000: /* "**End" */
	  /* printf("\nEnd**"); /**/
	  eof = 0;
	  break;
	}
    }
  /*
    if ( icharg == NULL ) {
    icharg = Int_alloca(nos);
    for (i = 0; i < nos; i++)
    *(icharg + i) = ZERO;
    }
  */
  fclose(fl0);
  
  return status;
}

/* Foram Comentados os FMolC, getmol, PMolC, getpmol, eeecp*/

