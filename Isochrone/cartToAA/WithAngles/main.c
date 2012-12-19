#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bodyAux.h"
#include "isochroneFunctions.h"
#include "io.h"
#include "types.h"
#include "vector.h"

//these hold the mass and scale radius of the potential

double M_ISO, B_ISO;



int main (int argc, const char * argv[]) {


  // reads in two text files containing the positions and velocities
  // outputs actions in a text file and optionally angles in a second file
	

  if (argc<6) {
    fprintf(stderr, "Usage: actionCalc M b posinfile velinfile outfile [angleoutfile]\n");
    fprintf(stderr, "\t M, b: potential parameters (mass in kpc^3/Myr^2 and scale radius in kpc)\n");
    fprintf(stderr, "\t [pos,vel]infile: text files of positions and velocities\n");
    fprintf(stderr, "\t outfile: name of file to which actions will be printed\n");
 fprintf(stderr, "\t angleoutfile: (optional) name of file to which angles will be printed\n");
    return 1;
  }
	
  //step one: read in file
	
  body **bs;
	
  size_t np;
	
  int stat = read_system_from_text_files(argv[3], argv[4], &bs, &np);
	
  if (stat) {
    return 1;
  }
	
  fprintf(stdout, "Read %ld particles from files %s and %s.\n", np, argv[3], argv[4]);
	
  //step two: calculate actions and so forth & write to file
	
  double J3, J2, J1, th3, th2, th1;
  double q[2],qdot[2];
 
	
  size_t i, n_write;
	
  FILE * outfile, *aoutfile;
  
  M_ISO = atof(argv[1]);
  B_ISO = atof(argv[2]);
  outfile = fopen(argv[5], "w");
  if(argc==7)
    aoutfile=fopen(argv[6],"w");
	
  n_write = 0;
  for (i=0; i<np; i++) {
		
    bodyStrip(bs[i], q, qdot);
		
    if(H(q,qdot) < 0) {
      //calculate actions & write to file
      J1 = bs[i]->q[0]*bs[i]->p[1] - bs[i]->q[1]*bs[i]->p[0];
      J2 = Lmag(q, qdot);
      J3 = Jr(q, qdot);
		
      fprintf(outfile, "%lg %lg %lg\n", J1, J2, J3);
      n_write++;

      //if desired, calculate angles and write to file
      if(argc==5) {
	double etaval = etaZero(q, qdot);
	double ph = phi(bs[i]->q);
	th1=theta1(J1, J2, bs[i]->q);
	th2=theta2(etaval,ph,q,qdot);
	th3=theta3(etaval,q,qdot);
	fprintf(aoutfile,"%lg %lg %lg\n", th1, th2, th3);
      }
    }
    else {
      //if orbit is unbound write nonsense values to that line
      fprintf(outfile, "-100 -100 -100\n");
      if(argc==5) {
	fprintf(aoutfile,"-100 -100 -100\n");
      }
    }
  }

  fprintf(stdout, "Wrote %ld valid actions to file %s; %ld objects had unbound orbits\n", n_write, argv[3], np-n_write);
	
  fclose(outfile);
	
}
