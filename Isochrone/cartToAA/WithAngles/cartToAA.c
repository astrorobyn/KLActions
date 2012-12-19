#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bodyAux.h"
#include "isochroneMap.h"
#include "io.h"
#include "types.h"
#include "vector.h"
#include "AminaVersion.h"

//these hold the mass and scale radius of the potential

double M_ISO, B_ISO;

double theta1(double J1, double J2, double q[3]) {
  // computes the angle associated with Lz, which is a constant
  // (the longitude of the ascending node, natch)

  double r = norm(q);
  double ph = phi(q);
  double z = q[2];
 
  double sinu = J1 * z / sqrt((J2*J2-J1*J1)*(r*r-z*z));

  return ph - asin(sinu); 

}






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
	
  double E;
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
	
    //calculate energy to check for bound orbit
    E = H(q,qdot);

    if(E < 0) {	
      //calculate actions & angles &  write to file
      
      ActionsFrequenciesAndDerivatives result;

      //auxiliary stuff

      double LZ = bs[i]->q[0]*bs[i]->p[1] - bs[i]->q[1]*bs[i]->p[0];
      double L = Lmag(q, qdot);

      result = computeActionsFrequenciesAndDerivatives(E, L, LZ, q[0], qdot[0], M_ISO, B_ISO);
      result.Thetatheta = theta2(etaval,ph,q,qdot);
      result.Thetaphi = theta1(LZ, L, bs[i]->q)
