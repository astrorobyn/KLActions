#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bodyAux.h"
#include "isochroneMap.h"
#include "io.h"
#include "types.h"
#include "vector.h"

int main (int argc, const char * argv[]) {


  // reads in two text files containing the positions and velocities
  // outputs actions in a text file
	

  if (argc!=4) {
    fprintf(stderr, "Usage: actionCalc posinfile velinfile outfile\n");
    fprintf(stderr, "\t [pos,vel]infile: text files of positions and velocities\n");
    fprintf(stderr, "\t outfile: name of file to which actions/frequencies/energies will be printed\n");
    return 1;
  }
	
  //step one: read in file
	
  body **bs;
	
  size_t np;
	
  int stat = read_system_from_text_files(argv[1], argv[2], &bs, &np);
	
  if (stat) {
    return 1;
  }
	
  fprintf(stdout, "Read %ld particles from files %s and %s.\n", np, argv[1], argv[2]);
	
  //step two: calculate actions and so forth & write to file
	
  double J3, J2, J1;
  double q[2],qdot[2];
	
  size_t i, n_write;
	
  FILE * outfile;
	
  outfile = fopen(argv[3], "w");
	
  n_write = 0;
  for (i=0; i<np; i++) {
		
    bodyStrip(bs[i], q, qdot);
		
    if(H(q,qdot) < 0) {
      J1 = bs[i].q[0]*bs[i].p[1] - bs[i].q[1]*bs[i].p[0];
      J2 = Lmag(q, qdot);
      J3 = Jr(q, qdot);
		
      fprintf(outfile, "%lg %lg %lg\n", J1, J2, J3);
      n_write++;
    }
    else 
      fprintf(outfile, "0 0 0\n");
  }

  fprintf(stdout, "Wrote %ld valid actions to file %s; %ld objects had unbound orbits\n", n_write, argv[3], np-n_write);
	
  fclose(outfile);
	
}
