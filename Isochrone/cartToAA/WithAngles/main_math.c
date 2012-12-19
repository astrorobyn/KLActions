#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mathlink.h>

#include "bodyAux.h"
#include "types.h"
#include "vector.h"

#define KM_S_TO_KPC_MYR 0.001022729843472133

static double M_ISO;
static double B_ISO;

//functions for energy, angular momentum, etc

double Phi(double r) {
	double p = -1.0 * M_ISO / (B_ISO + sqrt( r * r + B_ISO * B_ISO));
	return p;
}

double H(double q[2], double qdot[2]) {
  	double h = Phi(q[0]) + 0.5 * (qdot[0] * qdot[0] + q[0] * q[0] * qdot[1] * qdot[1]);

	return h;
}

double Lmag(double q[2], double qdot[2]) {
	return q[0] * q[0] * qdot[1];
}


double Jr(double q[2], double qdot[2]) {
	
	double ham = H(q,qdot);
	double J2 = Lmag(q, qdot);
	
	return M_ISO / sqrt(-2.0 * ham) - 0.5 * (J2 + sqrt(J2 * J2 + 4.0 * M_ISO * B_ISO));
	
}

//main function

void cart_to_aa(double Miso, double biso, double x[3], long xlen, double v[3], long vlen){
	
  body *b;
	
  double J[3];
  double q[2],qdot[2];
 
	
  size_t i;
	
  M_ISO = Miso;
  B_ISO = biso;

  //set up body with x and v

  for(i=0;i<3;i++)
    v[i]*=KM_S_TO_KPC_MYR;
  b=malloc(sizeof(body));
  memcpy(b->q,x,3*sizeof(double));
  memcpy(b->p,v,3*sizeof(double));
  b->m = 1.0;

  bodyStrip(b, q, qdot);
		
  if(H(q,qdot) > 0) {
    //    MLEvaluateString(stdlink,"Print[\"Error: Orbit is unbound\"]"); 
    for(i=0;i<3;i++) 
      J[i]=-100.00;
    MLPutReal64List(stdlink, J, 3);

  } else {

    //calculate actions 
    J[0] = b->q[0]*b->p[1] - b->q[1]*b->p[0];
    J[1] = Lmag(q, qdot);
    J[2] = Jr(q, qdot);
    MLPutReal64List(stdlink, J, 3);

  }

  free(b);
  return;
}

//mathlink wrapper

int main(int argc, char * argv[]) {
  return MLMain(argc,argv);
}
