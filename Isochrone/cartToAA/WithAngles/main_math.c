#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mathlink.h>

#include "isochroneFunctions.h"
#include "bodyAux.h"
#include "types.h"
#include "vector.h"

#define KM_S_TO_KPC_MYR 0.001022729843472133

//these hold the mass and scale radius of the potential

double M_ISO, B_ISO;

//main function

void cart_to_aa(double Miso, double biso, double x[3], long xlen, double v[3], long vlen){
	
  body *b;
	
  double Jth[6];
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
    for(i=0;i<6;i++) {
      Jth[i]=-100.00;
      }
    MLPutReal64List(stdlink, Jth, 6);
    

  } else {

    //calculate actions 
    Jth[0] = b->q[0]*b->p[1] - b->q[1]*b->p[0];
    Jth[1] = Lmag(q, qdot);
    Jth[2] = Jr(q, qdot);


    //calculate angles
    double etaval = etaZero(q, qdot);
    double ph = phi(b->q);
    Jth[3]=theta1(Jth[0], Jth[1], b->q);
    Jth[4]=theta2(etaval,ph,q,qdot);
    Jth[5]=theta3(etaval,q,qdot);
 
    MLPutReal64List(stdlink, Jth, 6);
  }

  free(b);
  return;
}

//mathlink wrapper

int main(int argc, char * argv[]) {
  return MLMain(argc,argv);
}
