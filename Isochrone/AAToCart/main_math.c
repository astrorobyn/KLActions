/*
 *  main.c
 *  
 *
 *  Created by Robyn Sanderson on 7/20/12.
 *
 */
#include <string.h>

// MathLink Libraries
#include <mathlink.h>

// Error codes

#define ERR_UNBOUND 1
#define ERR_UNPHYS 2
#define WARN_RADORB 3
#define SUCCESS 0

#include "aa_to_cart.h"
#include "eta_from_aa.h"

//Mathematica passes the lengths of the lists as well as the lists themselves so there are two extra args in the main call

void aa_to_cart(double J[3], long Jlen, double th[3], long thlen) {

  double  x[3], v[3], xp[3], vp[3], xpp[3], vpp[3];
  double xorb[2], vorb[2];    //pos & vel in polar coords of orbital pln
  double w[6]; //full phase space position in Cartesian
  
  double h, c, ecc, o2, o3, eta, inc; 

  int radorb;
  
   // get auxiliary variables
 
  h = H_from_J(J);

  // check that orbit is bound (E<0)
  if (h>0) {

   MLEvaluateString(stdlink,"Print[\"Error: Orbit is unbound\"]"); 
   MLPutRealList(stdlink, w, 6);
 
    return;
  }

  //check that orbit is physical (L>=Lz)
  if (J[1]<J[0]) {
   MLEvaluateString(stdlink,"Print[\"Error: Orbit has unphysical (L, Lz)\"]"); 
   MLPutRealList(stdlink, w, 6);
 
    return;
  }

 //check for radial orbit
  radorb=(fabs(J[1])<TINY);


  c = c_from_H(h);
  ecc = e_from_J(J,c);

  o3 = Omega3_from_J(J);
  o2 = Omega2_from_J(J,o3);
      
  //get eta
  eta=eta_from_aa(th[2],ecc,c);

  //get phi in the orbital plane
  xorb[1] = phiorb_from_aa(th,J,o2,o3, ecc, eta, c);

  //get r in the orbital plane
  xorb[0] = rorb_from_aa(eta, ecc, c);

  //angular momentum gives v_phi
  vorb[1] = J[1]/xorb[0];

  //consv of energy gives v_r
  vorb[0] = vrorb_from_econsv(h,J[1],xorb);

  //get the Cartesian position and velocity
  cartesianFromPolar(xorb,vorb,xpp,vpp);

  //two rotations: about y axis for the inclination
  
  inc = radorb ? 0 : acos(J[0]/J[1]);
  rotateY(xpp,vpp,xp,vp,inc);

  //then about z axis to put the line of nodes in the right place
  rotateZ(xp,vp,x,v,th[0]);

  //package the result in one vector for returning
  memcpy(w, x, 3*sizeof(double));
  memcpy(w + 3, v, 3*sizeof(double));

  //Print a warning if the orbit was radial
  if(radorb)
    MLEvaluateString(stdlink,"Print[\"Warning: Radial orbit, angle of orbital plane is undefined. Placing orbit in x-y plane.\"]");

  //This returns the Cartesian coos to mathematica
  MLPutRealList(stdlink, w, 6);

    return;
	
}

int main(int argc, char* argv[]) {
  return MLMain(argc, argv);
}
