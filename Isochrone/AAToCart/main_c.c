/*
 *  main.c
 *  
 *
 *  Created by Robyn Sanderson on 7/20/12.
 *
 */
#include <string.h>

#include "aa_to_cart.h"
#include "eta_from_aa.h"

int main(void) {


   double  x[3], v[3], xp[3], vp[3], xpp[3], vpp[3];
  double xorb[2], vorb[2];    //pos & vel in polar coords of orbital pln
  double w[6]; //full phase space position in Cartesian

  double h, c, ecc, o2, o3, eta, inc; 

  //assign debugging values to J and th
  double J[3] = {-5.42583, 8.25335, 1.98504};
  double th[3] = {1.00059, 2.00102, 6.00004};
 
  // get auxiliary variables
 
  h = H_from_J(J);

  //check that orbit is bound (E<0)

  if (h>0) {
    fprintf(stderr, "Error: Orbit is unbound, E = %.3lg\n", h);
    return 1;
  }

  //check that orbit is physical (L>=Lz)
  if (J[1]<J[0]) {
    fprintf(stderr, "Error: Orbit has unphysical (L, Lz) = (%.3lg, %.3lg)\n", J[1], J[0]);
    return 1;
  }

  //check for radial orbit
  if(fabs(J[1])<TINY) 
    fprintf(stderr, "Warning: Radial orbit, angle of orbital plane is undefined.\nPlacing orbit in x-y plane.\n");

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
  inc = (fabs(J[1])<TINY) ? 0 : acos(J[0]/J[1]);
  rotateY(xpp,vpp,xp,vp,inc);

  //then about z axis to put the line of nodes in the right place
  rotateZ(xp,vp,x,v,th[0]);

  //print the result
  fprintf(stdout, "(Lz, L, Jr) = (%0.5g, %0.5g, %0.5g)\n", J[0], J[1], J[2]);
  fprintf(stdout, "(th1, th2, th3) = (%0.5g, %0.5g, %0.5g)\n", th[0], th[1], th[2]);
  fprintf(stdout, "(x,y,z) = (%0.5g, %0.5g, %0.5g)\n", x[0], x[1], x[2]);
  fprintf(stdout, "(vx, vy, vz) = (%0.5g, %0.5g, %0.5g)\n", v[0], v[1], v[2]);

   return 0;
	
}

