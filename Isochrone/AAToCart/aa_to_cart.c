#include "aa_to_cart.h"

double H_from_J(double J[3]) {

  // decrement pointer so that formula reads with the right J-subscripts
  J--;

  double denomsqrt = (J[3] + 0.5 *(J[2] + sqrt(J[2] * J[2] + FourGMb)));

  //put the indexing back!
  J++;

  return -(GM_ISO * GM_ISO) / (2.0 * denomsqrt * denomsqrt);

}

double c_from_H(double h) {

  return GM_ISO /(-2.0 * h) - b_ISO;


}

double e_from_J(double J[3], double c) {
 
  // decrement pointer so that formula reads with the right J-subscripts
  J--;

  double eccsq;
  if(J[3]<TINY)        //this avoids truncation error problems with circular orbits
    eccsq=0;
  else
    eccsq = 1. - J[2]*J[2]/(GM_ISO * c) * (1. + b_ISO/c);

 //put the indexing back!
  J++;

  return sqrt(eccsq);

}

double Omega3_from_J(double J[3]) {

 // decrement pointer so that formula reads with the right J-subscripts
  J--;

  double denom13 = J[3] + 0.5 * (J[2] + sqrt(J[2] * J[2] + FourGMb));

 //put the indexing back!
  J++;

  return GM_ISO * GM_ISO / (denom13 * denom13 * denom13);

}

double Omega2_from_J(double J[3], double omega3) {

 // decrement pointer so that formula reads with the right J-subscripts
  J--;

  double pref = 1. + J[2] / sqrt(J[2] * J[2] + FourGMb);

 //put the indexing back!
  J++;

  return 0.5 * pref * omega3;
}



double phiorb_from_aa(double th[3], double J[3], double omega2, double omega3, double ecc, double eta, double c) {

// decrement pointers so that formula reads with the right subscripts
  J--;
  th--;

  double tanterm1, tanterm2, phiorb;

  if(fabs(1-ecc)<TINY)     //checks for circular orbits to avoid overflow
    tanterm1=M_PI/2.0;
  else
    tanterm1 = atan(sqrt((1.+ecc)/(1.-ecc)) * tan(0.5*eta));

  tanterm2 = sqrt((1. + ecc+2.* b_ISO/c)/(1. - ecc+2.* b_ISO/c)) * tan(0.5*eta);
  phiorb = th[2] - omega2/omega3 * th[3] + tanterm1 + atan(tanterm2)/sqrt(1. + FourGMb /(J[2] * J[2]));

 //put the indexing back!
  J++;
  th++;

  return phiorb;

}

double rorb_from_aa(double eta, double ecc, double c) {

  double sqterm = 1.0 + c/b_ISO *(1.0 - ecc * cos(eta));

  return b_ISO * sqrt(sqterm * sqterm - 1.0);

}


double vrorb_from_econsv(double energy, double angmom, double xorb[2]) {

  double pot = Phi(xorb[0]);
  double vrsq = 2.0 * (energy - pot) - angmom * angmom / (xorb[0] * xorb[0]);

  // the following test avoids truncation error problems with circular orbits
  double vrmag = (fabs(vrsq)<TINY) ? 0.0 : sqrt(vrsq);

  double sign = (sin(xorb[1]) < 1.0) ? 1.0 : -1.0;

  return sign * vrmag;


}

void cartesianFromPolar(double xorb[2], double vorb[2], double xpp[3], double vpp[3]) {

  double cosphi = cos(xorb[1]);
  double sinphi = sin(xorb[1]);

  xpp[0] = xorb[0] * cosphi;
  xpp[1] = xorb[0] * sinphi;

  vpp[0] = vorb[0] * cosphi - vorb[1] * sinphi;
  vpp[1] = vorb[0] * sinphi  + vorb[1] * cosphi;

}

void rotateY(double xpp[3], double vpp[3], double xp[3], double vp[3], double angle) {

  double cosi = cos(angle);
  double sini = sin(angle);

  xp[0] = xpp[0] * cosi + xpp[2] * sini;
  xp[1] = xpp[1];
  xp[2] = -xpp[0] * sini + xpp[2] * cosi;

  vp[0] = vpp[0] * cosi + vpp[2] * sini;
  vp[1] = vpp[1];
  vp[2] = -vpp[0] * sini + vpp[2] * cosi;

}


void rotateZ(double xp[3], double vp[3], double x[3], double v[3], double angle) {

  double costh1 = cos(angle);
  double sinth1 = sin(angle);
  
  x[0] = xp[0] * costh1 - xp[1] * sinth1;
  x[1] = xp[0] * sinth1 + xp[1] * costh1;
  x[2] = xp[2];

  v[0] = vp[0] * costh1 - vp[1] * sinth1;
  v[1] = vp[0] * sinth1 + vp[1] * costh1;
  v[2] = vp[2];

}


double Phi(double r) {
	double p = -1.0 * GM_ISO / (b_ISO + sqrt( r * r + b_ISO * b_ISO));
	return p;
}
