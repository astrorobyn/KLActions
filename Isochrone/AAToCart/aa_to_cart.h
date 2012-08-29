//functions for computing Cartesian (x,v) from actions (J1, J2, J3)=(Lz, L, Jr)
//and angles (th1, th2, th3)=(Omega, th_th, th_r)
//(Omega is the longitude of the ascending node)

#include <math.h>

//parameters for isochrone potential
#include "isoparams.h"

//threshold for overflow/underflow checks
#define TINY 1e-10

double H_from_J(double J[3]);
double c_from_H(double h);
double e_from_J(double J[3], double c);

double Omega3_from_J(double J[3]);
double Omega2_from_J(double J[3], double omega3);


double phiorb_from_aa(double th[3], double J[3], double omega2, double omega3, double ecc, double eta, double c);
double rorb_from_aa(double eta, double ecc, double c);

double vrorb_from_econsv(double energy, double angmom, double xorb[2]);

void cartesianFromPolar(double xorb[2], double vorb[2], double xpp[3], double vpp[3]);

void rotateY(double xpp[3], double vpp[3], double xp[3], double vp[3], double angle);

void rotateZ(double xp[3], double vp[3], double x[3], double v[3], double angle);


double Phi(double r);
