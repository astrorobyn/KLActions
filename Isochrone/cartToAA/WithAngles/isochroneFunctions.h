#ifndef __ISOFUNC_H__
#define __ISOFUNC_H__

#include <math.h>
#include <gsl/gsl_math.h>

#include "vector.h"

double Phi(double r);
double H(double q[2], double qdot[2]);
double Lmag(double q[2], double qdot[2]);
double Jr(double q[2], double qdot[2]);

double eccentricity(double q[2], double qdot[2]);
double c_aux(double q[2], double qdot[2]);
double etaZero(double q[2], double qdot[2]);
double Omega2(double q[2], double qdot[2]);
double Omega3(double q[2], double qdot[2]);

double theta3(double etaval, double q[2], double qdot[2]);
double theta2(double etaval, double phi, double q[2], double qdot[2]);
double theta1(double J1, double J2, double q[3]);

#endif  //__ISOFUNC_H__
