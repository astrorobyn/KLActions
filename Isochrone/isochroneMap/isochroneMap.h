/*
 *  isochroneMap.h
 *  
 *
 *  Created by Robyn Sanderson on 7/6/10.
 *
 */


#ifndef __ISOMAP_H__
#define __ISOMAP_H__

#include <assert.h>
#include <stdio.h>
#include <math.h>

// root finder to solve Kepler's equation needs
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

//type for calling function with parameters in root-finding algorithm
typedef struct {double A, theta;} keplerEquationParamsType;

//parameters for isochrone potential
#define M_ISO 12.2
#define B_ISO 8.0

//various auxiliary functions
double keplerEquation(double eta, void *p);
double Phi(double r);
double H(double q[2], double qdot[2]);
double Lmag(double q[2], double qdot[2]);
double eccentricity(double q[2], double qdot[2]);
double c_aux(double q[2], double qdot[2]);
double Omega2(double q[2], double qdot[2]);
double Omega3(double q[2], double qdot[2]);
double eta(double q[2], double qdot[2], double t, double eta0);
double etaZero(double q[2], double qdot[2]);
double theta3(double etaval, double q[2], double qdot[2]);
double Jr(double q[2], double qdot[2]);

void forwardMap(double q[2], double qdot[2], double t, double x[2], double v[2]);

#endif // __ISOMAP_H__
