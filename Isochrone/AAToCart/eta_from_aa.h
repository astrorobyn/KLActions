

#include <stdio.h>
#include <math.h>

// root finder to solve Kepler's equation needs
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

//isochrone parameters
#include "isoparams.h"

//type for calling function with parameters in root-finding algorithm
typedef struct {double A, theta;} keplerEquationParamsType;

double eta_from_aa(double th3, double ecc, double c);
double keplerEquation(double eta, void *p);
