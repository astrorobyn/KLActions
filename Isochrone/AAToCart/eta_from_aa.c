#include "eta_from_aa.h"

double eta_from_aa(double th3, double ecc, double c) {

 const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	gsl_root_fsolver *s;
	gsl_function keq;


	// iteration limits for root solver
	
	double epsabs=0.0;
	double epsrel = 1e-10;
	size_t max_iter=100;
	size_t iter=0;
	int status;
	double xlo,xhi;
	double r;


	// set up function parameters
	
	keplerEquationParamsType pkeq={ ecc * c / (c + b_ISO), th3 };

	//set up function for rootfinding

	keq.function = &keplerEquation;
	keq.params = &pkeq;

	//initialize rootfinder

	s=gsl_root_fsolver_alloc(T);
	xlo = 0.0;
	xhi = 2.0 * M_PI;

	gsl_root_fsolver_set(s, &keq, xlo, xhi);

	//iterate till root is solved
	
	do {
		iter++;
		status=gsl_root_fsolver_iterate(s);
		xlo=gsl_root_fsolver_x_lower(s);
		xhi=gsl_root_fsolver_x_upper(s);
		status=gsl_root_test_interval(xlo,xhi,epsabs,epsrel);
	} while (status==GSL_CONTINUE && iter<max_iter);
	
	r = gsl_root_fsolver_root(s);
	
	gsl_root_fsolver_free(s);
	
	return r;

						
}

//set up kepler equation for root finding

double keplerEquation(double eta, void *p) {
	keplerEquationParamsType * params = (keplerEquationParamsType *)p;
	double A = (params->A);
	double theta = (params->theta);
	
	return eta - A * sin(eta) - theta;
	
}
