/*
 *  isochroneMap.c
 *  
 *
 *  Created by Robyn Sanderson on 7/6/10.
 *
 */


#include "isochroneMap.h"

static double M_ISO=12.2;
static double B_ISO=8.0;

//set up kepler equation for root finding

double keplerEquation(double eta, void *p) {
	keplerEquationParamsType * params = (keplerEquationParamsType *)p;
	double A = (params->A);
	double theta = (params->theta);
	
	return eta - A * sin(eta) - theta;
	
}

// auxiliary functions

//***NOTE the definitions below:
//  q = (qr, qphi);
//  qdot = (qrdot, qphidot);
//	x = (r, phi);
//	v = (rdot, phidot);
//	G = 1;
//	energies/angular momenta are SPECIFIC (per mass);
//***



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


double eccentricity(double q[2], double qdot[2]) {
	double J2 = Lmag(q,qdot);
	double c = c_aux(q,qdot);
	
	return sqrt(1.0 - J2 * J2 / (M_ISO * c) * (1.0 + B_ISO/c));
}

double c_aux(double q[2], double qdot[2]) {

	return M_ISO / (-2.0 * H(q,qdot)) - B_ISO;

}

double Omega2(double q[2], double qdot[2]) {
	
	double J2 = Lmag(q,qdot);

	return Omega3(q,qdot) * (1 + J2 / sqrt(J2 * J2 + 4.0 * M_ISO * B_ISO));
	
}

double Omega3(double q[2], double qdot[2]) {
	
	double ham = H(q,qdot);
	
	return sqrt(-8.0 * ham * ham * ham) / M_ISO;
	
}

double eta(double q[2], double qdot[2], double t, double eta0) {

	// calculates the auxiliary variable eta by solving Kepler's equation
	// the equation is solved by searching for a root
	// using the brent's algorithm root bracketer in GSL
	
	// definitions for root solver
	
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
	
	//compute e, c, omega3, initial theta3
	
	double e = eccentricity(q, qdot);
	double c = c_aux(q, qdot);
	double omega3 = Omega3(q,qdot);
	double theta30 = theta3(eta0,q,qdot);
	
	// set up function parameters
	
	keplerEquationParamsType pkeq={ e * c / (c + B_ISO), omega3 * t + theta30};

	//set up function for rootfinding

	keq.function = &keplerEquation;
	keq.params = &pkeq;

	//initialize rootfinder

	s=gsl_root_fsolver_alloc(T);
	
	//set it to bracket root of kepler eqn.  The initial guess should be near eta=Omega3*t + theta30
	
	xlo = omega3 * t + theta30 - M_PI;
	xhi = omega3 * t + theta30 + M_PI;

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

double etaZero(double q[2], double qdot[2]) {

	double e = eccentricity(q, qdot);
	double c = c_aux(q, qdot);
	
	double coseta = (1.0 - B_ISO / c * (sqrt(1.0 + q[0]*q[0]/(B_ISO * B_ISO))-1.0))/e;

	return acos(coseta);
	
}

double theta3(double etaval, double q[2], double qdot[2]) {
	
	double e = eccentricity(q, qdot);
	double c = c_aux(q, qdot);
	
	return etaval - e * c / (c + B_ISO) * sin(etaval);
	
}

double theta2(double etaval, double phi, double q[2], double qdot[2]) {

	double e = eccentricity(q, qdot);
	double c = c_aux(q, qdot);
	double J2 = Lmag(q, qdot);
	double omega2 = Omega2(q,qdot);
	double omega3 = Omega3(q,qdot);
	double th3 = theta3(etaval, q, qdot);

	double tan1 = atan(sqrt((1.0 + e)/(1.0 - e)) * tan(etaval/2.0)) +  M_PI * floor((etaval+M_PI)/(2.0 * M_PI));
	double tan2 = (atan(sqrt((1.0 + e + 2.0 * B_ISO / c)/(1.0 - e + 2.0 * B_ISO / c)) * tan(etaval/2.0)) + M_PI * floor((etaval+M_PI)/(2.0 * M_PI)))/sqrt(1 + 4.0 * M_PI * B_ISO / (J2 * J2));
	
	return phi + omega2/omega3 * th3 - tan1 - tan2;
	
}


void forwardMap(double q[2], double qdot[2], double t, double x[2], double v[2]) {
	
	// maps from (q, qdot) -> (x(t), v(t)) (forward in time)
	// to go backward, let t -> -t
	
	// first step: find auxiliary variables
	
	// eta0 is between 0 and pi but extends to between 0 and 2pi based on the sign of the initial r-velocity


	double eta0 = etaZero(q,qdot);
	

	if (qdot[0]<0.0) {
		eta0 = 2.0 * M_PI - eta0;
	}
	
	// etabase is between 0 and pi
	double etabase = eta(q,qdot,t, eta0);

	double omega3 = Omega3(q,qdot);
	double omega2 = Omega2(q,qdot);
	double J2 = Lmag(q,qdot);
	double ham = H(q,qdot);
	double e = eccentricity(q, qdot);
	double c = c_aux(q, qdot);

	
	// compute number of pis to add for correct branch (theta2 and theta3 must increase continuously)
	
	double Tr = 2.0 * M_PI / omega3;
	double tperi = (theta3(0.0, q, qdot) - theta3(eta0,q, qdot))/omega3;
	double tapo = (theta3(M_PI, q, qdot) - theta3(eta0,q, qdot))/omega3;
	
	// both these should be larger than 0 (time of *next* pericenter and apocenter)
	
	while (tperi<0) {
		tperi += Tr;
	}
	while (tapo<0) {
		tapo += Tr;
	}

	// calculate number of half-periods (no. of peri- or apocenter passages)
	int nstar = ((int) floor((t - GSL_MIN(tapo, tperi))/(Tr/2.0))) + 1;
	
	
	// get r
	
	double base = (1.0 + c/B_ISO * (1.0 - e * cos(etabase)));
	x[0] = B_ISO * sqrt(base * base - 1.0);
	
	// get phi
	
	double tan1 = atan(sqrt((1.0 + e)/(1.0 - e)) * tan(etabase/2.0)) + M_PI * floor((etabase+M_PI)/(2.0 * M_PI));
	double tan2 = (atan(sqrt((1.0 + e + 2.0 * B_ISO / c)/(1.0 - e + 2.0 * B_ISO / c)) * tan(etabase/2.0)) + M_PI * floor((etabase+M_PI)/(2.0 * M_PI)))/sqrt(1 + 4.0 * M_PI * B_ISO / (J2 * J2));


	double th20 = theta2(eta0, q[1], q, qdot);
	double th3 = theta3(etabase, q, qdot);
	
//	fprintf(stdout, "%lg %lg %lg %lg\n", t, tan1, tan2,th3);
	
	x[1] = th20 + omega2 * t - omega2/omega3 * th3 + tan1 + tan2;
		
	// conservation of angular momentum gives phidot
	
	v[1] = J2 / (x[0] * x[0]);
	
	// conservation of energy gives magnitude of rdot
	
	v[0] = M_SQRT2 * sqrt(ham - J2 * J2 / (2.0 * x[0] * x[0]) - Phi(x[0]));
	
	// figure out the sign - if it started heading toward pericenter 
	// and has gone thru an even number of half-periods since then, vr<0
	// likewise if it started heading toward apo and has gone thru an odd nm
	// of half-periods, vr<0
	if((GSL_IS_EVEN(nstar) && (tperi<tapo)) || (GSL_IS_ODD(nstar) && (tapo<tperi)))
		v[0] *= -1.0;
	   
}

double Jr(double q[2], double qdot[2]) {
	
	double ham = H(q,qdot);
	double J2 = Lmag(q, qdot);
	
	return M_ISO / sqrt(-2.0 * ham) - 0.5 * (J2 + sqrt(J2 * J2 + 4.0 * M_ISO * B_ISO));
	
}
