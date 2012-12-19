//these hold the mass and scale radius of the potential

#include "isochroneFunctions.h"

double M_ISO, B_ISO;


//functions for energy, angular momentum, actions

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


double Jr(double q[2], double qdot[2]) {
	
	double ham = H(q,qdot);
	double J2 = Lmag(q, qdot);
	
	return M_ISO / sqrt(-2.0 * ham) - 0.5 * (J2 + sqrt(J2 * J2 + 4.0 * M_ISO * B_ISO));
	
}


//functions for frequencies

double Omega2(double q[2], double qdot[2]) {
	
	double J2 = Lmag(q,qdot);

	return Omega3(q,qdot) * (1 + J2 / sqrt(J2 * J2 + 4.0 * M_ISO * B_ISO));
	
}

double Omega3(double q[2], double qdot[2]) {
	
	double ham = H(q,qdot);
	
	return sqrt(-8.0 * ham * ham * ham) / M_ISO;
}

//auxiliary functions for the angles

double eccentricity(double q[2], double qdot[2]) {
	double J2 = Lmag(q,qdot);
	double c = c_aux(q,qdot);
	
	return sqrt(1.0 - J2 * J2 / (M_ISO * c) * (1.0 + B_ISO/c));
}

double c_aux(double q[2], double qdot[2]) {

	return M_ISO / (-2.0 * H(q,qdot)) - B_ISO;

}


double etaZero(double q[2], double qdot[2]) {

	double e = eccentricity(q, qdot);
	double c = c_aux(q, qdot);
	
	double coseta = (1.0 - B_ISO / c * (sqrt(1.0 + q[0]*q[0]/(B_ISO * B_ISO))-1.0))/e;

	double etamag = acos(coseta);

	if(qdot[0]<0)
	  return 2.0 * M_PI - etamag;
	else
	  return etamag;
	
}


//functions for the angles 

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

double theta1(double J1, double J2, double q[3]) {

  double r = norm(q);
  double ph = phi(q);
  double z = q[2];
 
  double sinu = J1 * z / sqrt((J2*J2-J1*J1)*(r*r-z*z));

  return ph - asin(sinu); 

}
