ActionsFrequenciesAndDerivatives computeActionsFrequenciesAndDerivatives(double E, double L, double LZ, double R, double VR, double GM, double bISO)
{
	ActionsFrequenciesAndDerivatives result;
    double term = L * L + 4.0 * GM * bISO;

	// Actions
    result.Jphi = LZ;
    result.Jtheta = L - abs(LZ);
    result.Jr = GM / sqrt(-2.0 * E) - 1.0 / 2.0 * (L + sqrt(term));

	// Stuff needed for angles
	double c = GM / (-2.0 * E) - bISO;
	double ee = sqrt(1.0 - L * L / (GM * c) * (1 + bISO / c));
	//double ee = sqrt(1.0 + 2.0 * E * L * L / pow(2.0 * bISO * E + GM, 2.0));
	double s = sqrt(1.0 + R * R / (bISO * bISO)) + 1.0;
	double coseta = 1.0 / ee * (1.0 - bISO / c * (s - 2.0));
	double eta = acos(coseta);

	// Angles
	if (VR < 0)
		eta = 2.0 * M_PI - eta;
		
    result.Thetar = eta - ee * c / (c + bISO) * sin(eta);
    result.Thetatheta = 0;
    result.Thetaphi = 0;

	// Frequencies
    result.Omegar = pow(-2.0 * E, 3.0 / 2.0 ) / GM;
    result.Omegatheta = 1.0 / 2.0 * result.Omegar * (1.0 + L / sqrt(term));
    result.Omegaphi = result.Omegatheta * copysign(1.0, result.Jphi);

	// Derivatives of frequencies w.r.t. actions
	result.dOmegardJr     	  = 3.0 / (2.0 * E) * result.Omegar * result.Omegar;
	result.dOmegardJtheta 	  = 3.0 / (2.0 * E) * result.Omegar * result.Omegatheta;
	result.dOmegardJphi   	  = 3.0 / (2.0 * E) * result.Omegar * result.Omegaphi;
	result.dOmegathetadJr 	  = 3.0 / (2.0 * E) * result.Omegar * result.Omegatheta;
	result.dOmegathetadJtheta = 3.0 / (2.0 * E) * result.Omegatheta * result.Omegatheta + result.Omegar * 2 * GM * bISO / pow(term, 3.0/2.0);
	result.dOmegathetadJphi   = result.dOmegathetadJtheta * copysign(1.0, result.Jphi);
	result.dOmegaphidJr       = 3.0 / (2.0 * E) * result.Omegar * result.Omegaphi;
	result.dOmegaphidJtheta   = result.dOmegathetadJtheta * copysign(1.0, result.Jphi);
	result.dOmegaphidJphi     = result.dOmegathetadJtheta;

	double tmp1 = 2 * bISO * E * GM + GM * GM + E * L * L;
	double tmp2 = sqrt(GM * GM * (pow(2 * bISO * E + GM, 2)+2 * E * L * L));
	result.peri = sqrt((tmp1 - tmp2) / (2 * E * E));
	result.apo  = sqrt((tmp1 + tmp2) / (2 * E * E));

    return result;

}

