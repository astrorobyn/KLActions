//header file for Amina's action calculator

struct ActionsFrequenciesAndDerivatives
{
	double Jr, Jtheta, Jphi;
	double Thetar, Thetatheta, Thetaphi;
	double Omegar, Omegatheta, Omegaphi;
	double dOmegardJr, dOmegardJtheta, dOmegardJphi;
	double dOmegathetadJr, dOmegathetadJtheta, dOmegathetadJphi;
	double dOmegaphidJr, dOmegaphidJtheta, dOmegaphidJphi;
	double apo, peri;
};

ActionsFrequenciesAndDerivatives computeActionsFrequenciesAndDerivatives(double E, double L, double LZ, double R, double VR, double GM, double bISO)
