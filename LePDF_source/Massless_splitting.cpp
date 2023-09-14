/*splitting functions. Notations: Pxy is the splitting y \to x, f = fermion, v = vector, h = scalar. If Pxyz, z is the kind of interaction, g = gauge, y = yukawa*/
double Pffg(double z){ 
	double a;
	a = (1+pow(z,2))/(1-z);
	return a;
}

double Pvpfl(double z){ 
	double a;
	a = pow(1-z,2)/z;
	return a;
}

double Pvmfl(double z){ 
	double a;
	a = 1/z;
	return a;
}

double Pvf(double z){ 
	double a;
	a = Pvpfl(z)+Pvmfl(z);
	return a;
}

double Pflvp(double z){ 
	double a;
	a = pow(1-z,2);
	return a;
}

double Pflvm(double z){ 
	double a;
	a = pow(z,2);
	return a;
}

double Pfv(double z){ 
	double a;
	a = Pflvp(z)+Pflvm(z);
	return a;
}

double Pvpvp(double z){ 
	double a;
	a = (1+pow(z,4))/(z*(1-z));
	return a;
}

double Pvpvm(double z){ 
	double a;
	a = pow(1-z,3)/z;
	return a;
}

double Pvv(double z){ 
	double a;
	a = Pvpvp(z)+Pvpvm(z);
	return a;
}

double Phh(double z){ 
	double a;
	a = 2*z/(1-z);
	return a;
}

double Pvh(double z){ 
	double a;
	a = (1-z)/z;
	return a;
}

double Phv(double z){ 
	double a;
	a = z*(1-z);
	return a;
}

double Pffy(double z){ 
	double a;
	a = (1-z)/2;
	return a;
}

double Phf(double z){ 
	double a;
	a = z/2;
	return a;
}

double Pfh(double z){ 
	double a;
	a = 0.5;
	return a;
}
