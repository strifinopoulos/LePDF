/*One loop running of the elctromagnetic coupling alpha=e^2/4pi in the first phase of the evolution*/
double alphaQED(double t){
	double alpha;
	alpha = aem0/(1-(b0QED*aem0*(t-t0)/2));
	return alpha;
}

/*One loop running of the U(1)_Y coupling alpha1=g1^2/4pi*/
double alpha1(double t){
	double alpha;
	alpha = a10/(1-(b01*a10*(t-t0r)/2));
	return alpha;
}

/*One loop running of the SU(2)_L coupling alpha2=g2^2/4pi*/
double alpha2(double t){
	double alpha;
	alpha = a20/(1+(b02*a20*(t-t0r)/2));
	return alpha;
}

/*This function linearly interpolates between the points in the positions i and j of the nx2 matrix A*/
double interpolation(double t, double **A, int i, int j, int n){
	double x;
	if(i>n-1){
		printf("The first point for interpolation is out of the array\n");
		return EXIT_FAILURE;
	}
	if(j>n-1){
		printf("The second point for interpolation is out of the array\n");
		return EXIT_FAILURE;
	}
	else{
		x = (A[j][1]-A[i][1])*t/(A[j][0]-A[i][0]) + (A[i][1]*A[j][0]-A[j][1]*A[i][0])/(A[j][0]-A[i][0]);
	}
	return x;
}

/*Returns alpha_s after the order 1 interpolation of the matrix alpha*/
double alphas(double t, double **alpha){
	if(t>alpha[Na-1][0]){ /*The values considered must be within the maximum scale reported in the file */
		printf("Not enough points in the alpha_s file. Try lower energies\n");
		return EXIT_FAILURE;
	}
	double a;
	if(t<tQCD){ /*Below tQCD we do not consider QCD and we set alpha_s = 0*/
		a= 0;
	}
	else{
		int ti;
		double tf;
		tf = (t-alpha[0][0])/dta;
		ti = (int)tf;
		a = interpolation(t,alpha,ti,ti+1,Na);
	}
	return a;
}

/*Running of the Yukawa coupling alphay=yt^2/4pi*/
double alphay(double t, double **alpha){
	double b;
	b = double(9)/(2*alphas(t,alpha))-(double(9)/(2*alphas(t0r,alpha))-1/ay0)*pow(alphas(t0r,alpha)/alphas(t,alpha),8/7);
	return 1/b;
}

/*cos(thetaW)*/
double cw(double t){
	return sqrt(alpha2(t)/(alpha1(t)+alpha2(t)));
}

/*sin(thetaW)*/
double sw(double t){
	return sqrt(alpha1(t)/(alpha1(t)+alpha2(t)));
}

/*tan(thetaW)*/
double tw(double t){
	return sw(t)/cw(t);
}

/*cos(2thetaW)*/
double c2w(double t){
	return cw(t)*cw(t)-sw(t)*sw(t);
}

/*electromagnetic coupling in the SM*/
double alphaem(double t){
	double alpha;
	alpha = alpha2(t)*sw(t)*sw(t);
	return alpha;
}

/*Z coupling*/
double alphaz(double t){
	double alpha;
	alpha = alpha2(t)/pow(cw(t),2);
	return alpha;
}

/*Z charges*/
double Qzul(double t, double Q){ /*for up quarks and neutrinos*/
	double qz;
	qz = 0.5-Q*sw(t)*sw(t);
	return qz;
}

double Qzdl(double t, double Q){ /*for down quarks and left leptons*/
	double qz;
	qz = -0.5-Q*sw(t)*sw(t);
	return qz;
}

double Qzr(double t, double Q){ /*for right fermions*/
	double qz;
	qz = -Q*sw(t)*sw(t);
	return qz;
}

/*Mixed coupling for the Z-photon inteference*/
double alpham(double t){
	return sqrt(alphaem(t)*alpha2(t));
}
