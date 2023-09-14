/*W and Z PDFs in effective vector approximation as a function of t and x. Notation for the Ws: the two letters following fw are the electric charge and the polarization*/
double fwmpeva(double t, double x){
	double f;
	f = alpha2(t)*(1-x)*(1-x)*(log((m0*m0*exp(t)+(1-x)*mw*mw)/(m0*m0+(1-x)*mw*mw))-m0*m0*exp(t)/(m0*m0*exp(t)+(1-x)*mw*mw))/(8*pi*x);
	return f;
}

double fwmmeva(double t, double x){
	double f;
	f = alpha2(t)*(log((m0*m0*exp(t)+(1-x)*mw*mw)/(m0*m0+(1-x)*mw*mw))-m0*m0*exp(t)/(m0*m0*exp(t)+(1-x)*mw*mw))/(8*pi*x);
	return f;
}

double fzpeva(double t, double x){
	double f;
	f = alpha2(t)*(pow(0.5-sw(t)*sw(t),2)*(1-x)*(1-x)/x+pow(sw(t),4)/x)*(log((m0*m0*exp(t)+(1-x)*mz*mz)/(m0*m0+(1-x)*mz*mz))-m0*m0*exp(t)/(m0*m0*exp(t)+(1-x)*mz*mz))/(4*pi*pow(cw(t),2));
	return f;
}

double fzmeva(double t, double x){
	double f;
	f = alpha2(t)*0.5*(pow(0.5-sw(t)*sw(t),2)/x+pow(sw(t),4)*(1-x)*(1-x)/x)*(log((m0*m0*exp(t)+(1-x)*mz*mz)/(m0*m0+(1-x)*mz*mz))-m0*m0*exp(t)/(m0*m0*exp(t)+(1-x)*mz*mz))/(4*pi*pow(cw(t),2));
	return f;
}

double fwmleva(double t, double x){
	double f;
	f = alpha2(t)*(1-x)*m0*m0*exp(t)/(4*pi*(m0*m0*exp(t)+(1-x)*mw*mw)*x);
	return f;
}

double fzleva(double t, double x){
	double f;
	f = alpha2(t)*(0.25-sw(t)*sw(t)+2*pow(sw(t),4))*(1-x)*m0*m0*exp(t)/(2*pi*(m0*m0*exp(t)+(1-x)*mz*mz)*x);
	return f;
}
