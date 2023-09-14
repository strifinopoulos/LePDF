/*correction to the usual splitting functions due to the mass*/
double prop2(double t, double ma, double mb, double mc, double z){
	double p;
	p = 1/pow(1+((1-z)*mb*mb+z*mc*mc-z*(1-z)*ma*ma)/(mu(t)*mu(t)),2);
	return p;
}

/*correction to the usual splitting functions due to the mass. It is the same as before but with no ^2: this has to be used for the mixed PDFs*/
double prop(double t, double ma, double mb, double mc, double z){
	double p;
	p = 1/(1+((1-z)*mb*mb+z*mc*mc-z*(1-z)*ma*ma)/(mu(t)*mu(t)));
	return p;
}

/*same but computed at z=1 (needed to regulate the 1/zb divergencies). Since only mc matters, this is valid also for the mixed pdfs*/
double regprop(double t, double mc){
	double p;
	p = 1/pow(1+mc*mc/(mu(t)*mu(t)),2);
	return p;
}
