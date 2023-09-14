/*Introduces the bound of integration in the DGLAP equations*/
double cut(double x, double z, double Q){ 
	double lambda;
	lambda=x*(1+(mew/Q)); /*lambda = lower bound of integration in the DGLAP equations*/
	if(z>=lambda){
		return 1; /*we consider only values of z above the cut...*/
	}
	else{
		return 0; /*...otherwise the function gives zero*/
	}
}

double theta(double x, double y){ /*theta function*/
	if(x<y){
		return 0;
	}
	else{
		return 1;
	}
}
