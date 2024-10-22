/*Discretized QED-QCD DGLAP equations for the independent PDFs fe (sea leptons), fmu (valence lepton), fu, fd, fb (quarks), fph (photon) and fg (gluon)*/
/*Since we use the RUnge-Kutta algorithm, the discretized equations here reported are already in the form used to compute the k coefficients*/
double ke(int j, int k, double *x, double *dx, double t, double dt, double fe, double fph){
	double b;
	b = dt*(alphaQED(t)/(2*pi))*dx[k]*(Pffg(x[j]/x[k])*fe+Pfv(x[j]/x[k])*fph)/x[k];
	return b;
}

double ktau(int j, int k, double *x, double *dx, double t, double dt, double ftau, double fph){
	double b;
	b = theta(t,ttau)*dt*(alphaQED(t)/(2*pi))*dx[k]*(Pffg(x[j]/x[k])*ftau+Pfv(x[j]/x[k])*fph)/x[k];
	return b;
}

double kmu(int j, int k, double *x, double *dx, double t, double dt, double fm, double fph){
	double b;
	b = dt*(alphaQED(t)/(2*pi))*dx[k]*(Pffg(x[j]/x[k])*fm+Pfv(x[j]/x[k])*fph)/x[k];
	return b;
}

double ku(int j, int k, double *x, double *dx, double t, double dt, double fu, double fph, double fg, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*((alphaQED(t)/(2*pi))*Qu*Qu*(Pffg(x[j]/x[k])*fu+3*Pfv(x[j]/x[k])*fph)+(alphas(t,alpha)/(2*pi))*(Cf*Pffg(x[j]/x[k])*fu+0.5*Pfv(x[j]/x[k])*fg));
	return b;
}

double kc(int j, int k, double *x, double *dx, double t, double dt, double fc, double fph, double fg, double **alpha){
	double b;
	b = theta(t,tc)*dt*(dx[k]/x[k])*((alphaQED(t)/(2*pi))*Qu*Qu*(Pffg(x[j]/x[k])*fc+3*Pfv(x[j]/x[k])*fph)+(alphas(t,alpha)/(2*pi))*(Cf*Pffg(x[j]/x[k])*fc+0.5*Pfv(x[j]/x[k])*fg));
	return b;
}

double kd(int j, int k, double *x, double *dx, double t, double dt, double fd, double fph, double fg, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*((alphaQED(t)/(2*pi))*Qd*Qd*(Pffg(x[j]/x[k])*fd+3*Pfv(x[j]/x[k])*fph)+(alphas(t,alpha)/(2*pi))*(Cf*Pffg(x[j]/x[k])*fd+0.5*Pfv(x[j]/x[k])*fg));
	return b;
}

double kb(int j, int k, double *x, double *dx, double t, double dt, double fb, double fph, double fg, double **alpha){
	double b;
	b = theta(t,tb)*dt*(dx[k]/x[k])*((alphaQED(t)/(2*pi))*Qd*Qd*(Pffg(x[j]/x[k])*fb+3*Pfv(x[j]/x[k])*fph)+(alphas(t,alpha)/(2*pi))*(Cf*Pffg(x[j]/x[k])*fb+0.5*Pfv(x[j]/x[k])*fg));
	return b;
}

double kph(int j, int k, double *x, double *dx, double t, double dt, double fe, double fm, double fu, double fd, double fb){
	double b;
	b = dt*(alphaQED(t)/(2*pi))*dx[k]*Pvf(x[j]/x[k])*(5*fe+fm+4*Qu*Qu*fu+4*Qd*Qd*fd+2*Qd*Qd*fb)/x[k];
	return b;
}

double kg(int j, int k, double *x, double *dx, double t, double dt, double fu, double fd, double fg, double fb, double **alpha){
	double b;
	b = dt*(alphas(t,alpha)/(2*pi))*dx[k]*(Pvf(x[j]/x[k])*Cf*(4*fu+4*fd+2*fb)+Ca*Pvv(x[j]/x[k])*fg)/x[k];
	return b;
}

/*Function to update the matrices with the evolution coefficients with x_k contribution*/
double **RK0(int i, int j, int k, int Nx, double t, double dt, double *grid, double *deltagrid, double fm, double fu, double fd, double fph, double fg, double fe, double fb, double fc, double ftau, double **K, double **alpha){
	K[i][j]+=kmu(j,k,grid,deltagrid,t,dt,fm,fph);
	K[i][j+Nx]+=ku(j,k,grid,deltagrid,t,dt,fu,fph,fg,alpha);
	K[i][j+2*Nx]+=kd(j,k,grid,deltagrid,t,dt,fd,fph,fg,alpha);
	K[i][j+3*Nx]+=kph(j,k,grid,deltagrid,t,dt,fe,fm,fu,fd,fb);
	K[i][j+4*Nx]+=kg(j,k,grid,deltagrid,t,dt,fu,fd,fg,fb,alpha);
	K[i][j+5*Nx]+=ke(j,k,grid,deltagrid,t,dt,fe,fph);
	K[i][j+6*Nx]+=kb(j,k,grid,deltagrid,t,dt,fb,fph,fg,alpha);
	K[i][j+7*Nx]+=kc(j,k,grid,deltagrid,t,dt,fc,fph,fg,alpha);
	K[i][j+8*Nx]+=ktau(j,k,grid,deltagrid,t,dt,ftau,fph);
	return K;
}

/*Function to compute the term X_j*/
double newterm(double *x, double *dx, int i, int n){ 
	int j;
	double s;
	s = 0;
	for(j = i+1;j<n;j++){
		s+=dx[j]*x[i]/(x[j]*(x[j]-x[i]));
	}
	return s;
}

/*Effective QED number of fermions. It depends on t due to the thresholds*/
double NfQED(double t){
	double nf;
	if(t<tc){
		nf = 4;
	}
	else if(t>tc && t<ttau){
		nf = 16/3;
	}
	else if(t>ttau && t<tb){
		nf = 19/3;
	}
	else{
		nf=20/3;
	}
	return nf;
}


/*Function to update the matrices with the evolution coefficients with x_j contribution*/
double **RKj0(int i, int j, int Nx, double t, double dt, double *grid, double *deltagrid, double fm, double fu, double fd, double fph, double fg, double fe, double fb, double fc, double ftau, double **K, double a, double **alpha){
	K[i][j] += (1.5+2*(log(1-grid[j])-a))*dt*(alphaQED(t)/(2*pi))*fm;
	K[i][j+Nx] += (1.5+2*(log(1-grid[j])-a))*dt*((alphaQED(t)/(2*pi))*Qu*Qu+Cf*(alphas(t,alpha)/(2*pi)))*fu;
	K[i][j+2*Nx] += (1.5+2*(log(1-grid[j])-a))*dt*((alphaQED(t)/(2*pi))*Qd*Qd+Cf*(alphas(t,alpha)/(2*pi)))*fd;
	K[i][j+3*Nx] += -2*NfQED(t)*dt*(alphaQED(t)/(2*pi))*fph/3; 
	K[i][j+4*Nx] += ((25-2*theta(t,tb))/6+6*(log(1-grid[j])-a))*dt*(alphas(t,alpha)/(2*pi))*fg;
	K[i][j+5*Nx] += (1.5+2*(log(1-grid[j])-a))*dt*(alphaQED(t)/(2*pi))*fe;
	K[i][j+6*Nx] += theta(t,tb)*(1.5+2*(log(1-grid[j])-a))*dt*((alphaQED(t)/(2*pi))*Qd*Qd+Cf*(alphas(t,alpha)/(2*pi)))*fb;
	K[i][j+7*Nx] += theta(t,tc)*(1.5+2*(log(1-grid[j])-a))*dt*((alphaQED(t)/(2*pi))*Qu*Qu+Cf*(alphas(t,alpha)/(2*pi)))*fc;
	K[i][j+8*Nx] += theta(t,ttau)*(1.5+2*(log(1-grid[j])-a))*dt*(alphaQED(t)/(2*pi))*ftau;
	return K;
}

/*Imposition of momentum conservation*/
double Lt0(double *x, double *dx, int n, int i, double **fe, double **fm, double **ftau, double **fu, double **fc, double **fd, double **fph, double **fg, double **fb){
	double L;
	L = 1;
	int j;
	for(j=0;j<n-1;j++){
		L-=x[j]*dx[j]*(3*fe[i][j]+2*ftau[i][j]+fm[i][j]+2*fu[i][j]+2*fc[i][j]+4*fd[i][j]+2*fb[i][j]+fph[i][j]+fg[i][j]);
	}
	return L;
}
