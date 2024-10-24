/*Evolution for neutrinos*/
double knu(int j, int k, double *x, double *dx, double t, double dt, double fn, double fe, double fzp, double fzm, double fwpp, double fwpm){
	double b;
	b = dt*(dx[k]/x[k])*(theta(t,tz)*alphaz(t)*0.25*Sffv(t,j,k,x,fn,fzp,fzm,1,1,mz,0)+alpha2(t)*0.5*Sffvcut(t,j,k,x,fe,fwpp,fwpm,1,0,0))/(2*pi);
	return b;
}

double knub(int j, int k, double *x, double *dx, double t, double dt, double fnb, double feb, double fzp, double fzm, double fwmp, double fwmm){
	double b;
	b = dt*(dx[k]/x[k])*(theta(t,tz)*alphaz(t)*0.25*Sffv(t,j,k,x,fnb,fzm,fzp,1,1,mz,0)+alpha2(t)*0.5*Sffvcut(t,j,k,x,feb,fwmm,fwmp,1,0,0))/(2*pi);
	return b;
}

/*Evolution for left handed leptons (same expression for el, ml and tal)*/
double kel(int j, int k, double *x, double *dx, double t, double dt, double fe, double fn, double fphp, double fphm, double fzp, double fzm, double fwmp, double fwmm, double fzphp, double fzphm){
	double b;
	b = dt*(dx[k]/x[k])*(alphaem(t)*Ql*Ql*Sffv(t,j,k,x,fe,fphp,fphm,1,1,0,0)+theta(t,tz)*alphaz(t)*pow(Qzdl(t,Ql),2)*Sffv(t,j,k,x,fe,fzp,fzm,1,1,mz,0)+alpha2(t)*0.5*Sffvcut(t,j,k,x,fn,fwmp,fwmm,1,0,0)+theta(t,tz)*(alpham(t)/cw(t))*Ql*Qzdl(t,Ql)*Sffzgamma(t,j,k,x,fzphp,fzphm,0))/(2*pi);
	return b;
}

double kelb(int j, int k, double *x, double *dx, double t, double dt, double feb, double fnb, double fphp, double fphm, double fzp, double fzm, double fwpp, double fwpm, double fzphp, double fzphm){
	double b;
	b = dt*(dx[k]/x[k])*(alphaem(t)*Ql*Ql*Sffv(t,j,k,x,feb,fphm,fphp,1,1,0,0)+theta(t,tz)*alphaz(t)*pow(Qzdl(t,Ql),2)*Sffv(t,j,k,x,feb,fzm,fzp,1,1,mz,0)+alpha2(t)*0.5*Sffvcut(t,j,k,x,fnb,fwpm,fwpp,1,0,0)+theta(t,tz)*(alpham(t)/cw(t))*Ql*Qzdl(t,Ql)*Sffzgamma(t,j,k,x,fzphm,fzphp,0))/(2*pi);
	return b;
}

/*Evolution for left handed up quarks (u,c)*/
double kul(int j, int k, double *x, double *dx, double t, double dt, double fu, double fd, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fwpp, double fwpm, double fzphp, double fzphm, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(alphas(t,alpha)*Sffv(t,j,k,x,fu,fgp,fgm,Cf,Tf,0,0)+alphaem(t)*Qu*Qu*Sffv(t,j,k,x,fu,fphp,fphm,1,3,0,0)+theta(t,tz)*alphaz(t)*pow(Qzul(t,Qu),2)*Sffv(t,j,k,x,fu,fzp,fzm,1,3,mz,0)+alpha2(t)*0.5*Sffvcut(t,j,k,x,fd,fwpp,fwpm,3,0,0)+3*theta(t,tz)*(alpham(t)/cw(t))*Qu*Qzul(t,Qu)*Sffzgamma(t,j,k,x,fzphp,fzphm,0))/(2*pi);
	return b;
}

double kulb(int j, int k, double *x, double *dx, double t, double dt, double fub, double fdb, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fwmp, double fwmm, double fzphp, double fzphm, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(alphas(t,alpha)*Sffv(t,j,k,x,fub,fgm,fgp,Cf,Tf,0,0)+alphaem(t)*Qu*Qu*Sffv(t,j,k,x,fub,fphm,fphp,1,3,0,0)+theta(t,tz)*alphaz(t)*pow(Qzul(t,Qu),2)*Sffv(t,j,k,x,fub,fzm,fzp,1,3,mz,0)+alpha2(t)*0.5*Sffvcut(t,j,k,x,fdb,fwmm,fwmp,3,0,0)+3*theta(t,tz)*(alpham(t)/cw(t))*Qu*Qzul(t,Qu)*Sffzgamma(t,j,k,x,fzphm,fzphp,0))/(2*pi);
	return b;
}

/*Evolution for left handed down quarks (d,s)*/
double kdl(int j, int k, double *x, double *dx, double t, double dt, double fu, double fd, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fwmp, double fwmm, double fzphp, double fzphm, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(alphas(t,alpha)*Sffv(t,j,k,x,fd,fgp,fgm,Cf,Tf,0,0)+alphaem(t)*Qd*Qd*Sffv(t,j,k,x,fd,fphp,fphm,1,3,0,0)+theta(t,tz)*alphaz(t)*pow(Qzdl(t,Qd),2)*Sffv(t,j,k,x,fd,fzp,fzm,1,3,mz,0)+alpha2(t)*0.5*Sffvcut(t,j,k,x,fu,fwmp,fwmm,3,0,0)+3*theta(t,tz)*(alpham(t)/cw(t))*Qd*Qzdl(t,Qd)*Sffzgamma(t,j,k,x,fzphp,fzphm,0))/(2*pi);
	return b;
}

double kdlb(int j, int k, double *x, double *dx, double t, double dt, double fub, double fdb, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fwpp, double fwpm, double fzphp, double fzphm, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(alphas(t,alpha)*Sffv(t,j,k,x,fdb,fgm,fgp,Cf,Tf,0,0)+alphaem(t)*Qd*Qd*Sffv(t,j,k,x,fdb,fphm,fphp,1,3,0,0)+theta(t,tz)*alphaz(t)*pow(Qzdl(t,Qd),2)*Sffv(t,j,k,x,fdb,fzm,fzp,1,3,mz,0)+alpha2(t)*0.5*Sffvcut(t,j,k,x,fub,fwpm,fwpp,3,0,0)+3*theta(t,tz)*(alpham(t)/cw(t))*Qd*Qzdl(t,Qd)*Sffzgamma(t,j,k,x,fzphm,fzphp,0))/(2*pi);
	return b;
}

/*Evolution for left handed top*/
double ktl(int j, int k, double *x, double *dx, double t, double dt, double ftl, double ftr, double fb, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fwpp, double fwpm, double fzphp, double fzphm, double fh, double fzl, double fhzl, double **alpha){
	double b;
	b = theta(t,tt)*dt*(dx[k]/x[k])*(alphay(t,alpha)*0.5*Stly(t,j,k,x,ftr,fh,fzl,fhzl)+alphas(t,alpha)*Sffv(t,j,k,x,ftl,fgp,fgm,Cf,Tf,0,mt)+alphaem(t)*Qu*Qu*Sffv(t,j,k,x,ftl,fphp,fphm,1,3,0,mt)+theta(t,tz)*alphaz(t)*pow(Qzul(t,Qu),2)*Sffv(t,j,k,x,ftl,fzp,fzm,1,3,mz,mt)+alpha2(t)*0.5*Sffvcut(t,j,k,x,fb,fwpp,fwpm,3,mt,0)+3*theta(t,tz)*(alpham(t)/cw(t))*Qu*Qzul(t,Qu)*Sffzgamma(t,j,k,x,fzphp,fzphm,mt))/(2*pi);
	return b;
}

double ktlb(int j, int k, double *x, double *dx, double t, double dt, double ftlb, double ftrb, double fbb, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fwmp, double fwmm, double fzphp, double fzphm, double fh, double fzl, double fhzl, double **alpha){
	double b;
	b = theta(t,tt)*dt*(dx[k]/x[k])*(alphay(t,alpha)*0.5*Stly(t,j,k,x,ftrb,fh,fzl,-fhzl)+alphas(t,alpha)*Sffv(t,j,k,x,ftlb,fgm,fgp,Cf,Tf,0,mt)+alphaem(t)*Qu*Qu*Sffv(t,j,k,x,ftlb,fphm,fphp,1,3,0,mt)+theta(t,tz)*alphaz(t)*pow(Qzul(t,Qu),2)*Sffv(t,j,k,x,ftlb,fzm,fzp,1,3,mz,mt)+alpha2(t)*0.5*Sffvcut(t,j,k,x,fbb,fwmm,fwmp,3,mt,0)+3*theta(t,tz)*(alpham(t)/cw(t))*Qu*Qzul(t,Qu)*Sffzgamma(t,j,k,x,fzphm,fzphp,mt))/(2*pi);
	return b;
}

/*Evolution for left handed bottom*/
double kbl(int j, int k, double *x, double *dx, double t, double dt, double ftl, double ftr, double fb, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fwmp, double fwmm, double fzphp, double fzphm, double fwml, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(theta(t,tt)*alphay(t,alpha)*Sbly(t,j,k,x,ftr,fwml)+alphas(t,alpha)*Sffv(t,j,k,x,fb,fgp,fgm,Cf,Tf,0,0)+alphaem(t)*Qd*Qd*Sffv(t,j,k,x,fb,fphp,fphm,1,3,0,0)+theta(t,tz)*alphaz(t)*pow(Qzdl(t,Qd),2)*Sffv(t,j,k,x,fb,fzp,fzm,1,3,mz,0)+theta(t,tt)*alpha2(t)*0.5*Sffvcut(t,j,k,x,ftl,fwmp,fwmm,3,0,mt)+3*theta(t,tz)*(alpham(t)/cw(t))*Qd*Qzdl(t,Qd)*Sffzgamma(t,j,k,x,fzphp,fzphm,0))/(2*pi);
	return b;
}

double kblb(int j, int k, double *x, double *dx, double t, double dt, double ftlb, double ftrb, double fbb, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fwpp, double fwpm, double fzphp, double fzphm, double fwpl, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(theta(t,tt)*alphay(t,alpha)*Sbly(t,j,k,x,ftrb,fwpl)+alphas(t,alpha)*Sffv(t,j,k,x,fbb,fgm,fgp,Cf,Tf,0,0)+alphaem(t)*Qd*Qd*Sffv(t,j,k,x,fbb,fphm,fphp,1,3,0,0)+theta(t,tz)*alphaz(t)*pow(Qzdl(t,Qd),2)*Sffv(t,j,k,x,fbb,fzm,fzp,1,3,mz,0)+theta(t,tt)*alpha2(t)*0.5*Sffvcut(t,j,k,x,ftlb,fwpm,fwpp,3,0,mt)+3*theta(t,tz)*(alpham(t)/cw(t))*Qd*Qzdl(t,Qd)*Sffzgamma(t,j,k,x,fzphm,fzphp,0))/(2*pi);
	return b;
}

/*Evolution for right handed leptons (same expression for er, mr and tar)*/
double ker(int j, int k, double *x, double *dx, double t, double dt, double fe, double fphp, double fphm, double fzp, double fzm, double fzphp, double fzphm){
	double b;
	b = dt*(dx[k]/x[k])*(alphaem(t)*Ql*Ql*Sffv(t,j,k,x,fe,fphm,fphp,1,1,0,0)+theta(t,tz)*alphaz(t)*pow(Qzr(t,Ql),2)*Sffv(t,j,k,x,fe,fzm,fzp,1,1,mz,0)+theta(t,tz)*(alpham(t)/cw(t))*Ql*Qzr(t,Ql)*Sffzgamma(t,j,k,x,fzphm,fzphp,0))/(2*pi);
	return b;
}

double kerb(int j, int k, double *x, double *dx, double t, double dt, double feb, double fphp, double fphm, double fzp, double fzm, double fzphp, double fzphm){
	double b;
	b = dt*(dx[k]/x[k])*(alphaem(t)*Ql*Ql*Sffv(t,j,k,x,feb,fphp,fphm,1,1,0,0)+theta(t,tz)*alphaz(t)*pow(Qzr(t,Ql),2)*Sffv(t,j,k,x,feb,fzp,fzm,1,1,mz,0)+theta(t,tz)*(alpham(t)/cw(t))*Ql*Qzr(t,Ql)*Sffzgamma(t,j,k,x,fzphp,fzphm,0))/(2*pi);
	return b;
}

/*Evolution for right handed up quarks (u,c)*/
double kur(int j, int k, double *x, double *dx, double t, double dt, double fu, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fzphp, double fzphm, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(alphas(t,alpha)*Sffv(t,j,k,x,fu,fgm,fgp,Cf,Tf,0,0)+alphaem(t)*Qu*Qu*Sffv(t,j,k,x,fu,fphm,fphp,1,3,0,0)+theta(t,tz)*alphaz(t)*pow(Qzr(t,Qu),2)*Sffv(t,j,k,x,fu,fzm,fzp,1,3,mz,0)+3*theta(t,tz)*(alpham(t)/cw(t))*Qu*Qzr(t,Qu)*Sffzgamma(t,j,k,x,fzphm,fzphp,0))/(2*pi);
	return b;
}

double kurb(int j, int k, double *x, double *dx, double t, double dt, double fub, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fzphp, double fzphm, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(alphas(t,alpha)*Sffv(t,j,k,x,fub,fgp,fgm,Cf,Tf,0,0)+alphaem(t)*Qu*Qu*Sffv(t,j,k,x,fub,fphp,fphm,1,3,0,0)+theta(t,tz)*alphaz(t)*pow(Qzr(t,Qu),2)*Sffv(t,j,k,x,fub,fzp,fzm,1,3,mz,0)+3*theta(t,tz)*(alpham(t)/cw(t))*Qu*Qzr(t,Qu)*Sffzgamma(t,j,k,x,fzphp,fzphm,0))/(2*pi);
	return b;
}

/*Evolution for right handed down quarks (same expression for d and b)*/
double kdr(int j, int k, double *x, double *dx, double t, double dt, double fd, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fzphp, double fzphm, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(alphas(t,alpha)*Sffv(t,j,k,x,fd,fgm,fgp,Cf,Tf,0,0)+alphaem(t)*Qd*Qd*Sffv(t,j,k,x,fd,fphm,fphp,1,3,0,0)+theta(t,tz)*alphaz(t)*pow(Qzr(t,Qd),2)*Sffv(t,j,k,x,fd,fzm,fzp,1,3,mz,0)+3*theta(t,tz)*(alpham(t)/cw(t))*Qd*Qzr(t,Qd)*Sffzgamma(t,j,k,x,fzphm,fzphp,0))/(2*pi);
	return b;
}

double kdrb(int j, int k, double *x, double *dx, double t, double dt, double fdb, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fzphp, double fzphm, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(alphas(t,alpha)*Sffv(t,j,k,x,fdb,fgp,fgm,Cf,Tf,0,0)+alphaem(t)*Qd*Qd*Sffv(t,j,k,x,fdb,fphp,fphm,1,3,0,0)+theta(t,tz)*alphaz(t)*pow(Qzr(t,Qd),2)*Sffv(t,j,k,x,fdb,fzp,fzm,1,3,mz,0)+3*theta(t,tz)*(alpham(t)/cw(t))*Qd*Qzr(t,Qd)*Sffzgamma(t,j,k,x,fzphp,fzphm,0))/(2*pi);
	return b;
}

/*Evolution for right handed top*/
double ktr(int j, int k, double *x, double *dx, double t, double dt, double ftl, double ftr, double fb, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fzphp, double fzphm, double fh, double fzl, double fhzl, double fwpl, double **alpha){
	double b;
	b = theta(t,tt)*dt*(dx[k]/x[k])*(alphay(t,alpha)*Stry(t,j,k,x,ftl,fb,fh,fzl,-fhzl,fwpl)+alphas(t,alpha)*Sffv(t,j,k,x,ftr,fgm,fgp,Cf,Tf,0,mt)+alphaem(t)*Qu*Qu*Sffv(t,j,k,x,ftr,fphm,fphp,1,3,0,mt)+theta(t,tz)*alphaz(t)*pow(Qzr(t,Qu),2)*Sffv(t,j,k,x,ftr,fzm,fzp,1,3,mz,mt)+3*theta(t,tz)*(alpham(t)/cw(t))*Qu*Qzr(t,Qu)*Sffzgamma(t,j,k,x,fzphm,fzphp,mt))/(2*pi);
	return b;
}

double ktrb(int j, int k, double *x, double *dx, double t, double dt, double ftlb, double ftrb, double fbb, double fgp, double fgm, double fphp, double fphm, double fzp, double fzm, double fzphp, double fzphm, double fh, double fzl, double fhzl, double fwml, double **alpha){
	double b;
	b = theta(t,tt)*dt*(dx[k]/x[k])*(alphay(t,alpha)*Stry(t,j,k,x,ftlb,fbb,fh,fzl,fhzl,fwml)+alphas(t,alpha)*Sffv(t,j,k,x,ftrb,fgp,fgm,Cf,Tf,0,mt)+alphaem(t)*Qu*Qu*Sffv(t,j,k,x,ftrb,fphp,fphm,1,3,0,mt)+theta(t,tz)*alphaz(t)*pow(Qzr(t,Qu),2)*Sffv(t,j,k,x,ftrb,fzp,fzm,1,3,mz,mt)+3*theta(t,tz)*(alpham(t)/cw(t))*Qu*Qzr(t,Qu)*Sffzgamma(t,j,k,x,fzphp,fzphm,mt))/(2*pi);
	return b;
}

/*Evolution for scalars (Higgs and longitudinal gauge bosons)*/
double kh(int j, int k, double *x, double *dx, double t, double dt, double ftl, double ftlb, double ftr, double ftrb, double fwpp, double fwpm, double fwmp, double fwmm, double fzp, double fzm, double fzl, double fwpl, double fwml, double **alpha){
	double b;
	b = theta(t,th)*dt*(dx[k]/x[k])*(alpha2(t)*0.25*Sswlwt(t,j,k,x,fwpl,fwml,fwpp,fwpm,fwmp,fwmm,mh)+alphaz(t)*0.25*Ssszt(t,j,k,x,fzl,fzp,fzm,mh,mz)+theta(t,tt)*alphay(t,alpha)*0.5*Sstt(t,j,k,x,ftl,ftlb,ftr,ftrb,mh))/(2*pi);
	return b;
}

double kzl(int j, int k, double *x, double *dx, double t, double dt, double ftl, double ftlb, double ftr, double ftrb, double fwpp, double fwpm, double fwmp, double fwmm, double fzp, double fzm, double fh, double fwpl, double fwml, double **alpha){
	double b;
	b = theta(t,tz)*dt*(dx[k]/x[k])*(alpha2(t)*0.25*Sswlwt(t,j,k,x,fwpl,fwml,fwpp,fwpm,fwmp,fwmm,mz)+alphaz(t)*0.25*Ssszt(t,j,k,x,fh,fzp,fzm,mz,mh)+theta(t,tt)*alphay(t,alpha)*0.5*Sstt(t,j,k,x,ftl,ftlb,ftr,ftrb,mz))/(2*pi);
	return b;
}

double khzl(int j, int k, double *x, double *dx, double t, double dt, double ftl, double ftlb, double ftr, double ftrb, double fwpp, double fwpm, double fwmp, double fwmm, double fwpl, double fwml, double **alpha){
	double b;
	b = theta(t,th)*dt*(dx[k]/x[k])*(alpha2(t)*0.25*Shzlwlwt(t,j,k,x,fwpl,fwml,fwpp,fwpm,fwmp,fwmm)+theta(t,tt)*alphay(t,alpha)*Shzltt(t,j,k,x,ftl,ftlb,ftr,ftrb))/(2*pi);
	return b;
}

double kwpl(int j, int k, double *x, double *dx, double t, double dt, double fwpp, double fwpm, double fwpl, double fh, double fzl, double fhzl, double fzp, double fzm, double fphp, double fphm, double fzphp, double fzphm, double ftr, double fblb, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(alpha2(t)*0.25*Swlswt(t,j,k,x,fh,fzl,-fhzl,fwpp,fwpm)+theta(t,tz)*alphaz(t)*pow(c2w(t),2)*0.25*Swlwlvt(t,j,k,x,fwpl,fzp,fzm,mz)+alphaem(t)*Swlwlvt(t,j,k,x,fwpl,fphp,fphm,0)+theta(t,tz)*alpham(t)*c2w(t)*0.5*Swlwlzph(t,j,k,x,fzphp,fzphm)/cw(t)+theta(t,tt)*alphay(t,alpha)*Swltb(t,j,k,x,ftr,fblb))/(2*pi);
	return b;
}

double kwml(int j, int k, double *x, double *dx, double t, double dt, double fwmp, double fwmm, double fwml, double fh, double fzl, double fhzl, double fzp, double fzm, double fphp, double fphm, double fzphp, double fzphm, double ftrb, double fbl, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*(alpha2(t)*0.25*Swlswt(t,j,k,x,fh,fzl,fhzl,fwmp,fwmm)+theta(t,tz)*alphaz(t)*pow(c2w(t),2)*0.25*Swlwlvt(t,j,k,x,fwml,fzp,fzm,mz)+alphaem(t)*Swlwlvt(t,j,k,x,fwml,fphp,fphm,0)+theta(t,tz)*alpham(t)*c2w(t)*0.5*Swlwlzph(t,j,k,x,fzphp,fzphm)/cw(t)+theta(t,tt)*alphay(t,alpha)*Swltb(t,j,k,x,ftrb,fbl))/(2*pi);
	return b;
}

/*Evolution for transverse gauge bosons*/
double kphp(int j, int k, double *x, double *dx, double t, double dt, double fwpp, double fwmp, double fwpm, double fwmm, double fwpl, double fwml, double fel, double fml, double ftal, double felb, double ftalb, double fer, double fmr, double ftar, double ferb, double ftarb, double ful, double fcl, double fdl, double fsl, double ftl, double fbl, double fulb, double fclb, double fdlb, double fslb, double ftlb, double fblb, double fur, double fcr, double fdr, double ftr, double fbr, double furb, double fcrb, double fdrb, double ftrb, double fbrb){
	double b;
	b = dt*(dx[k]/x[k])*alphaem(t)*(cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fwpp+fwmp,fwpm+fwmm,mw,0,mw)+Svtss(t,j,k,x,fwpl+fwml,mw,0,mw)+Sphpff(t,j,k,x,fel,fml,ftal,felb,ftalb,fer,fmr,ftar,ferb,ftarb,ful,fulb,fcl,fclb,ftl,ftlb,fdl,fdlb,fsl,fslb,fbl,fblb,fur,furb,fcr,fcrb,ftr,ftrb,fdr,fdrb,fbr,fbrb))/(2*pi);
	return b;
}

double kphm(int j, int k, double *x, double *dx, double t, double dt, double fwpp, double fwmp, double fwpm, double fwmm, double fwpl, double fwml, double fel, double fml, double ftal, double felb, double ftalb, double fer, double fmr, double ftar, double ferb, double ftarb, double ful, double fcl, double fdl, double fsl, double ftl, double fbl, double fulb, double fclb, double fdlb, double fslb, double ftlb, double fblb, double fur, double fcr, double fdr, double ftr, double fbr, double furb, double fcrb, double fdrb, double ftrb, double fbrb){
	double b;
	b = dt*(dx[k]/x[k])*alphaem(t)*(cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fwpm+fwmm,fwpp+fwmp,mw,0,mw)+Svtss(t,j,k,x,fwpl+fwml,mw,0,mw)+Sphmff(t,j,k,x,fel,fml,ftal,felb,ftalb,fer,fmr,ftar,ferb,ftarb,ful,fulb,fcl,fclb,ftl,ftlb,fdl,fdlb,fsl,fslb,fbl,fblb,fur,furb,fcr,fcrb,ftr,ftrb,fdr,fdrb,fbr,fbrb))/(2*pi);
	return b;
}

double kgp(int j, int k, double *x, double *dx, double t, double dt, double fgp, double fgm, double ful, double fcl, double fdl, double fsl, double ftl, double fbl, double fulb, double fclb, double fdlb, double fslb, double ftlb, double fblb, double fur, double fcr, double fdr, double ftr, double fbr, double furb, double fcrb, double fdrb, double ftrb, double fbrb, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*alphas(t,alpha)*(Ca*Svvv(t,j,k,x,fgp,fgm,0,0,0)+Cf*Sgpff(t,j,k,x,ful,fulb,fcl,fclb,ftl,ftlb,fdl,fdlb,fsl,fslb,fbl,fblb,fur,furb,fcr,fcrb,ftr,ftrb,fdr,fdrb,fbr,fbrb))/(2*pi);
	return b;
}

double kgm(int j, int k, double *x, double *dx, double t, double dt, double fgp, double fgm, double ful, double fcl, double fdl, double fsl, double ftl, double fbl, double fulb, double fclb, double fdlb, double fslb, double ftlb, double fblb, double fur, double fcr, double fdr, double ftr, double fbr, double furb, double fcrb, double fdrb, double ftrb, double fbrb, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*alphas(t,alpha)*(Ca*Svvv(t,j,k,x,fgm,fgp,0,0,0)+Cf*Sgmff(t,j,k,x,ful,fulb,fcl,fclb,ftl,ftlb,fdl,fdlb,fsl,fslb,fbl,fblb,fur,furb,fcr,fcrb,ftr,ftrb,fdr,fdrb,fbr,fbrb))/(2*pi);
	return b;
}

double kzp(int j, int k, double *x, double *dx, double t, double dt, double fwpp, double fwmp, double fwpm, double fwmm, double fwpl, double fwml, double fh, double fzl, double fnue, double fnum, double fnuta, double fnueb, double fnutab, double fel, double fml, double ftal, double felb, double ftalb, double fer, double fmr, double ftar, double ferb, double ftarb, double ful, double fcl, double fdl, double fsl, double ftl, double fbl, double fulb, double fclb, double fdlb, double fslb, double ftlb, double fblb, double fur, double fcr, double fdr, double ftr, double fbr, double furb, double fcrb, double fdrb, double ftrb, double fbrb){
	double b;
	b = dt*(dx[k]/x[k])*alphaz(t)*(pow(cw(t),4)*cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fwpp+fwmp,fwpm+fwmm,mw,mz,mw)+0.25*pow(c2w(t),2)*Svtss(t,j,k,x,fwpl+fwml,mw,mz,mw)+0.25*Svtss(t,j,k,x,fh,mh,mz,mz)+0.25*Svtss(t,j,k,x,fzl,mz,mz,mh)+Szpff(t,j,k,x,fnue,fnum,fnuta,fnueb,fnutab,fel,fml,ftal,felb,ftalb,fer,fmr,ftar,ferb,ftarb,ful,fulb,fcl,fclb,ftl,ftlb,fdl,fdlb,fsl,fslb,fbl,fblb,fur,furb,fcr,fcrb,ftr,ftrb,fdr,fdrb,fbr,fbrb))/(2*pi);
	return b;
}

double kzm(int j, int k, double *x, double *dx, double t, double dt, double fwpp, double fwmp, double fwpm, double fwmm, double fwpl, double fwml, double fh, double fzl, double fnue, double fnum, double fnuta, double fnueb, double fnutab, double fel, double fml, double ftal, double felb, double ftalb, double fer, double fmr, double ftar, double ferb, double ftarb, double ful, double fcl, double fdl, double fsl, double ftl, double fbl, double fulb, double fclb, double fdlb, double fslb, double ftlb, double fblb, double fur, double fcr, double fdr, double ftr, double fbr, double furb, double fcrb, double fdrb, double ftrb, double fbrb){
	double b;
	b = dt*(dx[k]/x[k])*alphaz(t)*(pow(cw(t),4)*cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fwpm+fwmm,fwpp+fwmp,mw,mz,mw)+0.25*pow(c2w(t),2)*Svtss(t,j,k,x,fwpl+fwml,mw,mz,mw)+0.25*Svtss(t,j,k,x,fh,mh,mz,mz)+0.25*Svtss(t,j,k,x,fzl,mz,mz,mh)+Szmff(t,j,k,x,fnue,fnum,fnuta,fnueb,fnutab,fel,fml,ftal,felb,ftalb,fer,fmr,ftar,ferb,ftarb,ful,fulb,fcl,fclb,ftl,ftlb,fdl,fdlb,fsl,fslb,fbl,fblb,fur,furb,fcr,fcrb,ftr,ftrb,fdr,fdrb,fbr,fbrb))/(2*pi);
	return b;
}

double kzphp(int j, int k, double *x, double *dx, double t, double dt, double fwpp, double fwmp, double fwpm, double fwmm, double fwpl, double fwml, double fel, double fml, double ftal, double felb, double ftalb, double fer, double fmr, double ftar, double ferb, double ftarb, double ful, double fcl, double fdl, double fsl, double ftl, double fbl, double fulb, double fclb, double fdlb, double fslb, double ftlb, double fblb, double fur, double fcr, double fdr, double ftr, double fbr, double furb, double fcrb, double fdrb, double ftrb, double fbrb){
	double b;
	b = dt*(dx[k]/x[k])*alpham(t)*(2*cw(t)*cut(x[j],x[k],mu(t))*Szphwtwt(t,j,k,x,fwpp+fwmp,fwpm+fwmm)+c2w(t)*Szphwlwl(t,j,k,x,fwpl+fwml)/cw(t)+2*Szphpff(t,j,k,x,fel,fml,ftal,felb,ftalb,fer,fmr,ftar,ferb,ftarb,ful,fulb,fcl,fclb,ftl,ftlb,fdl,fdlb,fsl,fslb,fbl,fblb,fur,furb,fcr,fcrb,ftr,ftrb,fdr,fdrb,fbr,fbrb)/cw(t))/(2*pi);
	return b;
}

double kzphm(int j, int k, double *x, double *dx, double t, double dt, double fwpp, double fwmp, double fwpm, double fwmm, double fwpl, double fwml, double fel, double fml, double ftal, double felb, double ftalb, double fer, double fmr, double ftar, double ferb, double ftarb, double ful, double fcl, double fdl, double fsl, double ftl, double fbl, double fulb, double fclb, double fdlb, double fslb, double ftlb, double fblb, double fur, double fcr, double fdr, double ftr, double fbr, double furb, double fcrb, double fdrb, double ftrb, double fbrb){
	double b;
	b = dt*(dx[k]/x[k])*alpham(t)*(2*cw(t)*cut(x[j],x[k],mu(t))*Szphwtwt(t,j,k,x,fwpm+fwmm,fwpp+fwmp)+c2w(t)*Szphwlwl(t,j,k,x,fwpl+fwml)/cw(t)+2*Szphmff(t,j,k,x,fel,fml,ftal,felb,ftalb,fer,fmr,ftar,ferb,ftarb,ful,fulb,fcl,fclb,ftl,ftlb,fdl,fdlb,fsl,fslb,fbl,fblb,fur,furb,fcr,fcrb,ftr,ftrb,fdr,fdrb,fbr,fbrb)/cw(t))/(2*pi);
	return b;
}

double kwpp(int j, int k, double *x, double *dx, double t, double dt, double fwpp, double fwpm, double fwpl, double fzp, double fzm, double fphp, double fphm, double fzphp, double fzphm, double fh, double fzl, double fhzl, double fnue, double fnum, double fnuta, double felb, double ftalb, double ful, double fcl, double ftl, double fdlb, double fslb, double fblb){
	double b;
	b = dt*(dx[k]/x[k])*(alpha2(t)*pow(cw(t),2)*(Svvv(t,j,k,x,fwpp,fwpm,mw,mw,mz)+cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fzp,fzm,mz,mw,mw))+alphaem(t)*(Svvv(t,j,k,x,fwpp,fwpm,mw,mw,0)+cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fphp,fphm,0,mw,mw))+cut(x[j],x[k],mu(t))*alpham(t)*cw(t)*Swtzphwt(t,j,k,x,fzphp,fzphm)+0.25*alpha2(t)*(Svtss(t,j,k,x,fh,mh,mw,mw)+Svtss(t,j,k,x,fzl,mz,mw,mw)+Swthzlwl(t,j,k,x,fhzl)+Svtss(t,j,k,x,fwpl,mw,mw,mh)+Svtss(t,j,k,x,fwpl,mw,mw,mz))+0.5*alpha2(t)*Swppff(t,j,k,x,fnue,fnum,fnuta,felb,ftalb,ful,fcl,ftl,fdlb,fslb,fblb))/(2*pi);
	return b;
}

double kwpm(int j, int k, double *x, double *dx, double t, double dt, double fwpp, double fwpm, double fwpl, double fzp, double fzm, double fphp, double fphm, double fzphp, double fzphm, double fh, double fzl, double fhzl, double fnue, double fnum, double fnuta, double felb, double ftalb, double ful, double fcl, double ftl, double fdlb, double fslb, double fblb){
	double b;
	b = dt*(dx[k]/x[k])*(alpha2(t)*pow(cw(t),2)*(Svvv(t,j,k,x,fwpm,fwpp,mw,mw,mz)+cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fzm,fzp,mz,mw,mw))+alphaem(t)*(Svvv(t,j,k,x,fwpm,fwpp,mw,mw,0)+cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fphm,fphp,0,mw,mw))+cut(x[j],x[k],mu(t))*alpham(t)*cw(t)*Swtzphwt(t,j,k,x,fzphm,fzphp)+0.25*alpha2(t)*(Svtss(t,j,k,x,fh,mh,mw,mw)+Svtss(t,j,k,x,fzl,mz,mw,mw)+Swthzlwl(t,j,k,x,fhzl)+Svtss(t,j,k,x,fwpl,mw,mw,mh)+Svtss(t,j,k,x,fwpl,mw,mw,mz))+0.5*alpha2(t)*Swpmff(t,j,k,x,fnue,fnum,fnuta,felb,ftalb,ful,fcl,ftl,fdlb,fslb,fblb))/(2*pi);
	return b;
}

double kwmp(int j, int k, double *x, double *dx, double t, double dt, double fwmp, double fwmm, double fwml, double fzp, double fzm, double fphp, double fphm, double fzphp, double fzphm, double fh, double fzl, double fhzl, double fnueb, double fnutab, double fel, double fml, double ftal, double fdl, double fsl, double fbl, double fulb, double fclb, double ftlb){
	double b;
	b = dt*(dx[k]/x[k])*(alpha2(t)*pow(cw(t),2)*(Svvv(t,j,k,x,fwmp,fwmm,mw,mw,mz)+cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fzp,fzm,mz,mw,mw))+alphaem(t)*(Svvv(t,j,k,x,fwmp,fwmm,mw,mw,0)+cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fphp,fphm,0,mw,mw))+cut(x[j],x[k],mu(t))*alpham(t)*cw(t)*Swtzphwt(t,j,k,x,fzphp,fzphm)+0.25*alpha2(t)*(Svtss(t,j,k,x,fh,mh,mw,mw)+Svtss(t,j,k,x,fzl,mz,mw,mw)-Swthzlwl(t,j,k,x,fhzl)+Svtss(t,j,k,x,fwml,mw,mw,mh)+Svtss(t,j,k,x,fwml,mw,mw,mz))+0.5*alpha2(t)*Swmpff(t,j,k,x,fnueb,fnutab,fel,fml,ftal,fulb,fclb,ftlb,fdl,fsl,fbl))/(2*pi);
	return b;
}

double kwmm(int j, int k, double *x, double *dx, double t, double dt, double fwmp, double fwmm, double fwml, double fzp, double fzm, double fphp, double fphm, double fzphp, double fzphm, double fh, double fzl, double fhzl, double fnueb, double fnutab, double fel, double fml, double ftal, double fdl, double fsl, double fbl, double fulb, double fclb, double ftlb){
	double b;
	b = dt*(dx[k]/x[k])*(alpha2(t)*pow(cw(t),2)*(Svvv(t,j,k,x,fwmm,fwmp,mw,mw,mz)+cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fzm,fzp,mz,mw,mw))+alphaem(t)*(Svvv(t,j,k,x,fwmm,fwmp,mw,mw,0)+cut(x[j],x[k],mu(t))*Svvv(t,j,k,x,fphm,fphp,0,mw,mw))+cut(x[j],x[k],mu(t))*alpham(t)*cw(t)*Swtzphwt(t,j,k,x,fzphm,fzphp)+0.25*alpha2(t)*(Svtss(t,j,k,x,fh,mh,mw,mw)+Svtss(t,j,k,x,fzl,mz,mw,mw)-Swthzlwl(t,j,k,x,fhzl)+Svtss(t,j,k,x,fwml,mw,mw,mh)+Svtss(t,j,k,x,fwml,mw,mw,mz))+0.5*alpha2(t)*Swmmff(t,j,k,x,fnueb,fnutab,fel,fml,ftal,fulb,fclb,ftlb,fdl,fsl,fbl))/(2*pi);
	return b;
}
