/*Ultra-collinear evolution for left-handed leptons (same expression for el, elb, ml)*/
double keluc(int j, int k, double *x, double *dx, double t, double dt, double fe, double fnu, double fwl, double fzl){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(cut(x[j],x[k],mu(t))*prop2(t,0,0,mw,x[j]/x[k])*Pucfl2fl1wl(x[j]/x[k],alpha2(t),0,0)*fnu+prop2(t,0,0,mz,x[j]/x[k])*theta(t,tz)*Pucflflzl(x[j]/x[k],alpha2(t),0,t,-0.5,-1)*fe+prop2(t,mw,0,0,x[j]/x[k])*Pucfl1wlflb2(x[j]/x[k],alpha2(t),0,0,1)*fwl+prop2(t,mz,0,0,x[j]/x[k])*theta(t,tz)*Pucflzlflb(x[j]/x[k],alpha2(t),0,t,-0.5,-1,1)*fzl);
	return b;
}

/*Ultra-collinear evolution for neutrinos (same expression for nue, nueb, num)*/
double knuuc(int j, int k, double *x, double *dx, double t, double dt, double fe, double fnu, double fwl, double fzl){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(cut(x[j],x[k],mu(t))*prop2(t,0,0,mw,x[j]/x[k])*Pucfl2fl1wl(x[j]/x[k],alpha2(t),0,0)*fe+prop2(t,0,0,mz,x[j]/x[k])*theta(t,tz)*Pucflflzl(x[j]/x[k],alpha2(t),0,t,0.5,0)*fnu+prop2(t,mw,0,0,x[j]/x[k])*Pucfl1wlflb2(x[j]/x[k],alpha2(t),0,0,1)*fwl+prop2(t,mz,0,0,x[j]/x[k])*theta(t,tz)*Pucflzlflb(x[j]/x[k],alpha2(t),0,t,0.5,0,1)*fzl);
	return b;

}

/*Ultra-collinear evolution for left-handed up quarks (same expression for ul, ulb)*/
double kuluc(int j, int k, double *x, double *dx, double t, double dt, double fu, double fd, double fwl, double fzl){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(cut(x[j],x[k],mu(t))*prop2(t,0,0,mw,x[j]/x[k])*Pucfl2fl1wl(x[j]/x[k],alpha2(t),0,0)*fd+prop2(t,0,0,mz,x[j]/x[k])*theta(t,tz)*Pucflflzl(x[j]/x[k],alpha2(t),0,t,0.5,Qu)*fu+prop2(t,mw,0,0,x[j]/x[k])*Pucfl1wlflb2(x[j]/x[k],alpha2(t),0,0,3)*fwl+prop2(t,mz,0,0,x[j]/x[k])*theta(t,tz)*Pucflzlflb(x[j]/x[k],alpha2(t),0,t,0.5,Qu,3)*fzl);
	return b;
}

/*Ultra-collinear evolution for left-handed top quark (same expression for tl, tlb)*/
double ktluc(int j, int k, double *x, double *dx, double t, double dt, double ftl, double ftr, double fb, double fwl, double fzl, double fh, double fhzl, double fph, double fz, double fzph, double fg, double **alpha){
	double b;
	b = theta(t,tt)*dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(cut(x[j],x[k],mu(t))*prop2(t,0,mt,mw,x[j]/x[k])*Pucfl2fl1wl(x[j]/x[k],alpha2(t),0,alphay(t,alpha))*fb+prop2(t,mt,mt,mz,x[j]/x[k])*Pucflflzl(x[j]/x[k],alpha2(t),alphay(t,alpha),t,0.5,Qu)*ftl+prop2(t,mw,mt,0,x[j]/x[k])*Pucfl1wlflb2(x[j]/x[k],alpha2(t),alphay(t,alpha),0,3)*fwl+prop2(t,mz,mt,mt,x[j]/x[k])*Pucflzlflb(x[j]/x[k],alpha2(t),alphay(t,alpha),t,0.5,Qu,3)*fzl+prop2(t,mh,mt,mt,x[j]/x[k])*Puctht(x[j]/x[k],alphay(t,alpha))*fh+prop2(t,mt,mt,mh,x[j]/x[k])*Puctth(x[j]/x[k],alphay(t,alpha))*ftl+prop(t,mh,mt,mt,x[j]/x[k])*prop(t,mz,mt,mt,x[j]/x[k])*Puctlhzltlb(x[j]/x[k],alpha2(t),alphay(t,alpha),t)*fhzl+prop2(t,0,mt,mt,x[j]/x[k])*Puctgt(x[j]/x[k],alphas(t,alpha),alphay(t,alpha))*fg+prop2(t,0,mt,mt,x[j]/x[k])*Puctpht(x[j]/x[k],alphaem(t),alphay(t,alpha))*fph+prop2(t,mz,mt,mt,x[j]/x[k])*Puctlzmtrb(x[j]/x[k],alpha2(t),alphay(t,alpha),t)*fz+prop2(t,mt,mt,0,x[j]/x[k])*Pucttg(x[j]/x[k],alphas(t,alpha),alphay(t,alpha))*ftr+prop2(t,mt,mt,0,x[j]/x[k])*Pucttph(x[j]/x[k],alphaem(t),alphay(t,alpha))*ftr+prop2(t,mt,mt,mz,x[j]/x[k])*Puctltrzp(x[j]/x[k],alpha2(t),alphay(t,alpha),t)*ftr+prop(t,mz,mt,mt,x[j]/x[k])*prop(t,0,mt,mt,x[j]/x[k])*Puctlzphmtrb(x[j]/x[k],alpham(t),alphay(t,alpha),t)*fzph);
	return b;
}

/*Ultra-collinear evolution for left-handed down quarks (same expression for dl, dlb)*/
double kdluc(int j, int k, double *x, double *dx, double t, double dt, double fd, double fu, double fwl, double fzl){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(cut(x[j],x[k],mu(t))*prop2(t,0,0,mw,x[j]/x[k])*Pucfl2fl1wl(x[j]/x[k],alpha2(t),0,0)*fu+prop2(t,0,0,mz,x[j]/x[k])*theta(t,tz)*Pucflflzl(x[j]/x[k],alpha2(t),0,t,-0.5,Qd)*fd+prop2(t,mw,0,0,x[j]/x[k])*Pucfl1wlflb2(x[j]/x[k],alpha2(t),0,0,3)*fwl+prop2(t,mz,0,0,x[j]/x[k])*theta(t,tz)*Pucflzlflb(x[j]/x[k],alpha2(t),0,t,-0.5,Qd,3)*fzl);
	return b;
}

/*Ultra-collinear evolution for left-handed bottom quark (same expression for bl, blb)*/
double kbluc(int j, int k, double *x, double *dx, double t, double dt, double fb, double ftl, double ftr, double fwl, double fzl, double fwt, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(cut(x[j],x[k],mu(t))*prop2(t,mt,0,mw,x[j]/x[k])*Pucfl2fl1wl(x[j]/x[k],alpha2(t),alphay(t,alpha),0)*ftl+prop2(t,0,0,mz,x[j]/x[k])*theta(t,tz)*Pucflflzl(x[j]/x[k],alpha2(t),0,t,-0.5,Qd)*fb+prop2(t,mw,0,mt,x[j]/x[k])*Pucfl1wlflb2(x[j]/x[k],alpha2(t),0,alphay(t,alpha),3)*fwl+prop2(t,mz,0,0,x[j]/x[k])*theta(t,tz)*Pucflzlflb(x[j]/x[k],alpha2(t),0,t,-0.5,Qd,3)*fzl+prop2(t,mt,0,mw,x[j]/x[k])*Pucbltrwpp(x[j]/x[k],alpha2(t),alphay(t,alpha))*ftr+prop2(t,mw,0,mt,x[j]/x[k])*Pucblwmmtrb(x[j]/x[k],alpha2(t),alphay(t,alpha))*fwt);
	return b;
}

/*Ultra-collinear evolution for right-handed leptons (same expression for er, erb, mr)*/
double keruc(int j, int k, double *x, double *dx, double t, double dt, double fe, double fzl){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*theta(t,tz)*(prop2(t,0,0,mz,x[j]/x[k])*Pucfrfrzl(x[j]/x[k],alpha2(t),0,t,0,Ql)*fe+prop2(t,mz,0,0,x[j]/x[k])*Pucfrzlfrb(x[j]/x[k],alpha2(t),0,t,0,Ql,1)*fzl);
	return b;
}

/*Ultra-collinear evolution for right-handed up quarks (same expression for ur, urb)*/
double kuruc(int j, int k, double *x, double *dx, double t, double dt, double fu, double fzl){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*theta(t,tz)*(prop2(t,0,0,mz,x[j]/x[k])*Pucfrfrzl(x[j]/x[k],alpha2(t),0,t,0,Qu)*fu+prop2(t,mz,0,0,x[j]/x[k])*Pucfrzlfrb(x[j]/x[k],alpha2(t),0,t,0,Qu,3)*fzl);
	return b;
}

/*Ultra-collinear evolution for right-handed top quark (same expression for tr, trb)*/
double ktruc(int j, int k, double *x, double *dx, double t, double dt, double ftr, double ftl, double fbl, double fzl, double fh, double fhzl, double fph, double fz, double fzph, double fg, double fwt, double **alpha){
	double b;
	b = theta(t,tt)*dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop2(t,mt,mt,mz,x[j]/x[k])*Pucfrfrzl(x[j]/x[k],alpha2(t),alphay(t,alpha),t,0.5,Qu)*ftr+prop2(t,mz,mt,mt,x[j]/x[k])*Pucfrzlfrb(x[j]/x[k],alpha2(t),alphay(t,alpha),t,0.5,Qu,3)*fzl+prop2(t,mh,mt,mt,x[j]/x[k])*Puctht(x[j]/x[k],alphay(t,alpha))*fh+prop2(t,mt,mt,mh,x[j]/x[k])*Puctth(x[j]/x[k],alphay(t,alpha))*ftr+prop(t,mh,mt,mt,x[j]/x[k])*prop(t,mz,mt,mt,x[j]/x[k])*Puctlhzltlb(x[j]/x[k],alpha2(t),alphay(t,alpha),t)*fhzl+prop2(t,0,mt,mt,x[j]/x[k])*Puctgt(x[j]/x[k],alphas(t,alpha),alphay(t,alpha))*fg+prop2(t,0,mt,mt,x[j]/x[k])*Puctpht(x[j]/x[k],alphaem(t),alphay(t,alpha))*fph+prop2(t,mz,mt,mt,x[j]/x[k])*Puctrzptlb(x[j]/x[k],alpha2(t),alphay(t,alpha),t)*fz+prop2(t,mt,mt,0,x[j]/x[k])*Pucttg(x[j]/x[k],alphas(t,alpha),alphay(t,alpha))*ftl+prop2(t,mt,mt,0,x[j]/x[k])*Pucttph(x[j]/x[k],alphaem(t),alphay(t,alpha))*ftl+prop2(t,mt,mt,mz,x[j]/x[k])*Puctrtlzm(x[j]/x[k],alpha2(t),alphay(t,alpha),t)*ftl+prop(t,mz,mt,mt,x[j]/x[k])*prop(t,0,mt,mt,x[j]/x[k])*Puctrzphptlb(x[j]/x[k],alpham(t),alphay(t,alpha),t)*fzph+prop2(t,0,mt,mw,x[j]/x[k])*Puctrblwmm(x[j]/x[k],alpha2(t),alphay(t,alpha))*fbl+prop2(t,mw,mt,0,x[j]/x[k])*Puctrwppblb(x[j]/x[k],alpha2(t),alphay(t,alpha))*fwt);
	return b;
}

/*Ultra-collinear evolution for right-handed down quarks (same expression for dr, drb, br, brb)*/
double kdruc(int j, int k, double *x, double *dx, double t, double dt, double fd, double fzl){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*theta(t,tz)*(prop2(t,0,0,mz,x[j]/x[k])*Pucfrfrzl(x[j]/x[k],alpha2(t),0,t,0,Qd)*fd+prop2(t,mz,0,0,x[j]/x[k])*Pucfrzlfrb(x[j]/x[k],alpha2(t),0,t,0,Qd,3)*fzl);
	return b;
}

/*Ultra-collinear evolution for scalars (Higgs and longitudinal gauge bosons)*/
double khuc(int j, int k, double *x, double *dx, double t, double dt, double fwpl, double fwml, double fwpp, double fwmp, double fwpm, double fwmm, double fzl, double fzp, double fzm, double fh, double ftl, double ftr, double ftlb, double ftrb,double **alpha){
	double b;
	b = theta(t,th)*dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop2(t,mh,mh,mh,x[j]/x[k])*Puchhh(x[j]/x[k])*fh+cut(x[j],x[k],mu(t))*prop2(t,mw,mh,mw,x[j]/x[k])*Puchwlwl(x[j]/x[k],alpha2(t))*(fwpl+fwml)+prop2(t,mz,mh,mz,x[j]/x[k])*cut(x[j],x[k],mu(t))*Puchzlzl(x[j]/x[k],alpha2(t),t)*fzl+prop2(t,mw,mh,mw,x[j]/x[k])*Puchwtwt(x[j]/x[k],alpha2(t))*(fwpp+fwmp+fwpm+fwmm)+prop2(t,mz,mh,mz,x[j]/x[k])*Puchztzt(x[j]/x[k],alpha2(t),t)*(fzp+fzm)+prop2(t,mt,mh,mt,x[j]/x[k])*theta(t,tt)*Puchtt(x[j]/x[k],alphay(t,alpha))*(ftl+ftlb+ftr+ftrb));
	return b;
}

double kwlpuc(int j, int k, double *x, double *dx, double t, double dt, double fwl, double fwp, double fwm, double fzl, double fzp, double fzm, double fphp, double fphm, double fzphp, double fzphm, double fh, double fhzl, double felb, double fnue, double fnum, double ful, double ftl, double fdlb, double fblb, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop2(t,0,mw,0,x[j]/x[k])*Pucwlfl1fl2(x[j]/x[k],alpha2(t),0,0)*(2*fnue+fnum+3*felb+2*ful+2*fdlb)+prop2(t,mt,mw,0,x[j]/x[k])*theta(t,tt)*Pucwlfl1fl2(x[j]/x[k],alpha2(t),alphay(t,alpha),0)*ftl+prop2(t,0,mw,mt,x[j]/x[k])*theta(t,tt)*Pucwlfl1fl2(x[j]/x[k],alpha2(t),0,alphay(t,alpha))*fblb+prop2(t,mw,mw,0,x[j]/x[k])*Pucwlwtph(x[j]/x[k],alphaem(t),alpha2(t))*(fwp+fwm)+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwlwtzt(x[j]/x[k],alpha2(t),t)*(fwp+fwm)+prop2(t,mz,mw,mw,x[j]/x[k])*theta(t,tz)*Pucwlztwt(x[j]/x[k],alpha2(t),t)*(fzp+fzm)+prop2(t,0,mw,mw,x[j]/x[k])*Pucwlphwt(x[j]/x[k],alphaem(t),alpha2(t))*(fphp+fphm)+prop(t,0,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*theta(t,tz)*Pucwlzphwt(x[j]/x[k],alpham(t),alpha2(t),t)*(fzphp+fzphm)+prop2(t,mh,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,th)*Pucwlhwl(x[j]/x[k],alpha2(t))*fh+prop2(t,mz,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,tz)*Pucwlzlwl(x[j]/x[k],alpha2(t),t)*fzl+prop2(t,mw,mw,mh,x[j]/x[k])*theta(t,th)*Pucwlwlh(x[j]/x[k],alpha2(t))*fwl+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwlwlzl(x[j]/x[k],alpha2(t),t)*fwl+prop(t,mh,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,th)*Pucwlhzlwl(x[j]/x[k],alpha2(t),t)*fhzl);
	return b;
}

double kwlmuc(int j, int k, double *x, double *dx, double t, double dt, double fwl, double fwp, double fwm, double fzl, double fzp, double fzm, double fphp, double fphm, double fzphp, double fzphm, double fh, double fhzl, double fel, double fml, double fnub, double fulb, double ftlb, double fdl, double fbl, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop2(t,0,mw,0,x[j]/x[k])*Pucwlfl1fl2(x[j]/x[k],alpha2(t),0,0)*(3*fnub+2*fel+fml+2*fulb+2*fdl)+prop2(t,mt,mw,0,x[j]/x[k])*theta(t,tt)*Pucwlfl1fl2(x[j]/x[k],alpha2(t),alphay(t,alpha),0)*ftlb+prop2(t,0,mw,mt,x[j]/x[k])*theta(t,tt)*Pucwlfl1fl2(x[j]/x[k],alpha2(t),0,alphay(t,alpha))*fbl+prop2(t,mw,mw,0,x[j]/x[k])*Pucwlwtph(x[j]/x[k],alphaem(t),alpha2(t))*(fwp+fwm)+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwlwtzt(x[j]/x[k],alpha2(t),t)*(fwp+fwm)+prop2(t,mz,mw,mw,x[j]/x[k])*theta(t,tz)*Pucwlztwt(x[j]/x[k],alpha2(t),t)*(fzp+fzm)+prop2(t,0,mw,mw,x[j]/x[k])*Pucwlphwt(x[j]/x[k],alphaem(t),alpha2(t))*(fphp+fphm)+prop(t,0,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*theta(t,tz)*Pucwlzphwt(x[j]/x[k],alpham(t),alpha2(t),t)*(fzphp+fzphm)+prop2(t,mh,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,th)*Pucwlhwl(x[j]/x[k],alpha2(t))*fh+prop2(t,mz,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,tz)*Pucwlzlwl(x[j]/x[k],alpha2(t),t)*fzl+prop2(t,mw,mw,mh,x[j]/x[k])*theta(t,th)*Pucwlwlh(x[j]/x[k],alpha2(t))*fwl+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwlwlzl(x[j]/x[k],alpha2(t),t)*fwl-prop(t,mh,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,th)*Pucwlhzlwl(x[j]/x[k],alpha2(t),t)*fhzl);
	return b;
}

double kzluc(int j, int k, double *x, double *dx, double t, double dt, double fwpl, double fwml, double fwpp, double fwmp, double fwpm, double fwmm, double fzl, double fh, double fel, double fml, double felb, double fnue, double fnum, double fnub, double fer, double fmr, double ferb, double ful, double fulb, double ftl, double ftlb, double fdl, double fdlb, double fbl, double fblb, double fur, double furb, double ftr, double ftrb, double fdr, double fdrb, double fbr, double fbrb, double **alpha){
	double b;
	b = theta(t,tz)*dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop2(t,0,mz,0,x[j]/x[k])*Puczlflfl(x[j]/x[k],alpha2(t),0,t,0.5,0)*(2*fnue+fnum+3*fnub)+prop2(t,0,mz,0,x[j]/x[k])*Puczlflfl(x[j]/x[k],alpha2(t),0,t,-0.5,Ql)*(2*fel+fml+3*felb)+prop2(t,0,mz,0,x[j]/x[k])*Puczlflfl(x[j]/x[k],alpha2(t),0,t,0.5,Qu)*(2*ful+2*fulb)+prop2(t,0,mz,0,x[j]/x[k])*Puczlflfl(x[j]/x[k],alpha2(t),0,t,-0.5,Qd)*(2*fdl+fbl+2*fdlb+fblb)+prop2(t,0,mz,0,x[j]/x[k])*Puczlfrfr(x[j]/x[k],alpha2(t),0,t,0,Ql)*(2*fer+fmr+3*ferb)+prop2(t,0,mz,0,x[j]/x[k])*Puczlfrfr(x[j]/x[k],alpha2(t),0,t,0,Qu)*(2*fur+2*furb)+prop2(t,0,mz,0,x[j]/x[k])*Puczlfrfr(x[j]/x[k],alpha2(t),0,t,0,Qd)*(2*fdr+fbr+2*fdrb+fbrb)+prop2(t,mt,mz,mt,x[j]/x[k])*theta(t,tt)*(Puczlflfl(x[j]/x[k],alpha2(t),alphay(t,alpha),t,0.5,Qu)*(ftl+ftlb)+Puczlfrfr(x[j]/x[k],alpha2(t),alphay(t,alpha),t,0.5,Qu)*(ftr+ftrb))+prop2(t,mw,mz,mw,x[j]/x[k])*Puczlwtwt(x[j]/x[k],alpha2(t))*(fwpp+fwmp+fwpm+fwmm)+prop2(t,mz,mz,mh,x[j]/x[k])*theta(t,th)*Puczlzlh(x[j]/x[k],alpha2(t),t)*fzl+prop2(t,mh,mz,mz,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,th)*Puczlhzl(x[j]/x[k],alpha2(t),t)*fh+prop2(t,mw,mz,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*Puczlwlwl(x[j]/x[k],alpha2(t),t)*(fwpl+fwml));
	return b;
}

double khzluc(int j, int k, double *x, double *dx, double t, double dt, double fwpp, double fwpm, double fwmp, double fwmm, double fwpl, double fwml, double ftl, double ftr, double ftlb, double ftrb, double **alpha){
	double b;
	b = theta(t,th)*dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop(t,mw,mh,mw,x[j]/x[k])*prop(t,mw,mz,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*Puchzlwlwl(x[j]/x[k],alpha2(t),t)*(fwpl-fwml)+prop(t,mw,mh,mw,x[j]/x[k])*prop(t,mw,mz,mw,x[j]/x[k])*Puchzlwtwt(x[j]/x[k],alpha2(t))*(fwpp+fwpm-fwmp-fwmm)+prop(t,mt,mh,mt,x[j]/x[k])*prop(t,mt,mz,mt,x[j]/x[k])*theta(t,tt)*(Puchzltltl(x[j]/x[k],alpha2(t),alphay(t,alpha),t)*(ftl-ftlb)+Puchzltrtr(x[j]/x[k],alpha2(t),alphay(t,alpha),t)*(ftr-ftrb)));
	return b;
}

/*Ultra-collinear evolution for the Ws*/
double kwppuc(int j, int k, double *x, double *dx, double t, double dt, double fwpl, double fwpp, double fphp, double fzp, double fzphp, double fh, double fzl, double fhzl, double ftr, double fblb, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop2(t,mw,mw,0,x[j]/x[k])*Pucwtwlph(x[j]/x[k],alphaem(t),alpha2(t))*fwpl+prop2(t,0,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*Pucwtphwl(x[j]/x[k],alphaem(t),alpha2(t))*fphp+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwtwlzt(x[j]/x[k],alpha2(t),t)*fwpl+prop2(t,mz,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,tz)*Pucwtztwl(x[j]/x[k],alpha2(t),t)*fzp+prop(t,mz,mw,mw,x[j]/x[k])*prop(t,0,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,tz)*Pucwtzphwl(x[j]/x[k],alpham(t),alpha2(t),t)*fzphp+(prop2(t,mw,mw,mh,x[j]/x[k])*theta(t,th)*Pucwtwth(x[j]/x[k],alpha2(t))+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwtwtzl(x[j]/x[k],alpha2(t)))*fwpp+prop2(t,mh,mw,mw,x[j]/x[k])*theta(t,th)*Pucwthwt(x[j]/x[k],alpha2(t))*fh+prop2(t,mz,mw,mw,x[j]/x[k])*theta(t,tz)*Pucwtzlwt(x[j]/x[k],alpha2(t))*fzl+prop(t,mh,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*theta(t,th)*Pucwthzlwt(x[j]/x[k],alpha2(t))*fhzl+theta(t,tt)*prop2(t,0,mw,mt,x[j]/x[k])*Pucwmmbltr(x[j]/x[k],alpha2(t),alphay(t,alpha))*fblb+theta(t,tt)*prop2(t,mt,mw,0,x[j]/x[k])*Pucwpptrbl(x[j]/x[k],alpha2(t),alphay(t,alpha))*ftr);
	return b;
}

double kwpmuc(int j, int k, double *x, double *dx, double t, double dt, double fwpl, double fwpm, double fphm, double fzm, double fzphm, double fh, double fzl, double fhzl){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop2(t,mw,mw,0,x[j]/x[k])*Pucwtwlph(x[j]/x[k],alphaem(t),alpha2(t))*fwpl+prop2(t,0,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*Pucwtphwl(x[j]/x[k],alphaem(t),alpha2(t))*fphm+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwtwlzt(x[j]/x[k],alpha2(t),t)*fwpl+prop2(t,mz,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,tz)*Pucwtztwl(x[j]/x[k],alpha2(t),t)*fzm+prop(t,mz,mw,mw,x[j]/x[k])*prop(t,0,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,tz)*Pucwtzphwl(x[j]/x[k],alpham(t),alpha2(t),t)*fzphm+(prop2(t,mw,mw,mh,x[j]/x[k])*theta(t,th)*Pucwtwth(x[j]/x[k],alpha2(t))+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwtwtzl(x[j]/x[k],alpha2(t)))*fwpm+prop2(t,mh,mw,mw,x[j]/x[k])*theta(t,th)*Pucwthwt(x[j]/x[k],alpha2(t))*fh+prop2(t,mz,mw,mw,x[j]/x[k])*theta(t,tz)*Pucwtzlwt(x[j]/x[k],alpha2(t))*fzl+prop(t,mh,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*theta(t,th)*Pucwthzlwt(x[j]/x[k],alpha2(t))*fhzl);
	return b;
}

double kwmpuc(int j, int k, double *x, double *dx, double t, double dt, double fwml, double fwmp, double fphp, double fzp, double fzphp, double fh, double fzl, double fhzl){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop2(t,mw,mw,0,x[j]/x[k])*Pucwtwlph(x[j]/x[k],alphaem(t),alpha2(t))*fwml+prop2(t,0,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*Pucwtphwl(x[j]/x[k],alphaem(t),alpha2(t))*fphp+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwtwlzt(x[j]/x[k],alpha2(t),t)*fwml+prop2(t,mz,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,tz)*Pucwtztwl(x[j]/x[k],alpha2(t),t)*fzp+prop(t,mz,mw,mw,x[j]/x[k])*prop(t,0,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,tz)*Pucwtzphwl(x[j]/x[k],alpham(t),alpha2(t),t)*fzphp+(prop2(t,mw,mw,mh,x[j]/x[k])*theta(t,th)*Pucwtwth(x[j]/x[k],alpha2(t))+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwtwtzl(x[j]/x[k],alpha2(t)))*fwmp+prop2(t,mh,mw,mw,x[j]/x[k])*theta(t,th)*Pucwthwt(x[j]/x[k],alpha2(t))*fh+prop2(t,mz,mw,mw,x[j]/x[k])*theta(t,tz)*Pucwtzlwt(x[j]/x[k],alpha2(t))*fzl-prop(t,mh,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*theta(t,th)*Pucwthzlwt(x[j]/x[k],alpha2(t))*fhzl);
	return b;
}

double kwmmuc(int j, int k, double *x, double *dx, double t, double dt, double fwml, double fwmm, double fphm, double fzm, double fzphm, double fh, double fzl, double fhzl, double ftrb, double fbl, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop2(t,mw,mw,0,x[j]/x[k])*Pucwtwlph(x[j]/x[k],alphaem(t),alpha2(t))*fwml+prop2(t,0,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*Pucwtphwl(x[j]/x[k],alphaem(t),alpha2(t))*fphm+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwtwlzt(x[j]/x[k],alpha2(t),t)*fwml+prop2(t,mz,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,tz)*Pucwtztwl(x[j]/x[k],alpha2(t),t)*fzm+prop(t,mz,mw,mw,x[j]/x[k])*prop(t,0,mw,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*theta(t,tz)*Pucwtzphwl(x[j]/x[k],alpham(t),alpha2(t),t)*fzphm+(prop2(t,mw,mw,mh,x[j]/x[k])*theta(t,th)*Pucwtwth(x[j]/x[k],alpha2(t))+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz)*Pucwtwtzl(x[j]/x[k],alpha2(t)))*fwmm+prop2(t,mh,mw,mw,x[j]/x[k])*theta(t,th)*Pucwthwt(x[j]/x[k],alpha2(t))*fh+prop2(t,mz,mw,mw,x[j]/x[k])*theta(t,tz)*Pucwtzlwt(x[j]/x[k],alpha2(t))*fzl-prop(t,mh,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*theta(t,th)*Pucwthzlwt(x[j]/x[k],alpha2(t))*fhzl+theta(t,tt)*prop2(t,0,mw,mt,x[j]/x[k])*Pucwmmbltr(x[j]/x[k],alpha2(t),alphay(t,alpha))*fbl+theta(t,tt)*prop2(t,mt,mw,0,x[j]/x[k])*Pucwpptrbl(x[j]/x[k],alpha2(t),alphay(t,alpha))*ftrb);
	return b;
}

/*Ultra-collinear evolution for neutral gauge bosons (same expression for the two transverse polarizations)*/
double kphuc(int j, int k, double *x, double *dx, double t, double dt, double fwp, double fwm, double fwpl, double fwml, double ftl, double ftr, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop2(t,mw,0,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*Pucphwtwl(x[j]/x[k],alphaem(t),alpha2(t))*(fwp+fwm)+prop2(t,mw,0,mw,x[j]/x[k])*Pucphwlwt(x[j]/x[k],alphaem(t),alpha2(t))*(fwpl+fwml)+prop2(t,mt,0,mt,x[j]/x[k])*theta(t,tt)*Pucphtt(x[j]/x[k],alphaem(t),alphay(t,alpha))*(ftl+ftr));
	return b;
}

double kzuc(int j, int k, double *x, double *dx, double t, double dt, double fwp, double fwm, double fwpl, double fwml, double fh, double fz, double ftl, double ftr, double **alpha){
	double b;
	b = theta(t,tz)*dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop2(t,mw,mz,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*Pucztwtwl(x[j]/x[k],alpha2(t),t)*(fwp+fwm)+prop2(t,mw,mz,mw,x[j]/x[k])*Pucztwlwt(x[j]/x[k],alpha2(t),t)*(fwpl+fwml)+prop2(t,mz,mz,mh,x[j]/x[k])*theta(t,th)*Pucztzth(x[j]/x[k],alpha2(t),t)*fz+prop2(t,mh,mz,mz,x[j]/x[k])*theta(t,th)*Puczthzt(x[j]/x[k],alpha2(t),t)*fh+prop2(t,mt,mz,mt,x[j]/x[k])*theta(t,tt)*Puczptrtl(x[j]/x[k],alpha2(t),alphay(t,alpha),t)*ftr+prop2(t,mt,mz,mt,x[j]/x[k])*theta(t,tt)*Puczmtltr(x[j]/x[k],alpha2(t),alphay(t,alpha),t)*ftl);
	return b;
}

double kzphuc(int j, int k, double *x, double *dx, double t, double dt, double fwp, double fwm, double fwpl, double fwml, double ftl, double ftr, double **alpha){
	double b;
	b = theta(t,tz)*dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*(prop(t,mw,0,mw,x[j]/x[k])*prop(t,mw,mz,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*Puczphwtwl(x[j]/x[k],alpham(t),alpha2(t),t)*(fwp+fwm)+prop(t,mw,0,mw,x[j]/x[k])*prop(t,mw,mz,mw,x[j]/x[k])*Puczphwlwt(x[j]/x[k],alpham(t),alpha2(t),t)*(fwpl+fwml)+prop(t,mt,0,mt,x[j]/x[k])*prop(t,mt,mz,mt,x[j]/x[k])*theta(t,tt)*(Puczphptrtl(x[j]/x[k],alpham(t),alphay(t,alpha),t)*ftr+Puczphmtltr(x[j]/x[k],alpham(t),alphay(t,alpha),t)*ftl));
	return b;
}

double kguc(int j, int k, double *x, double *dx, double t, double dt, double ftl, double ftr, double **alpha){
	double b;
	b = dt*(dx[k]/x[k])*pow(double(1)/mu(t),2)*theta(t,tt)*prop2(t,mt,0,mt,x[j]/x[k])*Pucgtt(x[j]/x[k],alphas(t,alpha),alphay(t,alpha))*(ftl+ftr);
	return b;
}
