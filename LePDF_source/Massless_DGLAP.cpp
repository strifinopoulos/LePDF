/*Terms contributing to DGLAP equations. Notation: S stands for splitting, followed by the parton B and a specification different for each one (see below in details)*/

double Sffv(double t, int j, int k, double *x, double ff, double fvp, double fvm, double cpffg, double cpfv, double mv, double mf){/*A and B = same fermion, C = transverse neutral gauge boson*/
	double s;
	s = prop2(t,mf,mf,mv,x[j]/x[k])*cpffg*Pffg(x[j]/x[k])*ff+cpfv*prop2(t,mv,mf,mf,x[j]/x[k])*(Pflvp(x[j]/x[k])*fvp+Pflvm(x[j]/x[k])*fvm);
	return s;
}

double Sffvcut(double t, int j, int k, double *x, double ff, double fvp, double fvm, double Nc, double mfb, double mf){/*A and B = two different fermions, C = transverse W*/
	double s;
	s = prop2(t,mf,mfb,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*Pffg(x[j]/x[k])*ff+prop2(t,mw,mfb,mf,x[j]/x[k])*Nc*(Pflvp(x[j]/x[k])*fvp+Pflvm(x[j]/x[k])*fvm);
	return s;
}

double Sffzgamma(double t, int j, int k, double *x, double fvp, double fvm, double mf){/*A and B = same fermion, C = transverse Z-photon*/
	double s;
	s = prop(t,0,mf,mf,x[j]/x[k])*prop(t,mz,mf,mf,x[j]/x[k])*(Pflvp(x[j]/x[k])*fvp+Pflvm(x[j]/x[k])*fvm);
	return s;
}

double Stly(double t, int j, int k, double *x, double ftr, double fh, double fzl, double fhzl){/*Yukawa interactions with A and B = t_L*/
	double s;
	s = Pffy(x[j]/x[k])*(prop2(t,mt,mt,mh,x[j]/x[k])+prop2(t,mt,mt,mz,x[j]/x[k]))*ftr+3*Pfh(x[j]/x[k])*(prop2(t,mh,mt,mt,x[j]/x[k])*fh+prop2(t,mz,mt,mt,x[j]/x[k])*fzl+prop(t,mh,mt,mt,x[j]/x[k])*prop(t,mz,mt,mt,x[j]/x[k])*fhzl);
	return s;
}

double Sbly(double t, int j, int k, double *x, double ftr, double fwl){/*Yukawa interactions with B = b_L*/
	double s;
	s = prop2(t,mt,0,mw,x[j]/x[k])*Pffy(x[j]/x[k])*ftr+3*prop2(t,mw,0,mt,x[j]/x[k])*Pfh(x[j]/x[k])*fwl;
	return s;
}

double Stry(double t, int j, int k, double *x, double ftl, double fbl, double fh, double fzl, double fhzl, double fwl){/*Yukawa interactions with A and B = t_R*/
	double s;
	s = Pffy(x[j]/x[k])*0.5*(prop2(t,mt,mt,mh,x[j]/x[k])+prop2(t,mt,mt,mz,x[j]/x[k]))*ftl+prop2(t,0,mt,mw,x[j]/x[k])*Pffy(x[j]/x[k])*fbl+3*Pfh(x[j]/x[k])*(0.5*prop2(t,mh,mt,mt,x[j]/x[k])*fh+0.5*prop2(t,mz,mt,mt,x[j]/x[k])*fzl+0.5*prop(t,mh,mt,mt,x[j]/x[k])*prop(t,mz,mt,mt,x[j]/x[k])*fhzl+prop2(t,mw,mt,0,x[j]/x[k])*fwl);
	return s;
}

double Sswlwt(double t, int j, int k, double *x, double fwpl, double fwml, double fwpp, double fwpm, double fwmp, double fwmm, double ms){/*B = scalar (h or Z_L) with a longitudinal and a transverse W*/
	double s;
	s = cut(x[j],x[k],mu(t))*prop2(t,mw,ms,mw,x[j]/x[k])*Phh(x[j]/x[k])*(fwpl+fwml)+prop2(t,mw,ms,mw,x[j]/x[k])*Phv(x[j]/x[k])*(fwpp+fwmp+fwpm+fwmm);
	return s;
}

double Shzlwlwt(double t, int j, int k, double *x, double fwpl, double fwml, double fwpp, double fwpm, double fwmp, double fwmm){/*B = mixed h-Z_L with a longitudinal and a transverse W*/
	double s;
	s = prop(t,mw,mh,mw,x[j]/x[k])*prop(t,mw,mz,mw,x[j]/x[k])*cut(x[j],x[k],mu(t))*Phh(x[j]/x[k])*(fwml-fwpl)+prop(t,mw,mh,mw,x[j]/x[k])*prop(t,mw,mz,mw,x[j]/x[k])*Phv(x[j]/x[k])*(fwpp-fwmp+fwpm-fwmm);
	return s;
}

double Ssszt(double t, int j, int k, double *x, double fs, double fzp, double fzm, double msb, double ms){/*B = scalar (h or Z_L) with itself and a transverse Z*/
	double s;
	s = cut(x[j],x[k],mu(t))*prop2(t,ms,msb,mz,x[j]/x[k])*Phh(x[j]/x[k])*fs+prop2(t,mz,msb,ms,x[j]/x[k])*Phv(x[j]/x[k])*(fzp+fzm);
	return s;
}

double Sstt(double t, int j, int k, double *x, double ftl, double ftlb, double ftr, double ftrb, double ms){/*B = scalar (h or Z_L) with two top quarks*/
	double s;
	s = prop2(t,mt,ms,mt,x[j]/x[k])*Phf(x[j]/x[k])*(ftl+ftlb+ftr+ftrb);
	return s;
}

double Shzltt(double t, int j, int k, double *x, double ftl, double ftlb, double ftr, double ftrb){/*B = mixed h-Z_L with two top quarks*/
	double s;
	s = prop(t,mt,mh,mt,x[j]/x[k])*prop(t,mt,mz,mt,x[j]/x[k])*Phf(x[j]/x[k])*(ftl-ftlb-ftr+ftrb);
	return s;
}

double Swlswt(double t, int j, int k, double *x, double fh, double fzl, double fhzl, double fwp, double fwm){/*B = W_L, with a scalar (h, Z_L or h-Z_L) and a transverse W*/
	double s;
	s = cut(x[j],x[k],mu(t))*Phh(x[j]/x[k])*(prop2(t,mh,mw,mw,x[j]/x[k])*theta(t,th)*fh+prop2(t,mz,mw,mw,x[j]/x[k])*theta(t,tz)*fzl+prop(t,mh,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*theta(t,th)*fhzl)+(prop2(t,mw,mw,mh,x[j]/x[k])*theta(t,th)+prop2(t,mw,mw,mz,x[j]/x[k])*theta(t,tz))*Phv(x[j]/x[k])*(fwp+fwm);
	return s;
}

double Swlwlvt(double t, int j, int k, double *x, double fwl, double fvp, double fvm, double mv){/*B = W_L with itself and a transverse vector (Z or photon)*/
	double s;
	s = prop2(t,mw,mw,mv,x[j]/x[k])*Phh(x[j]/x[k])*fwl+prop2(t,mv,mw,mw,x[j]/x[k])*Phv(x[j]/x[k])*(fvp+fvm);
	return s;
}

double Swlwlzph(double t, int j, int k, double *x, double fzphp, double fzphm){/*B = W_L with itself and the mixed Z-photon*/
	double s;
	s = prop(t,0,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*Phv(x[j]/x[k])*(fzphp+fzphm);
	return s;
}

double Swltb(double t, int j, int k, double *x, double ft, double fb){/*Yukawa interactions with B=W_L*/
	double s;
	s = Phf(x[j]/x[k])*(prop2(t,mt,mw,0,x[j]/x[k])*ft+prop2(t,0,mw,mt,x[j]/x[k])*fb);
	return s;
}

double Svvv(double t, int j, int k, double *x, double fvp, double fvm, double ma, double mb, double mc){ /*Three transverse gauge bosons*/
	double s;
	s = prop2(t,ma,mb,mc,x[j]/x[k])*(Pvpvp(x[j]/x[k])*fvp+Pvpvm(x[j]/x[k])*fvm);
	return s;
}

double Svtss(double t, int j, int k, double *x, double fs, double msa, double mv, double msc){/*B = transverse vector and two scalars*/
	double s;
	s = prop2(t,msa,mv,msc,x[j]/x[k])*Pvh(x[j]/x[k])*fs;
	return s;
}

double Swthzlwl(double t, int j, int k, double *x, double fhzl){/*B = W_T, A = h-Z_L, C = W_L*/
	double s;
	s = prop(t,mh,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*Pvh(x[j]/x[k])*fhzl;
	return s;
}


double Szphwtwt(double t, int j, int k, double *x, double fvp, double fvm){/*B = Z-photon, A and C = transverse W*/
	double s;
	s = prop(t,mw,0,mw,x[j]/x[k])*prop(t,mw,mz,mw,x[j]/x[k])*(Pvpvp(x[j]/x[k])*fvp+Pvpvm(x[j]/x[k])*fvm);
	return s;
}

double Swtzphwt(double t, int j, int k, double *x, double fvp, double fvm){/*A = Z-photon, B and C = transverse W*/
	double s;
	s = prop(t,0,mw,mw,x[j]/x[k])*prop(t,mz,mw,mw,x[j]/x[k])*(Pvpvp(x[j]/x[k])*fvp+Pvpvm(x[j]/x[k])*fvm);
	return s;
}

double Szphwlwl(double t, int j, int k, double *x, double fs){/*B = Z-photon, A and C = longitudinal W*/
	double s;
	s = prop(t,mw,0,mw,x[j]/x[k])*prop(t,mw,mz,mw,x[j]/x[k])*Pvh(x[j]/x[k])*fs;
	return s;
}

double Sphpff(double t, int j, int k, double *x, double fel, double fml, double felb, double fer, double fmr, double ferb, double ful, double fulb, double ftl, double ftlb, double fdl, double fdlb, double fbl, double fblb, double fur, double furb, double ftr, double ftrb, double fdr, double fdrb, double fbr, double fbrb){/*B = photon_+, A and C = fermions*/
	double s;
	s = Pvpfl(x[j]/x[k])*(Ql*Ql*(2*fel+fml+3*ferb)+Qu*Qu*(2*ful+2*furb+theta(t,tt)*prop2(t,mt,0,mt,x[j]/x[k])*ftl+prop2(t,mt,0,mt,x[j]/x[k])*theta(t,tt)*ftrb)+Qd*Qd*(2*fdl+fbl+2*fdrb+fbrb))+Pvmfl(x[j]/x[k])*(Ql*Ql*(3*felb+2*fer+fmr)+Qu*Qu*(2*fulb+2*fur+prop2(t,mt,0,mt,x[j]/x[k])*theta(t,tt)*ftlb+prop2(t,mt,0,mt,x[j]/x[k])*theta(t,tt)*ftr)+Qd*Qd*(2*fdlb+fblb+2*fdr+fbr));
	return s;
}

double Sphmff(double t, int j, int k, double *x, double fel, double fml, double felb, double fer, double fmr, double ferb, double ful, double fulb, double ftl, double ftlb, double fdl, double fdlb, double fbl, double fblb, double fur, double furb, double ftr, double ftrb, double fdr, double fdrb, double fbr, double fbrb){/*B = photon_-, A and C = fermions*/
	double s;
	s = Pvmfl(x[j]/x[k])*(Ql*Ql*(2*fel+fml+3*ferb)+Qu*Qu*(2*ful+2*furb+theta(t,tt)*prop2(t,mt,0,mt,x[j]/x[k])*ftl+prop2(t,mt,0,mt,x[j]/x[k])*theta(t,tt)*ftrb)+Qd*Qd*(2*fdl+fbl+2*fdrb+fbrb))+Pvpfl(x[j]/x[k])*(Ql*Ql*(3*felb+2*fer+fmr)+Qu*Qu*(2*fulb+2*fur+prop2(t,mt,0,mt,x[j]/x[k])*theta(t,tt)*ftlb+prop2(t,mt,0,mt,x[j]/x[k])*theta(t,tt)*ftr)+Qd*Qd*(2*fdlb+fblb+2*fdr+fbr));
	return s;
}

double Sgpff(double t, int j, int k, double *x, double ful, double fulb, double ftl, double ftlb, double fdl, double fdlb, double fbl, double fblb, double fur, double furb, double ftr, double ftrb, double fdr, double fdrb, double fbr, double fbrb){/*B = g_+, A and C = quarks*/
	double s;
	s = Pvpfl(x[j]/x[k])*(2*ful+2*fdl+theta(t,tt)*prop2(t,mt,0,mt,x[j]/x[k])*ftl+fbl+2*furb+2*fdrb+theta(t,tt)*prop2(t,mt,0,mt,x[j]/x[k])*ftrb+fbrb)+Pvmfl(x[j]/x[k])*(2*fulb+2*fdlb+theta(t,tt)*prop2(t,mt,0,mt,x[j]/x[k])*ftlb+fblb+2*fur+2*fdr+theta(t,tt)*prop2(t,mt,0,mt,x[j]/x[k])*ftr+fbr);
	return s;
}

double Sgmff(double t, int j, int k, double *x, double ful, double fulb, double ftl, double ftlb, double fdl, double fdlb, double fbl, double fblb, double fur, double furb, double ftr, double ftrb, double fdr, double fdrb, double fbr, double fbrb){/*B = g_-, A and C = quarks*/
	double s;
	s = Pvmfl(x[j]/x[k])*(2*ful+2*fdl+theta(t,tt)*prop2(t,mt,0,mt,x[j]/x[k])*ftl+fbl+2*furb+2*fdrb+theta(t,tt)*prop2(t,mt,0,mt,x[j]/x[k])*ftrb+fbrb)+Pvpfl(x[j]/x[k])*(2*fulb+2*fdlb+theta(t,tt)*prop2(t,mt,0,mt,x[j]/x[k])*ftlb+fblb+2*fur+2*fdr+theta(t,tt)*prop2(t,mt,0,mt,x[j]/x[k])*ftr+fbr);
	return s;
}

double Szpff(double t, int j, int k, double *x, double fnue, double fnum, double fnub, double fel, double fml, double felb, double fer, double fmr, double ferb, double ful, double fulb, double ftl, double ftlb, double fdl, double fdlb, double fbl, double fblb, double fur, double furb, double ftr, double ftrb, double fdr, double fdrb, double fbr, double fbrb){/*B = Z_+, A and C = fermions*/
	double s;
	s = Pvpfl(x[j]/x[k])*(prop2(t,0,mz,0,x[j]/x[k])*0.25*(2*fnue+fnum)+prop2(t,0,mz,0,x[j]/x[k])*Qzdl(t,Ql)*Qzdl(t,Ql)*(2*fel+fml)+prop2(t,0,mz,0,x[j]/x[k])*Qzr(t,Ql)*Qzr(t,Ql)*3*ferb+Qzul(t,Qu)*Qzul(t,Qu)*(prop2(t,0,mz,0,x[j]/x[k])*2*ful+theta(t,tt)*prop2(t,mt,mz,mt,x[j]/x[k])*ftl)+prop2(t,0,mz,0,x[j]/x[k])*Qzdl(t,Qd)*Qzdl(t,Qd)*(2*fdl+fbl)+Qzr(t,Qu)*Qzr(t,Qu)*(prop2(t,0,mz,0,x[j]/x[k])*2*furb+theta(t,tt)*prop2(t,mt,mz,mt,x[j]/x[k])*ftrb)+prop2(t,0,mz,0,x[j]/x[k])*Qzr(t,Qd)*Qzr(t,Qd)*(2*fdrb+fbrb))+Pvmfl(x[j]/x[k])*(prop2(t,0,mz,0,x[j]/x[k])*0.25*3*fnub+prop2(t,0,mz,0,x[j]/x[k])*Qzdl(t,Ql)*Qzdl(t,Ql)*3*felb+prop2(t,0,mz,0,x[j]/x[k])*Qzr(t,Ql)*Qzr(t,Ql)*(2*fer+fmr)+Qzul(t,Qu)*Qzul(t,Qu)*(prop2(t,0,mz,0,x[j]/x[k])*2*fulb+theta(t,tt)*prop2(t,mt,mz,mt,x[j]/x[k])*ftlb)+prop2(t,0,mz,0,x[j]/x[k])*Qzdl(t,Qd)*Qzdl(t,Qd)*(2*fdlb+fblb)+Qzr(t,Qu)*Qzr(t,Qu)*(prop2(t,0,mz,0,x[j]/x[k])*2*fur+theta(t,tt)*prop2(t,mt,mz,mt,x[j]/x[k])*ftr)+prop2(t,0,mz,0,x[j]/x[k])*Qzr(t,Qd)*Qzr(t,Qd)*(2*fdr+fbr));
	return s;
}

double Szmff(double t, int j, int k, double *x, double fnue, double fnum, double fnub, double fel, double fml, double felb, double fer, double fmr, double ferb, double ful, double fulb, double ftl, double ftlb, double fdl, double fdlb, double fbl, double fblb, double fur, double furb, double ftr, double ftrb, double fdr, double fdrb, double fbr, double fbrb){/*B = Z_-, A and C = fermions*/
	double s;
	s = Pvmfl(x[j]/x[k])*(prop2(t,0,mz,0,x[j]/x[k])*0.25*(2*fnue+fnum)+prop2(t,0,mz,0,x[j]/x[k])*Qzdl(t,Ql)*Qzdl(t,Ql)*(2*fel+fml)+prop2(t,0,mz,0,x[j]/x[k])*Qzr(t,Ql)*Qzr(t,Ql)*3*ferb+Qzul(t,Qu)*Qzul(t,Qu)*(prop2(t,0,mz,0,x[j]/x[k])*2*ful+theta(t,tt)*prop2(t,mt,mz,mt,x[j]/x[k])*ftl)+prop2(t,0,mz,0,x[j]/x[k])*Qzdl(t,Qd)*Qzdl(t,Qd)*(2*fdl+fbl)+Qzr(t,Qu)*Qzr(t,Qu)*(prop2(t,0,mz,0,x[j]/x[k])*2*furb+theta(t,tt)*prop2(t,mt,mz,mt,x[j]/x[k])*ftrb)+prop2(t,0,mz,0,x[j]/x[k])*Qzr(t,Qd)*Qzr(t,Qd)*(2*fdrb+fbrb))+Pvpfl(x[j]/x[k])*(prop2(t,0,mz,0,x[j]/x[k])*0.25*3*fnub+prop2(t,0,mz,0,x[j]/x[k])*Qzdl(t,Ql)*Qzdl(t,Ql)*3*felb+prop2(t,0,mz,0,x[j]/x[k])*Qzr(t,Ql)*Qzr(t,Ql)*(2*fer+fmr)+Qzul(t,Qu)*Qzul(t,Qu)*(prop2(t,0,mz,0,x[j]/x[k])*2*fulb+theta(t,tt)*prop2(t,mt,mz,mt,x[j]/x[k])*ftlb)+prop2(t,0,mz,0,x[j]/x[k])*Qzdl(t,Qd)*Qzdl(t,Qd)*(2*fdlb+fblb)+Qzr(t,Qu)*Qzr(t,Qu)*(prop2(t,0,mz,0,x[j]/x[k])*2*fur+theta(t,tt)*prop2(t,mt,mz,mt,x[j]/x[k])*ftr)+prop2(t,0,mz,0,x[j]/x[k])*Qzr(t,Qd)*Qzr(t,Qd)*(2*fdr+fbr));
	return s;
}

double Szphpff(double t, int j, int k, double *x, double fel, double fml, double felb, double fer, double fmr, double ferb, double ful, double fulb, double ftl, double ftlb, double fdl, double fdlb, double fbl, double fblb, double fur, double furb, double ftr, double ftrb, double fdr, double fdrb, double fbr, double fbrb){/*B = Z-photon_+, A and C = fermions*/
	double s;
	s = prop(t,0,mz,0,x[j]/x[k])*Pvpfl(x[j]/x[k])*(Ql*Qzdl(t,Ql)*(2*fel+fml)+Qu*Qzul(t,Qu)*2*ful+Qd*Qzdl(t,Qd)*(2*fdl+fbl)+Ql*Qzr(t,Ql)*3*ferb+Qu*Qzr(t,Qu)*2*furb+Qd*Qzr(t,Qd)*(2*fdrb+fbrb))+prop(t,0,mz,0,x[j]/x[k])*Pvmfl(x[j]/x[k])*(Ql*Qzdl(t,Ql)*3*felb+Qu*Qzul(t,Qu)*2*fulb+Qd*Qzdl(t,Qd)*(2*fdlb+fblb)+Ql*Qzr(t,Ql)*(2*fer+fmr)+Qu*Qzr(t,Qu)*2*fur+Qd*Qzr(t,Qd)*(2*fdr+fbr))+theta(t,tt)*prop(t,mt,mz,mt,x[j]/x[k])*prop(t,mt,0,mt,x[j]/x[k])*(Pvpfl(x[j]/x[k])*(Qu*Qzul(t,Qu)*ftl+Qu*Qzr(t,Qu)*ftrb)+Pvmfl(x[j]/x[k])*(Qu*Qzul(t,Qu)*ftlb+Qu*Qzr(t,Qu)*ftr));
	return s;
}

double Szphmff(double t, int j, int k, double *x, double fel, double fml, double felb, double fer, double fmr, double ferb, double ful, double fulb, double ftl, double ftlb, double fdl, double fdlb, double fbl, double fblb, double fur, double furb, double ftr, double ftrb, double fdr, double fdrb, double fbr, double fbrb){/*B = Z-photon_-, A and C = fermions*/
	double s;
	s = prop(t,0,mz,0,x[j]/x[k])*Pvmfl(x[j]/x[k])*(Ql*Qzdl(t,Ql)*(2*fel+fml)+Qu*Qzul(t,Qu)*2*ful+Qd*Qzdl(t,Qd)*(2*fdl+fbl)+Ql*Qzr(t,Ql)*3*ferb+Qu*Qzr(t,Qu)*2*furb+Qd*Qzr(t,Qd)*(2*fdrb+fbrb))+prop(t,0,mz,0,x[j]/x[k])*Pvpfl(x[j]/x[k])*(Ql*Qzdl(t,Ql)*3*felb+Qu*Qzul(t,Qu)*2*fulb+Qd*Qzdl(t,Qd)*(2*fdlb+fblb)+Ql*Qzr(t,Ql)*(2*fer+fmr)+Qu*Qzr(t,Qu)*2*fur+Qd*Qzr(t,Qd)*(2*fdr+fbr))+theta(t,tt)*prop(t,mt,mz,mt,x[j]/x[k])*prop(t,mt,0,mt,x[j]/x[k])*(Pvmfl(x[j]/x[k])*(Qu*Qzul(t,Qu)*ftl+Qu*Qzr(t,Qu)*ftrb)+Pvpfl(x[j]/x[k])*(Qu*Qzul(t,Qu)*ftlb+Qu*Qzr(t,Qu)*ftr));
	return s;
}

double Swppff(double t, int j, int k, double *x, double fnue, double fnum, double felb, double ful, double ftl, double fdlb, double fblb){/*B = W^+_+, A and C = fermions*/
	double s;
	s = Pvpfl(x[j]/x[k])*(prop2(t,0,mw,0,x[j]/x[k])*2*ful+prop2(t,mt,mw,0,x[j]/x[k])*theta(t,tt)*ftl+prop2(t,0,mw,0,x[j]/x[k])*2*fnue+prop2(t,0,mw,0,x[j]/x[k])*fnum)+Pvmfl(x[j]/x[k])*(prop2(t,0,mw,0,x[j]/x[k])*2*fdlb+prop2(t,0,mw,mt,x[j]/x[k])*fblb+prop2(t,0,mw,0,x[j]/x[k])*3*felb);
	return s;
}

double Swpmff(double t, int j, int k, double *x, double fnue, double fnum, double felb, double ful, double ftl, double fdlb, double fblb){/*B = W^+_-, A and C = fermions*/
	double s;
	s = Pvmfl(x[j]/x[k])*(prop2(t,0,mw,0,x[j]/x[k])*2*ful+prop2(t,mt,mw,0,x[j]/x[k])*theta(t,tt)*ftl+prop2(t,0,mw,0,x[j]/x[k])*2*fnue+prop2(t,0,mw,0,x[j]/x[k])*fnum)+Pvpfl(x[j]/x[k])*(prop2(t,0,mw,0,x[j]/x[k])*2*fdlb+prop2(t,0,mw,mt,x[j]/x[k])*fblb+prop2(t,0,mw,0,x[j]/x[k])*3*felb);
	return s;
}

double Swmpff(double t, int j, int k, double *x, double fnub, double fel, double fml, double fulb, double ftlb, double fdl, double fbl){/*B = W^-_+, A and C = fermions*/
	double s;
	s = Pvmfl(x[j]/x[k])*(prop2(t,0,mw,0,x[j]/x[k])*2*fulb+prop2(t,mt,mw,0,x[j]/x[k])*theta(t,tt)*ftlb+prop2(t,0,mw,0,x[j]/x[k])*3*fnub)+Pvpfl(x[j]/x[k])*(prop2(t,0,mw,0,x[j]/x[k])*2*fdl+prop2(t,0,mw,mt,x[j]/x[k])*fbl+prop2(t,0,mw,0,x[j]/x[k])*2*fel+prop2(t,0,mw,0,x[j]/x[k])*fml);
	return s;
}

double Swmmff(double t, int j, int k, double *x, double fnub, double fel, double fml, double fulb, double ftlb, double fdl, double fbl){/*B = W^-_-, A and C = fermions*/
	double s;
	s = Pvpfl(x[j]/x[k])*(prop2(t,0,mw,0,x[j]/x[k])*2*fulb+prop2(t,mt,mw,0,x[j]/x[k])*theta(t,tt)*ftlb+prop2(t,0,mw,0,x[j]/x[k])*3*fnub)+Pvmfl(x[j]/x[k])*(prop2(t,0,mw,0,x[j]/x[k])*2*fdl+prop2(t,0,mw,mt,x[j]/x[k])*fbl+prop2(t,0,mw,0,x[j]/x[k])*2*fel+prop2(t,0,mw,0,x[j]/x[k])*fml);
	return s;
}
