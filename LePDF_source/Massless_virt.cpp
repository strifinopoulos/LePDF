/*Contributions to the coefficients P_A^v. Notation: Iv (Integral virtual) followed by the partons A, B, C*/

double Ivffzt(double t){
	double I;
	I = -(pow(mz,2)*pow(mu(t),2)*(pow(mz,4)+pow(mu(t),4))+(-(pow(mz,4)*pow(mu(t),4))+2*pow(mz,2)*pow(mu(t),6)+pow(mu(t),8))*log(pow(mu(t),2)/(pow(mz,2)+pow(mu(t),2))))/(pow(mz,4)*pow(pow(mz,2)+pow(mu(t),2),2));
	return I;
}

double Ivffwt(double t){
	double I;
	if(t>22.9158){
		I = 4*mew/mu(t)-0.5*(3+4*log(mew/mu(t)))+(34*mw*mw-9*mew*mew+24*mw*mw*log(mew/mu(t)))/(6*mu(t)*mu(t));
	}
	else{
		I = -(pow(mu(t),6)*(3-(2*(2*mew+pow(mw,2)/mu(t)))/(mew+(mew*pow(mw,2))/pow(mu(t),2))+(mew-(mew*pow(mw,2))/pow(mu(t),2)+(2*pow(mw,2))/mu(t))/(mew+(mew*pow(mw,2)*(1-mew/mu(t)))/pow(mu(t),2))+pow(mw,4)/pow(mu(t),4)+(mew*pow(mw,2))/pow(mu(t),3)-(2*pow(mw,2))/pow(mu(t),2)-2*log(1+pow(mw,2)/pow(mu(t),2))+((2+pow(mw,6)/pow(mu(t),6)+(3*pow(mw,2))/pow(mu(t),2))*log(1+(pow(mw,2)*(1-mew/mu(t)))/pow(mu(t),2)))/pow(1+pow(mw,2)/pow(mu(t),2),2)))/pow(mw,6)-(2*log(mew/mu(t)))/pow(1+pow(mw,2)/pow(mu(t),2),2);
	}
	return I;
}

double Ivttv(double t, double mv){
	double I;
	I = -(pow(mu(t),4)*(-(((pow(-2*pow(mt,2)*mv + pow(mv,3),2) + 4*pow(mt,4)*pow(mu(t),2) + (2*pow(mt,2) + pow(mv,2))*pow(mu(t),4))*
            sqrt(-pow(mv,4) + 4*pow(mt,2)*(pow(mv,2) + pow(mu(t),2))))/(pow(mt,2) + pow(mu(t),2))) + 
       2*(2*pow(mt,2)*(pow(mv,2) - 2*pow(mu(t),2))*(pow(mv,2) + pow(mu(t),2)) + pow(mv,2)*pow(mu(t),2)*(2*pow(mv,2) + pow(mu(t),2)))*
        atan(pow(mv,2)/sqrt(-pow(mv,4) + 4*pow(mt,2)*(pow(mv,2) + pow(mu(t),2)))) + 
       2*(2*pow(mt,2)*(pow(mv,2) - 2*pow(mu(t),2))*(pow(mv,2) + pow(mu(t),2)) + pow(mv,2)*pow(mu(t),2)*(2*pow(mv,2) + pow(mu(t),2)))*
        atan((2*pow(mt,2) - pow(mv,2))/sqrt(-pow(mv,4) + 4*pow(mt,2)*(pow(mv,2) + pow(mu(t),2)))) + 
       pow(-pow(mv,4) + 4*pow(mt,2)*(pow(mv,2) + pow(mu(t),2)),1.5)*(-log(pow(mt,2) + pow(mu(t),2)) + log(pow(mv,2) + pow(mu(t),2)))))/
   (pow(pow(mv,2) + pow(mu(t),2),2)*pow(-pow(mv,4) + 4*pow(mt,2)*(pow(mv,2) + pow(mu(t),2)),1.5));
   return I;
}

double Ivtlbwt(double t){
	double I;
	I = -(-((pow(mu(t),4)*(pow(mw,2)/pow(mt,2) + (2*pow(mt,2))/(pow(mw,2) + pow(mu(t),2)) - (pow(mw,2) + 3*pow(mu(t),2))/(pow(mw,2) + pow(mu(t),2)) - 
          (4*(pow(mt,2) + pow(mu(t),2))*atan((-pow(mt,2) - pow(mw,2))/sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/
           sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/(pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))) + 
   (pow(mu(t),4)*(pow(mt,2) + pow(mu(t),2))*((-pow(mt,2) + pow(mw,2))/(pow(mt,2)*pow(mu(t),2)) - 
        (4*atan((pow(mt,2) - pow(mw,2))/sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/
         sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/(pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2))) - 
   (pow(mu(t),4)*(((pow(mw,2) + pow(mu(t),2))*(2*pow(mt,4) + pow(mw,2)*(pow(mw,2) + pow(mu(t),2)) - pow(mt,2)*(pow(mw,2) + 3*pow(mu(t),2))))/
         (pow(mt,2)*(pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))) - 
        (2*(pow(mt,6) + 3*pow(mw,2)*pow(mu(t),4) + 2*pow(mu(t),6) - 3*pow(mt,4)*(pow(mw,2) + 2*pow(mu(t),2)) + 
             pow(mt,2)*(2*pow(mw,4) + 4*pow(mw,2)*pow(mu(t),2) + 5*pow(mu(t),4)))*
           atan((-pow(mt,2) + pow(mw,2))/sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/
         pow(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)),1.5) + log(pow(mu(t),2))))/pow(pow(mw,2) + pow(mu(t),2),2) + 
   (pow(mu(t),4)*((mu(t)*(pow(mw,2) + pow(mu(t),2))*(mew*(-2*pow(mt,6) + pow(mw,4)*(pow(mw,2) + pow(mu(t),2)) + pow(mt,4)*(4*pow(mw,2) + 6*pow(mu(t),2)) - 
                pow(mt,2)*(3*pow(mw,4) + 5*pow(mw,2)*pow(mu(t),2) + 2*pow(mu(t),4))) + 
             mu(t)*(2*pow(mt,6) - pow(mw,2)*pow(pow(mw,2) + pow(mu(t),2),2) - 4*pow(mt,4)*(pow(mw,2) + 2*pow(mu(t),2)) + 
                pow(mt,2)*(3*pow(mw,4) + 6*pow(mw,2)*pow(mu(t),2) + 5*pow(mu(t),4)))))/
         (pow(mt,2)*(-(pow(mew,2)*pow(mt,2)) + mew*(pow(mt,2) + pow(mw,2))*mu(t) - pow(mu(t),2)*(pow(mw,2) + pow(mu(t),2)))*
           (pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))) - 
        (2*(pow(mt,6) + 3*pow(mw,2)*pow(mu(t),4) + 2*pow(mu(t),6) - 3*pow(mt,4)*(pow(mw,2) + 2*pow(mu(t),2)) + 
             pow(mt,2)*(2*pow(mw,4) + 4*pow(mw,2)*pow(mu(t),2) + 5*pow(mu(t),4)))*
           atan((-2*mew*pow(mt,2) + (pow(mt,2) + pow(mw,2))*mu(t))/(mu(t)*sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2))))))/
         pow(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)),1.5) - 2*log(mew/mu(t)) + 
        log((mew*pow(mt,2)*(mew - mu(t)) - pow(mw,2)*(mew - mu(t))*mu(t) + pow(mu(t),4))/pow(mu(t),2))))/pow(pow(mw,2) + pow(mu(t),2),2));
   return I;
}


double Ivtts(double t, double ms){
	double I;
	I = 0.5*((pow(ms,2) - 2*pow(mt,2))*pow(mu(t),4))/((pow(mt,2) + pow(mu(t),2))*(-pow(ms,4) + 4*pow(ms,2)*pow(mt,2) + 4*pow(mt,2)*pow(mu(t),2))) - 
   (pow(ms,2)*pow(mu(t),4)*(atan(pow(ms,2)/sqrt(-pow(ms,4) + 4*pow(ms,2)*pow(mt,2) + 4*pow(mt,2)*pow(mu(t),2))) - 
        atan((pow(ms,2) - 2*pow(mt,2))/sqrt(-pow(ms,4) + 4*pow(ms,2)*pow(mt,2) + 4*pow(mt,2)*pow(mu(t),2)))))/
    pow(-pow(ms,4) + 4*pow(ms,2)*pow(mt,2) + 4*pow(mt,2)*pow(mu(t),2),1.5);
   return I;
}

double Ivtrbwl(double t){
	double I;
	I = -((-pow(mt,2) + pow(mw,2))*pow(mu(t),2))/(2*pow(pow(mt,2) - pow(mw,2),2) - 8*pow(mt,2)*pow(mu(t),2)) -
   ((pow(mt,2) + pow(mw,2))*pow(mu(t),4)*(atan(((mt - mw)*(mt + mw))/sqrt(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mt,2)*pow(mu(t),2))) + 
        atan((pow(mt,2) + pow(mw,2))/sqrt(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mt,2)*pow(mu(t),2)))))/
    pow(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mt,2)*pow(mu(t),2),1.5);
   return I;
}

double Ivbltrwl(double t){
	double I;
	I = -(pow(mu(t),4)*((-pow(mt,2) + pow(mw,2))/pow(mu(t),2) + (1 + pow(mt,2)/pow(mu(t),2))*log(1 + pow(mt,2)/pow(mu(t),2)) - 
       (1 + pow(mt,2)/pow(mu(t),2))*log(1 + pow(mw,2)/pow(mu(t),2))))/(2.*pow(pow(mt,2) - pow(mw,2),2)*(1 + pow(mt,2)/pow(mu(t),2)));
   return I;
}

double Ivbltlwt(double t){
	double I;
	if(t>22.9158){
		I = -((4*mew)/mu(t) + (-3 - 4*log(mew/mu(t)))/2. + (-9*pow(mew,2) - 16*pow(mt,2) + 34*pow(mw,2) + 24*pow(mw,2)*log(mew/mu(t)))/(6.*pow(mu(t),2)));
	}
	else{
		I = -(((-1 + mew/mu(t))*pow(mu(t),4))/pow(pow(mt,2) - pow(mw,2),2) + ((1 + pow(mt,2)*((-2*pow(mw,2))/pow(mu(t),4) + 2/pow(mu(t),2)) + (2*pow(mt,4))/pow(mu(t),4) + 
        pow(mw,4)/pow(mu(t),4))*pow(mu(t),6))/(pow(pow(mt,2) - pow(mw,2),3)*(1 + pow(mw,2)/pow(mu(t),2))) - 
   ((1 + pow(mt,2)*((-2*pow(mw,2))/pow(mu(t),4) + 2/pow(mu(t),2)) + (2*pow(mt,4))/pow(mu(t),4) + pow(mw,4)/pow(mu(t),4))*(1 + pow(mt,2)/pow(mu(t),2))*pow(mu(t),6))/
    (pow(pow(mt,2) - pow(mw,2),3)*(1 + pow(mw,2)*(-(mew/pow(mu(t),3)) + pow(mu(t),-2)) + (mew*pow(mt,2))/pow(mu(t),3))*(1 + pow(mw,2)/pow(mu(t),2))) + 
   ((2 + pow(mt,2)*((-3*pow(mw,4))/pow(mu(t),6) + (6*pow(mw,2))/pow(mu(t),4) + 3/pow(mu(t),2)) - (2*pow(mt,6))/pow(mu(t),6) + (6*pow(mt,4)*pow(mw,2))/pow(mu(t),6) + 
        pow(mw,6)/pow(mu(t),6) + (3*pow(mw,2))/pow(mu(t),2))*pow(mu(t),6)*log(1 + pow(mt,2)/pow(mu(t),2)))/(pow(pow(mt,2) - pow(mw,2),3)*pow(1 + pow(mw,2)/pow(mu(t),2),2))
     - (pow(mu(t),6)*(1 + pow(mt,2)*(-(pow(mw,2)/pow(mu(t),4)) + pow(mu(t),-2)) + pow(mt,4)/pow(mu(t),4) + pow(mw,4)/pow(mu(t),4) + pow(mw,2)/pow(mu(t),2) + 
        2*pow(1 + pow(mt,2)/pow(mu(t),2),2)*log(1 + pow(mt,2)/pow(mu(t),2))))/(pow(pow(mt,2) - pow(mw,2),3)*(1 + pow(mt,2)/pow(mu(t),2))) - 
   (pow(mu(t),6)*((1 + pow(mt,2)*((-2*pow(mw,2))/pow(mu(t),4) + 2/pow(mu(t),2)) + (2*pow(mt,4))/pow(mu(t),4) + pow(mw,4)/pow(mu(t),4))/(1 + pow(mw,2)/pow(mu(t),2)) + 
        2*(1 + pow(mt,2)/pow(mu(t),2))*log(1 + pow(mw,2)/pow(mu(t),2))))/pow(-pow(mt,2) + pow(mw,2),3) + 
   ((-2 + (2*pow(mt,6))/pow(mu(t),6) - (6*pow(mt,4)*pow(mw,2))/pow(mu(t),6) - pow(mw,6)/pow(mu(t),6) - (3*pow(mw,2))/pow(mu(t),2) + 
        (3*pow(mt,2)*(-1 + pow(mw,4)/pow(mu(t),4) - (2*pow(mw,2))/pow(mu(t),2)))/pow(mu(t),2))*pow(mu(t),6)*
      log(1 + (mew*pow(mt,2))/pow(mu(t),3) + (pow(mw,2)*(1 - mew/mu(t)))/pow(mu(t),2)))/(pow(pow(mt,2) - pow(mw,2),3)*pow(1 + pow(mw,2)/pow(mu(t),2),2)) - 
   (2*log(mew/mu(t)))/pow(1 + pow(mw,2)/pow(mu(t),2),2));
	}
	return I;
}

double Ivstt(double t, double ms){
	double I;
	I = -0.5*(pow(mu(t),4)*(1/((pow(mt,2) + pow(mu(t),2))*(-pow(ms,2) + 4*(pow(mt,2) + pow(mu(t),2)))) + 
     (4*atan(ms/sqrt(-pow(ms,2) + 4*(pow(mt,2) + pow(mu(t),2)))))/(ms*pow(-pow(ms,2) + 4*(pow(mt,2) + pow(mu(t),2)),1.5))));
	return I;
}

double Ivzlhzt(double t){
	double I;
	I = -(pow(mu(t),4)*((4*(sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2)) + 
          pow(mh,2)*atan(pow(mh,2)/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2)))))/
      pow(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2),1.5) + 
     (2*((pow(mh,2) + 2*pow(mu(t),2))/(pow(mz,2) + pow(mu(t),2)) + (2*pow(mh,2)*
             atan((pow(mh,2) - 2*pow(mz,2))/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))))/
           sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))))/(pow(mh,4) - 4*pow(mh,2)*pow(mz,2) - 4*pow(mz,2)*pow(mu(t),2)) - 
     ((2*(pow(mz,2) + pow(mu(t),2))*(pow(mh,2) + 2*pow(mu(t),2)))/(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2)) - 
        (2*(pow(mh,6) - 6*pow(mh,4)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2)*(pow(mz,2) - pow(mu(t),2)) + 
             2*pow(mh,2)*(2*pow(mz,4) - 5*pow(mz,2)*pow(mu(t),2) - pow(mu(t),4)))*atan(pow(mh,2)/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2)))
           )/pow(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2),1.5) + log(pow(mh,2) + pow(mu(t),2)))/pow(pow(mz,2) + pow(mu(t),2),2) + 
     ((-2*mu(t)*(pow(mz,2) + pow(mu(t),2))*(pow(mh,4)*mu(t) + 2*mew*pow(mz,2)*pow(mu(t),2) - 2*pow(mz,2)*pow(mu(t),3) + 2*pow(mu(t),5) + 2*pow(mh,2)*mu(t)*(-pow(mz,2) + pow(mu(t),2)) + 
             mew*pow(mh,2)*(2*pow(mz,2) + pow(mu(t),2))))/
         ((pow(mh,4) - 4*pow(mh,2)*pow(mz,2) - 4*pow(mz,2)*pow(mu(t),2))*
           (pow(mew,2)*pow(mz,2) + mew*(pow(mh,2) - 2*pow(mz,2))*mu(t) + pow(mu(t),2)*(pow(mz,2) + pow(mu(t),2)))) + 
        (2*(pow(mh,6) - 6*pow(mh,4)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2)*(pow(mz,2) - pow(mu(t),2)) + 
             2*pow(mh,2)*(2*pow(mz,4) - 5*pow(mz,2)*pow(mu(t),2) - pow(mu(t),4)))*
           atan((-pow(mh,2) + 2*pow(mz,2)*(1 - mew/mu(t)))/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))))/
         pow(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2),1.5) - 2*log(mew/mu(t)) + 
        log((pow(mew,2)*pow(mz,2) + mew*pow(mh,2)*mu(t) - 2*mew*pow(mz,2)*mu(t) + pow(mz,2)*pow(mu(t),2) + pow(mu(t),4))/pow(mu(t),2)))/pow(pow(mz,2) + pow(mu(t),2),2)));
	return I;
}

double Ivsvtvl(double t, double ms, double mv){
	double I;
	I = -(pow(mu(t),4)*((4*(sqrt(-pow(ms,2) + 4*(pow(mv,2) + pow(mu(t),2))) + ms*atan(ms/sqrt(-pow(ms,2) + 4*(pow(mv,2) + pow(mu(t),2))))))/
      (pow(ms,2)*pow(-pow(ms,2) + 4*(pow(mv,2) + pow(mu(t),2)),1.5)) + 
     (2*((-pow(ms,2) + 2*(pow(mv,2) + pow(mu(t),2)))/((pow(mv,2) + pow(mu(t),2))*(pow(ms,2) - 4*(pow(mv,2) + pow(mu(t),2)))) + 
          (2*ms*atan(ms/sqrt(-pow(ms,2) + 4*(pow(mv,2) + pow(mu(t),2)))))/pow(-pow(ms,2) + 4*(pow(mv,2) + pow(mu(t),2)),1.5)))/pow(ms,2) - 
     ((2*(pow(mv,2) + pow(mu(t),2))*(pow(ms,2) - 2*(pow(mv,2) + pow(mu(t),2))))/(pow(ms,4) - 4*pow(ms,2)*(pow(mv,2) + pow(mu(t),2))) + 
        (2*(pow(ms,4) - 6*pow(ms,2)*(pow(mv,2) + pow(mu(t),2)) + 6*pow(pow(mv,2) + pow(mu(t),2),2))*atan(ms/sqrt(-pow(ms,2) + 4*(pow(mv,2) + pow(mu(t),2)))))/
         (ms*pow(-pow(ms,2) + 4*(pow(mv,2) + pow(mu(t),2)),1.5)) + log(pow(mv,2) + pow(mu(t),2)))/pow(pow(mv,2) + pow(mu(t),2),2) - 
     ((2*(pow(mv,2) + pow(mu(t),2))*(pow(ms,4)*(1 - mew/mu(t)) + (pow(ms,2)*(3*mew - 4*mu(t))*(pow(mv,2) + pow(mu(t),2)))/mu(t) + 2*pow(pow(mv,2) + pow(mu(t),2),2)))/
         (pow(ms,2)*(pow(mv,2) + (mew*pow(ms,2)*(mew - mu(t)))/pow(mu(t),2) + pow(mu(t),2))*(pow(ms,2) - 4*(pow(mv,2) + pow(mu(t),2)))) + 
        (2*(pow(ms,4) - 6*pow(ms,2)*(pow(mv,2) + pow(mu(t),2)) + 6*pow(pow(mv,2) + pow(mu(t),2),2))*
           atan((ms*(-2*mew + mu(t)))/(mu(t)*sqrt(-pow(ms,2) + 4*(pow(mv,2) + pow(mu(t),2))))))/(ms*pow(-pow(ms,2) + 4*(pow(mv,2) + pow(mu(t),2)),1.5)) + 2*log(mew/mu(t)) - 
        log(pow(mv,2) + (mew*pow(ms,2)*(mew - mu(t)))/pow(mu(t),2) + pow(mu(t),2)))/pow(pow(mv,2) + pow(mu(t),2),2)));
	return I;
}

double Ivwltb(double t){
	double I;
	I = -((pow(mu(t),2)*(pow(mt,4) - pow(mt,2)*pow(mw,2) - 2*pow(mw,2)*pow(mu(t),2)))/
    (2.*(pow(mt,2) + pow(mu(t),2))*(pow(pow(mt,2) - pow(mw,2),2) - 4*pow(mw,2)*pow(mu(t),2))) - 
   (2*pow(mw,2)*pow(mu(t),4)*(atan(((mt - mw)*(mt + mw))/sqrt(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mw,2)*pow(mu(t),2))) - 
        atan((pow(mt,2) + pow(mw,2))/sqrt(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mw,2)*pow(mu(t),2)))))/
    pow(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mw,2)*pow(mu(t),2),1.5));
	if(t<tt){
    	return 0;
	}
	else{
		return I;
	}
}

double Ivwlwts(double t, double ms){
	double I;
	I = -(pow(mu(t),4)*((4*(sqrt(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2)) + 
          pow(ms,2)*atan(pow(ms,2)/sqrt(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2)))))/
      pow(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5) + 
     (2*((pow(ms,2) + 2*pow(mu(t),2))/(pow(mw,2) + pow(mu(t),2)) + (2*pow(ms,2)*
             atan((pow(ms,2) - 2*pow(mw,2))/sqrt(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))/
           sqrt(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))/(pow(ms,4) - 4*pow(ms,2)*pow(mw,2) - 4*pow(mw,2)*pow(mu(t),2)) - 
     ((2*(pow(mw,2) + pow(mu(t),2))*(pow(ms,2) + 2*pow(mu(t),2)))/(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2)) - 
        (2*(pow(ms,6) - 6*pow(ms,4)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2)*(pow(mw,2) - pow(mu(t),2)) + 
             2*pow(ms,2)*(2*pow(mw,4) - 5*pow(mw,2)*pow(mu(t),2) - pow(mu(t),4)))*atan(pow(ms,2)/sqrt(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2)))
           )/pow(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5) + log(pow(ms,2) + pow(mu(t),2)))/pow(pow(mw,2) + pow(mu(t),2),2) + 
     ((-2*mu(t)*(pow(mw,2) + pow(mu(t),2))*(pow(ms,4)*mu(t) + 2*mew*pow(mw,2)*pow(mu(t),2) - 2*pow(mw,2)*pow(mu(t),3) + 2*pow(mu(t),5) + 2*pow(ms,2)*mu(t)*(-pow(mw,2) + pow(mu(t),2)) + 
             mew*pow(ms,2)*(2*pow(mw,2) + pow(mu(t),2))))/
         ((pow(ms,4) - 4*pow(ms,2)*pow(mw,2) - 4*pow(mw,2)*pow(mu(t),2))*
           (pow(mew,2)*pow(mw,2) + mew*(pow(ms,2) - 2*pow(mw,2))*mu(t) + pow(mu(t),2)*(pow(mw,2) + pow(mu(t),2)))) + 
        (2*(pow(ms,6) - 6*pow(ms,4)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2)*(pow(mw,2) - pow(mu(t),2)) + 
             2*pow(ms,2)*(2*pow(mw,4) - 5*pow(mw,2)*pow(mu(t),2) - pow(mu(t),4)))*
           atan((-pow(ms,2) + 2*pow(mw,2)*(1 - mew/mu(t)))/sqrt(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))/
         pow(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5) - 2*log(mew/mu(t)) + 
        log((pow(mew,2)*pow(mw,2) + mew*pow(ms,2)*mu(t) - 2*mew*pow(mw,2)*mu(t) + pow(mw,2)*pow(mu(t),2) + pow(mu(t),4))/pow(mu(t),2)))/pow(pow(mw,2) + pow(mu(t),2),2)));
	return I;
}

double Ivwlwlvt(double t, double mv){
	double I;
	I = ((pow(mu(t),4)*(2*(pow(mv,6) - 2*pow(mv,4)*pow(mw,2) + 2*pow(mv,2)*pow(mw,2)*pow(mu(t),2) + 4*pow(mw,2)*pow(mu(t),4))*
          atan(pow(mv,2)/sqrt(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))) - 
         2*(pow(mv,6) - 2*pow(mv,4)*pow(mw,2) + 2*pow(mv,2)*pow(mw,2)*pow(mu(t),2) + 4*pow(mw,2)*pow(mu(t),4))*
          atan((pow(mv,2) - 2*pow(mw,2))/sqrt(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))) + 
         sqrt(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))*
          (-2*pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2) + 
            (pow(mv,4) - 4*pow(mv,2)*pow(mw,2) - 4*pow(mw,2)*pow(mu(t),2))*(log(pow(mv,2) + pow(mu(t),2)) - log(pow(mw,2) + pow(mu(t),2))))))/
     (pow(pow(mv,2) + pow(mu(t),2),2)*pow(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5)));
	return I;
}

double Ivvtff0(double t, double mv){
	double I;
	I = -((pow(mu(t),2)*(-((pow(mv,3) - 2*mv*pow(mu(t),2))/(pow(mv,2) - 4*pow(mu(t),2))) + 
       (8*pow(mu(t),4)*atan(mv/sqrt(-pow(mv,2) + 4*pow(mu(t),2))))/pow(-pow(mv,2) + 4*pow(mu(t),2),1.5)))/pow(mv,3));
	return I;
}

double Ivzttt(double t){
	double I;
	I = -((pow(mu(t),4)*(mz*(-2*pow(mt,2) + pow(mz,2) - 2*pow(mu(t),2))*sqrt(4*pow(mt,2) - pow(mz,2) + 4*pow(mu(t),2)) + 
       8*pow(pow(mt,2) + pow(mu(t),2),2)*atan(mz/sqrt(4*pow(mt,2) - pow(mz,2) + 4*pow(mu(t),2)))))/
   (pow(mz,3)*(pow(mt,2) + pow(mu(t),2))*pow(4*pow(mt,2) - pow(mz,2) + 4*pow(mu(t),2),1.5)));
	return I;
}

double Ivzthzl(double t){
	double I;
	I = -((2*pow(mu(t),4)*(sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2)) - 
       (pow(mh,2) + 2*pow(mu(t),2))*(atan(pow(mh,2)/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))) - 
          atan((pow(mh,2) - 2*pow(mz,2))/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))))))/
   pow(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2),1.5));
	return I;
}

double Ivztwlwl(double t){
	double I;
	I = (pow(mu(t),4)*(-(mz*sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2))) + (4*pow(mw,2) - 2*pow(mz,2) + 4*pow(mu(t),2))*
        atan(mz/sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2)))))/(pow(mz,3)*pow(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2),1.5));
	return I;
}

double Ivztwtwt(double t){
	double I;
	I = -((pow(mu(t),4)*((2*(-2*pow(mw,2) - 2*pow(mu(t),2))*pow(pow(mw,2) - pow(mz,2) + pow(mu(t),2),2))/(pow(mz,4)*(-4*pow(mw,2) + pow(mz,2) - 4*pow(mu(t),2))) - 
        (2*(2*pow(mw,6) + pow(mz,6) - 6*pow(mz,4)*pow(mu(t),2) + 3*pow(mz,2)*pow(mu(t),4) + 2*pow(mu(t),6) + 3*pow(mw,4)*(pow(mz,2) + 2*pow(mu(t),2)) + 
             6*pow(mw,2)*(-pow(mz,4) + pow(mz,2)*pow(mu(t),2) + pow(mu(t),4)))*atan(mz/sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2))))/
         (pow(mz,3)*pow(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2),1.5)) + 
        ((pow(mw,4) - pow(mz,4) + 2*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4))*log(pow(mw,2) + pow(mu(t),2)))/pow(mz,4)))/pow(pow(mw,2) + pow(mu(t),2),2) - 
   (pow(mu(t),4)*((2*(-2*pow(mw,2) + pow(mz,2)*(1 - mew/mu(t)) - 2*pow(mu(t),2))*(pow(mw,2) + pow(mu(t),2))*pow(pow(mw,2) - pow(mz,2) + pow(mu(t),2),2))/
         (pow(mz,4)*(-4*pow(mw,2) + pow(mz,2) - 4*pow(mu(t),2))*(pow(mw,2) - (mew*pow(mz,2)*(1 - mew/mu(t)))/mu(t) + pow(mu(t),2))) + 
        (2*(2*pow(mw,6) + pow(mz,6) - 6*pow(mz,4)*pow(mu(t),2) + 3*pow(mz,2)*pow(mu(t),4) + 2*pow(mu(t),6) + 3*pow(mw,4)*(pow(mz,2) + 2*pow(mu(t),2)) + 
             6*pow(mw,2)*(-pow(mz,4) + pow(mz,2)*pow(mu(t),2) + pow(mu(t),4)))*atan((mz*(-1 + 2*(1 - mew/mu(t))))/sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2))))/
         (pow(mz,3)*pow(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2),1.5)) + 2*log(mew/mu(t)) + 
        ((pow(mw,4) - pow(mz,4) + 2*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4))*log(pow(mw,2) - (mew*pow(mz,2)*(1 - mew/mu(t)))/mu(t) + pow(mu(t),2)))/pow(mz,4)))/
    pow(pow(mw,2) + pow(mu(t),2),2));
	return I;
}

double Ivwptb(double t){
	double I;
	I = -((pow(mu(t),4)*((2*pow(mw,4)*(pow(mt,2) - pow(mw,2) + 2*pow(mu(t),2)))/((pow(mt,2) + pow(mu(t),2))*(pow(pow(mt,2) - pow(mw,2),2) - 4*pow(mw,2)*pow(mu(t),2))) - 
       (8*pow(mw,4)*pow(mu(t),2)*(atan(((mt - mw)*(mt + mw))/sqrt(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mw,2)*pow(mu(t),2))) - 
            atan((pow(mt,2) + pow(mw,2))/sqrt(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mw,2)*pow(mu(t),2)))))/
        pow(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mw,2)*pow(mu(t),2),1.5)))/(2.*pow(mw,4)));
    if(t<tt){
    	return 0;
	}
	else{
		return I;
	}
}

double Ivwmtb(double t){
	double I;
	I = (2*pow(mw,4)*pow(mu(t),2)*(pow(mt,2) - pow(mw,2) + 2*pow(mu(t),2))*sqrt(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mw,2)*pow(mu(t),2)) + 
     8*pow(mw,4)*pow(mu(t),4)*(pow(mt,2) + pow(mu(t),2))*atan(((mt - mw)*(mt + mw))/sqrt(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mw,2)*pow(mu(t),2))) - 
     8*pow(mw,4)*pow(mu(t),4)*(pow(mt,2) + pow(mu(t),2))*atan((pow(mt,2) + pow(mw,2))/sqrt(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mw,2)*pow(mu(t),2))))/
   (2.*pow(mw,4)*pow(-pow(pow(mt,2) - pow(mw,2),2) + 4*pow(mw,2)*pow(mu(t),2),1.5));
	if(t<tt){
    	return 0;
	}
	else{
		return I;
	}
}

double Ivwtwls(double t, double ms){
	double I;
	I = -((2*pow(mu(t),4)*(sqrt(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2)) - 
       (pow(ms,2) + 2*pow(mu(t),2))*(atan(pow(ms,2)/sqrt(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))) - 
          atan((pow(ms,2) - 2*pow(mw,2))/sqrt(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))))/
   pow(-pow(ms,4) + 4*pow(ms,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5));
	return I;
}

double Ivwtwtvt(double t, double mv){
	double I;
	I = -(pow(mu(t),4)*((-4*(pow(mv,6)*(pow(mw,2) + pow(mu(t),2)) + 2*pow(mw,2)*pow(mu(t),2)*(pow(mw,4) + pow(mw,2)*pow(mu(t),2) - pow(mu(t),4)) + 
          pow(mv,4)*(-2*pow(mw,4) - 2*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) + pow(mv,2)*(pow(mw,6) - pow(mw,4)*pow(mu(t),2) - 5*pow(mw,2)*pow(mu(t),4))))/
      (pow(mw,4)*(pow(mv,2) + pow(mu(t),2))*(pow(mw,2) + pow(mu(t),2))*(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))) - 
     (2*(pow(mv,8) + pow(mv,6)*(-4*pow(mw,2) + 2*pow(mu(t),2)) + pow(mv,4)*(5*pow(mw,4) - 10*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) - 
          2*pow(mw,2)*pow(mu(t),2)*(pow(mw,4) - 3*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) - 2*pow(mv,2)*(pow(mw,6) - 5*pow(mw,4)*pow(mu(t),2) + 4*pow(mw,2)*pow(mu(t),4))))/
      (pow(mw,4)*pow(pow(mv,2) + pow(mu(t),2),2)*(pow(mv,4) - 4*pow(mv,2)*pow(mw,2) - 4*pow(mw,2)*pow(mu(t),2))) + 
     (2*mu(t)*(mew*(pow(mv,6)*(pow(mw,2) + pow(mu(t),2)) - 4*pow(mv,4)*pow(mw,2)*(pow(mw,2) + pow(mu(t),2)) + 
             2*pow(mw,4)*(-pow(mw,4) + pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) + pow(mv,2)*(5*pow(mw,6) + pow(mw,4)*pow(mu(t),2) - 3*pow(mw,2)*pow(mu(t),4))) + 
          mu(t)*(-2*pow(mv,2)*pow(mw,2)*(2*pow(mw,4) + pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) + pow(mv,4)*(2*pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) + 
             2*pow(mw,2)*(pow(mw,6) - pow(mu(t),6)))))/
      (pow(mw,4)*(pow(mw,2) + pow(mu(t),2))*(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))*
        (pow(mew,2)*pow(mw,2) + mew*(pow(mv,2) - 2*pow(mw,2))*mu(t) + pow(mu(t),2)*(pow(mw,2) + pow(mu(t),2)))) - 
     (2*(pow(mv,10) + 8*pow(mw,4)*pow(mu(t),4)*(pow(mw,2) + pow(mu(t),2)) + pow(mv,8)*(-6*pow(mw,2) + 2*pow(mu(t),2)) + 
          pow(mv,6)*(3*pow(mw,4) - 18*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) + 2*pow(mv,4)*(pow(mw,6) + 6*pow(mw,4)*pow(mu(t),2) - 9*pow(mw,2)*pow(mu(t),4)) + 
          2*pow(mv,2)*(5*pow(mw,6)*pow(mu(t),2) + 9*pow(mw,4)*pow(mu(t),4) - 3*pow(mw,2)*pow(mu(t),6)))*
        atan(pow(mv,2)/sqrt(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))/
      (pow(mw,4)*pow(pow(mv,2) + pow(mu(t),2),2)*pow(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5)) + 
     (2*(pow(mv,6)*(2*pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) - 6*pow(mv,4)*pow(mw,2)*(2*pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) + 
          4*pow(mw,4)*(pow(mw,6) + 4*pow(mw,4)*pow(mu(t),2) + 2*pow(mw,2)*pow(mu(t),4) + pow(mu(t),6)) + 
          6*pow(mv,2)*(pow(mw,8) - 2*pow(mw,6)*pow(mu(t),2) - 2*pow(mw,4)*pow(mu(t),4) - pow(mw,2)*pow(mu(t),6)))*
        atan(pow(mv,2)/sqrt(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))/
      (pow(mw,4)*pow(pow(mw,2) + pow(mu(t),2),2)*pow(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5)) + 
     (2*(pow(mv,10) + 8*pow(mw,4)*pow(mu(t),4)*(pow(mw,2) + pow(mu(t),2)) + pow(mv,8)*(-6*pow(mw,2) + 2*pow(mu(t),2)) + 
          pow(mv,6)*(3*pow(mw,4) - 18*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) + 2*pow(mv,4)*(pow(mw,6) + 6*pow(mw,4)*pow(mu(t),2) - 9*pow(mw,2)*pow(mu(t),4)) + 
          2*pow(mv,2)*(5*pow(mw,6)*pow(mu(t),2) + 9*pow(mw,4)*pow(mu(t),4) - 3*pow(mw,2)*pow(mu(t),6)))*
        atan((pow(mv,2) - 2*pow(mw,2))/sqrt(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))/
      (pow(mw,4)*pow(pow(mv,2) + pow(mu(t),2),2)*pow(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5)) + 
     (2*(pow(mv,6)*(2*pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) - 6*pow(mv,4)*pow(mw,2)*(2*pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) + 
          4*pow(mw,4)*(pow(mw,6) + 4*pow(mw,4)*pow(mu(t),2) + 2*pow(mw,2)*pow(mu(t),4) + pow(mu(t),6)) + 
          6*pow(mv,2)*(pow(mw,8) - 2*pow(mw,6)*pow(mu(t),2) - 2*pow(mw,4)*pow(mu(t),4) - pow(mw,2)*pow(mu(t),6)))*
        atan((-pow(mv,2) + 2*pow(mw,2)*(1 - mew/mu(t)))/sqrt(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))/
      (pow(mw,4)*pow(pow(mw,2) + pow(mu(t),2),2)*pow(-pow(mv,4) + 4*pow(mv,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5)) - 
     (2*log(mew/mu(t)))/pow(pow(mw,2) + pow(mu(t),2),2) - ((pow(mv,4) - pow(mw,4) + 2*pow(mv,2)*pow(mu(t),2) + pow(mu(t),4))*log(pow(mv,2) + pow(mu(t),2)))/
      (pow(mw,4)*pow(pow(mv,2) + pow(mu(t),2),2)) + ((2*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4))*log(pow(mv,2) + pow(mu(t),2)))/
      (pow(mw,4)*pow(pow(mw,2) + pow(mu(t),2),2)) + ((pow(mv,4) - pow(mw,4) + 2*pow(mv,2)*pow(mu(t),2) + pow(mu(t),4))*log(pow(mw,2) + pow(mu(t),2)))/
      (pow(mw,4)*pow(pow(mv,2) + pow(mu(t),2),2)) - ((2*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4))*
        log((pow(mew,2)*pow(mw,2) + mew*pow(mv,2)*mu(t) - 2*mew*pow(mw,2)*mu(t) + pow(mw,2)*pow(mu(t),2) + pow(mu(t),4))/pow(mu(t),2)))/
      (pow(mw,4)*pow(pow(mw,2) + pow(mu(t),2),2))));
	return I;
}

double Ivphwtwt(double t){
	double I;
	I = 2*((pow(mew - mu(t),2)*(3*pow(mew,2) - 2*mew*mu(t) + 11*pow(mu(t),2)))/(12.*pow(mu(t),4)) + log(mew/mu(t)));
	return I;
}
