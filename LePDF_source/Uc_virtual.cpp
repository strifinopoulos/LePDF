/*Ultra-collinear contributions to the coefficients P_A^v. Same notation as the massless ones*/

double Ivffzluc(double t){
	double I;
	I =-(pow(mu(t),6)*(pow(mz,2) - (2*pow(mz,2) + pow(mu(t),2))*log((pow(mz,2) + pow(mu(t),2))/pow(mu(t),2))))/(pow(mz,4)*pow(pow(mz,2) + pow(mu(t),2),2));
	return I;
}

double Ivffwluc(double t){
	double I;
	I =-((pow(mu(t),3)*(mew*pow(mw,2) + pow(mu(t),3) - (2*pow(mu(t),5))/(pow(mw,2) + pow(mu(t),2)) + 
       pow(mu(t),8)/((pow(mw,2) + pow(mu(t),2))*(-(mew*pow(mw,2)) + pow(mw,2)*mu(t) + pow(mu(t),3))) + 
       (-(pow(mw,6)*mu(t)*log(mew/mu(t))) + pow(mu(t),3)*(pow(mw,2)*(2*pow(mw,2) + pow(mu(t),2))*log(pow(mu(t),2)) - 2*pow(pow(mw,2) + pow(mu(t),2),2)*log(pow(mw,2) + pow(mu(t),2)) + 
             pow(mu(t),2)*(3*pow(mw,2) + 2*pow(mu(t),2))*log(pow(mw,2)*(1 - mew/mu(t)) + pow(mu(t),2))))/pow(pow(mw,2) + pow(mu(t),2),2)))/pow(mw,6));
	return I;
}

double Ivtltlzluc(double t, double alphaz, double alphay){
	double I;
	I =-((pow(mu(t),4)*((-4*(pow(alphay,2)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,4) - 3*pow(mt,2)*(pow(mz,2) + pow(mu(t),2))) + 
            4*alphay*alphaz*pow(mt,2)*Qzul(t,Qu)*(pow(mz,2) + pow(mu(t),2))*(-pow(mz,4) + pow(mt,2)*(3*pow(mz,2) + 2*pow(mu(t),2))) + 
            4*pow(alphaz,2)*pow(mt,4)*pow(Qzul(t,Qu),2)*(pow(mz,2)*(pow(mz,2) + pow(mu(t),2)) - pow(mt,2)*(3*pow(mz,2) + 4*pow(mu(t),2)))))/
        ((pow(mz,2) + pow(mu(t),2))*(pow(mz,4) - 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))) + 
       (-8*alphay*alphaz*pow(mt,2)*Qzul(t,Qu)*pow(pow(mz,2) + pow(mu(t),2),2)*(2*pow(mt,4) + pow(mz,4) - 2*pow(mt,2)*(2*pow(mz,2) + pow(mu(t),2))) + 
          2*pow(alphay,2)*pow(pow(mz,2) + pow(mu(t),2),2)*(pow(mz,6) + 2*pow(mt,4)*(pow(mz,2) + pow(mu(t),2)) - 
             pow(mt,2)*(4*pow(mz,4) + 3*pow(mz,2)*pow(mu(t),2))) + 
          8*pow(alphaz,2)*pow(mt,4)*pow(Qzul(t,Qu),2)*(2*pow(mt,4)*(pow(mz,2) + pow(mu(t),2)) + pow(mz,2)*pow(pow(mz,2) + pow(mu(t),2),2) - 
             pow(mt,2)*(4*pow(mz,4) + 9*pow(mz,2)*pow(mu(t),2) + 6*pow(mu(t),4))))/(pow(pow(mz,2) + pow(mu(t),2),2)*(pow(mz,4) - 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2))))
         + (8*pow(alphaz,2)*pow(mt,4)*(2*pow(mt,2) - pow(mz,2))*pow(Qzul(t,Qu),2)*(pow(mt,2) + pow(mu(t),2)) - 
          8*alphay*alphaz*pow(mt,2)*Qzul(t,Qu)*(pow(mt,2) + pow(mu(t),2))*(-pow(mz,4) + 2*pow(mt,2)*(pow(mz,2) + pow(mu(t),2))) + 
          pow(alphay,2)*(-2*pow(mz,6)*pow(mu(t),2) + 4*pow(mt,4)*pow(pow(mz,2) + pow(mu(t),2),2) + 
             pow(mt,2)*(-2*pow(mz,6) + 4*pow(mz,4)*pow(mu(t),2) + 6*pow(mz,2)*pow(mu(t),4))))/
        ((pow(mt,2) + pow(mu(t),2))*(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))) + 
       (4*(8*pow(alphaz,2)*pow(mt,6)*pow(Qzul(t,Qu),2)*(pow(mt,2) + pow(mu(t),2)) + 
            2*alphay*alphaz*pow(mt,2)*Qzul(t,Qu)*(pow(mz,6) + 4*pow(mt,4)*(pow(mz,2) + pow(mu(t),2)) - 6*pow(mt,2)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2))) - 
            pow(alphay,2)*(pow(mz,8) - 6*pow(mt,2)*pow(mz,4)*(pow(mz,2) + pow(mu(t),2)) + 6*pow(mt,4)*pow(pow(mz,2) + pow(mu(t),2),2)))*
          atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5) + 
       (2*(-4*pow(alphaz,2)*pow(mt,6)*pow(Qzul(t,Qu),2)*(pow(mz,2) + 2*pow(mu(t),2))*
             (-pow(mz,4) + 2*pow(mz,2)*pow(mu(t),2) + 2*pow(mu(t),4) + 6*pow(mt,2)*(pow(mz,2) + pow(mu(t),2))) - 
            4*alphay*alphaz*pow(mt,2)*Qzul(t,Qu)*pow(pow(mz,2) + pow(mu(t),2),2)*
             (pow(mz,6) - 6*pow(mt,2)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2)) + pow(mt,4)*(6*pow(mz,2) + 8*pow(mu(t),2))) + 
            pow(alphay,2)*pow(pow(mz,2) + pow(mu(t),2),2)*(2*pow(mz,8) - pow(mt,2)*(13*pow(mz,6) + 12*pow(mz,4)*pow(mu(t),2)) + 
               6*pow(mt,4)*(3*pow(mz,4) + 5*pow(mz,2)*pow(mu(t),2) + 2*pow(mu(t),4))))*atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/
        (pow(pow(mz,2) + pow(mu(t),2),2)*pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5)) + 
       (4*(8*pow(alphaz,2)*pow(mt,6)*pow(Qzul(t,Qu),2)*(pow(mt,2) + pow(mu(t),2)) + 
            2*alphay*alphaz*pow(mt,2)*Qzul(t,Qu)*(pow(mz,6) + 4*pow(mt,4)*(pow(mz,2) + pow(mu(t),2)) - 6*pow(mt,2)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2))) - 
            pow(alphay,2)*(pow(mz,8) - 6*pow(mt,2)*pow(mz,4)*(pow(mz,2) + pow(mu(t),2)) + 6*pow(mt,4)*pow(pow(mz,2) + pow(mu(t),2),2)))*
          atan((2*pow(mt,2) - pow(mz,2))/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5)
         - (2*(-4*pow(alphaz,2)*pow(mt,6)*pow(Qzul(t,Qu),2)*(pow(mz,2) + 2*pow(mu(t),2))*
             (-pow(mz,4) + 2*pow(mz,2)*pow(mu(t),2) + 2*pow(mu(t),4) + 6*pow(mt,2)*(pow(mz,2) + pow(mu(t),2))) - 
            4*alphay*alphaz*pow(mt,2)*Qzul(t,Qu)*pow(pow(mz,2) + pow(mu(t),2),2)*
             (pow(mz,6) - 6*pow(mt,2)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2)) + pow(mt,4)*(6*pow(mz,2) + 8*pow(mu(t),2))) + 
            pow(alphay,2)*pow(pow(mz,2) + pow(mu(t),2),2)*(2*pow(mz,8) - pow(mt,2)*(13*pow(mz,6) + 12*pow(mz,4)*pow(mu(t),2)) + 
               6*pow(mt,4)*(3*pow(mz,4) + 5*pow(mz,2)*pow(mu(t),2) + 2*pow(mu(t),4))))*
          atan((-2*pow(mt,2) + pow(mz,2))/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/
        (pow(pow(mz,2) + pow(mu(t),2),2)*pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5)) + 
       2*alphay*(alphay*pow(mz,2) - 2*alphaz*pow(mt,2)*Qzul(t,Qu))*log(pow(mt,2) + pow(mu(t),2)) + 
       (pow(alphay,2)*(pow(mt,2) - 2*pow(mz,2)) + 4*alphay*alphaz*pow(mt,2)*Qzul(t,Qu) - (4*pow(alphaz,2)*pow(mt,6)*pow(Qzul(t,Qu),2))/pow(pow(mz,2) + pow(mu(t),2),2))*
        log(pow(mt,2) + pow(mu(t),2)) - 2*alphay*(alphay*pow(mz,2) - 2*alphaz*pow(mt,2)*Qzul(t,Qu))*log(pow(mz,2) + pow(mu(t),2)) + 
       (-(pow(alphay,2)*(pow(mt,2) - 2*pow(mz,2))) - 4*alphay*alphaz*pow(mt,2)*Qzul(t,Qu) + (4*pow(alphaz,2)*pow(mt,6)*pow(Qzul(t,Qu),2))/pow(pow(mz,2) + pow(mu(t),2),2))*
        log(pow(mz,2) + pow(mu(t),2))))/(8.*pow(mt,6)));
	return I;
}

double Ivtlblwluc(double t, double alpha2, double alphay){
	double I;
	I =-(pow(mu(t),4)*(-(pow(alphay,2)/pow(mt,4)) + (pow(alphay,2)*mew)/(pow(mt,4)*mu(t)) - 
     (-2*alpha2*alphay*pow(mt,2)*(pow(mw,2) + pow(mu(t),2))*(pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + pow(mu(t),2))) + 
        pow(alpha2,2)*(pow(mt,8) + pow(mt,4)*pow(mw,2)*(pow(mw,2) + pow(mu(t),2)) - pow(mt,6)*(2*pow(mw,2) + 3*pow(mu(t),2))) + 
        pow(alphay,2)*(pow(mw,2) + pow(mu(t),2))*(pow(mw,6) + pow(mt,4)*(pow(mw,2) + pow(mu(t),2)) - pow(mt,2)*(2*pow(mw,4) + 3*pow(mw,2)*pow(mu(t),2))))/
      (pow(mt,6)*(pow(mw,2) + pow(mu(t),2))*(pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))) + 
     (mu(t)*(2*alpha2*alphay*pow(mt,2)*(pow(mw,2) + pow(mu(t),2))*(mew*(pow(mt,2) - pow(mw,2))*(pow(mt,4) + pow(mw,4) - pow(mt,2)*(2*pow(mw,2) + 3*pow(mu(t),2))) + 
             mu(t)*(-pow(mt,6) + pow(mw,4)*(pow(mw,2) + pow(mu(t),2)) + pow(mt,4)*(3*pow(mw,2) + 4*pow(mu(t),2)) - 
                pow(mt,2)*(3*pow(mw,4) + 5*pow(mw,2)*pow(mu(t),2) + 2*pow(mu(t),4)))) - 
          pow(alphay,2)*(pow(mw,2) + pow(mu(t),2))*(-((pow(mt,2) - pow(mw,2))*mu(t)*(pow(mw,2) + pow(mu(t),2))*
                (pow(mt,4) + pow(mw,4) - pow(mt,2)*(2*pow(mw,2) + 3*pow(mu(t),2)))) + 
             mew*(-pow(mw,8) + pow(mt,6)*(pow(mw,2) + pow(mu(t),2)) + pow(mt,2)*(3*pow(mw,6) + 4*pow(mw,4)*pow(mu(t),2)) - 
                pow(mt,4)*(3*pow(mw,4) + 5*pow(mw,2)*pow(mu(t),2) + 2*pow(mu(t),4)))) + 
          pow(alpha2,2)*pow(mt,4)*(mew*(-pow(mt,6) + pow(mw,4)*(pow(mw,2) + pow(mu(t),2)) + pow(mt,4)*(3*pow(mw,2) + 4*pow(mu(t),2)) - 
                pow(mt,2)*(3*pow(mw,4) + 5*pow(mw,2)*pow(mu(t),2) + 2*pow(mu(t),4))) + 
             mu(t)*(pow(mt,6) - pow(mw,2)*pow(pow(mw,2) + pow(mu(t),2),2) - pow(mt,4)*(3*pow(mw,2) + 5*pow(mu(t),2)) + 
                pow(mt,2)*(3*pow(mw,4) + 7*pow(mw,2)*pow(mu(t),2) + 5*pow(mu(t),4))))))/
      (pow(mt,6)*(pow(mw,2) + pow(mu(t),2))*(-(pow(mew,2)*pow(mt,2)) + mew*(pow(mt,2) + pow(mw,2))*mu(t) - pow(mu(t),2)*(pow(mw,2) + pow(mu(t),2)))*
        (pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))) + 
     ((-2*alpha2*alphay*pow(mt,2)*(pow(mt,2) - pow(mw,2))*pow(pow(mw,2) + pow(mu(t),2),2)*(pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 3*pow(mu(t),2))) + 
          pow(alpha2,2)*pow(mt,6)*(pow(mt,6) - pow(mw,6) + 6*pow(mw,2)*pow(mu(t),4) + 4*pow(mu(t),6) - 3*pow(mt,4)*(pow(mw,2) + 2*pow(mu(t),2)) + 
             3*pow(mt,2)*(pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2) + 2*pow(mu(t),4))) - 
          pow(alphay,2)*pow(pow(mw,2) + pow(mu(t),2),2)*(pow(mt,8) + 2*pow(mw,8) - pow(mt,6)*(5*pow(mw,2) + 6*pow(mu(t),2)) - 
             pow(mt,2)*(7*pow(mw,6) + 12*pow(mw,4)*pow(mu(t),2)) + 3*pow(mt,4)*(3*pow(mw,4) + 6*pow(mw,2)*pow(mu(t),2) + 4*pow(mu(t),4))))*
        atan((-pow(mt,2) + pow(mw,2))/sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/
      (pow(mt,6)*pow(pow(mw,2) + pow(mu(t),2),2)*pow(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)),1.5)) - 
     ((-2*alpha2*alphay*pow(mt,2)*(pow(mt,2) - pow(mw,2))*pow(pow(mw,2) + pow(mu(t),2),2)*(pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 3*pow(mu(t),2))) + 
          pow(alpha2,2)*pow(mt,6)*(pow(mt,6) - pow(mw,6) + 6*pow(mw,2)*pow(mu(t),4) + 4*pow(mu(t),6) - 3*pow(mt,4)*(pow(mw,2) + 2*pow(mu(t),2)) + 
             3*pow(mt,2)*(pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2) + 2*pow(mu(t),4))) - 
          pow(alphay,2)*pow(pow(mw,2) + pow(mu(t),2),2)*(pow(mt,8) + 2*pow(mw,8) - pow(mt,6)*(5*pow(mw,2) + 6*pow(mu(t),2)) - 
             pow(mt,2)*(7*pow(mw,6) + 12*pow(mw,4)*pow(mu(t),2)) + 3*pow(mt,4)*(3*pow(mw,4) + 6*pow(mw,2)*pow(mu(t),2) + 4*pow(mu(t),4))))*
        atan((-2*mew*pow(mt,2) + (pow(mt,2) + pow(mw,2))*mu(t))/(mu(t)*sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2))))))/
      (pow(mt,6)*pow(pow(mw,2) + pow(mu(t),2),2)*pow(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)),1.5)) - 
     (pow(alpha2,2)*log(mew/mu(t)))/pow(pow(mw,2) + pow(mu(t),2),2) + (((2*alpha2*alphay)/pow(mt,4) + (pow(alphay,2)*(pow(mt,2) - 2*pow(mw,2)))/pow(mt,6) - 
          pow(alpha2,2)/pow(pow(mw,2) + pow(mu(t),2),2))*log(pow(mu(t),2)))/2. + 
     ((pow(alpha2,2)*(-pow(mt,6) + pow(mt,4)*pow(mw,2)) + 2*alpha2*alphay*pow(mt,2)*(-pow(mw,4) + pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2))) - 
           pow(alphay,2)*(-pow(mw,6) + pow(mt,4)*pow(mu(t),2) + pow(mt,2)*(pow(mw,4) + 3*pow(mw,2)*pow(mu(t),2))))/
         (pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2))) + 
        (2*(2*pow(alpha2,2)*pow(mt,6)*pow(mu(t),2) + alpha2*alphay*pow(mt,2)*
              (-pow(mt,6) + pow(mw,6) + pow(mt,4)*(3*pow(mw,2) + 2*pow(mu(t),2)) - 3*pow(mt,2)*(pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2))) + 
             pow(alphay,2)*(-pow(mw,8) + pow(mt,6)*(pow(mw,2) + 2*pow(mu(t),2)) + 3*pow(mt,2)*(pow(mw,6) + 2*pow(mw,4)*pow(mu(t),2)) - 
                3*pow(mt,4)*(pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2) + 2*pow(mu(t),4))))*
           atan((pow(mt,2) - pow(mw,2))/sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/
         pow(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)),1.5) + alphay*(-(alpha2*pow(mt,2)) + alphay*pow(mw,2))*log(pow(mu(t),2)))/pow(mt,6) - 
     (-(pow(alphay,2)*pow(mt,2)) + (-2*alpha2*alphay*pow(mt,2)*(pow(mw,2) + pow(mu(t),2))*(pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + pow(mu(t),2))) + 
           pow(alpha2,2)*(pow(mt,8) + pow(mt,4)*pow(mw,2)*(pow(mw,2) + pow(mu(t),2)) - pow(mt,6)*(2*pow(mw,2) + 3*pow(mu(t),2))) + 
           pow(alphay,2)*(pow(mw,2) + pow(mu(t),2))*(pow(mw,6) + pow(mt,4)*(pow(mw,2) + pow(mu(t),2)) - pow(mt,2)*(2*pow(mw,4) + 3*pow(mw,2)*pow(mu(t),2))))/
         ((pow(mw,2) + pow(mu(t),2))*(pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))) + 
        (2*(2*pow(alpha2,2)*pow(mt,6)*pow(mu(t),2) + alpha2*alphay*pow(mt,2)*
              (-pow(mt,6) + pow(mw,6) + pow(mt,4)*(3*pow(mw,2) + 2*pow(mu(t),2)) - 3*pow(mt,2)*(pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2))) + 
             pow(alphay,2)*(-pow(mw,8) + pow(mt,6)*(pow(mw,2) + 2*pow(mu(t),2)) + 3*pow(mt,2)*(pow(mw,6) + 2*pow(mw,4)*pow(mu(t),2)) - 
                3*pow(mt,4)*(pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2) + 2*pow(mu(t),4))))*
           atan((-pow(mt,2) - pow(mw,2))/sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/
         pow(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)),1.5) + alphay*(-(alpha2*pow(mt,2)) + alphay*pow(mw,2))*log(pow(mw,2) + pow(mu(t),2)))/
      pow(mt,6) + ((pow(alpha2,2)*pow(mt,6) - 2*alpha2*alphay*pow(mt,2)*pow(pow(mw,2) + pow(mu(t),2),2) - 
          pow(alphay,2)*(pow(mt,2) - 2*pow(mw,2))*pow(pow(mw,2) + pow(mu(t),2),2))*log((mew*pow(mt,2)*(mew - mu(t)) - pow(mw,2)*(mew - mu(t))*mu(t) + pow(mu(t),4))/pow(mu(t),2)))/
      (2.*pow(mt,6)*pow(pow(mw,2) + pow(mu(t),2),2))));
	return I;
}

double Ivtthuc(double t){
	double I;
	I =-((pow(mu(t),4)*((2*(-8*pow(mt,6) - 6*pow(mt,4)*pow(mu(t),2) + pow(mh,2)*pow(mt,2)*(2*pow(mt,2) + pow(mu(t),2))))/
        ((pow(mt,2) + pow(mu(t),2))*(pow(mh,4) - 4*pow(mh,2)*pow(mt,2) - 4*pow(mt,2)*pow(mu(t),2))) - 
       (2*(pow(mh,6) - 6*pow(mh,4)*pow(mt,2) + 16*pow(mt,4)*pow(mu(t),2) + pow(mh,2)*(8*pow(mt,4) - 6*pow(mt,2)*pow(mu(t),2)))*
          (atan(pow(mh,2)/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mt,2) + 4*pow(mt,2)*pow(mu(t),2))) - 
            atan((pow(mh,2) - 2*pow(mt,2))/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mt,2) + 4*pow(mt,2)*pow(mu(t),2)))))/
        pow(-pow(mh,4) + 4*pow(mh,2)*pow(mt,2) + 4*pow(mt,2)*pow(mu(t),2),1.5) - log(pow(mh,2) + pow(mu(t),2)) + log(pow(mt,2) + pow(mu(t),2))))/(2.*pow(mt,4)));
	return I;
}

double Ivtltrztuc(double t){
	double I;
	I =-((pow(mu(t),4)*((18*pow(mt,8) + pow(mt,6)*pow(mz,2)*(-9 + 24*pow(sw(t),2) - 32*pow(sw(t),4)) - 
          16*pow(mt,4)*pow(sw(t),2)*(pow(mz,2)*pow(sw(t),2) + pow(mt,2)*(-3 + 2*pow(sw(t),2)))*pow(mu(t),2))/
        ((pow(mt,2) + pow(mu(t),2))*(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))) + 
       (2*pow(mt,2)*((-8*pow(mz,6)*pow(sw(t),4) + 48*pow(mt,2)*pow(mz,2)*pow(sw(t),4)*(pow(mz,2) + pow(mu(t),2)) + 
               pow(mt,4)*(pow(mz,2)*(9 - 48*pow(sw(t),2)) - 48*pow(sw(t),2)*pow(mu(t),2)))*atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))) + 
            (-8*pow(mz,6)*pow(sw(t),4) + 48*pow(mt,2)*pow(mz,2)*pow(sw(t),4)*(pow(mz,2) + pow(mu(t),2)) + 
               pow(mt,4)*(pow(mz,2)*(9 - 48*pow(sw(t),2)) - 48*pow(sw(t),2)*pow(mu(t),2)))*
             atan((2*pow(mt,2) - pow(mz,2))/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))) + 
            4*pow(sw(t),4)*pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5)*(log(pow(mt,2) + pow(mu(t),2)) - log(pow(mz,2) + pow(mu(t),2)))))/
        pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5)))/(36.*pow(mt,6)));
	return I;
}

double Ivttvtuc(double t){
	double I;
	I =-((pow(mu(t),4)*(-(pow(mt,2)/(pow(mt,2) + pow(mu(t),2))) - log(pow(mu(t),2)) + log(pow(mt,2) + pow(mu(t),2))))/(2.*pow(mt,4)));
	return I;
}

double Ivtrtrzluc(double t, double alphaz, double alphay){
	double I;
	I =-((pow(mu(t),4)*((36*pow(alphay,2)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2))*(-pow(mz,4) + 3*pow(mt,2)*(pow(mz,2) + pow(mu(t),2))) - 
          96*alphay*alphaz*pow(mt,2)*pow(sw(t),2)*(pow(mz,2) + pow(mu(t),2))*(-pow(mz,4) + pow(mt,2)*(3*pow(mz,2) + 2*pow(mu(t),2))) + 
          64*pow(alphaz,2)*pow(mt,4)*pow(sw(t),4)*(-(pow(mz,2)*(pow(mz,2) + pow(mu(t),2))) + pow(mt,2)*(3*pow(mz,2) + 4*pow(mu(t),2))))/
        ((pow(mz,2) + pow(mu(t),2))*(pow(mz,4) - 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))) + 
       (-48*alphay*alphaz*pow(mt,2)*pow(sw(t),2)*pow(pow(mz,2) + pow(mu(t),2),2)*(2*pow(mt,4) + pow(mz,4) - 2*pow(mt,2)*(2*pow(mz,2) + pow(mu(t),2))) + 
          18*pow(alphay,2)*pow(pow(mz,2) + pow(mu(t),2),2)*(pow(mz,6) + 2*pow(mt,4)*(pow(mz,2) + pow(mu(t),2)) - 
             pow(mt,2)*(4*pow(mz,4) + 3*pow(mz,2)*pow(mu(t),2))) + 32*pow(alphaz,2)*pow(mt,4)*pow(sw(t),4)*
           (2*pow(mt,4)*(pow(mz,2) + pow(mu(t),2)) + pow(mz,2)*pow(pow(mz,2) + pow(mu(t),2),2) - pow(mt,2)*(4*pow(mz,4) + 9*pow(mz,2)*pow(mu(t),2) + 6*pow(mu(t),4))))/
        (pow(pow(mz,2) + pow(mu(t),2),2)*(pow(mz,4) - 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))) + 
       (32*pow(alphaz,2)*pow(mt,4)*(2*pow(mt,2) - pow(mz,2))*pow(sw(t),4)*(pow(mt,2) + pow(mu(t),2)) - 
          48*alphay*alphaz*pow(mt,2)*pow(sw(t),2)*(pow(mt,2) + pow(mu(t),2))*(-pow(mz,4) + 2*pow(mt,2)*(pow(mz,2) + pow(mu(t),2))) + 
          18*pow(alphay,2)*(-(pow(mz,6)*pow(mu(t),2)) + 2*pow(mt,4)*pow(pow(mz,2) + pow(mu(t),2),2) + 
             pow(mt,2)*(-pow(mz,6) + 2*pow(mz,4)*pow(mu(t),2) + 3*pow(mz,2)*pow(mu(t),4))))/
        ((pow(mt,2) + pow(mu(t),2))*(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))) + 
       (4*(32*pow(alphaz,2)*pow(mt,6)*pow(sw(t),4)*(pow(mt,2) + pow(mu(t),2)) + 
            12*alphay*alphaz*pow(mt,2)*pow(sw(t),2)*(pow(mz,6) + 4*pow(mt,4)*(pow(mz,2) + pow(mu(t),2)) - 6*pow(mt,2)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2))) - 
            9*pow(alphay,2)*(pow(mz,8) - 6*pow(mt,2)*pow(mz,4)*(pow(mz,2) + pow(mu(t),2)) + 6*pow(mt,4)*pow(pow(mz,2) + pow(mu(t),2),2)))*
          atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5) + 
       (2*(-16*pow(alphaz,2)*pow(mt,6)*pow(sw(t),4)*(pow(mz,2) + 2*pow(mu(t),2))*
             (-pow(mz,4) + 2*pow(mz,2)*pow(mu(t),2) + 2*pow(mu(t),4) + 6*pow(mt,2)*(pow(mz,2) + pow(mu(t),2))) - 
            24*alphay*alphaz*pow(mt,2)*pow(sw(t),2)*pow(pow(mz,2) + pow(mu(t),2),2)*
             (pow(mz,6) - 6*pow(mt,2)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2)) + pow(mt,4)*(6*pow(mz,2) + 8*pow(mu(t),2))) + 
            9*pow(alphay,2)*pow(pow(mz,2) + pow(mu(t),2),2)*(2*pow(mz,8) - pow(mt,2)*(13*pow(mz,6) + 12*pow(mz,4)*pow(mu(t),2)) + 
               6*pow(mt,4)*(3*pow(mz,4) + 5*pow(mz,2)*pow(mu(t),2) + 2*pow(mu(t),4))))*atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/
        (pow(pow(mz,2) + pow(mu(t),2),2)*pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5)) + 
       (4*(32*pow(alphaz,2)*pow(mt,6)*pow(sw(t),4)*(pow(mt,2) + pow(mu(t),2)) + 
            12*alphay*alphaz*pow(mt,2)*pow(sw(t),2)*(pow(mz,6) + 4*pow(mt,4)*(pow(mz,2) + pow(mu(t),2)) - 6*pow(mt,2)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2))) - 
            9*pow(alphay,2)*(pow(mz,8) - 6*pow(mt,2)*pow(mz,4)*(pow(mz,2) + pow(mu(t),2)) + 6*pow(mt,4)*pow(pow(mz,2) + pow(mu(t),2),2)))*
          atan((2*pow(mt,2) - pow(mz,2))/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5)\
        - (2*(-16*pow(alphaz,2)*pow(mt,6)*pow(sw(t),4)*(pow(mz,2) + 2*pow(mu(t),2))*
             (-pow(mz,4) + 2*pow(mz,2)*pow(mu(t),2) + 2*pow(mu(t),4) + 6*pow(mt,2)*(pow(mz,2) + pow(mu(t),2))) - 
            24*alphay*alphaz*pow(mt,2)*pow(sw(t),2)*pow(pow(mz,2) + pow(mu(t),2),2)*
             (pow(mz,6) - 6*pow(mt,2)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2)) + pow(mt,4)*(6*pow(mz,2) + 8*pow(mu(t),2))) + 
            9*pow(alphay,2)*pow(pow(mz,2) + pow(mu(t),2),2)*(2*pow(mz,8) - pow(mt,2)*(13*pow(mz,6) + 12*pow(mz,4)*pow(mu(t),2)) + 
               6*pow(mt,4)*(3*pow(mz,4) + 5*pow(mz,2)*pow(mu(t),2) + 2*pow(mu(t),4))))*
          atan((-2*pow(mt,2) + pow(mz,2))/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/
        (pow(pow(mz,2) + pow(mu(t),2),2)*pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5)) + 
       6*alphay*(3*alphay*pow(mz,2) - 4*alphaz*pow(mt,2)*pow(sw(t),2))*log(pow(mt,2) + pow(mu(t),2)) + 
       (9*pow(alphay,2)*(pow(mt,2) - 2*pow(mz,2)) + 24*alphay*alphaz*pow(mt,2)*pow(sw(t),2) - 
          (16*pow(alphaz,2)*pow(mt,6)*pow(sw(t),4))/pow(pow(mz,2) + pow(mu(t),2),2))*log(pow(mt,2) + pow(mu(t),2)) - 
       6*alphay*(3*alphay*pow(mz,2) - 4*alphaz*pow(mt,2)*pow(sw(t),2))*log(pow(mz,2) + pow(mu(t),2)) + 
       (-9*pow(alphay,2)*(pow(mt,2) - 2*pow(mz,2)) - 24*alphay*alphaz*pow(mt,2)*pow(sw(t),2) + 
          (16*pow(alphaz,2)*pow(mt,6)*pow(sw(t),4))/pow(pow(mz,2) + pow(mu(t),2),2))*log(pow(mz,2) + pow(mu(t),2))))/(72.*pow(mt,6)));
	return I;
}

double Ivtrtlztuc(double t){
	double I;
	I =-((pow(mu(t),4)*((2*(-2*pow(mz,6)*pow(3 - 4*pow(sw(t),2),2) + 6*pow(mt,2)*pow(mz,2)*(-3 + 4*pow(sw(t),2))*
             (pow(mz,2)*(-5 + 4*pow(sw(t),2)) + (-3 + 4*pow(sw(t),2))*pow(mu(t),2)) + pow(mt,4)*(pow(mz,2)*(-90 + 96*pow(sw(t),2)) + 24*(-3 + 4*pow(sw(t),2))*pow(mu(t),2))))/
        (pow(mz,4) - 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2))) + (2*(-18*pow(mt,6) + pow(mz,6)*pow(3 - 4*pow(sw(t),2),2) - 
            pow(mt,2)*pow(mz,2)*(-3 + 4*pow(sw(t),2))*(2*pow(mz,2)*(-9 + 8*pow(sw(t),2)) + 3*(-3 + 4*pow(sw(t),2))*pow(mu(t),2)) + 
            pow(mt,4)*(pow(mz,2)*(81 - 120*pow(sw(t),2) + 32*pow(sw(t),4)) + 2*(27 - 48*pow(sw(t),2) + 16*pow(sw(t),4))*pow(mu(t),2))))/
        (pow(mz,4) - 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2))) - (2*(pow(mz,6)*pow(3 - 4*pow(sw(t),2),2)*pow(mu(t),2) + 18*pow(mt,6)*(pow(mz,2) + pow(mu(t),2)) + 
            pow(mt,2)*pow(mz,2)*(-3 + 4*pow(sw(t),2))*(pow(mz,4)*(-3 + 4*pow(sw(t),2)) + 4*pow(mz,2)*(3 - 2*pow(sw(t),2))*pow(mu(t),2) + 3*(3 - 4*pow(sw(t),2))*pow(mu(t),4)) - 
            pow(mt,4)*(4*pow(mz,4)*(9 - 18*pow(sw(t),2) + 8*pow(sw(t),4)) + pow(mz,2)*(9 - 72*pow(sw(t),2) + 64*pow(sw(t),4))*pow(mu(t),2) + 2*(-9 + 16*pow(sw(t),4))*pow(mu(t),4))))
         /((pow(mt,2) + pow(mu(t),2))*(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))) + 
       (4*(-(pow(mz,8)*pow(3 - 4*pow(sw(t),2),2)) + 18*pow(mt,6)*(pow(mz,2) + pow(mu(t),2)) - 
            6*pow(mt,4)*(-3 + 4*pow(sw(t),2))*(pow(mz,2) + pow(mu(t),2))*(pow(mz,2)*(-6 + 4*pow(sw(t),2)) + (-3 + 4*pow(sw(t),2))*pow(mu(t),2)) + 
            3*pow(mt,2)*pow(mz,4)*(-3 + 4*pow(sw(t),2))*(pow(mz,2)*(-7 + 8*pow(sw(t),2)) + 2*(-3 + 4*pow(sw(t),2))*pow(mu(t),2)))*
          atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5) + 
       (2*(2*pow(mz,8)*pow(3 - 4*pow(sw(t),2),2) + 6*pow(mt,4)*(-3 + 4*pow(sw(t),2))*(pow(mz,2) + pow(mu(t),2))*
             (3*pow(mz,2)*(-5 + 4*pow(sw(t),2)) + 2*(-3 + 4*pow(sw(t),2))*pow(mu(t),2)) - 
            pow(mt,2)*pow(mz,4)*(-3 + 4*pow(sw(t),2))*(pow(mz,2)*(-45 + 52*pow(sw(t),2)) + 12*(-3 + 4*pow(sw(t),2))*pow(mu(t),2)) + 
            6*pow(mt,6)*(pow(mz,2)*(-15 + 16*pow(sw(t),2)) + 2*(-9 + 8*pow(sw(t),2))*pow(mu(t),2)))*
          atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5) + 
       (4*(-(pow(mz,8)*pow(3 - 4*pow(sw(t),2),2)) + 18*pow(mt,6)*(pow(mz,2) + pow(mu(t),2)) - 
            6*pow(mt,4)*(-3 + 4*pow(sw(t),2))*(pow(mz,2) + pow(mu(t),2))*(pow(mz,2)*(-6 + 4*pow(sw(t),2)) + (-3 + 4*pow(sw(t),2))*pow(mu(t),2)) + 
            3*pow(mt,2)*pow(mz,4)*(-3 + 4*pow(sw(t),2))*(pow(mz,2)*(-7 + 8*pow(sw(t),2)) + 2*(-3 + 4*pow(sw(t),2))*pow(mu(t),2)))*
          atan((2*pow(mt,2) - pow(mz,2))/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5)\
        - (2*(2*pow(mz,8)*pow(3 - 4*pow(sw(t),2),2) + 6*pow(mt,4)*(-3 + 4*pow(sw(t),2))*(pow(mz,2) + pow(mu(t),2))*
             (3*pow(mz,2)*(-5 + 4*pow(sw(t),2)) + 2*(-3 + 4*pow(sw(t),2))*pow(mu(t),2)) - 
            pow(mt,2)*pow(mz,4)*(-3 + 4*pow(sw(t),2))*(pow(mz,2)*(-45 + 52*pow(sw(t),2)) + 12*(-3 + 4*pow(sw(t),2))*pow(mu(t),2)) + 
            6*pow(mt,6)*(pow(mz,2)*(-15 + 16*pow(sw(t),2)) + 2*(-9 + 8*pow(sw(t),2))*pow(mu(t),2)))*
          atan((-2*pow(mt,2) + pow(mz,2))/sqrt(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)))))/pow(-pow(mz,4) + 4*pow(mt,2)*(pow(mz,2) + pow(mu(t),2)),1.5)
         + (-3 + 4*pow(sw(t),2))*(2*pow(mz,2)*(3 - 4*pow(sw(t),2)) + pow(mt,2)*(-9 + 4*pow(sw(t),2)))*log(pow(mt,2) + pow(mu(t),2)) + 
       2*(-3 + 4*pow(sw(t),2))*(3*pow(mt,2) + pow(mz,2)*(-3 + 4*pow(sw(t),2)))*log(pow(mt,2) + pow(mu(t),2)) - 
       (-3 + 4*pow(sw(t),2))*(2*pow(mz,2)*(3 - 4*pow(sw(t),2)) + pow(mt,2)*(-9 + 4*pow(sw(t),2)))*log(pow(mz,2) + pow(mu(t),2)) - 
       2*(-3 + 4*pow(sw(t),2))*(3*pow(mt,2) + pow(mz,2)*(-3 + 4*pow(sw(t),2)))*log(pow(mz,2) + pow(mu(t),2))))/(72.*pow(mt,6)));
	return I;
}

double Ivtrblwtuc(double t){
	double I;
	I =-((pow(mu(t),4)*((-2*(pow(mt,2) - pow(mw,2))*(pow(mt,4) + pow(mw,4) - pow(mt,2)*(2*pow(mw,2) + 3*pow(mu(t),2))))/
        (pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2))) - 
       (2*(-pow(mw,6) + pow(mt,4)*pow(mu(t),2) + pow(mt,2)*(pow(mw,4) + 3*pow(mw,2)*pow(mu(t),2))))/
        (pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2))) - 
       (4*(pow(mw,6) + pow(mt,4)*(pow(mw,2) + pow(mu(t),2)) - pow(mt,2)*(2*pow(mw,4) + 3*pow(mw,2)*pow(mu(t),2))))/
        (pow(mt,4) + pow(mw,4) - 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2))) - 
       (4*(-pow(mw,8) + pow(mt,6)*(pow(mw,2) + 2*pow(mu(t),2)) + 3*pow(mt,2)*(pow(mw,6) + 2*pow(mw,4)*pow(mu(t),2)) - 
            3*pow(mt,4)*(pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2) + 2*pow(mu(t),4)))*
          atan((-pow(mt,2) - pow(mw,2))/sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/
        pow(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)),1.5) + 
       (4*(-pow(mw,8) + pow(mt,6)*(pow(mw,2) + 2*pow(mu(t),2)) + 3*pow(mt,2)*(pow(mw,6) + 2*pow(mw,4)*pow(mu(t),2)) - 
            3*pow(mt,4)*(pow(mw,4) + 2*pow(mw,2)*pow(mu(t),2) + 2*pow(mu(t),4)))*
          atan((pow(mt,2) - pow(mw,2))/sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/
        pow(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)),1.5) - 
       (2*(pow(mt,8) + 2*pow(mw,8) - pow(mt,6)*(5*pow(mw,2) + 6*pow(mu(t),2)) - pow(mt,2)*(7*pow(mw,6) + 12*pow(mw,4)*pow(mu(t),2)) + 
            3*pow(mt,4)*(3*pow(mw,4) + 6*pow(mw,2)*pow(mu(t),2) + 4*pow(mu(t),4)))*
          atan((-pow(mt,2) + pow(mw,2))/sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/
        pow(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)),1.5) + 
       (2*(pow(mt,8) + 2*pow(mw,8) - pow(mt,6)*(5*pow(mw,2) + 6*pow(mu(t),2)) - pow(mt,2)*(7*pow(mw,6) + 12*pow(mw,4)*pow(mu(t),2)) + 
            3*pow(mt,4)*(3*pow(mw,4) + 6*pow(mw,2)*pow(mu(t),2) + 4*pow(mu(t),4)))*
          atan((pow(mt,2) + pow(mw,2))/sqrt(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)))))/
        pow(-pow(mt,4) - pow(mw,4) + 2*pow(mt,2)*(pow(mw,2) + 2*pow(mu(t),2)),1.5) + 2*pow(mw,2)*log(pow(mu(t),2)) + (pow(mt,2) - 2*pow(mw,2))*log(pow(mu(t),2)) - 
       2*pow(mw,2)*log(pow(mw,2) + pow(mu(t),2)) - (pow(mt,2) - 2*pow(mw,2))*log(pow(mw,2) + pow(mu(t),2))))/(2.*pow(mt,6)));
	return I;
}

double Ivbltlwluc(double t, double alpha2, double alphay){
	double I;
	I =-(pow(mu(t),4)*((pow(alpha2 - alphay,2)*(mew - mu(t)))/(pow(pow(mt,2) - pow(mw,2),2)*mu(t)) + 
     pow(alpha2*(pow(mt,2) + pow(mu(t),2)) - alphay*(pow(mw,2) + pow(mu(t),2)),2)/(pow(pow(mt,2) - pow(mw,2),3)*(pow(mw,2) + pow(mu(t),2))) - 
     (mu(t)*(pow(mt,2) + pow(mu(t),2))*pow(alpha2*(pow(mt,2) + pow(mu(t),2)) - alphay*(pow(mw,2) + pow(mu(t),2)),2))/
      (pow(pow(mt,2) - pow(mw,2),3)*(pow(mw,2) + pow(mu(t),2))*(mew*pow(mt,2) - mew*pow(mw,2) + pow(mw,2)*mu(t) + pow(mu(t),3))) - 
     (pow(alpha2,2)*log(mew/mu(t)))/pow(pow(mw,2) + pow(mu(t),2),2) - ((pow(alpha2,2)*(pow(mt,2) - 3*pow(mw,2) - 2*pow(mu(t),2))*pow(pow(mt,2) + pow(mu(t),2),2) + 
          4*alpha2*alphay*(pow(mt,2) + pow(mu(t),2))*pow(pow(mw,2) + pow(mu(t),2),2) - 
          pow(alphay,2)*pow(pow(mw,2) + pow(mu(t),2),2)*(pow(mt,2) + pow(mw,2) + 2*pow(mu(t),2)))*log(pow(mt,2) + pow(mu(t),2)))/
      (pow(pow(mt,2) - pow(mw,2),3)*pow(pow(mw,2) + pow(mu(t),2),2)) + 
     (pow(alpha2 - alphay,2)*(-pow(mt,2) + pow(mw,2)) + pow(alpha2*(pow(mt,2) + pow(mu(t),2)) - alphay*(pow(mw,2) + pow(mu(t),2)),2)/(pow(mt,2) + pow(mu(t),2)) + 
        2*(pow(alpha2,2)*(pow(mt,2) + pow(mu(t),2)) + pow(alphay,2)*(pow(mw,2) + pow(mu(t),2)) - alpha2*alphay*(pow(mt,2) + pow(mw,2) + 2*pow(mu(t),2)))*
         log(pow(mt,2) + pow(mu(t),2)))/pow(-pow(mt,2) + pow(mw,2),3) - 
     (pow(alpha2*(pow(mt,2) + pow(mu(t),2)) - alphay*(pow(mw,2) + pow(mu(t),2)),2)/(pow(mw,2) + pow(mu(t),2)) + 
        2*(pow(alpha2,2)*(pow(mt,2) + pow(mu(t),2)) + pow(alphay,2)*(pow(mw,2) + pow(mu(t),2)) - alpha2*alphay*(pow(mt,2) + pow(mw,2) + 2*pow(mu(t),2)))*
         log(pow(mw,2) + pow(mu(t),2)))/pow(-pow(mt,2) + pow(mw,2),3) + 
     ((pow(alpha2,2)*(pow(mt,2) - 3*pow(mw,2) - 2*pow(mu(t),2))*pow(pow(mt,2) + pow(mu(t),2),2) + 
          4*alpha2*alphay*(pow(mt,2) + pow(mu(t),2))*pow(pow(mw,2) + pow(mu(t),2),2) - 
          pow(alphay,2)*pow(pow(mw,2) + pow(mu(t),2),2)*(pow(mt,2) + pow(mw,2) + 2*pow(mu(t),2)))*log((mew*pow(mt,2) - mew*pow(mw,2) + pow(mw,2)*mu(t) + pow(mu(t),3))/mu(t)))/
      (pow(pow(mt,2) - pow(mw,2),3)*pow(pow(mw,2) + pow(mu(t),2),2))));
	return I;
}

double Ivbltrwtuc(double t){
	double I;
	I =-((pow(mu(t),4)*(pow(mt,2) - 2*pow(mw,2) - pow(mu(t),2) + pow(pow(mw,2) + pow(mu(t),2),2)/(pow(mt,2) + pow(mu(t),2)) - 
       (pow(mt,2) + pow(mw,2) + 2*pow(mu(t),2))*log(pow(mt,2) + pow(mu(t),2)) + 2*(pow(mw,2) + pow(mu(t),2))*log((pow(mt,2) + pow(mu(t),2))/(pow(mw,2) + pow(mu(t),2))) + 
       (pow(mt,2) + pow(mw,2) + 2*pow(mu(t),2))*log(pow(mw,2) + pow(mu(t),2))))/pow(-pow(mt,2) + pow(mw,2),3));
	return I;
}

double Ivzlwlwluc(double t){
	double I;
	I =-(pow(mu(t),4)*((-2*pow(-1 + pow(tw(t),2),2))/pow(mz,4) - (2*pow(mew,2)*pow(-1 + pow(tw(t),2),2))/(pow(mz,4)*pow(mu(t),2)) + 
     (4*mew*pow(-1 + pow(tw(t),2),2))/(pow(mz,4)*mu(t)) - (2*pow(-2*pow(mz,2) + pow(mw,2)*(-1 + pow(tw(t),2)) + (-1 + pow(tw(t),2))*pow(mu(t),2),2))/
      (pow(mz,6)*(pow(mw,2) + pow(mu(t),2))) + (mu(t)*pow(-2*pow(mz,2) + pow(mw,2)*(-1 + pow(tw(t),2)) + (-1 + pow(tw(t),2))*pow(mu(t),2),2)*
        (mew*pow(mz,2) + 2*pow(mw,2)*mu(t) - pow(mz,2)*mu(t) + 2*pow(mu(t),3)))/
      (pow(mz,6)*(pow(mw,2) + pow(mu(t),2))*(pow(mew,2)*pow(mz,2) - mew*pow(mz,2)*mu(t) + pow(mw,2)*pow(mu(t),2) + pow(mu(t),4))) + 
     ((4*pow(mz,6) + 6*pow(mw,6)*pow(-1 + pow(tw(t),2),2) - 8*pow(mz,4)*pow(mu(t),2) - pow(mz,2)*(-7 + 6*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 
          6*pow(-1 + pow(tw(t),2),2)*pow(mu(t),6) - pow(mw,4)*(-1 + pow(tw(t),2))*(pow(mz,2)*(7 + pow(tw(t),2)) - 18*(-1 + pow(tw(t),2))*pow(mu(t),2)) - 
          2*pow(mw,2)*(4*pow(mz,4) + pow(mz,2)*(-7 + 6*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) - 9*pow(-1 + pow(tw(t),2),2)*pow(mu(t),4)))*
        atan(mz/sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2))))/(pow(mz,5)*pow(pow(mw,2) + pow(mu(t),2),2)*sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2))) + 
     ((4*pow(mz,6) + 6*pow(mw,6)*pow(-1 + pow(tw(t),2),2) - 8*pow(mz,4)*pow(mu(t),2) - pow(mz,2)*(-7 + 6*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 
          6*pow(-1 + pow(tw(t),2),2)*pow(mu(t),6) - pow(mw,4)*(-1 + pow(tw(t),2))*(pow(mz,2)*(7 + pow(tw(t),2)) - 18*(-1 + pow(tw(t),2))*pow(mu(t),2)) - 
          2*pow(mw,2)*(4*pow(mz,4) + pow(mz,2)*(-7 + 6*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) - 9*pow(-1 + pow(tw(t),2),2)*pow(mu(t),4)))*
        atan((mz*(-2*mew + mu(t)))/(mu(t)*sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2)))))/
      (pow(mz,5)*pow(pow(mw,2) + pow(mu(t),2),2)*sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2))) - (4*log(mew/mu(t)))/pow(pow(mw,2) + pow(mu(t),2),2) - 
     ((4*pow(mz,6) + 8*pow(mw,6)*pow(-1 + pow(tw(t),2),2) - pow(mz,2)*(-15 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 8*pow(-1 + pow(tw(t),2),2)*pow(mu(t),6) - 
          pow(mw,4)*(-1 + pow(tw(t),2))*(pow(mz,2)*(15 + pow(tw(t),2)) - 24*(-1 + pow(tw(t),2))*pow(mu(t),2)) - 
          2*pow(mw,2)*(-1 + pow(tw(t),2))*pow(mu(t),2)*(pow(mz,2)*(15 + pow(tw(t),2)) - 12*(-1 + pow(tw(t),2))*pow(mu(t),2)))*log(pow(mw,2) + pow(mu(t),2)))/
      (2.*pow(mz,6)*pow(pow(mw,2) + pow(mu(t),2),2)) + ((4*pow(mz,6) + 8*pow(mw,6)*pow(-1 + pow(tw(t),2),2) - 
          pow(mz,2)*(-15 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 8*pow(-1 + pow(tw(t),2),2)*pow(mu(t),6) - 
          pow(mw,4)*(-1 + pow(tw(t),2))*(pow(mz,2)*(15 + pow(tw(t),2)) - 24*(-1 + pow(tw(t),2))*pow(mu(t),2)) - 
          2*pow(mw,2)*(-1 + pow(tw(t),2))*pow(mu(t),2)*(pow(mz,2)*(15 + pow(tw(t),2)) - 12*(-1 + pow(tw(t),2))*pow(mu(t),2)))*
        log(pow(mw,2) + (mew*pow(mz,2)*(mew - mu(t)))/pow(mu(t),2) + pow(mu(t),2)))/(2.*pow(mz,6)*pow(pow(mw,2) + pow(mu(t),2),2))));
	return I;
}

double Ivzlhzluc(double t, double alphaz){
	double I;
	I =-((pow(mu(t),4)*((8*alphaz*pi*pow(v,2)*(alphaz*pow(mz,2)*pi*pow(v,2) + pow(mh,2)*(pow(mz,2) - 2*alphaz*pi*pow(v,2)))*(1 - mew/mu(t)))/pow(mz,6) - 
       (4*pow(alphaz,2)*pow(pi,2)*pow(v,4)*pow(mew - mu(t),2))/(pow(mz,4)*pow(mu(t),2)) + 
       (2*mu(t)*(mu(t)*(pow(mh,8)*pow(pow(mz,2) - 2*alphaz*pi*pow(v,2),2)*pow(pow(mz,2) + pow(mu(t),2),2) - 
               2*pow(mh,6)*pow(mz,2)*(pow(mz,4) - 6*alphaz*pow(mz,2)*pi*pow(v,2) + 8*pow(alphaz,2)*pow(pi,2)*pow(v,4))*pow(pow(mz,2) + pow(mu(t),2),2) + 
               8*pow(alphaz,2)*pow(mz,4)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(-pow(mz,6) + pow(mu(t),6)) - 
               8*alphaz*pow(mh,2)*pow(mz,4)*pi*pow(v,2)*(pow(mz,4)*pow(mu(t),2)*(pow(mz,2) + pow(mu(t),2)) + 
                  alphaz*pi*pow(v,2)*(pow(mz,6) - 2*pow(mz,4)*pow(mu(t),2) - 3*pow(mz,2)*pow(mu(t),4) - 2*pow(mu(t),6))) - 
               2*pow(mh,4)*pow(mz,2)*(pow(mz,4)*pow(mu(t),2)*pow(pow(mz,2) + pow(mu(t),2),2) + 
                  2*alphaz*pow(mz,2)*pi*pow(v,2)*(2*pow(mz,6) - 5*pow(mz,2)*pow(mu(t),4) - 3*pow(mu(t),6)) + 
                  2*pow(alphaz,2)*pow(pi,2)*pow(v,4)*(-5*pow(mz,6) - 4*pow(mz,4)*pow(mu(t),2) + 4*pow(mz,2)*pow(mu(t),4) + 4*pow(mu(t),6)))) + 
            mew*(pow(mh,10)*pow(pow(mz,2) - 2*alphaz*pi*pow(v,2),2)*(pow(mz,2) + pow(mu(t),2)) - 
               4*pow(mh,8)*pow(mz,2)*(pow(mz,4) - 5*alphaz*pow(mz,2)*pi*pow(v,2) + 6*pow(alphaz,2)*pow(pi,2)*pow(v,4))*(pow(mz,2) + pow(mu(t),2)) + 
               8*pow(alphaz,2)*pow(mz,6)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mz,4) - pow(mz,2)*pow(mu(t),2) - pow(mu(t),4)) + 
               2*pow(mh,4)*pow(mz,4)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,4)*pow(mu(t),2) + 2*alphaz*pow(mz,2)*pi*pow(v,2)*(2*pow(mz,2) - 7*pow(mu(t),2)) - 
                  4*pow(alphaz,2)*pow(pi,2)*pow(v,4)*(4*pow(mz,2) - 7*pow(mu(t),2))) + 
               pow(mh,6)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2))*(2*pow(mz,6) - 20*pow(alphaz,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) - 
                  3*pow(mz,4)*(8*alphaz*pi*pow(v,2) + pow(mu(t),2)) + 4*alphaz*pow(mz,2)*pi*pow(v,2)*(11*alphaz*pi*pow(v,2) + 4*pow(mu(t),2))) + 
               4*alphaz*pow(mh,2)*pow(mz,4)*pi*pow(v,2)*(2*pow(mz,2)*pow(mu(t),2)*(pow(mz,4) - pow(mu(t),4)) + 
                  alphaz*pi*pow(v,2)*(2*pow(mz,6) - 9*pow(mz,4)*pow(mu(t),2) - 5*pow(mz,2)*pow(mu(t),4) + 5*pow(mu(t),6))))))/
        (pow(mz,8)*(pow(mz,2) + pow(mu(t),2))*(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))*
          (pow(mew,2)*pow(mz,2) + mew*(pow(mh,2) - 2*pow(mz,2))*mu(t) + pow(mu(t),2)*(pow(mz,2) + pow(mu(t),2)))) + 
       (2*(pow(mh,10)*(pow(mz,4) - 8*alphaz*pow(mz,2)*pi*pow(v,2) + 12*pow(alphaz,2)*pow(pi,2)*pow(v,4))*pow(pow(mz,2) + pow(mu(t),2),2) - 
            2*pow(mh,8)*pow(mz,2)*(3*pow(mz,4) - 26*alphaz*pow(mz,2)*pi*pow(v,2) + 44*pow(alphaz,2)*pow(pi,2)*pow(v,4))*pow(pow(mz,2) + pow(mu(t),2),2) + 
            16*pow(alphaz,2)*pow(mz,6)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mz,6) - 4*pow(mz,4)*pow(mu(t),2) - 6*pow(mz,2)*pow(mu(t),4) - 3*pow(mu(t),6)) + 
            2*pow(mh,6)*pow(mz,2)*(2*pow(mz,10) - 40*pow(alphaz,2)*pow(pi,2)*pow(v,4)*pow(mu(t),6) + pow(mz,8)*(-36*alphaz*pi*pow(v,2) + pow(mu(t),2)) + 
               8*alphaz*pow(mz,2)*pi*pow(v,2)*pow(mu(t),4)*(alphaz*pi*pow(v,2) + 3*pow(mu(t),2)) + 
               pow(mz,6)*(90*pow(alphaz,2)*pow(pi,2)*pow(v,4) - 48*alphaz*pi*pow(v,2)*pow(mu(t),2) - 4*pow(mu(t),4)) + 
               pow(mz,4)*(136*pow(alphaz,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) + 12*alphaz*pi*pow(v,2)*pow(mu(t),4) - 3*pow(mu(t),6))) + 
            4*pow(mh,4)*pow(mz,4)*(pow(mz,4)*pow(mu(t),2)*pow(pow(mz,2) + pow(mu(t),2),2) + 
               2*alphaz*pow(mz,2)*pi*pow(v,2)*(2*pow(mz,2) - 15*pow(mu(t),2))*pow(pow(mz,2) + pow(mu(t),2),2) + 
               6*pow(alphaz,2)*pow(pi,2)*pow(v,4)*(-5*pow(mz,6) + 4*pow(mz,4)*pow(mu(t),2) + 20*pow(mz,2)*pow(mu(t),4) + 12*pow(mu(t),6))) + 
            8*alphaz*pow(mh,2)*pow(mz,4)*pi*pow(v,2)*(2*pow(mz,2)*pow(mu(t),2)*(pow(mz,2) - 3*pow(mu(t),2))*pow(pow(mz,2) + pow(mu(t),2),2) + 
               alphaz*pi*pow(v,2)*(2*pow(mz,8) - 23*pow(mz,6)*pow(mu(t),2) - 22*pow(mz,4)*pow(mu(t),4) + 12*pow(mz,2)*pow(mu(t),6) + 15*pow(mu(t),8))))*
          atan((-pow(mh,2) + 2*pow(mz,2)*(1 - mew/mu(t)))/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))))/
        (pow(mz,8)*pow(pow(mz,2) + pow(mu(t),2),2)*pow(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2),1.5)) - 
       (8*pow(alphaz,2)*pow(pi,2)*pow(v,4)*log(mew/mu(t)))/pow(pow(mz,2) + pow(mu(t),2),2) - 
       (4*pow(alphaz,2)*pow(mz,4)*pow(pi,2)*pow(v,4) + 8*alphaz*pow(mz,2)*pi*pow(v,2)*
           (alphaz*pow(mz,2)*pi*pow(v,2) + pow(mh,2)*(pow(mz,2) - 2*alphaz*pi*pow(v,2))) + 
          (2*(pow(mh,8)*pow(pow(mz,2) - 2*alphaz*pi*pow(v,2),2) - 
               4*pow(mh,6)*(pow(mz,6) - 5*alphaz*pow(mz,4)*pi*pow(v,2) + 6*pow(alphaz,2)*pow(mz,2)*pow(pi,2)*pow(v,4)) + 
               8*pow(alphaz,2)*pow(mz,4)*pow(pi,2)*pow(v,4)*(pow(mz,4) - 3*pow(mz,2)*pow(mu(t),2) + pow(mu(t),4)) - 
               8*alphaz*pow(mh,2)*pow(mz,4)*pi*pow(v,2)*(-pow(mz,4) - 5*alphaz*pi*pow(v,2)*pow(mu(t),2) + 2*pow(mz,2)*(2*alphaz*pi*pow(v,2) + pow(mu(t),2))) + 
               2*pow(mh,4)*(pow(mz,8) - 8*pow(alphaz,2)*pow(mz,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) - pow(mz,6)*(12*alphaz*pi*pow(v,2) + pow(mu(t),2)) + 
                  2*alphaz*pow(mz,4)*pi*pow(v,2)*(11*alphaz*pi*pow(v,2) + 3*pow(mu(t),2)))))/(pow(mh,4) - 4*pow(mh,2)*pow(mz,2) - 4*pow(mz,2)*pow(mu(t),2)) + 
          (2*(pow(mh,10)*(pow(mz,4) - 8*alphaz*pow(mz,2)*pi*pow(v,2) + 12*pow(alphaz,2)*pow(pi,2)*pow(v,4)) + 
               pow(mh,8)*(-6*pow(mz,6) + 56*alphaz*pow(mz,4)*pi*pow(v,2) - 96*pow(alphaz,2)*pow(mz,2)*pow(pi,2)*pow(v,4)) + 
               32*pow(alphaz,2)*pow(mz,6)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mz,2) - 3*pow(mu(t),2)) - 
               8*pow(mh,4)*(-(pow(mz,8)*pow(mu(t),2)) - 3*alphaz*pow(mz,6)*pi*pow(v,2)*(pow(mz,2) - 6*pow(mu(t),2)) + 
                  21*pow(alphaz,2)*pow(mz,4)*pow(pi,2)*pow(v,4)*(pow(mz,2) - 2*pow(mu(t),2))) + 
               pow(mh,6)*(6*pow(mz,8) - 80*pow(alphaz,2)*pow(mz,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) - 6*pow(mz,6)*(16*alphaz*pi*pow(v,2) + pow(mu(t),2)) + 
                  12*alphaz*pow(mz,4)*pi*pow(v,2)*(19*alphaz*pi*pow(v,2) + 4*pow(mu(t),2))) + 
               8*alphaz*pow(mh,2)*pow(mz,4)*pi*pow(v,2)*(4*pow(mz,4)*pow(mu(t),2) - 6*pow(mz,2)*pow(mu(t),4) + 
                  3*alphaz*pi*pow(v,2)*(pow(mz,4) - 11*pow(mz,2)*pow(mu(t),2) + 5*pow(mu(t),4))))*
             atan(pow(mh,2)/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))))/
           pow(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2),1.5) + 
          (8*alphaz*pow(mh,2)*pow(mz,2)*pi*pow(v,2)*(pow(mz,2) - 3*alphaz*pi*pow(v,2)) + 
             pow(mh,4)*(pow(mz,4) - 8*alphaz*pow(mz,2)*pi*pow(v,2) + 12*pow(alphaz,2)*pow(pi,2)*pow(v,4)) + 
             4*pow(alphaz,2)*pow(mz,2)*pow(pi,2)*pow(v,4)*(3*pow(mz,2) - 2*pow(mu(t),2)))*log(pow(mh,2) + pow(mu(t),2)))/pow(mz,8) + 
       ((-2*(pow(mz,2) + pow(mu(t),2))*(-(pow(mh,8)*pow(pow(mz,2) - 2*alphaz*pi*pow(v,2),2)*(pow(mz,2) + pow(mu(t),2))) + 
               pow(mh,6)*pow(mz,2)*(3*pow(mz,4) - 16*alphaz*pow(mz,2)*pi*pow(v,2) + 20*pow(alphaz,2)*pow(pi,2)*pow(v,4))*(pow(mz,2) + pow(mu(t),2)) + 
               8*pow(alphaz,2)*pow(mz,4)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mz,4) + pow(mz,2)*pow(mu(t),2) - pow(mu(t),4)) + 
               2*pow(mh,4)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,4)*pow(mu(t),2) + 6*alphaz*pow(mz,2)*pi*pow(v,2)*(pow(mz,2) - pow(mu(t),2)) + 
                  2*pow(alphaz,2)*pow(pi,2)*pow(v,4)*(-7*pow(mz,2) + 4*pow(mu(t),2))) + 
               4*alphaz*pow(mh,2)*pow(mz,4)*pi*pow(v,2)*(2*pow(mz,2)*pow(mu(t),2)*(pow(mz,2) + pow(mu(t),2)) + 
                  alphaz*pi*pow(v,2)*(3*pow(mz,4) - 3*pow(mz,2)*pow(mu(t),2) - 7*pow(mu(t),4)))))/(pow(mh,4) - 4*pow(mh,2)*pow(mz,2) - 4*pow(mz,2)*pow(mu(t),2)) + 
          (2*(pow(mh,10)*(pow(mz,4) - 8*alphaz*pow(mz,2)*pi*pow(v,2) + 12*pow(alphaz,2)*pow(pi,2)*pow(v,4))*pow(pow(mz,2) + pow(mu(t),2),2) - 
               2*pow(mh,8)*pow(mz,2)*(3*pow(mz,4) - 26*alphaz*pow(mz,2)*pi*pow(v,2) + 44*pow(alphaz,2)*pow(pi,2)*pow(v,4))*pow(pow(mz,2) + pow(mu(t),2),2) + 
               16*pow(alphaz,2)*pow(mz,6)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mz,6) - 4*pow(mz,4)*pow(mu(t),2) - 6*pow(mz,2)*pow(mu(t),4) - 3*pow(mu(t),6)) + 
               2*pow(mh,6)*pow(mz,2)*(2*pow(mz,10) - 40*pow(alphaz,2)*pow(pi,2)*pow(v,4)*pow(mu(t),6) + pow(mz,8)*(-36*alphaz*pi*pow(v,2) + pow(mu(t),2)) + 
                  8*alphaz*pow(mz,2)*pi*pow(v,2)*pow(mu(t),4)*(alphaz*pi*pow(v,2) + 3*pow(mu(t),2)) + 
                  pow(mz,6)*(90*pow(alphaz,2)*pow(pi,2)*pow(v,4) - 48*alphaz*pi*pow(v,2)*pow(mu(t),2) - 4*pow(mu(t),4)) + 
                  pow(mz,4)*(136*pow(alphaz,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) + 12*alphaz*pi*pow(v,2)*pow(mu(t),4) - 3*pow(mu(t),6))) + 
               4*pow(mh,4)*pow(mz,4)*(pow(mz,4)*pow(mu(t),2)*pow(pow(mz,2) + pow(mu(t),2),2) + 
                  2*alphaz*pow(mz,2)*pi*pow(v,2)*(2*pow(mz,2) - 15*pow(mu(t),2))*pow(pow(mz,2) + pow(mu(t),2),2) + 
                  6*pow(alphaz,2)*pow(pi,2)*pow(v,4)*(-5*pow(mz,6) + 4*pow(mz,4)*pow(mu(t),2) + 20*pow(mz,2)*pow(mu(t),4) + 12*pow(mu(t),6))) + 
               8*alphaz*pow(mh,2)*pow(mz,4)*pi*pow(v,2)*(2*pow(mz,2)*pow(mu(t),2)*(pow(mz,2) - 3*pow(mu(t),2))*pow(pow(mz,2) + pow(mu(t),2),2) + 
                  alphaz*pi*pow(v,2)*(2*pow(mz,8) - 23*pow(mz,6)*pow(mu(t),2) - 22*pow(mz,4)*pow(mu(t),4) + 12*pow(mz,2)*pow(mu(t),6) + 15*pow(mu(t),8))))*
             atan(pow(mh,2)/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))))/
           pow(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2),1.5) + 
          (-4*alphaz*pow(mh,2)*pow(mz,2)*pi*pow(v,2)*(-pow(mz,2) + 4*alphaz*pi*pow(v,2))*pow(pow(mz,2) + pow(mu(t),2),2) + 
             pow(mh,4)*(pow(mz,4) - 8*alphaz*pow(mz,2)*pi*pow(v,2) + 12*pow(alphaz,2)*pow(pi,2)*pow(v,4))*pow(pow(mz,2) + pow(mu(t),2),2) + 
             4*pow(alphaz,2)*pow(mz,2)*pow(pi,2)*pow(v,4)*(pow(mz,6) + 2*pow(mz,4)*pow(mu(t),2) - 2*pow(mz,2)*pow(mu(t),4) - 2*pow(mu(t),6)))*
           log(pow(mh,2) + pow(mu(t),2)))/(pow(mz,8)*pow(pow(mz,2) + pow(mu(t),2),2)) + 
       ((2*(-(pow(mh,8)*pow(pow(mz,2) - 2*alphaz*pi*pow(v,2),2)*(pow(mz,2) + pow(mu(t),2))) + 
               pow(mh,6)*pow(mz,2)*(3*pow(mz,4) - 16*alphaz*pow(mz,2)*pi*pow(v,2) + 20*pow(alphaz,2)*pow(pi,2)*pow(v,4))*(pow(mz,2) + pow(mu(t),2)) + 
               8*pow(alphaz,2)*pow(mz,4)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mz,4) + pow(mz,2)*pow(mu(t),2) - pow(mu(t),4)) + 
               2*pow(mh,4)*pow(mz,2)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,4)*pow(mu(t),2) + 6*alphaz*pow(mz,2)*pi*pow(v,2)*(pow(mz,2) - pow(mu(t),2)) + 
                  2*pow(alphaz,2)*pow(pi,2)*pow(v,4)*(-7*pow(mz,2) + 4*pow(mu(t),2))) + 
               4*alphaz*pow(mh,2)*pow(mz,4)*pi*pow(v,2)*(2*pow(mz,2)*pow(mu(t),2)*(pow(mz,2) + pow(mu(t),2)) + 
                  alphaz*pi*pow(v,2)*(3*pow(mz,4) - 3*pow(mz,2)*pow(mu(t),2) - 7*pow(mu(t),4)))))/
           ((pow(mz,2) + pow(mu(t),2))*(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))) + 
          (2*(pow(mh,10)*(pow(mz,4) - 8*alphaz*pow(mz,2)*pi*pow(v,2) + 12*pow(alphaz,2)*pow(pi,2)*pow(v,4)) + 
               pow(mh,8)*(-6*pow(mz,6) + 56*alphaz*pow(mz,4)*pi*pow(v,2) - 96*pow(alphaz,2)*pow(mz,2)*pow(pi,2)*pow(v,4)) + 
               32*pow(alphaz,2)*pow(mz,6)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mz,2) - 3*pow(mu(t),2)) - 
               8*pow(mh,4)*(-(pow(mz,8)*pow(mu(t),2)) - 3*alphaz*pow(mz,6)*pi*pow(v,2)*(pow(mz,2) - 6*pow(mu(t),2)) + 
                  21*pow(alphaz,2)*pow(mz,4)*pow(pi,2)*pow(v,4)*(pow(mz,2) - 2*pow(mu(t),2))) + 
               pow(mh,6)*(6*pow(mz,8) - 80*pow(alphaz,2)*pow(mz,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) - 6*pow(mz,6)*(16*alphaz*pi*pow(v,2) + pow(mu(t),2)) + 
                  12*alphaz*pow(mz,4)*pi*pow(v,2)*(19*alphaz*pi*pow(v,2) + 4*pow(mu(t),2))) + 
               8*alphaz*pow(mh,2)*pow(mz,4)*pi*pow(v,2)*(4*pow(mz,4)*pow(mu(t),2) - 6*pow(mz,2)*pow(mu(t),4) + 
                  3*alphaz*pi*pow(v,2)*(pow(mz,4) - 11*pow(mz,2)*pow(mu(t),2) + 5*pow(mu(t),4))))*
             atan((pow(mh,2) - 2*pow(mz,2))/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))))/
           pow(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2),1.5) + 
          (8*alphaz*pow(mh,2)*pow(mz,2)*pi*pow(v,2)*(pow(mz,2) - 3*alphaz*pi*pow(v,2)) + 
             pow(mh,4)*(pow(mz,4) - 8*alphaz*pow(mz,2)*pi*pow(v,2) + 12*pow(alphaz,2)*pow(pi,2)*pow(v,4)) + 
             4*pow(alphaz,2)*pow(mz,2)*pow(pi,2)*pow(v,4)*(3*pow(mz,2) - 2*pow(mu(t),2)))*log(pow(mz,2) + pow(mu(t),2)))/pow(mz,8) - 
       ((-4*alphaz*pow(mh,2)*pow(mz,2)*pi*pow(v,2)*(-pow(mz,2) + 4*alphaz*pi*pow(v,2))*pow(pow(mz,2) + pow(mu(t),2),2) + 
            pow(mh,4)*(pow(mz,4) - 8*alphaz*pow(mz,2)*pi*pow(v,2) + 12*pow(alphaz,2)*pow(pi,2)*pow(v,4))*pow(pow(mz,2) + pow(mu(t),2),2) + 
            4*pow(alphaz,2)*pow(mz,2)*pow(pi,2)*pow(v,4)*(pow(mz,6) + 2*pow(mz,4)*pow(mu(t),2) - 2*pow(mz,2)*pow(mu(t),4) - 2*pow(mu(t),6)))*
          log((pow(mz,2)*pow(mew - mu(t),2) + mew*pow(mh,2)*mu(t) + pow(mu(t),4))/pow(mu(t),2)))/(pow(mz,8)*pow(pow(mz,2) + pow(mu(t),2),2))))/(8.*pow(pi,2)*pow(v,4)));
	return I;
}

double Ivzlwtwtuc(double t){
	double I;
	I =-((pow(mu(t),4)*(-3*mz + (2*(6*pow(mw,2) - pow(mz,2) + 6*pow(mu(t),2))*atan(mz/sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2))))/
        sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2))))/pow(mz,5));
	return I;
}

double Ivvlffuc(double t, double mv){
	double I;
	I =-(-0.5*(mv*pow(mu(t),4)*(pow(mv,2) - 6*pow(mu(t),2))*sqrt(-pow(mv,2) + 4*pow(mu(t),2)) - 8*pow(mu(t),6)*(pow(mv,2) - 3*pow(mu(t),2))*atan(mv/sqrt(-pow(mv,2) + 4*pow(mu(t),2))))/
    (pow(mv,5)*pow(-pow(mv,2) + 4*pow(mu(t),2),1.5)));
	return I;
}

double Ivzlttuc(double t, double alphaz, double alphay){
	double I;
	I =-((pow(mu(t),4)*(pow(alphaz,2)*pow(mz,2)*(9 - 24*pow(sw(t),2) + 32*pow(sw(t),4)) + 
       (4*(18*pow(alphay,2)*pow(mz,4) + 6*alphay*alphaz*pow(mz,2)*(-3 + 8*pow(sw(t),2))*(pow(mt,2) + pow(mu(t),2)) + 
            pow(alphaz,2)*(9 - 24*pow(sw(t),2) + 32*pow(sw(t),4))*pow(pow(mt,2) + pow(mu(t),2),2)))/(4*pow(mt,2) - pow(mz,2) + 4*pow(mu(t),2)) - 
       (2*(2*pow(mt,2) - pow(mz,2) + 2*pow(mu(t),2))*(18*pow(alphay,2)*pow(mz,4) + 6*alphay*alphaz*pow(mz,2)*(-3 + 8*pow(sw(t),2))*(pow(mt,2) + pow(mu(t),2)) + 
            pow(alphaz,2)*(9 - 24*pow(sw(t),2) + 32*pow(sw(t),4))*pow(pow(mt,2) + pow(mu(t),2),2)))/((pow(mt,2) + pow(mu(t),2))*(4*pow(mt,2) - pow(mz,2) + 4*pow(mu(t),2))) + 
       (8*mz*(18*pow(alphay,2)*pow(mz,4) + 3*alphay*alphaz*pow(mz,2)*(-3 + 8*pow(sw(t),2))*(-2*pow(mt,2) + pow(mz,2) - 2*pow(mu(t),2)) - 
            pow(alphaz,2)*(9 - 24*pow(sw(t),2) + 32*pow(sw(t),4))*(pow(mt,2) + pow(mu(t),2))*(3*pow(mt,2) - pow(mz,2) + 3*pow(mu(t),2)))*
          atan(mz/sqrt(4*pow(mt,2) - pow(mz,2) + 4*pow(mu(t),2))))/pow(4*pow(mt,2) - pow(mz,2) + 4*pow(mu(t),2),1.5)))/(72.*pow(mz,6)));
	return I;
}

double Ivwlwlzluc(double t){
	double I;
	I =-(-0.5*(pow(mu(t),4)*(-4*pow(mw,4) + 8*pow(mw,2)*(-2*pow(mz,2) + pow(mw,2)*(3 + pow(tw(t),2))) + 
         (8*pow(mz,8)*pow(mu(t),2)*(pow(mz,2) + pow(mu(t),2)) + 8*pow(mw,10)*pow(-1 + pow(tw(t),2),2)*(3*pow(mz,2) + 4*pow(mu(t),2)) + 
            8*pow(mw,2)*pow(mz,4)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,4) - pow(mz,2)*(6 + pow(tw(t),2))*pow(mu(t),2) - 4*pow(mu(t),4)) - 
            2*pow(mw,8)*(-1 + pow(tw(t),2))*(pow(mz,4)*(41 + 7*pow(tw(t),2)) + pow(mz,2)*(83 - 3*pow(tw(t),2))*pow(mu(t),2) + 2*(23 - 7*pow(tw(t),2))*pow(mu(t),4)) + 
            2*pow(mw,6)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,4)*(21 + 26*pow(tw(t),2) + pow(tw(t),4)) + pow(mz,2)*(77 - 22*pow(tw(t),2) - 7*pow(tw(t),4))*pow(mu(t),2) - 
               2*(-15 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4)) - 2*pow(mw,4)*(pow(mz,2) + pow(mu(t),2))*
             (4*pow(mz,6)*(6 + pow(tw(t),2)) - pow(mz,4)*(5 + 26*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) - 4*pow(mz,2)*(10 + 3*pow(tw(t),2))*pow(mu(t),4) - 8*pow(mu(t),6)))/
          ((pow(mw,2) + pow(mu(t),2))*(pow(mz,2) + pow(mu(t),2))*(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))) - 
         (2*(12*pow(mz,10)*pow(pow(mz,2) + pow(mu(t),2),2) - 8*pow(mw,2)*pow(mz,6)*pow(pow(mz,2) + pow(mu(t),2),2)*(pow(mz,2)*(13 + pow(tw(t),2)) + 10*pow(mu(t),2)) - 
              6*pow(mw,6)*pow(pow(mz,2) + pow(mu(t),2),3)*(pow(mz,2)*(21 + 26*pow(tw(t),2) + pow(tw(t),4)) + 8*(3 + pow(tw(t),2))*pow(mu(t),2)) + 
              24*pow(mw,10)*pow(-1 + pow(tw(t),2),2)*(pow(mz,4) + 3*pow(mz,2)*pow(mu(t),2) + 2*pow(mu(t),4)) + 
              pow(mw,4)*pow(mz,2)*pow(pow(mz,2) + pow(mu(t),2),2)*
               (pow(mz,4)*(261 + 66*pow(tw(t),2) + pow(tw(t),4)) + 48*pow(mz,2)*(8 + pow(tw(t),2))*pow(mu(t),2) + 120*pow(mu(t),4)) + 
              2*pow(mw,8)*(pow(mz,6)*(-47 + 46*pow(tw(t),2) + pow(tw(t),4)) + 10*pow(mz,4)*(-15 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 
                 pow(mz,2)*(-153 + 130*pow(tw(t),2) + 23*pow(tw(t),4))*pow(mu(t),4) + 4*(-13 + 10*pow(tw(t),2) + 3*pow(tw(t),4))*pow(mu(t),6)))*
            atan((-2*pow(mw,2) + pow(mz,2))/sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))))/
          (pow(pow(mz,2) + pow(mu(t),2),2)*pow(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)),1.5)) + 
         (-12*pow(mz,4) - pow(mw,4)*(-3 + 18*pow(tw(t),2) + pow(tw(t),4)) + (4*pow(mw,8)*pow(-1 + pow(tw(t),2),2))/pow(pow(mz,2) + pow(mu(t),2),2) + 
            8*pow(mw,2)*(pow(mz,2)*(4 + pow(tw(t),2)) + pow(mu(t),2)))*log(pow(mw,2) + pow(mu(t),2))))/pow(mw,8) + 
   (pow(mu(t),4)*((16*pow(mw,10)*pow(-1 + pow(tw(t),2),2)*(pow(mz,2) + pow(mu(t),2)) - 8*pow(mz,8)*pow(pow(mz,2) + pow(mu(t),2),2) + 
           8*pow(mw,2)*pow(mz,4)*pow(pow(mz,2) + pow(mu(t),2),2)*(pow(mz,2)*(7 + pow(tw(t),2)) + 4*pow(mu(t),2)) - 
           2*pow(mw,4)*pow(pow(mz,2) + pow(mu(t),2),2)*(pow(mz,4)*(41 + 30*pow(tw(t),2) + pow(tw(t),4)) + 4*pow(mz,2)*(13 + 3*pow(tw(t),2))*pow(mu(t),2) + 8*pow(mu(t),4)) - 
           2*pow(mw,8)*(-1 + pow(tw(t),2))*(2*pow(mz,4)*(7 + 9*pow(tw(t),2)) + 8*pow(mz,2)*(3 + 5*pow(tw(t),2))*pow(mu(t),2) + 2*(3 + 13*pow(tw(t),2))*pow(mu(t),4)) + 
           2*pow(mw,6)*(pow(mz,2) + pow(mu(t),2))*(8*pow(mz,4)*(-4 + 7*pow(tw(t),2) + pow(tw(t),4)) + 2*pow(mz,2)*(-19 + 46*pow(tw(t),2) + 5*pow(tw(t),4))*pow(mu(t),2) + 
              2*(-3 + 18*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4)))/(pow(mz,4) - 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2))) - 
        (2*(12*pow(mz,10)*pow(pow(mz,2) + pow(mu(t),2),2) - 8*pow(mw,2)*pow(mz,6)*pow(pow(mz,2) + pow(mu(t),2),2)*(pow(mz,2)*(13 + pow(tw(t),2)) + 10*pow(mu(t),2)) - 
             6*pow(mw,6)*pow(pow(mz,2) + pow(mu(t),2),3)*(pow(mz,2)*(21 + 26*pow(tw(t),2) + pow(tw(t),4)) + 8*(3 + pow(tw(t),2))*pow(mu(t),2)) + 
             24*pow(mw,10)*pow(-1 + pow(tw(t),2),2)*(pow(mz,4) + 3*pow(mz,2)*pow(mu(t),2) + 2*pow(mu(t),4)) + 
             pow(mw,4)*pow(mz,2)*pow(pow(mz,2) + pow(mu(t),2),2)*(pow(mz,4)*(261 + 66*pow(tw(t),2) + pow(tw(t),4)) + 48*pow(mz,2)*(8 + pow(tw(t),2))*pow(mu(t),2) + 
                120*pow(mu(t),4)) + 2*pow(mw,8)*(pow(mz,6)*(-47 + 46*pow(tw(t),2) + pow(tw(t),4)) + 10*pow(mz,4)*(-15 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 
                pow(mz,2)*(-153 + 130*pow(tw(t),2) + 23*pow(tw(t),4))*pow(mu(t),4) + 4*(-13 + 10*pow(tw(t),2) + 3*pow(tw(t),4))*pow(mu(t),6)))*
           atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))))/pow(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)),1.5) + 
        (4*pow(mw,8)*pow(-1 + pow(tw(t),2),2) - 12*pow(mz,4)*pow(pow(mz,2) + pow(mu(t),2),2) - 
           pow(mw,4)*(-3 + 18*pow(tw(t),2) + pow(tw(t),4))*pow(pow(mz,2) + pow(mu(t),2),2) + 
           8*pow(mw,2)*pow(pow(mz,2) + pow(mu(t),2),2)*(pow(mz,2)*(4 + pow(tw(t),2)) + pow(mu(t),2)))*log(pow(mz,2) + pow(mu(t),2))))/
    (2.*pow(mw,8)*pow(pow(mz,2) + pow(mu(t),2),2)) + pow(mu(t),4)*(2/pow(mw,4) + (4*(-2*pow(mz,2) + pow(mw,2)*(1 + pow(tw(t),2))))/pow(mw,6) + 
      (-4*pow(mz,8)*pow(mu(t),2)*(pow(mz,2) + pow(mu(t),2)) - 4*pow(mw,10)*pow(-1 + pow(tw(t),2),2)*(3*pow(mz,2) + 4*pow(mu(t),2)) - 
         4*pow(mw,2)*pow(mz,4)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,4) - pow(mz,2)*(6 + pow(tw(t),2))*pow(mu(t),2) - 4*pow(mu(t),4)) + 
         pow(mw,8)*(-1 + pow(tw(t),2))*(pow(mz,4)*(41 + 7*pow(tw(t),2)) + pow(mz,2)*(83 - 3*pow(tw(t),2))*pow(mu(t),2) + 2*(23 - 7*pow(tw(t),2))*pow(mu(t),4)) - 
         pow(mw,6)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,4)*(21 + 26*pow(tw(t),2) + pow(tw(t),4)) + pow(mz,2)*(77 - 22*pow(tw(t),2) - 7*pow(tw(t),4))*pow(mu(t),2) - 
            2*(-15 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4)) + pow(mw,4)*(pow(mz,2) + pow(mu(t),2))*
          (4*pow(mz,6)*(6 + pow(tw(t),2)) - pow(mz,4)*(5 + 26*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) - 4*pow(mz,2)*(10 + 3*pow(tw(t),2))*pow(mu(t),4) - 8*pow(mu(t),6)))/
       (pow(mw,8)*(pow(mw,2) + pow(mu(t),2))*(pow(mz,2) + pow(mu(t),2))*(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))) + 
      ((16*pow(mw,14)*pow(-1 + pow(tw(t),2),2) + 12*pow(mz,10)*pow(mu(t),4) + 
           4*pow(mw,12)*(-1 + pow(tw(t),2))*(pow(mz,2)*(15 + pow(tw(t),2)) + (3 + 13*pow(tw(t),2))*pow(mu(t),2)) + 
           8*pow(mw,2)*pow(mz,6)*pow(mu(t),2)*(3*pow(mz,4) - pow(mz,2)*(12 + pow(tw(t),2))*pow(mu(t),2) - 10*pow(mu(t),4)) - 
           2*pow(mw,10)*(3*pow(mz,4)*(9 + 22*pow(tw(t),2) + pow(tw(t),4)) + pow(mz,2)*(159 + 34*pow(tw(t),2) - pow(tw(t),4))*pow(mu(t),2) + 
              4*(27 + 4*pow(tw(t),2) - 7*pow(tw(t),4))*pow(mu(t),4)) + pow(mw,8)*
            (pow(mz,6)*(209 + 62*pow(tw(t),2) + pow(tw(t),4)) - 12*pow(mz,4)*(-23 + 18*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) - 
              8*pow(mz,2)*(27 + 38*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 4*(-63 - 18*pow(tw(t),2) + 5*pow(tw(t),4))*pow(mu(t),6)) + 
           pow(mw,4)*pow(mz,2)*(12*pow(mz,8) - 16*pow(mz,6)*(12 + pow(tw(t),2))*pow(mu(t),2) + pow(mz,4)*(45 + 62*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 
              48*pow(mz,2)*(7 + pow(tw(t),2))*pow(mu(t),6) + 120*pow(mu(t),8)) - 
           2*pow(mw,6)*(4*pow(mz,8)*(12 + pow(tw(t),2)) - pow(mz,6)*(165 + 62*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 
              3*pow(mz,4)*(-107 + 6*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 3*pow(mz,2)*(-19 + 30*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),6) + 24*(2 + pow(tw(t),2))*pow(mu(t),8)))*
         atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))))/
       (pow(mw,8)*pow(pow(mw,2) + pow(mu(t),2),2)*pow(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)),1.5)) + 
      ((pow(mw,8)*(-15 + 14*pow(tw(t),2) + pow(tw(t),4)) + 12*pow(mz,4)*pow(mu(t),4) + 
           pow(mw,6)*(-8*pow(mz,2)*(3 + pow(tw(t),2)) + 2*(-15 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2)) + 
           pow(mw,4)*(12*pow(mz,4) - 16*pow(mz,2)*(3 + pow(tw(t),2))*pow(mu(t),2) + (-27 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4)) + 
           8*pow(mw,2)*(3*pow(mz,4)*pow(mu(t),2) - pow(mz,2)*(3 + pow(tw(t),2))*pow(mu(t),4) - pow(mu(t),6)))*log(pow(mz,2) + pow(mu(t),2)))/
       (2.*pow(mw,8)*pow(pow(mw,2) + pow(mu(t),2),2))) - pow(mu(t),4)*
    ((2*pow(mew,2))/(pow(mw,4)*pow(mu(t),2)) + (4*mew*(-2*pow(mz,2) + pow(mw,2)*(1 + pow(tw(t),2))))/(pow(mw,6)*mu(t)) - 
      (mu(t)*(mew*(-8*pow(mw,12)*pow(-1 + pow(tw(t),2),2) + 4*pow(mz,10)*pow(mu(t),2) + 
              2*pow(mw,10)*(-1 + pow(tw(t),2))*(pow(mz,2)*(7 + 9*pow(tw(t),2)) + (15 + pow(tw(t),2))*pow(mu(t),2)) + 
              4*pow(mw,2)*pow(mz,6)*(pow(mz,4) - pow(mz,2)*(7 + pow(tw(t),2))*pow(mu(t),2) - 5*pow(mu(t),4)) + 
              pow(mw,8)*(-8*pow(mz,4)*(-4 + 7*pow(tw(t),2) + pow(tw(t),4)) + pow(mz,2)*(-17 - 62*pow(tw(t),2) + 15*pow(tw(t),4))*pow(mu(t),2) + 
                 2*(-23 + 2*pow(tw(t),2) + 5*pow(tw(t),4))*pow(mu(t),4)) + pow(mw,4)*
               (-4*pow(mz,8)*(7 + pow(tw(t),2)) + pow(mz,6)*(21 + 30*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 8*pow(mz,4)*(9 + 2*pow(tw(t),2))*pow(mu(t),4) + 
                 20*pow(mz,2)*pow(mu(t),6)) + pow(mw,6)*(pow(mz,6)*(41 + 30*pow(tw(t),2) + pow(tw(t),4)) - 8*pow(mz,4)*(-13 + 5*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 
                 pow(mz,2)*(21 - 58*pow(tw(t),2) - 3*pow(tw(t),4))*pow(mu(t),4) - 8*(2 + pow(tw(t),2))*pow(mu(t),6))) + 
           mu(t)*(8*pow(mw,12)*pow(-1 + pow(tw(t),2),2) + 4*pow(mz,8)*pow(mu(t),4) - 
              2*pow(mw,10)*(-1 + pow(tw(t),2))*(pow(mz,2)*(13 + 3*pow(tw(t),2)) + (23 - 7*pow(tw(t),2))*pow(mu(t),2)) + 
              pow(mw,8)*(pow(mz,4)*(9 + 22*pow(tw(t),2) + pow(tw(t),4)) - 4*pow(mz,2)*(-25 + 6*pow(tw(t),2) + 3*pow(tw(t),4))*pow(mu(t),2) + 
                 4*(23 - 16*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4)) + 4*pow(mw,2)*
               (2*pow(mz,8)*pow(mu(t),2) - pow(mz,6)*(5 + pow(tw(t),2))*pow(mu(t),4) - 4*pow(mz,4)*pow(mu(t),6)) - 
              2*pow(mw,6)*(2*pow(mz,6)*(5 + pow(tw(t),2)) - pow(mz,4)*(-3 + 22*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 
                 pow(mz,2)*(-43 - 4*pow(tw(t),2) + 3*pow(tw(t),4))*pow(mu(t),4) + (-23 + 10*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),6)) + 
              pow(mw,4)*(4*pow(mz,8) - 8*pow(mz,6)*(5 + pow(tw(t),2))*pow(mu(t),2) + pow(mz,4)*(-27 + 22*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 
                 4*pow(mz,2)*(7 + 3*pow(tw(t),2))*pow(mu(t),6) + 8*pow(mu(t),8)))))/
       (pow(mw,8)*(pow(mw,2) + pow(mu(t),2))*(pow(mew,2)*pow(mw,2) + mew*(-2*pow(mw,2) + pow(mz,2))*mu(t) + pow(mu(t),2)*(pow(mw,2) + pow(mu(t),2)))*
         (-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))) - ((16*pow(mw,14)*pow(-1 + pow(tw(t),2),2) + 12*pow(mz,10)*pow(mu(t),4) + 
           4*pow(mw,12)*(-1 + pow(tw(t),2))*(pow(mz,2)*(15 + pow(tw(t),2)) + (3 + 13*pow(tw(t),2))*pow(mu(t),2)) + 
           8*pow(mw,2)*pow(mz,6)*pow(mu(t),2)*(3*pow(mz,4) - pow(mz,2)*(12 + pow(tw(t),2))*pow(mu(t),2) - 10*pow(mu(t),4)) - 
           2*pow(mw,10)*(3*pow(mz,4)*(9 + 22*pow(tw(t),2) + pow(tw(t),4)) + pow(mz,2)*(159 + 34*pow(tw(t),2) - pow(tw(t),4))*pow(mu(t),2) + 
              4*(27 + 4*pow(tw(t),2) - 7*pow(tw(t),4))*pow(mu(t),4)) + pow(mw,8)*
            (pow(mz,6)*(209 + 62*pow(tw(t),2) + pow(tw(t),4)) - 12*pow(mz,4)*(-23 + 18*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) - 
              8*pow(mz,2)*(27 + 38*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 4*(-63 - 18*pow(tw(t),2) + 5*pow(tw(t),4))*pow(mu(t),6)) + 
           pow(mw,4)*pow(mz,2)*(12*pow(mz,8) - 16*pow(mz,6)*(12 + pow(tw(t),2))*pow(mu(t),2) + pow(mz,4)*(45 + 62*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 
              48*pow(mz,2)*(7 + pow(tw(t),2))*pow(mu(t),6) + 120*pow(mu(t),8)) - 
           2*pow(mw,6)*(4*pow(mz,8)*(12 + pow(tw(t),2)) - pow(mz,6)*(165 + 62*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 
              3*pow(mz,4)*(-107 + 6*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 3*pow(mz,2)*(-19 + 30*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),6) + 24*(2 + pow(tw(t),2))*pow(mu(t),8)))*
         atan((-pow(mz,2) + pow(mw,2)*(2 - (2*mew)/mu(t)))/sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))))/
       (pow(mw,8)*pow(pow(mw,2) + pow(mu(t),2),2)*pow(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)),1.5)) + (4*log(mew/mu(t)))/pow(pow(mw,2) + pow(mu(t),2),2) + 
      ((pow(mw,8)*(-15 + 14*pow(tw(t),2) + pow(tw(t),4)) + 12*pow(mz,4)*pow(mu(t),4) + 
           pow(mw,6)*(-8*pow(mz,2)*(3 + pow(tw(t),2)) + 2*(-15 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2)) + 
           pow(mw,4)*(12*pow(mz,4) - 16*pow(mz,2)*(3 + pow(tw(t),2))*pow(mu(t),2) + (-27 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4)) + 
           8*pow(mw,2)*(3*pow(mz,4)*pow(mu(t),2) - pow(mz,2)*(3 + pow(tw(t),2))*pow(mu(t),4) - pow(mu(t),6)))*
         log((pow(mew,2)*pow(mw,2) - 2*mew*pow(mw,2)*mu(t) + mew*pow(mz,2)*mu(t) + pow(mw,2)*pow(mu(t),2) + pow(mu(t),4))/pow(mu(t),2)))/
       (2.*pow(mw,8)*pow(pow(mw,2) + pow(mu(t),2),2))));
	return I;
}

double Ivwlhwluc(double t, double alpha2){
	double I;
	I =-(-0.125*(pow(mu(t),4)*(4*pow(alpha2,2)*pow(mw,4)*pow(pi,2)*pow(v,4) + 
         8*alpha2*pow(mw,2)*pi*pow(v,2)*(alpha2*pow(mw,2)*pi*pow(v,2) + pow(mh,2)*(pow(mw,2) - 2*alpha2*pi*pow(v,2))) - 
         (2*(pow(mh,10)*pow(pow(mw,2) - 2*alpha2*pi*pow(v,2),2) + 
              8*pow(alpha2,2)*pow(mw,4)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mw,4) - 3*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) - 
              pow(mh,8)*(pow(mw,2) - 2*alpha2*pi*pow(v,2))*(4*pow(mw,4) + 2*alpha2*pi*pow(v,2)*pow(mu(t),2) - pow(mw,2)*(12*alpha2*pi*pow(v,2) + pow(mu(t),2))) + 
              pow(mh,6)*(2*pow(mw,8) - 40*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) - pow(mw,6)*(24*alpha2*pi*pow(v,2) + 6*pow(mu(t),2)) + 
                 4*alpha2*pow(mw,4)*pi*pow(v,2)*(11*alpha2*pi*pow(v,2) + 8*pow(mu(t),2))) - 
              2*pow(mh,4)*(pow(mw,6)*pow(mu(t),2)*(-pow(mw,2) + pow(mu(t),2)) + 
                 2*alpha2*pow(mw,4)*pi*pow(v,2)*(-2*pow(mw,4) + 10*pow(mw,2)*pow(mu(t),2) - 3*pow(mu(t),4)) + 
                 2*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*(8*pow(mw,4) - 21*pow(mw,2)*pow(mu(t),2) + 4*pow(mu(t),4))) + 
              4*alpha2*pow(mh,2)*pow(mw,4)*pi*pow(v,2)*(2*pow(mw,4)*pow(mu(t),2) - 4*pow(mw,2)*pow(mu(t),4) + 
                 alpha2*pi*pow(v,2)*(2*pow(mw,4) - 14*pow(mw,2)*pow(mu(t),2) + 12*pow(mu(t),4)))))/
          ((pow(mh,2) + pow(mu(t),2))*(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))) + 
         (2*(pow(mh,10)*(pow(mw,4) - 8*alpha2*pow(mw,2)*pi*pow(v,2) + 12*pow(alpha2,2)*pow(pi,2)*pow(v,4)) + 
              pow(mh,8)*(-6*pow(mw,6) + 56*alpha2*pow(mw,4)*pi*pow(v,2) - 96*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)) + 
              32*pow(alpha2,2)*pow(mw,6)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mw,2) - 3*pow(mu(t),2)) - 
              8*pow(mh,4)*(-(pow(mw,8)*pow(mu(t),2)) - 3*alpha2*pow(mw,6)*pi*pow(v,2)*(pow(mw,2) - 6*pow(mu(t),2)) + 
                 21*pow(alpha2,2)*pow(mw,4)*pow(pi,2)*pow(v,4)*(pow(mw,2) - 2*pow(mu(t),2))) + 
              pow(mh,6)*(6*pow(mw,8) - 80*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) - 6*pow(mw,6)*(16*alpha2*pi*pow(v,2) + pow(mu(t),2)) + 
                 12*alpha2*pow(mw,4)*pi*pow(v,2)*(19*alpha2*pi*pow(v,2) + 4*pow(mu(t),2))) + 
              8*alpha2*pow(mh,2)*pow(mw,4)*pi*pow(v,2)*(4*pow(mw,4)*pow(mu(t),2) - 6*pow(mw,2)*pow(mu(t),4) + 
                 3*alpha2*pi*pow(v,2)*(pow(mw,4) - 11*pow(mw,2)*pow(mu(t),2) + 5*pow(mu(t),4))))*
            atan(pow(mh,2)/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))/
          pow(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5) + 
         (8*alpha2*pow(mh,2)*pow(mw,2)*pi*pow(v,2)*(pow(mw,2) - 3*alpha2*pi*pow(v,2)) + 
            pow(mh,4)*(pow(mw,4) - 8*alpha2*pow(mw,2)*pi*pow(v,2) + 12*pow(alpha2,2)*pow(pi,2)*pow(v,4)) + 
            4*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*(3*pow(mw,2) - 2*pow(mu(t),2)))*log(pow(mh,2) + pow(mu(t),2))))/(pow(mw,8)*pow(pi,2)*pow(v,4)) + 
   (pow(mu(t),4)*((-(pow(mh,10)*pow(pow(mw,2) - 2*alpha2*pi*pow(v,2),2)*(pow(mw,2) + pow(mu(t),2))) + 
           8*pow(alpha2,2)*pow(mw,4)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mw,4)*pow(mu(t),2) + pow(mw,2)*pow(mu(t),4) - pow(mu(t),6)) - 
           pow(mh,8)*(pow(mw,2) - 2*alpha2*pi*pow(v,2))*(pow(mw,2) + pow(mu(t),2))*
            (-3*pow(mw,4) - 2*alpha2*pi*pow(v,2)*pow(mu(t),2) + pow(mw,2)*(10*alpha2*pi*pow(v,2) + pow(mu(t),2))) + 
           pow(mh,6)*pow(mw,2)*(pow(mw,2) + pow(mu(t),2))*(5*pow(mw,4)*pow(mu(t),2) - 4*alpha2*pow(mw,2)*pi*pow(v,2)*(-3*pow(mw,2) + 7*pow(mu(t),2)) + 
              4*pow(alpha2,2)*pow(pi,2)*pow(v,4)*(-7*pow(mw,2) + 9*pow(mu(t),2))) + 
           pow(mh,4)*(2*pow(mw,6)*pow(mu(t),4)*(pow(mw,2) + pow(mu(t),2)) + 
              4*alpha2*pow(mw,4)*pi*pow(v,2)*(pow(mw,2) + pow(mu(t),2))*(5*pow(mw,2)*pow(mu(t),2) - 3*pow(mu(t),4)) + 
              4*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*(3*pow(mw,6) - 10*pow(mw,4)*pow(mu(t),2) - 10*pow(mw,2)*pow(mu(t),4) + 4*pow(mu(t),6))) + 
           4*alpha2*pow(mh,2)*pow(mw,4)*pi*pow(v,2)*(alpha2*pi*pow(v,2)*(5*pow(mw,4)*pow(mu(t),2) - pow(mw,2)*pow(mu(t),4) - 9*pow(mu(t),6)) + 
              2*(pow(mw,4)*pow(mu(t),4) + pow(mw,2)*pow(mu(t),6))))/
         (pow(mw,8)*(pow(mh,2) + pow(mu(t),2))*(pow(mw,2) + pow(mu(t),2))*(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))) + 
        ((pow(mh,10)*(pow(mw,4) - 8*alpha2*pow(mw,2)*pi*pow(v,2) + 12*pow(alpha2,2)*pow(pi,2)*pow(v,4))*pow(pow(mw,2) + pow(mu(t),2),2) - 
             2*pow(mh,8)*pow(mw,2)*(3*pow(mw,4) - 26*alpha2*pow(mw,2)*pi*pow(v,2) + 44*pow(alpha2,2)*pow(pi,2)*pow(v,4))*pow(pow(mw,2) + pow(mu(t),2),2) + 
             16*pow(alpha2,2)*pow(mw,6)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mw,6) - 4*pow(mw,4)*pow(mu(t),2) - 6*pow(mw,2)*pow(mu(t),4) - 3*pow(mu(t),6)) + 
             2*pow(mh,6)*pow(mw,2)*(2*pow(mw,10) - 40*pow(alpha2,2)*pow(pi,2)*pow(v,4)*pow(mu(t),6) + pow(mw,8)*(-36*alpha2*pi*pow(v,2) + pow(mu(t),2)) + 
                8*alpha2*pow(mw,2)*pi*pow(v,2)*pow(mu(t),4)*(alpha2*pi*pow(v,2) + 3*pow(mu(t),2)) + 
                pow(mw,6)*(90*pow(alpha2,2)*pow(pi,2)*pow(v,4) - 48*alpha2*pi*pow(v,2)*pow(mu(t),2) - 4*pow(mu(t),4)) + 
                pow(mw,4)*(136*pow(alpha2,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) + 12*alpha2*pi*pow(v,2)*pow(mu(t),4) - 3*pow(mu(t),6))) + 
             4*pow(mh,4)*pow(mw,4)*(pow(mw,4)*pow(mu(t),2)*pow(pow(mw,2) + pow(mu(t),2),2) + 
                2*alpha2*pow(mw,2)*pi*pow(v,2)*(2*pow(mw,2) - 15*pow(mu(t),2))*pow(pow(mw,2) + pow(mu(t),2),2) + 
                6*pow(alpha2,2)*pow(pi,2)*pow(v,4)*(-5*pow(mw,6) + 4*pow(mw,4)*pow(mu(t),2) + 20*pow(mw,2)*pow(mu(t),4) + 12*pow(mu(t),6))) + 
             8*alpha2*pow(mh,2)*pow(mw,4)*pi*pow(v,2)*(2*pow(mw,2)*pow(mu(t),2)*(pow(mw,2) - 3*pow(mu(t),2))*pow(pow(mw,2) + pow(mu(t),2),2) + 
                alpha2*pi*pow(v,2)*(2*pow(mw,8) - 23*pow(mw,6)*pow(mu(t),2) - 22*pow(mw,4)*pow(mu(t),4) + 12*pow(mw,2)*pow(mu(t),6) + 15*pow(mu(t),8))))*
           atan(pow(mh,2)/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))/
         (pow(mw,8)*pow(pow(mw,2) + pow(mu(t),2),2)*pow(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5)) + 
        ((-4*alpha2*pow(mh,2)*pow(mw,2)*pi*pow(v,2)*(-pow(mw,2) + 4*alpha2*pi*pow(v,2))*pow(pow(mw,2) + pow(mu(t),2),2) + 
             pow(mh,4)*(pow(mw,4) - 8*alpha2*pow(mw,2)*pi*pow(v,2) + 12*pow(alpha2,2)*pow(pi,2)*pow(v,4))*pow(pow(mw,2) + pow(mu(t),2),2) + 
             4*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*(pow(mw,6) + 2*pow(mw,4)*pow(mu(t),2) - 2*pow(mw,2)*pow(mu(t),4) - 2*pow(mu(t),6)))*
           log(pow(mh,2) + pow(mu(t),2)))/(2.*pow(mw,8)*pow(pow(mw,2) + pow(mu(t),2),2))))/(4.*pow(pi,2)*pow(v,4)) + 
   (pow(mu(t),4)*((-2*(8*pow(alpha2,2)*pow(mw,4)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(-pow(mw,4) - pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)) - 
             pow(mh,8)*(pow(mw,2) - 2*alpha2*pi*pow(v,2))*(-pow(mw,4) + 2*alpha2*pi*pow(v,2)*pow(mu(t),2) - pow(mw,2)*(-2*alpha2*pi*pow(v,2) + pow(mu(t),2))) + 
             pow(mh,6)*(-3*pow(mw,8) - 20*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) - pow(mw,6)*(-16*alpha2*pi*pow(v,2) + 3*pow(mu(t),2)) + 
                4*alpha2*pow(mw,4)*pi*pow(v,2)*(-5*alpha2*pi*pow(v,2) + 4*pow(mu(t),2))) - 
             2*pow(mh,4)*(pow(mw,6)*pow(mu(t),2)*(pow(mw,2) + pow(mu(t),2)) + 2*alpha2*pow(mw,4)*pi*pow(v,2)*(3*pow(mw,4) - 3*pow(mu(t),4)) + 
                2*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*(-7*pow(mw,4) - 3*pow(mw,2)*pow(mu(t),2) + 4*pow(mu(t),4))) + 
             4*alpha2*pow(mh,2)*pow(mw,4)*pi*pow(v,2)*(-2*pow(mw,4)*pow(mu(t),2) - 2*pow(mw,2)*pow(mu(t),4) + 
                alpha2*pi*pow(v,2)*(-3*pow(mw,4) + 3*pow(mw,2)*pow(mu(t),2) + 7*pow(mu(t),4)))))/
         ((pow(mw,2) + pow(mu(t),2))*(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))) + 
        (2*(pow(mh,10)*(pow(mw,4) - 8*alpha2*pow(mw,2)*pi*pow(v,2) + 12*pow(alpha2,2)*pow(pi,2)*pow(v,4)) + 
             pow(mh,8)*(-6*pow(mw,6) + 56*alpha2*pow(mw,4)*pi*pow(v,2) - 96*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)) + 
             32*pow(alpha2,2)*pow(mw,6)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mw,2) - 3*pow(mu(t),2)) - 
             8*pow(mh,4)*(-(pow(mw,8)*pow(mu(t),2)) - 3*alpha2*pow(mw,6)*pi*pow(v,2)*(pow(mw,2) - 6*pow(mu(t),2)) + 
                21*pow(alpha2,2)*pow(mw,4)*pow(pi,2)*pow(v,4)*(pow(mw,2) - 2*pow(mu(t),2))) + 
             pow(mh,6)*(6*pow(mw,8) - 80*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) - 6*pow(mw,6)*(16*alpha2*pi*pow(v,2) + pow(mu(t),2)) + 
                12*alpha2*pow(mw,4)*pi*pow(v,2)*(19*alpha2*pi*pow(v,2) + 4*pow(mu(t),2))) + 
             8*alpha2*pow(mh,2)*pow(mw,4)*pi*pow(v,2)*(4*pow(mw,4)*pow(mu(t),2) - 6*pow(mw,2)*pow(mu(t),4) + 
                3*alpha2*pi*pow(v,2)*(pow(mw,4) - 11*pow(mw,2)*pow(mu(t),2) + 5*pow(mu(t),4))))*
           atan((pow(mh,2) - 2*pow(mw,2))/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))/
         pow(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5) + 
        (8*alpha2*pow(mh,2)*pow(mw,2)*pi*pow(v,2)*(pow(mw,2) - 3*alpha2*pi*pow(v,2)) + 
           pow(mh,4)*(pow(mw,4) - 8*alpha2*pow(mw,2)*pi*pow(v,2) + 12*pow(alpha2,2)*pow(pi,2)*pow(v,4)) + 
           4*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*(3*pow(mw,2) - 2*pow(mu(t),2)))*log(pow(mw,2) + pow(mu(t),2))))/(8.*pow(mw,8)*pow(pi,2)*pow(v,4)) - 
   (pow(mu(t),4)*((-4*alpha2*pi*pow(v,2)*(alpha2*pow(mw,2)*pi*pow(v,2) + pow(mh,2)*(pow(mw,2) - 2*alpha2*pi*pow(v,2)))*(1 - mew/mu(t)))/pow(mw,6) + 
        (2*pow(alpha2,2)*pow(pi,2)*pow(v,4)*pow(1 - mew/mu(t),2))/pow(mw,4) + 
        (-((mew*pow(mh,10)*pow(pow(mw,2) - 2*alpha2*pi*pow(v,2),2)*(pow(mw,2) + pow(mu(t),2)))/mu(t)) + 
           8*pow(alpha2,2)*pow(mw,4)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mw,6)*(1 - mew/mu(t)) + mew*pow(mw,4)*mu(t) + mew*pow(mw,2)*pow(mu(t),3) - pow(mu(t),6)) - 
           pow(mh,8)*(pow(mw,2) - 2*alpha2*pi*pow(v,2))*(pow(mw,2) + pow(mu(t),2))*
            (pow(mw,4)*(-3 + 4*(1 - mew/mu(t))) - 2*alpha2*pi*pow(v,2)*pow(mu(t),2) + pow(mw,2)*(-2*alpha2*pi*pow(v,2)*(-5 + 6*(1 - mew/mu(t))) + pow(mu(t),2))) + 
           pow(mh,6)*pow(mw,2)*(pow(mw,2) + pow(mu(t),2))*(2*pow(mw,6)*(1 - mew/mu(t)) + pow(mw,4)*(5 - 3*(1 - mew/mu(t)))*pow(mu(t),2) + 
              4*pow(alpha2,2)*pow(pi,2)*pow(v,4)*(pow(mw,2)*(-7 + 11*(1 - mew/mu(t))) + (9 - 5*(1 - mew/mu(t)))*pow(mu(t),2)) - 
              4*alpha2*pow(mw,2)*pi*pow(v,2)*(pow(mw,2)*(-3 + 6*(1 - mew/mu(t))) + (7 - 4*(1 - mew/mu(t)))*pow(mu(t),2))) + 
           pow(mh,4)*(2*pow(mw,6)*pow(mu(t),2)*(pow(mw,2) + pow(mu(t),2))*(pow(mw,2)*(1 - mew/mu(t)) + pow(mu(t),2)) + 
              4*alpha2*pow(mw,4)*pi*pow(v,2)*(pow(mw,2) + pow(mu(t),2))*(2*pow(mw,4)*(1 - mew/mu(t)) + pow(mw,2)*(5 - 7*(1 - mew/mu(t)))*pow(mu(t),2) - 3*pow(mu(t),4)) + 
              4*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*(pow(mw,6)*(3 - 8*(1 - mew/mu(t))) + 2*pow(mw,4)*(-5 + 3*(1 - mew/mu(t)))*pow(mu(t),2) + 
                 2*pow(mw,2)*(-5 + 7*(1 - mew/mu(t)))*pow(mu(t),4) + 4*pow(mu(t),6))) + 
           4*alpha2*pow(mh,2)*pow(mw,4)*pi*pow(v,2)*(2*(pow(mw,6)*(1 - mew/mu(t))*pow(mu(t),2) + pow(mw,4)*pow(mu(t),4) + mew*pow(mw,2)*pow(mu(t),5)) + 
              alpha2*pi*pow(v,2)*(2*pow(mw,6)*(1 - mew/mu(t)) + pow(mw,4)*(5 - 9*(1 - mew/mu(t)))*pow(mu(t),2) - pow(mw,2)*(1 + 5*(1 - mew/mu(t)))*pow(mu(t),4) + 
                 (-9 + 5*(1 - mew/mu(t)))*pow(mu(t),6))))/
         (pow(mw,8)*(pow(mw,2) + pow(mu(t),2))*(pow(mw,2)*pow(1 - mew/mu(t),2) + (mew*pow(mh,2))/mu(t) + pow(mu(t),2))*
           (-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))) - 
        ((pow(mh,10)*(pow(mw,4) - 8*alpha2*pow(mw,2)*pi*pow(v,2) + 12*pow(alpha2,2)*pow(pi,2)*pow(v,4))*pow(pow(mw,2) + pow(mu(t),2),2) - 
             2*pow(mh,8)*pow(mw,2)*(3*pow(mw,4) - 26*alpha2*pow(mw,2)*pi*pow(v,2) + 44*pow(alpha2,2)*pow(pi,2)*pow(v,4))*pow(pow(mw,2) + pow(mu(t),2),2) + 
             16*pow(alpha2,2)*pow(mw,6)*pow(pi,2)*pow(v,4)*pow(mu(t),2)*(pow(mw,6) - 4*pow(mw,4)*pow(mu(t),2) - 6*pow(mw,2)*pow(mu(t),4) - 3*pow(mu(t),6)) + 
             2*pow(mh,6)*pow(mw,2)*(2*pow(mw,10) - 40*pow(alpha2,2)*pow(pi,2)*pow(v,4)*pow(mu(t),6) + pow(mw,8)*(-36*alpha2*pi*pow(v,2) + pow(mu(t),2)) + 
                8*alpha2*pow(mw,2)*pi*pow(v,2)*pow(mu(t),4)*(alpha2*pi*pow(v,2) + 3*pow(mu(t),2)) + 
                pow(mw,6)*(90*pow(alpha2,2)*pow(pi,2)*pow(v,4) - 48*alpha2*pi*pow(v,2)*pow(mu(t),2) - 4*pow(mu(t),4)) + 
                pow(mw,4)*(136*pow(alpha2,2)*pow(pi,2)*pow(v,4)*pow(mu(t),2) + 12*alpha2*pi*pow(v,2)*pow(mu(t),4) - 3*pow(mu(t),6))) + 
             4*pow(mh,4)*pow(mw,4)*(pow(mw,4)*pow(mu(t),2)*pow(pow(mw,2) + pow(mu(t),2),2) + 
                2*alpha2*pow(mw,2)*pi*pow(v,2)*(2*pow(mw,2) - 15*pow(mu(t),2))*pow(pow(mw,2) + pow(mu(t),2),2) + 
                6*pow(alpha2,2)*pow(pi,2)*pow(v,4)*(-5*pow(mw,6) + 4*pow(mw,4)*pow(mu(t),2) + 20*pow(mw,2)*pow(mu(t),4) + 12*pow(mu(t),6))) + 
             8*alpha2*pow(mh,2)*pow(mw,4)*pi*pow(v,2)*(2*pow(mw,2)*pow(mu(t),2)*(pow(mw,2) - 3*pow(mu(t),2))*pow(pow(mw,2) + pow(mu(t),2),2) + 
                alpha2*pi*pow(v,2)*(2*pow(mw,8) - 23*pow(mw,6)*pow(mu(t),2) - 22*pow(mw,4)*pow(mu(t),4) + 12*pow(mw,2)*pow(mu(t),6) + 15*pow(mu(t),8))))*
           atan((-pow(mh,2) + 2*pow(mw,2)*(1 - mew/mu(t)))/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))/
         (pow(mw,8)*pow(pow(mw,2) + pow(mu(t),2),2)*pow(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5)) + 
        (4*pow(alpha2,2)*pow(pi,2)*pow(v,4)*log(mew/mu(t)))/pow(pow(mw,2) + pow(mu(t),2),2) + 
        ((-4*alpha2*pow(mh,2)*pow(mw,2)*pi*pow(v,2)*(-pow(mw,2) + 4*alpha2*pi*pow(v,2))*pow(pow(mw,2) + pow(mu(t),2),2) + 
             pow(mh,4)*(pow(mw,4) - 8*alpha2*pow(mw,2)*pi*pow(v,2) + 12*pow(alpha2,2)*pow(pi,2)*pow(v,4))*pow(pow(mw,2) + pow(mu(t),2),2) + 
             4*pow(alpha2,2)*pow(mw,2)*pow(pi,2)*pow(v,4)*(pow(mw,6) + 2*pow(mw,4)*pow(mu(t),2) - 2*pow(mw,2)*pow(mu(t),4) - 2*pow(mu(t),6)))*
           log(pow(mw,2)*pow(1 - mew/mu(t),2) + (mew*pow(mh,2))/mu(t) + pow(mu(t),2)))/(2.*pow(mw,8)*pow(pow(mw,2) + pow(mu(t),2),2))))/(4.*pow(pi,2)*pow(v,4)));
	return I;
}

double Ivwlphwtuc(double t){
	double I;
	I =-((pow(mu(t),4)*(3*mu(t)*atan(mw/mu(t)) + mw*(-3 - log(pow(mu(t),2)) + log(pow(mw,2) + pow(mu(t),2)))))/(2.*pow(mw,5)));
	return I;
}

double Ivwlztwtuc(double t){
	double I;
	I =-((2*pow(mu(t),4)*(-((-4*pow(mz,8) - 12*pow(mw,4)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,2)*(4 + pow(tw(t),2)) + 2*pow(mu(t),2)) + 
            2*pow(mw,2)*pow(mz,4)*(pow(mz,2)*(14 + pow(tw(t),2)) + 12*pow(mu(t),2)) + 
            pow(mw,6)*(1 + pow(tw(t),2))*(pow(mz,2)*(9 + pow(tw(t),2)) + 2*(5 + pow(tw(t),2))*pow(mu(t),2)))*
          atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2))))) - 
       (-4*pow(mz,8) - 12*pow(mw,4)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,2)*(4 + pow(tw(t),2)) + 2*pow(mu(t),2)) + 
          2*pow(mw,2)*pow(mz,4)*(pow(mz,2)*(14 + pow(tw(t),2)) + 12*pow(mu(t),2)) + 
          pow(mw,6)*(1 + pow(tw(t),2))*(pow(mz,2)*(9 + pow(tw(t),2)) + 2*(5 + pow(tw(t),2))*pow(mu(t),2)))*
        atan((2*pow(mw,2) - pow(mz,2))/sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))) + 
       sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))*(4*pow(mw,2)*pow(mz,4) + pow(mw,6)*pow(1 + pow(tw(t),2),2) - 
          2*pow(mw,4)*(pow(mz,2)*(7 + pow(tw(t),2)) + 6*pow(mu(t),2)) + 
          (-2*pow(mz,2) + pow(mw,2)*(2 + pow(tw(t),2)))*(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))*
           (log(pow(mw,2) + pow(mu(t),2)) - log(pow(mz,2) + pow(mu(t),2))))))/(pow(mw,6)*pow(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)),1.5)));
	return I;
}

double Ivwltbuc(double t, double alpha2, double alphay){
	double I;
	I =-((pow(mu(t),4)*(2*pow(alpha2,2)*pow(mw,4) + (2*(pow(alphay,2)*pow(mw,4)*(pow(mt,4) - pow(mt,2)*pow(mw,2) - 2*pow(mw,2)*pow(mu(t),2)) + 
            2*alpha2*alphay*pow(mw,2)*(-pow(mt,6) + pow(mt,4)*pow(mw,2) + 3*pow(mt,2)*pow(mw,2)*pow(mu(t),2) + pow(mw,4)*pow(mu(t),2)) + 
            pow(alpha2,2)*(pow(mt,8) - pow(mt,6)*pow(mw,2) - 4*pow(mt,4)*pow(mw,2)*pow(mu(t),2) - pow(mt,2)*pow(mw,4)*pow(mu(t),2) - pow(mw,6)*pow(mu(t),2) + 
               2*pow(mw,4)*pow(mu(t),4))))/(pow(mt,4) - 2*pow(mt,2)*pow(mw,2) + pow(mw,4) - 4*pow(mw,2)*pow(mu(t),2)) - 
       (4*(pow(alphay,2)*pow(mw,4)*(pow(mt,4) - 2*pow(mt,2)*pow(mw,2) + pow(mw,4) - 2*pow(mw,2)*pow(mu(t),2)) + 
            pow(alpha2,2)*(pow(mt,8) - 2*pow(mt,6)*pow(mw,2) + 2*pow(mt,2)*pow(mw,4)*pow(mu(t),2) + 2*pow(mw,4)*pow(mu(t),4) + 
               pow(mt,4)*(pow(mw,4) - 4*pow(mw,2)*pow(mu(t),2))) - 2*alpha2*alphay*pow(mw,2)*
             (pow(mt,6) - 2*pow(mt,4)*pow(mw,2) + pow(mw,4)*pow(mu(t),2) + pow(mt,2)*(pow(mw,4) - 3*pow(mw,2)*pow(mu(t),2)))))/
        (pow(mt,4) - 2*pow(mt,2)*pow(mw,2) + pow(mw,4) - 4*pow(mw,2)*pow(mu(t),2)) + 
       (2*(2*alpha2*alphay*pow(mw,2)*(-pow(mt,2) + pow(mw,2))*(pow(mt,2) + pow(mu(t),2))*
             (pow(mt,4) - 2*pow(mt,2)*pow(mw,2) + pow(mw,4) - 3*pow(mw,2)*pow(mu(t),2)) + 
            pow(alpha2,2)*(pow(mt,2) + pow(mu(t),2))*(pow(mt,8) - 3*pow(mt,6)*pow(mw,2) - pow(mw,6)*pow(mu(t),2) + 2*pow(mw,4)*pow(mu(t),4) + 
               pow(mt,4)*(3*pow(mw,4) - 4*pow(mw,2)*pow(mu(t),2)) - pow(mt,2)*(pow(mw,6) - 5*pow(mw,4)*pow(mu(t),2))) + 
            pow(alphay,2)*pow(mw,4)*(pow(mt,6) + pow(mt,4)*(-3*pow(mw,2) + pow(mu(t),2)) + pow(mt,2)*(3*pow(mw,4) - 5*pow(mw,2)*pow(mu(t),2)) - 
               pow(mw,2)*(pow(mw,4) - 4*pow(mw,2)*pow(mu(t),2) + 2*pow(mu(t),4)))))/
        ((pow(mt,2) + pow(mu(t),2))*(pow(mt,4) - 2*pow(mt,2)*pow(mw,2) + pow(mw,4) - 4*pow(mw,2)*pow(mu(t),2))) + 
       (2*(pow(alphay,2)*pow(mw,4)*(pow(mt,6) - 3*pow(mt,4)*pow(mw,2) - pow(mw,6) + 2*pow(mw,4)*pow(mu(t),2) + 
               3*pow(mt,2)*(pow(mw,4) - 2*pow(mw,2)*pow(mu(t),2))) + 
            pow(alpha2,2)*(3*pow(mt,10) - 9*pow(mt,8)*pow(mw,2) - 2*pow(mw,8)*pow(mu(t),2) + 6*pow(mw,6)*pow(mu(t),4) - 
               6*pow(mt,2)*pow(mw,4)*pow(mu(t),2)*(pow(mw,2) - 5*pow(mu(t),2)) + pow(mt,6)*(9*pow(mw,4) - 20*pow(mw,2)*pow(mu(t),2)) - 
               3*pow(mt,4)*(pow(mw,6) - 8*pow(mw,4)*pow(mu(t),2))) - 
            4*alpha2*alphay*pow(mw,2)*(pow(mt,8) - 3*pow(mt,6)*pow(mw,2) - 2*pow(mw,6)*pow(mu(t),2) + 6*pow(mw,4)*pow(mu(t),4) + 
               3*pow(mt,4)*(pow(mw,4) - 2*pow(mw,2)*pow(mu(t),2)) - pow(mt,2)*(pow(mw,6) - 6*pow(mw,4)*pow(mu(t),2))))*
          atan((-pow(mt,2) - pow(mw,2))/sqrt(-pow(mt,4) + 2*pow(mt,2)*pow(mw,2) - pow(mw,4) + 4*pow(mw,2)*pow(mu(t),2))))/
        pow(-pow(mt,4) + 2*pow(mt,2)*pow(mw,2) - pow(mw,4) + 4*pow(mw,2)*pow(mu(t),2),1.5) - 
       (2*(pow(alphay,2)*pow(mw,4)*(pow(mt,2) - pow(mw,2))*(pow(mt,4) - 2*pow(mt,2)*pow(mw,2) + pow(mw,4) - 6*pow(mw,2)*pow(mu(t),2)) - 
            2*alpha2*alphay*pow(mw,2)*(2*pow(mt,8) - 7*pow(mt,6)*pow(mw,2) + pow(mw,8) - 6*pow(mw,6)*pow(mu(t),2) + 12*pow(mw,4)*pow(mu(t),4) + 
               3*pow(mt,4)*(3*pow(mw,4) - 4*pow(mw,2)*pow(mu(t),2)) + pow(mt,2)*(-5*pow(mw,6) + 18*pow(mw,4)*pow(mu(t),2))) + 
            pow(alpha2,2)*(3*pow(mt,10) - 11*pow(mt,8)*pow(mw,2) + 2*pow(mw,6)*pow(mu(t),2)*(pow(mw,2) - 3*pow(mu(t),2)) + 
               5*pow(mt,6)*(3*pow(mw,4) - 4*pow(mw,2)*pow(mu(t),2)) - 9*pow(mt,4)*(pow(mw,6) - 4*pow(mw,4)*pow(mu(t),2)) + 
               2*pow(mt,2)*(pow(mw,8) - 9*pow(mw,6)*pow(mu(t),2) + 15*pow(mw,4)*pow(mu(t),4))))*
          atan((pow(mt,2) - pow(mw,2))/sqrt(-pow(mt,4) + 2*pow(mt,2)*pow(mw,2) - pow(mw,4) + 4*pow(mw,2)*pow(mu(t),2))))/
        pow(-pow(mt,4) + 2*pow(mt,2)*pow(mw,2) - pow(mw,4) + 4*pow(mw,2)*pow(mu(t),2),1.5) - 
       (2*(pow(alphay,2)*pow(mw,4)*(pow(mt,6) - 3*pow(mt,4)*pow(mw,2) - pow(mw,6) + 2*pow(mw,4)*pow(mu(t),2) + 
               3*pow(mt,2)*(pow(mw,4) - 2*pow(mw,2)*pow(mu(t),2))) + 
            pow(alpha2,2)*(3*pow(mt,10) - 9*pow(mt,8)*pow(mw,2) - 2*pow(mw,8)*pow(mu(t),2) + 6*pow(mw,6)*pow(mu(t),4) - 
               6*pow(mt,2)*pow(mw,4)*pow(mu(t),2)*(pow(mw,2) - 5*pow(mu(t),2)) + pow(mt,6)*(9*pow(mw,4) - 20*pow(mw,2)*pow(mu(t),2)) - 
               3*pow(mt,4)*(pow(mw,6) - 8*pow(mw,4)*pow(mu(t),2))) - 
            4*alpha2*alphay*pow(mw,2)*(pow(mt,8) - 3*pow(mt,6)*pow(mw,2) - 2*pow(mw,6)*pow(mu(t),2) + 6*pow(mw,4)*pow(mu(t),4) + 
               3*pow(mt,4)*(pow(mw,4) - 2*pow(mw,2)*pow(mu(t),2)) - pow(mt,2)*(pow(mw,6) - 6*pow(mw,4)*pow(mu(t),2))))*
          atan((-pow(mt,2) + pow(mw,2))/sqrt(-pow(mt,4) + 2*pow(mt,2)*pow(mw,2) - pow(mw,4) + 4*pow(mw,2)*pow(mu(t),2))))/
        pow(-pow(mt,4) + 2*pow(mt,2)*pow(mw,2) - pow(mw,4) + 4*pow(mw,2)*pow(mu(t),2),1.5) + 
       (2*(pow(alphay,2)*pow(mw,4)*(pow(mt,2) - pow(mw,2))*(pow(mt,4) - 2*pow(mt,2)*pow(mw,2) + pow(mw,4) - 6*pow(mw,2)*pow(mu(t),2)) - 
            2*alpha2*alphay*pow(mw,2)*(2*pow(mt,8) - 7*pow(mt,6)*pow(mw,2) + pow(mw,8) - 6*pow(mw,6)*pow(mu(t),2) + 12*pow(mw,4)*pow(mu(t),4) + 
               3*pow(mt,4)*(3*pow(mw,4) - 4*pow(mw,2)*pow(mu(t),2)) + pow(mt,2)*(-5*pow(mw,6) + 18*pow(mw,4)*pow(mu(t),2))) + 
            pow(alpha2,2)*(3*pow(mt,10) - 11*pow(mt,8)*pow(mw,2) + 2*pow(mw,6)*pow(mu(t),2)*(pow(mw,2) - 3*pow(mu(t),2)) + 
               5*pow(mt,6)*(3*pow(mw,4) - 4*pow(mw,2)*pow(mu(t),2)) - 9*pow(mt,4)*(pow(mw,6) - 4*pow(mw,4)*pow(mu(t),2)) + 
               2*pow(mt,2)*(pow(mw,8) - 9*pow(mw,6)*pow(mu(t),2) + 15*pow(mw,4)*pow(mu(t),4))))*
          atan((pow(mt,2) + pow(mw,2))/sqrt(-pow(mt,4) + 2*pow(mt,2)*pow(mw,2) - pow(mw,4) + 4*pow(mw,2)*pow(mu(t),2))))/
        pow(-pow(mt,4) + 2*pow(mt,2)*pow(mw,2) - pow(mw,4) + 4*pow(mw,2)*pow(mu(t),2),1.5) + 
       (-4*alpha2*alphay*pow(mt,2)*pow(mw,2) + pow(alphay,2)*pow(mw,4) + pow(alpha2,2)*(3*pow(mt,4) - 2*pow(mw,2)*pow(mu(t),2)))*log(pow(mu(t),2)) - 
       (pow(alphay,2)*pow(mw,4) + 2*alpha2*alphay*pow(mw,2)*(-2*pow(mt,2) + pow(mw,2)) + 
          pow(alpha2,2)*(3*pow(mt,4) - 2*pow(mt,2)*pow(mw,2) - 2*pow(mw,2)*pow(mu(t),2)))*log(pow(mu(t),2)) - 
       (-4*alpha2*alphay*pow(mt,2)*pow(mw,2) + pow(alphay,2)*pow(mw,4) + pow(alpha2,2)*(3*pow(mt,4) - 2*pow(mw,2)*pow(mu(t),2)))*log(pow(mt,2) + pow(mu(t),2)) + 
       (pow(alphay,2)*pow(mw,4) + 2*alpha2*alphay*pow(mw,2)*(-2*pow(mt,2) + pow(mw,2)) + 
          pow(alpha2,2)*(3*pow(mt,4) - 2*pow(mt,2)*pow(mw,2) - 2*pow(mw,2)*pow(mu(t),2)))*log(pow(mt,2) + pow(mu(t),2))))/(2.*pow(mw,8)));
	if(t<tt){
    	return 0;
	}
	else{
		return I;
	}
}

double Ivhvlvluc(double t, double alpha, double mv){
	double I;
	I =-((pow(mu(t),4)*((-4*(pow(mv,2) + pow(mu(t),2))*pow(2*pi*pow(v,2)*alpha*(pow(mv,2) + pow(mu(t),2)) + pow(mh,2)*(pow(mv,2) - 2*pi*pow(v,2)*alpha + pow(mu(t),2)),2))/
        (pow(mh,6) - 4*pow(mh,4)*(pow(mv,2) + pow(mu(t),2))) + (2*mu(t)*(pow(mv,2) + pow(mu(t),2))*(mew*pow(mh,2) - pow(mh,2)*mu(t) + 2*mu(t)*(pow(mv,2) + pow(mu(t),2)))*
          pow(2*pi*pow(v,2)*alpha*(pow(mv,2) + pow(mu(t),2)) + pow(mh,2)*(pow(mv,2) - 2*pi*pow(v,2)*alpha + pow(mu(t),2)),2))/
        (pow(mh,4)*(pow(mew,2)*pow(mh,2) - mew*pow(mh,2)*mu(t) + pow(mv,2)*pow(mu(t),2) + pow(mu(t),4))*(pow(mh,2) - 4*(pow(mv,2) + pow(mu(t),2)))) + 
       (2*(-8*pow(pi,2)*pow(v,4)*pow(alpha,2)*pow(pow(mv,2) + pow(mu(t),2),3) - 
            4*pow(mh,2)*pi*pow(v,2)*alpha*pow(pow(mv,2) + pow(mu(t),2),2)*(2*pow(mv,2) + 3*pi*pow(v,2)*alpha + 2*pow(mu(t),2)) + 
            pow(mh,6)*(pow(mv,4) - 4*pow(pi,2)*pow(v,4)*pow(alpha,2) + 2*pow(mv,2)*pow(mu(t),2) + pow(mu(t),4)) - 
            2*pow(mh,4)*(pow(mv,2) + pow(mu(t),2))*(pow(mv,4) - 12*pow(pi,2)*pow(v,4)*pow(alpha,2) + 2*pi*pow(v,2)*alpha*pow(mu(t),2) + pow(mu(t),4) + 
               2*pow(mv,2)*(pi*pow(v,2)*alpha + pow(mu(t),2))))*atan(mh/sqrt(-pow(mh,2) + 4*(pow(mv,2) + pow(mu(t),2)))))/
        (pow(mh,3)*pow(-pow(mh,2) + 4*(pow(mv,2) + pow(mu(t),2)),1.5)) + 
       (2*(-8*pow(pi,2)*pow(v,4)*pow(alpha,2)*pow(pow(mv,2) + pow(mu(t),2),3) - 
            4*pow(mh,2)*pi*pow(v,2)*alpha*pow(pow(mv,2) + pow(mu(t),2),2)*(2*pow(mv,2) + 3*pi*pow(v,2)*alpha + 2*pow(mu(t),2)) + 
            pow(mh,6)*(pow(mv,4) - 4*pow(pi,2)*pow(v,4)*pow(alpha,2) + 2*pow(mv,2)*pow(mu(t),2) + pow(mu(t),4)) - 
            2*pow(mh,4)*(pow(mv,2) + pow(mu(t),2))*(pow(mv,4) - 12*pow(pi,2)*pow(v,4)*pow(alpha,2) + 2*pi*pow(v,2)*alpha*pow(mu(t),2) + pow(mu(t),4) + 
               2*pow(mv,2)*(pi*pow(v,2)*alpha + pow(mu(t),2))))*atan((mh*(-2*mew + mu(t)))/(mu(t)*sqrt(-pow(mh,2) + 4*(pow(mv,2) + pow(mu(t),2))))))/
        (pow(mh,3)*pow(-pow(mh,2) + 4*(pow(mv,2) + pow(mu(t),2)),1.5)) - 8*pow(pi,2)*pow(v,4)*pow(alpha,2)*log(mew/mu(t)) + 
       ((4*pow(mh,2)*pi*pow(v,2)*alpha*pow(pow(mv,2) + pow(mu(t),2),2) + 4*pow(pi,2)*pow(v,4)*pow(alpha,2)*pow(pow(mv,2) + pow(mu(t),2),2) + 
            pow(mh,4)*(pow(mv,4) - 4*pow(pi,2)*pow(v,4)*pow(alpha,2) + 2*pow(mv,2)*pow(mu(t),2) + pow(mu(t),4)))*log(pow(mv,2) + pow(mu(t),2)))/pow(mh,4) - 
       ((4*pow(mh,2)*pi*pow(v,2)*alpha*pow(pow(mv,2) + pow(mu(t),2),2) + 4*pow(pi,2)*pow(v,4)*pow(alpha,2)*pow(pow(mv,2) + pow(mu(t),2),2) + 
            pow(mh,4)*(pow(mv,4) - 4*pow(pi,2)*pow(v,4)*pow(alpha,2) + 2*pow(mv,2)*pow(mu(t),2) + pow(mu(t),4)))*
          log(pow(mv,2) + (mew*pow(mh,2)*(mew - mu(t)))/pow(mu(t),2) + pow(mu(t),2)))/pow(mh,4)))/(8.*pow(pi,2)*pow(v,4)*pow(pow(mv,2) + pow(mu(t),2),2)));
	return I;
}

double Ivhhhuc(double t){
	double I;
	I =-((pow(mu(t),4)*(mh*sqrt(3*pow(mh,2) + 4*pow(mu(t),2)) - 2*(pow(mh,2) + 2*pow(mu(t),2))*atan(mh/sqrt(3*pow(mh,2) + 4*pow(mu(t),2)))))/
   (pow(mh,3)*pow(3*pow(mh,2) + 4*pow(mu(t),2),1.5)));
	return I;
}

double Ivhvtvtuc(double t, double mv){
	double I;
	I = (pow(mu(t),4)*(-(mh*sqrt(-pow(mh,2) + 4*(pow(mv,2) + pow(mu(t),2)))) + (-2*pow(mh,2) + 4*(pow(mv,2) + pow(mu(t),2)))*
        atan(mh/sqrt(-pow(mh,2) + 4*(pow(mv,2) + pow(mu(t),2))))))/(pow(mh,3)*pow(-pow(mh,2) + 4*(pow(mv,2) + pow(mu(t),2)),1.5));
	return I;
}

double Ivhttuc(double t){
	double I;
	I =-((pow(mu(t),4)*(-(mh/(pow(mt,2) + pow(mu(t),2))) + (4*atan(mh/sqrt(-pow(mh,2) + 4*(pow(mt,2) + pow(mu(t),2)))))/sqrt(-pow(mh,2) + 4*(pow(mt,2) + pow(mu(t),2)))))/pow(mh,3));
	return I;
}

double Ivphwwuc(double t){
	double I;
	I =-(0.5 - ((mew - 3*mu(t))*(mew - mu(t)))/(2.*pow(mu(t),2)) - log(mew/mu(t)));
	return I;
}

double Ivztwlwtuc(double t){
	double I;
	I =-((pow(mu(t),4)*((4*pow(mw,2)*pow(1 + pow(tw(t),2),2) - 2*pow(mz,2)*(5 - 2*pow(tw(t),2) + pow(tw(t),4)) + 4*pow(1 + pow(tw(t),2),2)*pow(mu(t),2))/
        (pow(mz,4)*(-4*pow(mw,2) + pow(mz,2) - 4*pow(mu(t),2))) - (8*
          (2*pow(mz,4) + pow(mw,4)*pow(1 + pow(tw(t),2),2) - 2*pow(mz,2)*(3 + pow(tw(t),2))*pow(mu(t),2) + pow(1 + pow(tw(t),2),2)*pow(mu(t),4) + 
            pow(mw,2)*(-2*pow(mz,2)*(3 + pow(tw(t),2)) + 2*pow(1 + pow(tw(t),2),2)*pow(mu(t),2))))/
        (pow(mz,4)*(-4*pow(mw,2) + pow(mz,2) - 4*pow(mu(t),2))*(pow(mw,2) + pow(mu(t),2))) - 
       (2*mu(t)*(mew*pow(mz,2)*(4*pow(mz,4) + pow(mw,4)*(9 + 10*pow(tw(t),2) + pow(tw(t),4)) - 4*pow(mz,2)*(4 + pow(tw(t),2))*pow(mu(t),2) + 
               (9 + 10*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + pow(mw,2)*(-4*pow(mz,2)*(4 + pow(tw(t),2)) + 2*(9 + 10*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2))) + 
            mu(t)*(-4*pow(mz,6) + 2*pow(mw,6)*pow(1 + pow(tw(t),2),2) + 4*pow(mz,4)*(5 + pow(tw(t),2))*pow(mu(t),2) - 
               pow(mz,2)*(21 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 2*pow(1 + pow(tw(t),2),2)*pow(mu(t),6) + 
               pow(mw,4)*(-(pow(mz,2)*(21 + 14*pow(tw(t),2) + pow(tw(t),4))) + 6*pow(1 + pow(tw(t),2),2)*pow(mu(t),2)) + 
               2*pow(mw,2)*(2*pow(mz,4)*(5 + pow(tw(t),2)) - pow(mz,2)*(21 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 3*pow(1 + pow(tw(t),2),2)*pow(mu(t),4)))))/
        (pow(mz,4)*(pow(mw,2) + pow(mu(t),2))*(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2))*(pow(mew,2)*pow(mz,2) - mew*pow(mz,2)*mu(t) + pow(mw,2)*pow(mu(t),2) + pow(mu(t),4)))
         + (4*(-2*pow(mw,2)*(-7 - 6*pow(tw(t),2) + pow(tw(t),4)) + pow(mz,2)*(1 - 6*pow(tw(t),2) + pow(tw(t),4)) - 2*(-7 - 6*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2))*
          atan(mz/sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2))))/(pow(mz,3)*pow(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2),1.5)) - 
       (2*(4*pow(mz,6) + 2*pow(mw,6)*(9 + 10*pow(tw(t),2) + pow(tw(t),4)) - 24*pow(mz,4)*pow(mu(t),2) - pow(mz,2)*(-23 + 2*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 
            2*(9 + 10*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),6) + pow(mw,4)*
             (-(pow(mz,2)*(-23 + 2*pow(tw(t),2) + pow(tw(t),4))) + 6*(9 + 10*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2)) + 
            pow(mw,2)*(-24*pow(mz,4) - 2*pow(mz,2)*(-23 + 2*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 6*(9 + 10*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4)))*
          atan(mz/sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2))))/(pow(mz,3)*pow(pow(mw,2) + pow(mu(t),2),2)*pow(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2),1.5)) - 
       (2*(4*pow(mz,6) + 2*pow(mw,6)*(9 + 10*pow(tw(t),2) + pow(tw(t),4)) - 24*pow(mz,4)*pow(mu(t),2) - pow(mz,2)*(-23 + 2*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4) + 
            2*(9 + 10*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),6) + pow(mw,4)*
             (-(pow(mz,2)*(-23 + 2*pow(tw(t),2) + pow(tw(t),4))) + 6*(9 + 10*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2)) + 
            pow(mw,2)*(-24*pow(mz,4) - 2*pow(mz,2)*(-23 + 2*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 6*(9 + 10*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4)))*
          atan((mz*(-2*mew + mu(t)))/(mu(t)*sqrt(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2)))))/
        (pow(mz,3)*pow(pow(mw,2) + pow(mu(t),2),2)*pow(4*pow(mw,2) - pow(mz,2) + 4*pow(mu(t),2),1.5)) - (8*log(mew/mu(t)))/pow(pow(mw,2) + pow(mu(t),2),2) + 
       ((-4*pow(mz,4) + pow(mw,4)*pow(1 + pow(tw(t),2),2) + 2*pow(mw,2)*pow(1 + pow(tw(t),2),2)*pow(mu(t),2) + pow(1 + pow(tw(t),2),2)*pow(mu(t),4))*
          log(pow(mw,2) + pow(mu(t),2)))/(pow(mz,4)*pow(pow(mw,2) + pow(mu(t),2),2)) - 
       ((-4*pow(mz,4) + pow(mw,4)*pow(1 + pow(tw(t),2),2) + 2*pow(mw,2)*pow(1 + pow(tw(t),2),2)*pow(mu(t),2) + pow(1 + pow(tw(t),2),2)*pow(mu(t),4))*
          log(pow(mw,2) + (mew*pow(mz,2)*(mew - mu(t)))/pow(mu(t),2) + pow(mu(t),2)))/(pow(mz,4)*pow(pow(mw,2) + pow(mu(t),2),2))))/2.);
	return I;
}

double Ivzthztuc(double t){
	double I;
	I = -((2*pow(mu(t),4)*(sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2)) - 
       (pow(mh,2) + 2*pow(mu(t),2))*(atan(pow(mh,2)/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))) - 
          atan((pow(mh,2) - 2*pow(mz,2))/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2))))))/
   pow(-pow(mh,4) + 4*pow(mh,2)*pow(mz,2) + 4*pow(mz,2)*pow(mu(t),2),1.5));
	return I;
}

double Ivztttuc(double t){
	double I;
	I =-((pow(mu(t),4)*((mz*(-18*pow(mt,2) + pow(mz,2)*(9 - 24*pow(sw(t),2) + 32*pow(sw(t),4)) - 18*pow(mu(t),2)))/
        ((pow(mt,2) + pow(mu(t),2))*(4*pow(mt,2) - pow(mz,2) + 4*pow(mu(t),2))) + 
       (8*(9*pow(mt,2) + 4*pow(mz,2)*pow(sw(t),2)*(-3 + 4*pow(sw(t),2)) + 9*pow(mu(t),2))*atan(mz/sqrt(4*pow(mt,2) - pow(mz,2) + 4*pow(mu(t),2))))/
        pow(4*pow(mt,2) - pow(mz,2) + 4*pow(mu(t),2),1.5)))/(36.*pow(mz,3)));
	return I;
}

double Ivwtphwluc(double t){
	double I;
	I =-((pow(mu(t),4)*(-((1 + log(pow(mu(t),2)))/pow(mw,4)) + (pow(mw,2)*pow(mu(t),2) + pow(mu(t),4) + (2*pow(mw,2)*pow(mu(t),2) + pow(mu(t),4))*log(pow(mu(t),2)))/
        (pow(mw,4)*pow(pow(mw,2) + pow(mu(t),2),2)) + (pow(mu(t),2)/(pow(mw,2) + pow(mu(t),2)) + log(pow(mw,2) + pow(mu(t),2)))/pow(mw,4) - 
       ((pow(mu(t),2)*(pow(mw,2) + pow(mu(t),2))*(pow(mw,2)*(-1 + mew/mu(t)) + pow(mu(t),2)))/(pow(mw,4)*((pow(mw,2)*pow(mew - mu(t),2))/pow(mu(t),2) + pow(mu(t),2))) + 
          ((3*pow(mw,2)*mu(t) + pow(mu(t),3))*atan((mw*(-mew + mu(t)))/pow(mu(t),2)))/pow(mw,3) + 2*log(mew/mu(t)) + 
          (pow(mu(t),2)*(2*pow(mw,2) + pow(mu(t),2))*log((pow(mw,2)*pow(mew - mu(t),2))/pow(mu(t),2) + pow(mu(t),2)))/pow(mw,4))/pow(pow(mw,2) + pow(mu(t),2),2)))/2.);
	return I;
}

double Ivwtztwluc(double t){
	double I;
	I =-((pow(mu(t),4)*((-2*(pow(mz,4)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),2) - pow(mw,4)*(1 + pow(tw(t),2))*(pow(mz,2)*(-5 + 3*pow(tw(t),2)) + 2*(-3 + pow(tw(t),2))*pow(mu(t),2)) + 
             pow(mw,2)*(-1 + pow(tw(t),2))*(pow(mz,4)*(-1 + pow(tw(t),2)) - pow(mz,2)*(1 + 3*pow(tw(t),2))*pow(mu(t),2) - 2*(-1 + pow(tw(t),2))*pow(mu(t),4))))/
         ((pow(mw,2) + pow(mu(t),2))*(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))) + 
        (2*(pow(mz,6)*pow(-1 + pow(tw(t),2),2) - 6*pow(mw,2)*pow(mz,2)*pow(-1 + pow(tw(t),2),2)*(pow(mz,2) + pow(mu(t),2)) + 
             2*pow(mw,4)*(1 + pow(tw(t),2))*(pow(mz,2)*(-5 + 3*pow(tw(t),2)) + 4*(-1 + pow(tw(t),2))*pow(mu(t),2)))*
           atan((-2*pow(mw,2) + pow(mz,2))/sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))))/
         pow(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)),1.5) + pow(-1 + pow(tw(t),2),2)*log(pow(mw,2) + pow(mu(t),2))))/(2.*pow(mw,4)) - 
   (pow(mu(t),4)*((2*pow(mz,4)*pow(-1 + pow(tw(t),2),2) + 4*pow(mw,4)*pow(1 + pow(tw(t),2),2) - 
           4*pow(mw,2)*(-1 + pow(tw(t),2))*(2*pow(mz,2)*pow(tw(t),2) + (-1 + pow(tw(t),2))*pow(mu(t),2)))/(pow(mz,4) - 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2))) + 
        (2*(pow(mz,6)*pow(-1 + pow(tw(t),2),2) - 6*pow(mw,2)*pow(mz,2)*pow(-1 + pow(tw(t),2),2)*(pow(mz,2) + pow(mu(t),2)) + 
             2*pow(mw,4)*(1 + pow(tw(t),2))*(pow(mz,2)*(-5 + 3*pow(tw(t),2)) + 4*(-1 + pow(tw(t),2))*pow(mu(t),2)))*
           atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))))/pow(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)),1.5) + 
        pow(-1 + pow(tw(t),2),2)*log(pow(mz,2) + pow(mu(t),2))))/(2.*pow(mw,4)) + 
   (pow(mu(t),4)*((2*(pow(mw,2) + pow(mu(t),2))*(-(pow(mz,4)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),2)) + 
             pow(mw,4)*(1 + pow(tw(t),2))*(pow(mz,2)*(-5 + 3*pow(tw(t),2)) + 2*(-3 + pow(tw(t),2))*pow(mu(t),2)) - 
             pow(mw,2)*(-1 + pow(tw(t),2))*(pow(mz,4)*(-1 + pow(tw(t),2)) - pow(mz,2)*(1 + 3*pow(tw(t),2))*pow(mu(t),2) - 2*(-1 + pow(tw(t),2))*pow(mu(t),4))))/
         (-(pow(mw,4)*pow(mz,4)) + 4*pow(mw,6)*(pow(mz,2) + pow(mu(t),2))) + 
        (2*(pow(mz,6)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),4) + 4*pow(mw,8)*pow(1 + pow(tw(t),2),2)*(pow(mz,2) + pow(mu(t),2)) + 
             2*pow(mw,2)*pow(mz,2)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),2)*(pow(mz,4) - 3*pow(mz,2)*pow(mu(t),2) - 3*pow(mu(t),4)) + 
             pow(mw,6)*(-6*pow(mz,4)*(5 - 2*pow(tw(t),2) + pow(tw(t),4)) + 2*pow(mz,2)*(-35 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 
                8*(-5 + 2*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4)) + pow(mw,4)*
              (pow(mz,6)*(5 - 2*pow(tw(t),2) + pow(tw(t),4)) - 12*pow(mz,4)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),2) - 8*pow(mz,2)*pow(-2 + pow(tw(t),2),2)*pow(mu(t),4) + 
                4*(-3 + 2*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),6)))*atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))))/
         (pow(mw,4)*pow(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)),1.5)) + 
        ((pow(mw,4)*(-3 - 2*pow(tw(t),2) + pow(tw(t),4)) + 2*pow(mw,2)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),2) + pow(-1 + pow(tw(t),2),2)*pow(mu(t),4))*
           log(pow(mz,2) + pow(mu(t),2)))/pow(mw,4)))/(2.*pow(pow(mw,2) + pow(mu(t),2),2)) - 
   (pow(mu(t),4)*((2*mu(t)*(pow(mw,2) + pow(mu(t),2))*(mu(t)*(-(pow(mz,4)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),4)) + 2*pow(mw,6)*pow(1 + pow(tw(t),2),2)*(pow(mz,2) + pow(mu(t),2)) + 
                2*pow(mw,2)*(-1 + pow(tw(t),2))*pow(mu(t),2)*(-(pow(mz,4)*(-1 + pow(tw(t),2))) + pow(mz,2)*(1 + pow(tw(t),2))*pow(mu(t),2) + (-1 + pow(tw(t),2))*pow(mu(t),4)) + 
                pow(mw,4)*(-(pow(mz,4)*(5 - 2*pow(tw(t),2) + pow(tw(t),4))) + 4*pow(mz,2)*(-4 + pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 4*(-3 + pow(tw(t),4))*pow(mu(t),4))) - 
             mew*(pow(mz,6)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),2) + 2*pow(mw,6)*pow(1 + pow(tw(t),2),2)*(pow(mz,2) + pow(mu(t),2)) + 
                pow(mw,2)*pow(mz,2)*(-1 + pow(tw(t),2))*(pow(mz,4)*(-1 + pow(tw(t),2)) - 4*pow(mz,2)*pow(tw(t),2)*pow(mu(t),2) - 3*(-1 + pow(tw(t),2))*pow(mu(t),4)) + 
                pow(mw,4)*(-4*pow(mz,4)*pow(tw(t),2)*(-1 + pow(tw(t),2)) - pow(mz,2)*(5 - 10*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 
                   2*(-3 + 2*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4)))))/
         (pow(mw,4)*(pow(mew,2)*pow(mw,2) + mew*(-2*pow(mw,2) + pow(mz,2))*mu(t) + pow(mu(t),2)*(pow(mw,2) + pow(mu(t),2)))*
           (-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))) - 
        (2*(pow(mz,6)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),4) + 4*pow(mw,8)*pow(1 + pow(tw(t),2),2)*(pow(mz,2) + pow(mu(t),2)) + 
             2*pow(mw,2)*pow(mz,2)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),2)*(pow(mz,4) - 3*pow(mz,2)*pow(mu(t),2) - 3*pow(mu(t),4)) + 
             pow(mw,6)*(-6*pow(mz,4)*(5 - 2*pow(tw(t),2) + pow(tw(t),4)) + 2*pow(mz,2)*(-35 + 14*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),2) + 
                8*(-5 + 2*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),4)) + pow(mw,4)*
              (pow(mz,6)*(5 - 2*pow(tw(t),2) + pow(tw(t),4)) - 12*pow(mz,4)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),2) - 8*pow(mz,2)*pow(-2 + pow(tw(t),2),2)*pow(mu(t),4) + 
                4*(-3 + 2*pow(tw(t),2) + pow(tw(t),4))*pow(mu(t),6)))*atan((-pow(mz,2) + pow(mw,2)*(2 - (2*mew)/mu(t)))/
             sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))))/(pow(mw,4)*pow(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)),1.5)) + 8*log(mew/mu(t)) + 
        ((pow(mw,4)*(-3 - 2*pow(tw(t),2) + pow(tw(t),4)) + 2*pow(mw,2)*pow(-1 + pow(tw(t),2),2)*pow(mu(t),2) + pow(-1 + pow(tw(t),2),2)*pow(mu(t),4))*
           log((pow(mew,2)*pow(mw,2) - 2*mew*pow(mw,2)*mu(t) + mew*pow(mz,2)*mu(t) + pow(mw,2)*pow(mu(t),2) + pow(mu(t),4))/pow(mu(t),2)))/pow(mw,4)))/
    (2.*pow(pow(mw,2) + pow(mu(t),2),2)));
	return I;
}

double Ivwtwtzluc(double t){
	double I;
	I = -((-2*pow(mu(t),4)*((-pow(mz,6) - 4*pow(mz,4)*pow(mu(t),2) + pow(mz,2)*pow(mu(t),4) + 2*pow(mu(t),6) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,2) + 4*pow(mu(t),2)))*
        atan(pow(mz,2)/sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))) + 
       (-pow(mz,6) - 4*pow(mz,4)*pow(mu(t),2) + pow(mz,2)*pow(mu(t),4) + 2*pow(mu(t),6) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2))*(pow(mz,2) + 4*pow(mu(t),2)))*
        atan((2*pow(mw,2) - pow(mz,2))/sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))) + 
       sqrt(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))*(-pow(mz,4) - pow(mu(t),4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)) + 
          (-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)))*(log(pow(mw,2) + pow(mu(t),2)) - log(pow(mz,2) + pow(mu(t),2))))))/
   (pow(pow(mz,2) + pow(mu(t),2),2)*pow(-pow(mz,4) + 4*pow(mw,2)*(pow(mz,2) + pow(mu(t),2)),1.5)));
	return I;
}

double Ivwthwtuc(double t){
	double I;
	I =-((2*pow(mu(t),4)*(sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2)) - 
       (pow(mh,2) + 2*pow(mu(t),2))*(atan(pow(mh,2)/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))) - 
          atan((pow(mh,2) - 2*pow(mw,2))/sqrt(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2))))))/
   pow(-pow(mh,4) + 4*pow(mh,2)*pow(mw,2) + 4*pow(mw,2)*pow(mu(t),2),1.5));
	return I;
}

double Ivwtbtuc(double t){
	double I;
	I =-((pow(mu(t),4)*((2*pow(mw,4)*(pow(mt,2) - pow(mw,2) + 2*pow(mu(t),2)))/((pow(mt,2) + pow(mu(t),2))*(pow(pow(mt,2) - pow(mw,2),2) - 4*pow(mw,2)*pow(mu(t),2))) - 
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

