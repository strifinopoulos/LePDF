/*ultra-collinear splitting functions. The names are Puc followed by the particles B, A and C. The factors 4pi^2 are there because the alphas are defined divided by 2pi*/

/*two fermions and a transverse vector. p and m are the vector transverse polarizations, l and r the helicities of the fermions*/
/*A and B = fermions and C = transverse vector*/
double Puctrblwmm(double z, double alpha2, double alphay){   
	double P;
	P = v*v*0.5*alphay*alpha2*(1-z); /*4Pi*Pi due to the definition of alphas and a 1/2 from Qfw^2, total 2Pi*Pi*/
	return P;
}

double Pucwmmbltr(double z, double alpha2, double alphay){   
	double P;
	P = v*v*0.5*alphay*alpha2*z; 
	return P;
}

double Pucbltrwpp(double z, double alpha2, double alphay){   
	double P;
	P = v*v*0.5*alphay*alpha2*(1-z)*z*z; 
	return P;
}

double Pucwpptrbl(double z, double alpha2, double alphay){   
	double P;
	P = v*v*0.5*alphay*alpha2*(1-z)*(1-z)*z; 
	return P;
}

double Puctrtlzm(double z, double alpha2, double alphay, double t){  
	double P;
	P = v*v*alpha2*alphay*(1-z)*pow(0.5-Qu*sw(t)*sw(t)*(1-z),2)/pow(cw(t),2);
	return P;
}

double Puczmtltr(double z, double alpha2, double alphay, double t){  
	double P;
	P = Puctrtlzm(1-z,alpha2,alphay,t);
	return P;
}

double Puctltrzp(double z, double alpha2, double alphay, double t){  
	double P;
	P = v*v*alpha2*alphay*(1-z)*pow(0.5*z+Qu*sw(t)*sw(t)*(1-z),2)/pow(cw(t),2);
	return P;
}

double Puczptrtl(double z, double alpha2, double alphay, double t){  
	double P;
	P = Puctltrzp(1-z,alpha2,alphay,t);
	return P;
}

double Pucttph(double z, double alphaem, double alphay){  /*for the photon the LRM and RLP splittings are the same*/
	double P;
	P = v*v*Qu*Qu*alphaem*alphay*pow(1-z,3);
	return P;
}

double Pucphtt(double z, double alphaem, double alphay){  
	double P;
	P = Pucttph(1-z,alphaem,alphay);
	return P;
}

double Puczphmtltr(double z, double alpham, double alphay, double t){  
	double P;
	P = Qu*v*v*2*alpham*alphay*z*z*(0.5-Qu*sw(t)*sw(t)*z)/cw(t);
	return P;
}

double Puczphptrtl(double z, double alpham, double alphay, double t){  
	double P;
	P = -Qu*v*v*2*alpham*alphay*z*z*(0.5*(1-z)+Qu*sw(t)*sw(t)*z)/cw(t);
	return P;
}

double Pucttg(double z, double alpha3, double alphay){  /*for the gluon the LRM and RLP splittings are the same*/
	double P;
	P = v*v*Cf*alpha3*alphay*pow(1-z,3);
	return P;
}

double Pucgtt(double z, double alpha3, double alphay){  
	double P;
	P = Pucttg(1-z,alpha3,alphay);
	return P;
}

/*A = transverse vector, B = fermion, C = antifermion and B<->C*/
double Pucblwmmtrb(double z, double alpha2, double alphay){ 
	double P;
	P = 3*v*v*0.5*alpha2*alphay*z*z;
	return P;
}

double Puctrwppblb(double z, double alpha2, double alphay){ 
	double P;
	P = 3*v*v*0.5*alpha2*alphay*(1-z)*(1-z);
	return P;
}

double Puctlzmtrb(double z, double alpha2, double alphay, double t){ 
	double P;
	P = 3*v*v*alpha2*alphay*pow(0.5*z-Qu*sw(t)*sw(t),2)/pow(cw(t),2);
	return P;
}

double Puctrzptlb(double z, double alpha2, double alphay, double t){ 
	double P;
	P = 3*v*v*alpha2*alphay*pow(0.5*(1-z)-Qu*sw(t)*sw(t),2)/pow(cw(t),2);
	return P;
}

double Puctpht(double z, double alphaem, double alphay){ 
	double P;
	P = 3*v*v*Qu*Qu*alphaem*alphay;
	return P;
}

double Puctlzphmtrb(double z, double alpham, double alphay, double t){  
	double P;
	P = 3*Qu*v*v*alpham*alphay*(0.5*z-Qu*sw(t)*sw(t))/cw(t);
	return P;
}

double Puctrzphptlb(double z, double alpham, double alphay, double t){  
	double P;
	P = 3*Qu*v*v*alpham*alphay*(0.5*(1-z)-Qu*sw(t)*sw(t))/cw(t);
	return P;
}

double Puctgt(double z, double alpha3, double alphay){ 
	double P;
	P = v*v*0.5*alpha3*alphay;
	return P;
}

/*two fermions and a longitudinal vector*/
/*A and B = fermions and C = longitudinal vector*/ 
double Pucfl2fl1wl(double z, double alpha2, double alphay1, double alphay2){   
	double P;
	P = v*v*0.5*pow(alphay1*z*(1-z)-alphay2*(1-z)-alpha2*z,2)/(1-z);
	return P;
}

double Pucwlfl1fl2(double z, double alpha2, double alphay1, double alphay2){  
	double P;
	P = Pucfl2fl1wl(1-z, alpha2, alphay1, alphay2);
	return P;
}

double Pucflflzl(double z, double alpha2, double alphay, double t, double T3, double Q){  
	double P;
	P = v*v*pow(T3*alphay*(1-z)*(1-z)+alpha2*z*(T3-Q*sw(t)*sw(t))/pow(cw(t),2),2)/(1-z);
	return P;
}

double Puczlflfl(double z, double alpha2, double alphay, double t, double T3, double Q){  
	double P;
	P = Pucflflzl(1-z, alpha2, alphay, t, T3, Q);
	return P;
}

double Pucfrfrzl(double z, double alpha2, double alphay, double t, double T3, double Q){  /*T3 = 0 for all except the top, for which it is 0.5*/
	double P;
	P = v*v*pow(T3*alphay*(1-z)*(1-z)+alpha2*z*Q*sw(t)*sw(t)/pow(cw(t),2),2)/(1-z);
	return P;
}

double Puczlfrfr(double z, double alpha2, double alphay, double t, double T3, double Q){  
	double P;
	P = Pucfrfrzl(1-z, alpha2, alphay, t, T3, Q);
	return P;
}

/*A = longitudinal vector, B = fermion, C = antifermion*/
double Pucfl1wlflb2(double z, double alpha2, double alphay1, double alphay2, double Ncf){ 
	double P;
	P = Ncf*0.5*v*v*pow(alphay1*(1-z)+alphay2*z-alpha2*z*(1-z),2);
	return P;
}

double Pucflzlflb(double z, double alpha2, double alphay, double t, double T3, double Q, double Ncf){ 
	double P;
	P = Ncf*v*v*pow(T3*alphay-(T3-Q*sw(t)*sw(t))*z*(1-z)*alpha2/pow(cw(t),2),2);
	return P;
}

double Pucfrzlfrb(double z, double alpha2, double alphay, double t, double T3, double Q, double Ncf){ 
	double P;
	P = Ncf*v*v*pow(T3*alphay+Q*sw(t)*sw(t)*z*(1-z)*alpha2/pow(cw(t),2),2);
	return P;
}

/*two tops and h*/
double Puctth(double z, double alphay){ /*A and B = fermion, C = h*/
	double P;
	P = 0.25*v*v*(1-z)*(1+z)*(1+z)*alphay*alphay;
	return P;
}

double Puchtt(double z, double alphay){ /*A and C = fermion, B = h*/
	double P;
	P = Puctth(1-z,alphay);
	return P;
}

double Puctht(double z, double alphay){ /*A = h, B and C = (anti)fermions*/
	double P;
	P = 3*0.25*v*v*alphay*alphay*(1-2*z)*(1-2*z);
	return P;
}

/*two transverse vectors and a longitudinal one: cases ABC = TTL and ABC = TLT. The transverse vectors have the SAME POLARIZATION*/
double Pucphwtwl(double z, double alphaem, double alpha2){ /*A = WT, B = photon, C = WL*/
	double P;
	P = alphaem*alpha2*z*z*z*v*v/(1-z);
	return P;
}

double Pucwlwtph(double z, double alphaem, double alpha2){ /*A = WT, C = photon, B = WL*/
	double P;
	P = Pucphwtwl(1-z,alphaem,alpha2);
	return P;
}

double Pucztwtwl(double z, double alpha2, double t){ /*A = WT, B = ZT, C = WL*/
	double P;
	P = cw(t)*cw(t)*0.25*alpha2*alpha2*v*v*z*pow(1+z+tw(t)*tw(t)*(1-z),2)/(1-z);
	return P;
}

double Pucwlwtzt(double z, double alpha2, double t){ /*A = WT, C = ZT, B = WL*/
	double P;
	P = Pucztwtwl(1-z,alpha2,t);
	return P;
}

double Puczphwtwl(double z, double alpham, double alpha2, double t){ /*A = WT, B = ZT, C = WL*/
	double P;
	P = cw(t)*alpham*alpha2*v*v*z*z*(1+z+tw(t)*tw(t)*(1-z))/(1-z);
	return P;
}

double Pucwtwtzl(double z, double alpha2){ /*A,B = WT, C = ZL*/
	double P;
	P = 0.25*alpha2*alpha2*v*v*z*(1+z)*(1+z)/(1-z);
	return P;
}

double Puczlwtwt(double z, double alpha2){ /*A,C = WT, B = ZL*/
	double P;
	P = Pucwtwtzl(1-z,alpha2);
	return P;
}

double Pucwtphwl(double z, double alphaem, double alpha2){ /*A = photon, B = WT and C = WL*/
	double P;
	P = alphaem*alpha2*v*v*z/(1-z);
	return P;
}

double Pucwlphwt(double z, double alphaem, double alpha2){ /*A = photon, C = WT and B = WL*/
	double P;
	P = Pucwtphwl(1-z,alphaem,alpha2);
	return P;
}

double Pucwtztwl(double z, double alpha2, double t){ /*A = ZT, B = WT and C = WL*/
	double P;
	P = cw(t)*cw(t)*0.25*alpha2*alpha2*v*v*z*pow(1+z-tw(t)*tw(t)*(1-z),2)/(1-z);
	return P;
}

double Pucwlztwt(double z, double alpha2, double t){ /*A = ZT, C = WT and B = WL*/
	double P;
	P = Pucwtztwl(1-z,alpha2,t);
	return P;
}

double Pucwtzphwl(double z, double alpham, double alpha2, double t){ /*A = Z-photon, B = WT and C = WL*/
	double P;
	P = cw(t)*0.5*alpham*alpha2*v*v*z*(1+z-tw(t)*tw(t)*(1-z))/(1-z);
	return P;
}

double Pucwlzphwt(double z, double alpham, double alpha2, double t){ /*A = Z-photon, C = WT and B = WL*/
	double P;
	P = Pucwtzphwl(1-z,alpham,alpha2,t);
	return P;
}

/*Longitudinal vector splitting into transverse, that is A =VL, B=C=VT. The vectors have OPPOSITE POLARIZATIONS*/
double Pucphwlwt(double z, double alphaem, double alpha2){ /*A = WL, B = photon and C = WT*/
	double P;
	P = alphaem*alpha2*v*v*z*z*z*(1-z);
	return P;
}

double Pucwtwlph(double z, double alphaem, double alpha2){ /*A = WL, C = photon and B = WT*/
	double P;
	P = Pucphwlwt(1-z,alphaem,alpha2);
	return P;
}

double Pucztwlwt(double z, double alpha2, double t){ /*A = WL, B = ZT and C = WT*/
	double P;
	P = cw(t)*cw(t)*0.25*alpha2*alpha2*v*v*z*(1-z)*pow(1-2*z+tw(t)*tw(t),2);
	return P;
}

double Pucwtwlzt(double z, double alpha2, double t){ /*A = WL, C = ZT and B = WT*/
	double P;
	P = Pucztwlwt(1-z,alpha2,t);
	return P;
}

double Puczphwlwt(double z, double alpham, double alpha2, double t){ /*A = WL, B = ZT and C = WT*/
	double P;
	P = -cw(t)*alpham*alpha2*v*v*z*z*(1-z)*(1-2*z+tw(t)*tw(t));
	return P;
}

double Pucwtzlwt(double z, double alpha2){ /*A = ZL, B and C = WT*/
	double P;
	P = 0.25*alpha2*alpha2*v*v*z*(1-z)*(1-2*z)*(1-2*z);
	return P;
}

/*two transverse vectors and the Higgs. Due to the simmetry in z<->1-z, ABC = VT VT h and ABC = VT h VT are the same. The vectors have the SAME POLARIZATION*/
double Pucwtwth(double z, double alpha2){ /*A = WT, B,C = WT and h*/
	double P;
	P = 0.25*alpha2*alpha2*v*v*z*(1-z);
	return P;
}

double Puchwtwt(double z, double alpha2){ /*A = h, B,C = WT*/
	double P;
	P = 0.25*alpha2*alpha2*v*v*z*(1-z);
	return P;
}

double Pucztzth(double z, double alpha2, double t){ /*A = ZT, B,C = ZT and h*/
	double P;
	P = 0.25*alpha2*alpha2*v*v*z*(1-z)/(pow(cw(t),4));
	return P;
}

double Puchztzt(double z, double alpha2, double t){ /*A = h, B,C = ZT*/
	double P;
	P = 0.25*alpha2*alpha2*v*v*z*(1-z)/(pow(cw(t),4));
	return P;
}

/*Higgs splitting into transverse vectors, that is A = h. The vectors have OPPOSITE POLARIZATIONS*/
double Puczthzt(double z, double alpha2, double t){ /*A = h, B and C = ZT*/
	double P;
	P = 0.125*alpha2*alpha2*v*v*z*(1-z)/(pow(cw(t),4));
	return P;
}

double Pucwthwt(double z, double alpha2){ /*A = h, B and C = WT*/
	double P;
	P = 0.25*alpha2*alpha2*v*v*z*(1-z);
	return P;
}

/*Three longitudinal gauge bosons*/
double Pucwlzlwl(double z, double alpha2, double t){ /*A = ZL, B and C = WL*/
	double P;
	P = 0.0625*alpha2*alpha2*v*v*(1-2*z)*(1-2*z)*pow(2+(z-z*z)*(1-tw(t)*tw(t)),2)/(z*(1-z));
	return P;
}

double Puczlwlwl(double z, double alpha2, double t){ /*A = WL, B = ZL and C = WL*/
	double P;
	P = 0.0625*alpha2*alpha2*v*v*pow((1-2*z)*(2+z-z*z)-tw(t)*tw(t)*(1-z)*(2-z),2)/(z*(1-z));
	return P;
}

double Pucwlwlzl(double z, double alpha2, double t){ /*A = WL, C = ZL and B = WL*/
	double P;
	P = Puczlwlwl(1-z,alpha2,t);
	return P;
}

/*Higgs into longitudinal gauge bosons*/
double Puczlhzl(double z, double alpha2, double t){ /*A = h, B and C = ZL*/
	double P;
	P = v*v*pow(alpha2*(1-z+z*z)/(cw(t)*cw(t))-mh*mh*z*(1-z)/(2*pi*v*v),2)/(8*z*(1-z));
	return P;
}

double Pucwlhwl(double z, double alpha2){ /*A = h, B and C = WL*/
	double P;
	P = v*v*pow(alpha2*(1-z+z*z)-mh*mh*z*(1-z)/(2*pi*v*v),2)/(4*z*(1-z));
	return P;
}

/*Longitudinal vector in itself and Higgs*/
double Puchwlwl(double z, double alpha2){ /*A = WL, B = h and C = WL*/
	double P;
	P = v*v*z*pow(alpha2*(1-z+z*z)+mh*mh*(1-z)/(2*pi*v*v),2)/(4*(1-z));
	return P;
}

double Pucwlwlh(double z, double alpha2){ /*A = WL, C = h and B = WL*/
	double P;
	P = Puchwlwl(1-z,alpha2);
	return P;
}

double Puchzlzl(double z, double alpha2, double t){ /*A = ZL, B = h and C = ZL*/
	double P;
	P = v*v*z*pow(alpha2*(1-z+z*z)/(cw(t)*cw(t))+mh*mh*(1-z)/(2*pi*v*v),2)/(4*(1-z));
	return P;
}

double Puczlzlh(double z, double alpha2, double t){ /*A = ZL, C = h and B = ZL*/
	double P;
	P = Puchzlzl(1-z,alpha2,t);
	return P;
}

/*Three Higgses*/
double Puchhh(double z){ 
	double P;
	P = 9*pow(mh,4)*z*(1-z)/(32*pi*pi*v*v);
	return P;
}

/*U.C. splitting functions involving the mixed h-Zl pdf as A. These ones are for B = Wp: notice that both functions are antisymmetric in z<->1-z, so for Wm the usual exchange brings just a minus sign*/
double Pucwthzlwt(double z, double alpha2){
	double P;
	P = -v*v*0.25*alpha2*alpha2*(1-z)*z*(1-2*z);
	return P;
}

double Pucwlhzlwl(double z, double alpha2, double t){
	double P;
	P = v*v*0.125*alpha2*(alpha2*(1-z+z*z)-mh*mh*z*(1-z)/(2*pi*v*v))*(1-2*z)*(2+(1-tw(t)*tw(t))*(z-z*z))/(z*(1-z));
	return P;
}

double Puctlhzltlb(double z, double alpha2, double alphay, double t){
	double P;
	P = -v*v*0.5*alphay*(1-2*z)*(0.5*alphay-alpha2*(0.5-Qu*sw(t)*sw(t))*z*(1-z)/pow(cw(t),2));
	return P;
}

double Puctrhzltrb(double z, double alpha2, double alphay, double t){
	double P;
	P = v*v*0.5*alphay*(1-2*z)*(0.5*alphay+alpha2*Qu*sw(t)*sw(t)*z*(1-z)/pow(cw(t),2));
	return P;
}

/*U.C. splitting functions involving the mixed h-Zl pdf as B.*/
double Puchzlwlwl(double z, double alpha2, double t){
	double P;
	P = v*v*0.25*alpha2*(alpha2*(1-z+z*z)+mh*mh*(1-z)/(2*pi*v*v))*((1-2*z)*(2+z-z*z)-tw(t)*tw(t)*(1-z)*(2-z))/(1-z);
	return P;
}

double Puchzlwtwt(double z, double alpha2){
	double P;
	P = v*v*0.5*alpha2*alpha2*(1-z)*(2-z);
	return P;
}

double Puchzltltl(double z, double alpha2, double alphay, double t){
	double P;
	P = v*v*alphay*(2-z)*(0.5*alphay*z*z+alpha2*(0.5-Qu*sw(t)*sw(t))*(1-z)/pow(cw(t),2));
	return P;
}

double Puchzltrtr(double z, double alpha2, double alphay, double t){
	double P;
	P = -v*v*alphay*(1-2*z)*(0.5*alphay*z*z+alpha2*Qu*sw(t)*sw(t)*(1-z)/pow(cw(t),2));
	return P;
}

/*Regularized u.c. splitting functions: they are the divergent ones, multiplied by zb*/

double zbPucflflzl(double z, double alpha2, double alphay, double t, double T3, double Q){
	double P;
	P = v*v*pow(T3*alphay*(1-z)*(1-z)+alpha2*z*(T3-Q*sw(t)*sw(t))/pow(cw(t),2),2);
	return P;
}

double zbPucfrfrzl(double z, double alpha2, double alphay, double t, double T3, double Q){
	double P;
	P = v*v*pow(T3*alphay*(1-z)*(1-z)+alpha2*z*Q*sw(t)*sw(t)/pow(cw(t),2),2);
	return P;
}

double zbPucfl2fl1wl(double z, double alpha2, double alphay1, double alphay2){
	double P;
	P = v*v*0.5*pow(alphay1*z*(1-z)-alphay2*(1-z)-alpha2*z,2);
	return P;
}

double zbPucphwtwl(double z, double alphaem, double alpha2){
	double P;
	P = alphaem*alpha2*z*z*z*v*v;
	return P;
}

double zbPucztwtwl(double z, double alpha2, double t){
	double P;
	P = cw(t)*cw(t)*0.25*alpha2*alpha2*v*v*z*pow(1+z+tw(t)*tw(t)*(1-z),2);
	return P;
}

double zbPuczphwtwl(double z, double alpham, double alpha2, double t){
	double P;
	P = cw(t)*alpham*alpha2*v*v*z*z*(1+z+tw(t)*tw(t)*(1-z));
	return P;
}

double zbPucwtwtzl(double z, double alpha2){
	double P;
	P = 0.25*alpha2*alpha2*v*v*z*(1+z)*(1+z);
	return P;
}

double zbPucwtphwl(double z, double alphaem, double alpha2){
	double P;
	P = alphaem*alpha2*v*v*z;
	return P;
}

double zbPucwtztwl(double z, double alpha2, double t){
	double P;
	P = cw(t)*cw(t)*0.25*alpha2*alpha2*v*v*z*pow(1+z-tw(t)*tw(t)*(1-z),2);
	return P;
}

double zbPucwtzphwl(double z, double alpham, double alpha2, double t){
	double P;
	P = cw(t)*0.5*alpham*alpha2*v*v*z*(1+z-tw(t)*tw(t)*(1-z));
	return P;
}

double zbPuczlwlwl(double z, double alpha2, double t){
	double P;
	P = 0.0625*alpha2*alpha2*v*v*pow((1-2*z)*(2+z-z*z)-tw(t)*tw(t)*(1-z)*(2-z),2)/z;
	return P;
}

double zbPucwlwlzl(double z, double alpha2, double t){
	double P;
	P = 0.0625*alpha2*alpha2*v*v*pow((2*z-1)*(2+z-z*z)-tw(t)*tw(t)*z*(1+z),2)/z;
	return P;
}

double zbPucwlzlwl(double z, double alpha2, double t){
	double P;
	P = 0.0625*alpha2*alpha2*v*v*(1-2*z)*(1-2*z)*pow(2+(z-z*z)*(1-tw(t)*tw(t)),2)/z;
	return P;
}

double zbPuczlhzl(double z, double alpha2, double t){
	double P;
	P = v*v*pow(alpha2*(1-z+z*z)/(cw(t)*cw(t))-mh*mh*z*(1-z)/(2*pi*v*v),2)/(8*z);
	return P;
}

double zbPucwlhwl(double z, double alpha2){
	double P;
	P = v*v*pow(alpha2*(1-z+z*z)-mh*mh*z*(1-z)/(2*pi*v*v),2)/(4*z);
	return P;
}

double zbPuchzlzl(double z, double alpha2, double t){
	double P;
	P = v*v*z*pow(alpha2*(1-z+z*z)/(cw(t)*cw(t))+mh*mh*(1-z)/(2*pi*v*v),2)/4;
	return P;
}

double zbPuchwlwl(double z, double alpha2){
	double P;
	P = v*v*z*pow(alpha2*(1-z+z*z)+mh*mh*(1-z)/(2*pi*v*v),2)/4;
	return P;
}

double zbPuchzlwlwl(double z, double alpha2, double t){
	double P;
	P = v*v*0.25*alpha2*(alpha2*(1-z+z*z)+mh*mh*(1-z)/(2*pi*v*v))*((1-2*z)*(2+z-z*z)-tw(t)*tw(t)*(1-z)*(2-z));
	return P;
}

double zbPucwlhzlwl(double z, double alpha2, double t){
	double P;
	P = v*v*0.125*alpha2*(alpha2*(1-z+z*z)-mh*mh*z*(1-z)/(2*pi*v*v))*(1-2*z)*(2+(1-tw(t)*tw(t))*(z-z*z))/z;
	return P;
}
