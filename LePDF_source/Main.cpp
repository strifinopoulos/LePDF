/* -------------------------------------------------------------- 
						LePDF v.1.1

						20/02/2024
						
	Francesco Garosi,	David Marzocca, Sokratis Trifinopoulos
	
	Reference: arXiv:2303.16964
	
  -------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>


/* ----------------------- Input Parameters ----------------------- */

/* CHOOSE THE LEPTON TYPE */
int lepton = 1; /* Defines the type of lepton beam: 0 = electron, any other integer (for instance 1) = muon */
/*for electron:  int lepton=0; */

/* CHOOSE THE FLAVOUR SCHEME */
int FS = 5; /* Defines the flavour scheme: 5 = without top, any other integer = with top */

/* CHOOSE THE NUMBER OF x-GRID POINTS */
int Nx = 1001; /* Number of x-grid points. Notice that in our paper we defined Nx to be the number of bins, so here Nx is just shifted by 1 */

/* CHOOSE HOW MANY OF THE x-GRID POINTS SAVE TO FILE */
/* The default is to save one every 5 x-steps, plus all the last 100 x-steps to describe better the muon PDF peak close to x=1 */
int ixStep = 5;
int NixLast = 100;

/* ----------------------- Definitions ----------------------- */

/* Definition of the constants used in the evolution as global variables */

double Cf = double(4)/3; /*SU(3) quadratic Casimir in the fundamental representation*/
double Ca = 3; /*SU(3) quadratic Casimir in the adjoint representation*/
double Tf = 0.5; /*SU(3) Dynkin index in the fundamental representation*/
double Qu = double(2)/3; /*Up quark electric charge*/
double Qd = -double(1)/3; /*Down quark electric charge*/
double Ql = -1; /*Leptons electric charge*/

double v = 246.22; /* Higgs vev */
double mmu = 0.1056584; /* Muon mass */
double me = 0.000510998; /* Electron mass */
double mz = 91.1876; /* Z mass */
double mw = 80.379; /* W mass */
double mew = 80.379; /* Matching scale at which we include the full EW evolution. We choose the W mass */
double mh = 125.25; /* Higgs mass */
double mt = 163; /* Top quark mass */
double mb = 4.18; /* b quark mass (only for threshold) */
double mtau = 1.78; /* tau mass (only for threshold) */
double mc = 1.29; /* c quark mass (only for threshold) */


double pi = 3.1415926536;


/* ----------------------- Preliminary setup ----------------------- */

/*Function to fill the grid in x*/
double xgrid(int i){
	double expi;
	double xi;
	expi = -6*pow(double(Nx-1-i)/(Nx-1),2.5);
	xi = pow(10,expi);
	return xi; 
}

double leptonmass(int l){
	double m;
	if(l==0){
		m = me;
	}
	else{
		m = mmu;
	}
	return m;
}

double m0 = leptonmass(lepton); /*Mass of the valence lepton. It is set to the electron/muon mass for the computation of the electron/muon PDFs. */

/*Beta function coefficients for U(1)_em according to the c, tau and b thresholds*/
double b0QED0 = 8/(3*pi); /*below the c threshold*/
double b0QEDc = 32/(9*pi); /*between the c and tau thresholds*/
double b0QEDtau = 38/(9*pi); /*between the tau and b thresholds*/
double b0QEDb = 40/(9*pi); /*above the b threshold*/

/*Beta function coefficients for U(1)_Y and SU(2)_L*/
double b01 = 41/(12*pi);
double b02 = 19/(12*pi);

/*Value of electromagnetic coupling at the W mass, used as initial condition for the running in the QED+QCD phase*/
double aem0 = 0.00783581;

/*Values of the U(1)_Y, SU(2)_L and Yukawa couplings at 200 GeV, used as initial conditions for the running in the EW phase*/
double a10 = 0.0102475;
double a20 = 0.0332946;
double ay0 = 0.0679086;

double tFromGeV(double Q){ /*Function to get t values from GeV scales*/
	double tt;
	if(Q>0){
		tt = 2*log(Q/m0);
	}
	else{
		tt=0;
		return EXIT_FAILURE;
	}
	return tt;
}

double mu(double t){ /*Inverse relation from t to the scale Q in GeV*/
	double m;
	m = m0*exp(0.5*t);
	return m;
}

/* ----------------------- t discretization ----------------------- */

int Nt0choice(int l){
	int N;
	if(lepton==0){
		N = 351;
	}
	else{
		N = 201;
	}
	return N;
}

int Ntchoice(int l){
	int N;
	if(lepton==0){
		N = 150;
	}
	else{
		N = 200;
	}
	return N;
}

int Nt0 = Nt0choice(lepton); /*Number of t values at which we compute the PDFs in the first phase of the evolution, including t=0*/
int Nt = Ntchoice(lepton); /*Number of t values at which we compute the PDFs in the second phase of the evolution. This means that in the second phase there are Nt-1 integration steps*/
double t0 = tFromGeV(mew); /*t at the matching scale. It is also the t at which we define the initial condition for the electromagnetic running coupling*/

/*Amplitude of the t intervals used for the evolution in the two phases. We consider in both phases uniform intervals*/
double dt0 = t0/(Nt0-1); /*In the first phase we take as amplitude the EW scale divided by the number of integration steps, t0/(Nt0-1)*/
double dt = dt0; /*In the second phase dt can be different from dt0, but we choose them to be equal*/
/*Notice that in this way the last PDF is computed at t = t0+(Nt-1)*dt = t0*(1+(Nt-1)/N0). */

/*Other relevant energy scales and thresholds (in t)*/
double t0r = tFromGeV(200); /*t at which we define the initial conditions for U(1)_Y and SU(2)_L running couplings, here 200 GeV*/  
double tW = tFromGeV(mw); /*t0 = 2log(Mw/Mmu), i.e. t at W mass*/
double tz = tFromGeV(mz); /*tz = 2log(Mz/Mmu), i.e. t at Z mass*/
double th = tFromGeV(mh); /*th = 2log(Mh/Mmu), i.e. t at H mass*/
double tb = tFromGeV(mb); /*tb = 2log(Mb/Mmu), i.e. t at b mass*/
double ttau = tFromGeV(mtau); /*ttau = 2log(Mtau/Mmu), i.e. t at tau mass*/
double tc = tFromGeV(mc); /*tb = 2log(Mc/Mmu), i.e. t at c mass*/
double tQCD = tFromGeV(0.7); /*t at which QCD starts contributing*/

double ttchoice(int s){
	double ttop;
	if(s==5){
		ttop = 100; /*In the 5FS the top is absent: we keep the same code, shifting the top threshold to a value above the maximum t reached during the evolution*/
	}
	else{
		ttop = tFromGeV(mt); /*tt = 2log(Mt/Mmu), i.e. t at top mass*/
	}
	return ttop;
}

double tt = ttchoice(FS); /* t at top mass */


/* Variables to read alpha_s from file */
int Na = 246; /*Number of alpha_s points in the alpha_s file*/
double dta = 0.1; /*dt spacing in the alpha_s file*/

double twrite = tFromGeV(10); /*t at which we start writing the PDFs in the LHAPDF file*/

/* ----------------------- Load Code ----------------------- */

#include "Matrix_alloc.cpp" /*To allocate and free matrices*/
#include "Couplings_RG.cpp" /*Runnings couplings for both phases*/
#include "Cuts.cpp" /*Cuts due to z_max (upper bound of integration in the DGLAP equations) and to the mass thresholds*/
#include "EVA.cpp" /*PDFs in EVA used for the matching at the EW scale*/
#include "Prop_corr.cpp" /*Corrections to the splitting functions coming from propagators due to the masses of the particles*/
#include "Massless_splitting.cpp" /*Splitting functions non vanishing with massless particles*/
#include "QED_QCD.cpp" /*Evolution equations in the QED-QCD phase*/
#include "Uc_splitting.cpp"
#include "Massless_DGLAP.cpp" /*Contributions to the DGLAP equations from the "massless" splittings*/
#include "Massless_eqs.cpp" /*DGLAP equations from the "massless" splittings*/
#include "Uc_eqs.cpp" /*DGLAP equations from the ultra-collinear splittings*/
#include "RK_sum.cpp" /*k coefficients of the Runge-Kutta algorithm (z-contributions)*/
#include "Massless_virt.cpp" /*Formulas for the virtual contributions P^v (from massless terms)*/
#include "Uc_virtual.cpp" /*Formulas for the virtual contributions P^v (from ultra-collinear terms)*/
#include "RK_virtual.cpp" /*k coefficients of the Runge-Kutta algorithm (x-contributions from massless terms)*/
#include "RK_virtual_uc.cpp" /*k coefficients of the Runge-Kutta algorithm (x-contributions from ultra-collinear terms)*/
#include "Writing.cpp" /*Contains the functions to write the results in LHAPDF format*/


/* ------------------------ Main Program ------------------------ */

int main(){

	if(lepton==0){
		printf("Computation of electron PDFs: m0 = %f, t0 = %f, dt0 = %f, dt = %f, Nx = %d, N0 = %d, Nt = %d, FS=%d.\n",m0,t0,dt0,dt,Nx,Nt0,Nt,FS);
	}
	else{
		printf("Computation of muon PDFs: m0 = %f, t0 = %f, dt0 = %f, dt = %f, Nx = %d, N0 = %d, Nt = %d, FS=%d.\n",m0,t0,dt0,dt,Nx,Nt0,Nt,FS);
	}
	
	int i,j; /*Variables to be used in the for cycles*/
	int iw; /*Index corresponding to the time at which we start writing the PDFs in the LHAPDF file. It is the integer part of twrite/dt0*/
	iw = int(twrite/dt0);
	int N0 = 9; /*number of PDFs in QED+QCD*/
	int N = 54; /*number of PDFs in the full EW evolution*/
	int Nttot = Nt0+Nt-1; /*total number of values of t at which we compute the PDFs (the -1 is due to the fact that the last value of the first phase and the first of the second are the same*/
	

	/* ----------------------- Preparing t and Q grids ----------------------- */
	/*Grids in the energy scale, both in t and Q. They contain the Nttot values at which we compute the PDFs*/
	double *tgrid, *qgrid;
	tgrid = (double*)malloc(Nttot*sizeof(double));
	qgrid = (double*)malloc(Nttot*sizeof(double));
	for (i=0;i<Nt0-1;i++){ /*For the t grid, we take equidistant points, with distance dt0 and dt for the two phases respectively*/
		tgrid[i] = i*dt0;
	}
	for (i=0;i<Nt;i++){
		tgrid[i+Nt0-1] = t0+i*dt;
	}
	for (i=0;i<Nttot;i++){ /*Then the Q-grid is computed from the one in t*/
		qgrid[i] = m0*exp(tgrid[i]/2);
	}
	printf("t and Q grids are ready!\n");
	
	/* ----------------------- Memory Allocation ----------------------- */
	
	/*Declaration and allocation of the matrices with the PDFs. The first line in the following is for QED+QCD, the second one for the full SM. We initially set all the entries to zero*/
	double **fm,**fu,**fc,**fd,**fph,**fg,**fe,**ftau,**fb;
	double **fwpp,**fwpm,**fwmp,**fwmm,**fwpl,**fwml,**fzp,**fzm,**fzl,**fphp,**fphm,**fzphp,**fzphm,**fgp,**fgm,**fh,**fhzl,**fel,**felb,**fml,**ftal,**ftalb,**fer,**ferb,**fmr,**ftar,**ftarb,**fnue,**fnueb,**fnum,**fnuta,**fnutab,**ful,**fulb,**fcl,**fclb,**ftl,**ftlb,**fdl,**fdlb,**fsl,**fslb,**fbl,**fblb,**fur,**furb,**fcr,**fcrb,**ftr,**ftrb,**fdr,**fdrb,**fbr,**fbrb;
	fm =zeros(Nt0,Nx);
	fu =zeros(Nt0,Nx);
	fc =zeros(Nt0,Nx);
	fd =zeros(Nt0,Nx);
	fph =zeros(Nt0,Nx);
	fg =zeros(Nt0,Nx);
	fe =zeros(Nt0,Nx);
	ftau =zeros(Nt0,Nx);
	fb =zeros(Nt0,Nx);
	fwpp =zeros(Nt,Nx);
	fwpm =zeros(Nt,Nx);
	fwmp =zeros(Nt,Nx);
	fwmm =zeros(Nt,Nx);
	fwpl =zeros(Nt,Nx);
	fwml =zeros(Nt,Nx);
	fzp =zeros(Nt,Nx);
	fzm =zeros(Nt,Nx);
	fzl =zeros(Nt,Nx);
	fphp =zeros(Nt,Nx);
	fphm =zeros(Nt,Nx);
	fzphp =zeros(Nt,Nx);
	fzphm =zeros(Nt,Nx);
	fhzl =zeros(Nt,Nx);
	fgp =zeros(Nt,Nx);
	fgm =zeros(Nt,Nx);
	fh =zeros(Nt,Nx);
	fel =zeros(Nt,Nx);
	felb =zeros(Nt,Nx);
	fml =zeros(Nt,Nx);
	ftal =zeros(Nt,Nx);
	ftalb =zeros(Nt,Nx);
	fer =zeros(Nt,Nx);
	ferb =zeros(Nt,Nx);
	fmr =zeros(Nt,Nx);
	ftar =zeros(Nt,Nx);
	ftarb =zeros(Nt,Nx);
	fnue =zeros(Nt,Nx);
	fnum =zeros(Nt,Nx);
	fnuta =zeros(Nt,Nx);
	fnueb =zeros(Nt,Nx);
	fnutab =zeros(Nt,Nx);
	ful =zeros(Nt,Nx);
	fulb =zeros(Nt,Nx);
	fur =zeros(Nt,Nx);
	furb =zeros(Nt,Nx);
	fcl =zeros(Nt,Nx);
	fclb =zeros(Nt,Nx);
	fcr =zeros(Nt,Nx);
	fcrb =zeros(Nt,Nx);
	ftl =zeros(Nt,Nx);
	ftlb =zeros(Nt,Nx);
	ftr =zeros(Nt,Nx);
	ftrb =zeros(Nt,Nx);
	fdl =zeros(Nt,Nx);
	fdlb =zeros(Nt,Nx);
	fdr =zeros(Nt,Nx);
	fdrb=zeros(Nt,Nx);
	fsl =zeros(Nt,Nx);
	fslb =zeros(Nt,Nx);
	fbl =zeros(Nt,Nx);
	fblb =zeros(Nt,Nx);
	fbr =zeros(Nt,Nx);
	fbrb =zeros(Nt,Nx);
	printf("Matrices for PDFs allocated!\n");
	
	/*Declaration and allocation of the matrices with the Runge-Kutta coefficients matrices k1, k2, k3 and k4. The additional 0 refers to the first phase. As before we initially set all the entries to zero*/
	double **k01,**k02,**k03,**k04,**k1,**k2,**k3,**k4;
	k01 =zeros(Nt0-1,N0*Nx);
	k02 =zeros(Nt0-1,N0*Nx);
	k03 =zeros(Nt0-1,N0*Nx);
	k04 =zeros(Nt0-1,N0*Nx);
	k1 =zeros(Nt-1,N*Nx);
	k2 =zeros(Nt-1,N*Nx);
	k3 =zeros(Nt-1,N*Nx);
	k4 =zeros(Nt-1,N*Nx);
	printf("Matrices for Runge-Kutta coefficients allocated!\n");
	
	/* ----------------------- Opening the output files ----------------------- */
	
	FILE *ia,*lhapdf,*lhapdfbar; 
	ia = fopen ("alphas_NNLO.dat", "r"); /*This is to read the file with alpha_s and to write the results in LHAPDF format*/
	if(ia==NULL){
		perror("Error in opening the file with alpha_s");
		return EXIT_FAILURE;
	}
	/*Now open the files with the output, depending on the choice of the valence lepton and the flavour scheme*/
	if(lepton==0){
		if(FS==5){
			lhapdf = fopen ("LePDF_e_5FS_0000.dat", "w"); /*In this file we save the PDFs for the electron 5FS in LHAPDF format*/
			lhapdfbar = fopen ("LePDF_eb_5FS_0000.dat", "w"); /*In this file we save the PDFs for the positron 5FS in LHAPDF format*/
		}
		else{
			lhapdf = fopen ("LePDF_e_6FS_0000.dat", "w"); /*In this file we save the PDFs for the electron 6FS in LHAPDF format*/
			lhapdfbar = fopen ("LePDF_eb_6FS_0000.dat", "w"); /*In this file we save the PDFs for the positron 6FS in LHAPDF format*/
		}
	}
	else{
		if(FS==5){
			lhapdf = fopen ("LePDF_mu_5FS_0000.dat", "w"); /*In this file we save the PDFs for the muon 5FS in LHAPDF format*/
			lhapdfbar = fopen ("LePDF_mub_5FS_0000.dat", "w"); /*In this file we save the PDFs for the antimuon 5FS in LHAPDF format*/
		}
		else{
			lhapdf = fopen ("LePDF_mu_6FS_0000.dat", "w"); /*In this file we save the PDFs for the muon 6FS in LHAPDF format*/
			lhapdfbar = fopen ("LePDF_mub_6FS_0000.dat", "w"); /*In this file we save the PDFs for the antimuon 6FS in LHAPDF format*/
		}
	}
	
	if(lhapdf==NULL){
		perror("Error in opening the file with LHAPDF");
		return EXIT_FAILURE;
	}
	if(lhapdfbar==NULL){
		perror("Error in opening the file with LHAPDF for the antiparticle");
		return EXIT_FAILURE;
	}
	
	/* ----------------------- Loading alpha strong ----------------------- */
	/*Loading the file with alpha_s and reporting it in the matrix talpha*/
	double **talpha;
	talpha = zeros(Na,2);
	for(i=0;i<Na;i++){
		for(j=0;j<2;j++){
		fscanf(ia,"%lf", &talpha[i][j]); /* The file loaded has {Q, alpha} */
		}
		talpha[i][0] = tFromGeV(talpha[i][0]); /* Here we change it to {t, alpha} */
	}
	printf("alpha_s ready\n"); /*Notice: the function alphas in "Running couplings.cpp" reads this matrix and gives the value of alpha_s at any scale through inear interpolation*/
	
	/* ----------------------- Preparing x-Grid ----------------------- */
	/*Declaration of the arrays containing the grid in x, the corresponding spacings and the indices of the entries to be written*/
	double *grid,*deltagrid;
	int *rgrid;
	
	int NixSteps = ((Nx-1) - NixLast) / ixStep;
	int lr = NixSteps + NixLast + 1; /*Number of elements of the grid we want to write in the output files. It is the length of rgrid*/
	grid = (double*)malloc(Nx*sizeof(double));
	deltagrid = (double*)malloc(Nx*sizeof(double));
	rgrid = (int*)malloc(lr*sizeof(int));
	for(i=0; i <= NixSteps; i++){
		rgrid[i] = ixStep*i;
	}
	for(i=0; i < NixLast; i++){
		rgrid[NixSteps + 1 +i] = Nx - NixLast + i;
	}
	
	/*We directly use the formula reported in appendix E of our paper to fill the grid. Of course the formula can be changed or the grid can be read directly from a file: we leave both options, commenting the latter*/
	/*
	FILE *ig;
	ig = fopen ("FILENAME", "r");
	if(ig==NULL){
		perror("Error in opening the file with grid");
		return EXIT_FAILURE;
	}
	for(i=0;i<Nx;i++){
		fscanf(ig,"%lf", &grid[i]);
	}
	*/
	for(i=0;i<Nx;i++){
		grid[i] = xgrid(i);
	}
	printf("The x-grid is ready!\n");
	deltagrid[0] = grid[0]; 
	for(i=1;i<Nx;i++){
		deltagrid[i] = grid[i]-grid[i-1];
	}
	printf("x spacings computed!\n");
	
	/* ----------------------- Setting initial conditions ----------------------- */
	/*Now we start the first phase of the evolution*/
	fm[0][Nx-1] = 1/deltagrid[Nx-1]; /*The initial conditions at t=0 correspond to fmu(x) = delta(1-x), while all the other PDFs vanish: then we just need to change the last entry of fmu*/
	printf("Initial condition set\n"); 
	
	/* ----------------------- QED+QCD DGLAP evolution ----------------------- */
	int k; /*Variable used in the for cycle to sum over the x_k contributions to the integrals*/
	printf("Starting the evolution for QED+QCD!\n");
	for(i=0;i<Nt0-1;i++){ /*Cycle over t*/
		printf("QED+QCD evolution: step %d of %d\n",i+1,Nt0-1);
		for(j=0;j<Nx-1;j++){ /*Cycle over x: we compute the four Runge-Kutta coefficients corresponding to each f(x_j,t_i). RKj0 contains the contributions from x_j, RK0 from x_k*/
			k01 = RKj0(i,j,Nx,i*dt0,dt0,grid,deltagrid,fm[i][j],fu[i][j],fd[i][j],fph[i][j],fg[i][j],fe[i][j],fb[i][j],fc[i][j],ftau[i][j],k01,newterm(grid,deltagrid,j,Nx),talpha);
			for(k=j+1;k<Nx;k++){
				k01 = RK0(i,j,k,Nx,i*dt0,dt0,grid,deltagrid,fm[i][k],fu[i][k],fd[i][k],fph[i][k],fg[i][k],fe[i][k],fb[i][k],fc[i][k],ftau[i][k],k01,talpha);
			}
		}
		for(j=0;j<Nx-1;j++){
			k02 = RKj0(i,j,Nx,(i+0.5)*dt0,dt0,grid,deltagrid,fm[i][j]+0.5*k01[i][j],fu[i][j]+0.5*k01[i][j+Nx],fd[i][j]+0.5*k01[i][j+2*Nx],fph[i][j]+0.5*k01[i][j+3*Nx],fg[i][j]+0.5*k01[i][j+4*Nx],fe[i][j]+0.5*k01[i][j+5*Nx],fb[i][j]+0.5*k01[i][j+6*Nx],fc[i][j]+0.5*k01[i][j+7*Nx],ftau[i][j]+0.5*k01[i][j+8*Nx],k02,newterm(grid,deltagrid,j,Nx),talpha);
			for(k=j+1;k<Nx;k++){
				k02 = RK0(i,j,k,Nx,(i+0.5)*dt0,dt0,grid,deltagrid,fm[i][k]+0.5*k01[i][k],fu[i][k]+0.5*k01[i][k+Nx],fd[i][k]+0.5*k01[i][k+2*Nx],fph[i][k]+0.5*k01[i][k+3*Nx],fg[i][k]+0.5*k01[i][k+4*Nx],fe[i][k]+0.5*k01[i][k+5*Nx],fb[i][k]+0.5*k01[i][k+6*Nx],fc[i][k]+0.5*k01[i][k+7*Nx],ftau[i][k]+0.5*k01[i][k+8*Nx],k02,talpha);
			}
		}
		for(j=0;j<Nx-1;j++){
			k03 = RKj0(i,j,Nx,(i+0.5)*dt0,dt0,grid,deltagrid,fm[i][j]+0.5*k02[i][j],fu[i][j]+0.5*k02[i][j+Nx],fd[i][j]+0.5*k02[i][j+2*Nx],fph[i][j]+0.5*k02[i][j+3*Nx],fg[i][j]+0.5*k02[i][j+4*Nx],fe[i][j]+0.5*k02[i][j+5*Nx],fb[i][j]+0.5*k02[i][j+6*Nx],fc[i][j]+0.5*k02[i][j+7*Nx],ftau[i][j]+0.5*k02[i][j+8*Nx],k03,newterm(grid,deltagrid,j,Nx),talpha);
			for(k=j+1;k<Nx;k++){
				k03 = RK0(i,j,k,Nx,(i+0.5)*dt0,dt0,grid,deltagrid,fm[i][k]+0.5*k02[i][k],fu[i][k]+0.5*k02[i][k+Nx],fd[i][k]+0.5*k02[i][k+2*Nx],fph[i][k]+0.5*k02[i][k+3*Nx],fg[i][k]+0.5*k02[i][k+4*Nx],fe[i][k]+0.5*k02[i][k+5*Nx],fb[i][k]+0.5*k02[i][k+6*Nx],fc[i][k]+0.5*k02[i][k+7*Nx],ftau[i][k]+0.5*k02[i][k+8*Nx],k03,talpha);
			}
		}
		for(j=0;j<Nx-1;j++){
			k04 = RKj0(i,j,Nx,(i+1)*dt0,dt0,grid,deltagrid,fm[i][j]+k03[i][j],fu[i][j]+k03[i][j+Nx],fd[i][j]+k03[i][j+2*Nx],fph[i][j]+k03[i][j+3*Nx],fg[i][j]+k03[i][j+4*Nx],fe[i][j]+k03[i][j+5*Nx],fb[i][j]+k03[i][j+6*Nx],fc[i][j]+k03[i][j+7*Nx],ftau[i][j]+k03[i][j+8*Nx],k04,newterm(grid,deltagrid,j,Nx),talpha);
			for(k=j+1;k<Nx;k++){
				k04 = RK0(i,j,k,Nx,(i+1)*dt0,dt0,grid,deltagrid,fm[i][k]+k03[i][k],fu[i][k]+k03[i][k+Nx],fd[i][k]+k03[i][k+2*Nx],fph[i][k]+k03[i][k+3*Nx],fg[i][k]+k03[i][k+4*Nx],fe[i][k]+k03[i][k+5*Nx],fb[i][k]+k03[i][k+6*Nx],fc[i][k]+k03[i][k+7*Nx],ftau[i][k]+k03[i][k+8*Nx],k04,talpha);
			}
		}
		for(j=0;j<Nx-1;j++){
			fm[i+1][j] = fm[i][j] + k01[i][j]/6 + k02[i][j]/3 + k03[i][j]/3 + k04[i][j]/6;
			fu[i+1][j] = fu[i][j] + k01[i][j+Nx]/6 + k02[i][j+Nx]/3 + k03[i][j+Nx]/3 + k04[i][j+Nx]/6;
			fd[i+1][j] = fd[i][j] + k01[i][j+2*Nx]/6 + k02[i][j+2*Nx]/3 + k03[i][j+2*Nx]/3 + k04[i][j+2*Nx]/6;
			fph[i+1][j] = fph[i][j] + k01[i][j+3*Nx]/6 + k02[i][j+3*Nx]/3 + k03[i][j+3*Nx]/3 + k04[i][j+3*Nx]/6;
			fg[i+1][j] = fg[i][j] + k01[i][j+4*Nx]/6 + k02[i][j+4*Nx]/3 + k03[i][j+4*Nx]/3 + k04[i][j+4*Nx]/6;
			fe[i+1][j] = fe[i][j] + k01[i][j+5*Nx]/6 + k02[i][j+5*Nx]/3 + k03[i][j+5*Nx]/3 + k04[i][j+5*Nx]/6;
			fb[i+1][j] = fb[i][j] + k01[i][j+6*Nx]/6 + k02[i][j+6*Nx]/3 + k03[i][j+6*Nx]/3 + k04[i][j+6*Nx]/6;
			fc[i+1][j] = fc[i][j] + k01[i][j+7*Nx]/6 + k02[i][j+7*Nx]/3 + k03[i][j+7*Nx]/3 + k04[i][j+7*Nx]/6;
			ftau[i+1][j] = ftau[i][j] + k01[i][j+8*Nx]/6 + k02[i][j+8*Nx]/3 + k03[i][j+8*Nx]/3 + k04[i][j+8*Nx]/6;
		}
		fm[i+1][Nx-1] = Lt0(grid,deltagrid,Nx,i+1,fe,fm,ftau,fu,fc,fd,fph,fg,fb)/deltagrid[Nx-1];
	}
	printf("QED+QCD evolution completed!\n");
	
	/* ----------------------- Matching at the EW scale ----------------------- */
	/*Here we start with the second phase of the evolution filling the first row of the PDFs matrices by matching with the first phase: the PDFs equally split between the two chiralities (fermions) or the two transverse polarizations (gauge bosons)*/
	for(i=0;i<Nx;i++){
		fel[0][i] = fe[Nt0-1][i]/2;
		felb[0][i] = fel[0][i];
		fml[0][i] = fm[Nt0-1][i]/2;
		ftal[0][i] = ftau[Nt0-1][i]/2;
		ftalb[0][i] = ftal[0][i];
		fer[0][i] = fel[0][i];
		ferb[0][i] = fel[0][i];
		ftar[0][i] = ftal[0][i];
		ftarb[0][i] = ftal[0][i];
		fmr[0][i] = fml[0][i];
		ful[0][i] = fu[Nt0-1][i]/2;
		fulb[0][i] = ful[0][i];
		fur[0][i] = ful[0][i];
		furb[0][i] = ful[0][i];
		fcl[0][i] = fc[Nt0-1][i]/2;
		fclb[0][i] = fcl[0][i];
		fcr[0][i] = fcl[0][i];
		fcrb[0][i] = fcl[0][i];
		fdl[0][i] = fd[Nt0-1][i]/2;
		fdlb[0][i] = fdl[0][i];
		fdr[0][i] = fdl[0][i];
		fdrb[0][i] = fdl[0][i];
		fsl[0][i] = fdl[0][i];
		fslb[0][i] = fdl[0][i];
		fbl[0][i] = fb[Nt0-1][i]/2;
		fblb[0][i] = fbl[0][i];
		fbr[0][i] = fbl[0][i];
		fbrb[0][i] = fbl[0][i];
		fphp[0][i] = fph[Nt0-1][i]/2;
		fphm[0][i] = fphp[0][i];
		fgp[0][i] = fg[Nt0-1][i]/2;
		fgm[0][i] = fgp[0][i];
	}
	/*Matching with the EVA for the EW gauge bosons*/
	for(i=0;i<Nx-1;i++){
		fwml[0][i] = fwmleva(t0,grid[i]);
		fwmp[0][i] = fwmpeva(t0,grid[i]);
		fwmm[0][i] = fwmmeva(t0,grid[i]);
		fzl[0][i] = fzleva(t0,grid[i]);
		fzp[0][i] = fzpeva(t0,grid[i]);
		fzm[0][i] = fzmeva(t0,grid[i]);
	}
	printf("Matching at EW scale done. Starting the evolution for the SM.\n");
	
	/* ----------------------- SM DGLAP evolution ----------------------- */
	/*Runge-Kutta for the second phase: as before we compute the four Runge-Kutta coefficients corresponding to each f(x_j,t_i). RKj and RKjuc contains the contributions from x_j (massless and ultra-collinear contributions respectively), RK0 from x_k*/
	for(i=0;i<Nt-1;i++){
		printf("SM evolution: step %d of %d\n",i+1,Nt-1);
		for(j=0;j<Nx-1;j++){ 
			k1 = RKj(i,j,Nx,t0+i*dt,dt,grid,deltagrid,fel[i][j],felb[i][j],fml[i][j],fnue[i][j],fnueb[i][j],fnum[i][j],fer[i][j],ferb[i][j],fmr[i][j],ful[i][j],fulb[i][j],ftl[i][j],ftlb[i][j],fdl[i][j],fdlb[i][j],fbl[i][j],fblb[i][j],fur[i][j],furb[i][j],ftr[i][j],ftrb[i][j],fdr[i][j],fdrb[i][j],fbr[i][j],fbrb[i][j],fh[i][j],fwpl[i][j],fwml[i][j],fzl[i][j],fwpp[i][j],fwpm[i][j],fwmp[i][j],fwmm[i][j],fzp[i][j],fzm[i][j],fphp[i][j],fphm[i][j],fzphp[i][j],fzphm[i][j],fgp[i][j],fgm[i][j],fhzl[i][j],ftal[i][j],ftalb[i][j],ftar[i][j],ftarb[i][j],fnuta[i][j],fnutab[i][j],fcl[i][j],fclb[i][j],fcr[i][j],fcrb[i][j],fsl[i][j],fslb[i][j],k1,newterm(grid,deltagrid,j,Nx),talpha);
			k1 = RKjuc(i,j,Nx,t0+i*dt,dt,grid,deltagrid,fel[i][j],felb[i][j],fml[i][j],fnue[i][j],fnueb[i][j],fnum[i][j],fer[i][j],ferb[i][j],fmr[i][j],ful[i][j],fulb[i][j],ftl[i][j],ftlb[i][j],fdl[i][j],fdlb[i][j],fbl[i][j],fblb[i][j],fur[i][j],furb[i][j],ftr[i][j],ftrb[i][j],fdr[i][j],fdrb[i][j],fbr[i][j],fbrb[i][j],fh[i][j],fwpl[i][j],fwml[i][j],fzl[i][j],fwpp[i][j],fwpm[i][j],fwmp[i][j],fwmm[i][j],fzp[i][j],fzm[i][j],fphp[i][j],fphm[i][j],fzphp[i][j],fzphm[i][j],fgp[i][j],fgm[i][j],fhzl[i][j],ftal[i][j],ftalb[i][j],ftar[i][j],ftarb[i][j],fnuta[i][j],fnutab[i][j],fcl[i][j],fclb[i][j],fcr[i][j],fcrb[i][j],fsl[i][j],fslb[i][j],k1,newterm(grid,deltagrid,j,Nx),talpha);
			for(k=j+1;k<Nx;k++){
				k1 = RK(i,j,k,Nx,t0+i*dt,dt,grid,deltagrid,fel[i][k],felb[i][k],fml[i][k],fnue[i][k],fnueb[i][k],fnum[i][k],fer[i][k],ferb[i][k],fmr[i][k],ful[i][k],fulb[i][k],ftl[i][k],ftlb[i][k],fdl[i][k],fdlb[i][k],fbl[i][k],fblb[i][k],fur[i][k],furb[i][k],ftr[i][k],ftrb[i][k],fdr[i][k],fdrb[i][k],fbr[i][k],fbrb[i][k],fh[i][k],fwpl[i][k],fwml[i][k],fzl[i][k],fwpp[i][k],fwpm[i][k],fwmp[i][k],fwmm[i][k],fzp[i][k],fzm[i][k],fphp[i][k],fphm[i][k],fzphp[i][k],fzphm[i][k],fgp[i][k],fgm[i][k],fhzl[i][k],ftal[i][k],ftalb[i][k],ftar[i][k],ftarb[i][k],fnuta[i][k],fnutab[i][k],fcl[i][k],fclb[i][k],fcr[i][k],fcrb[i][k],fsl[i][k],fslb[i][k],k1,talpha);
			}
		}
		for(j=0;j<Nx-1;j++){
			k2 = RKj(i,j,Nx,t0+(i+0.5)*dt,dt,grid,deltagrid,fel[i][j]+0.5*k1[i][j],felb[i][j]+0.5*k1[i][j+Nx],fml[i][j]+0.5*k1[i][j+2*Nx],fnue[i][j]+0.5*k1[i][j+3*Nx],fnueb[i][j]+0.5*k1[i][j+4*Nx],fnum[i][j]+0.5*k1[i][j+5*Nx],fer[i][j]+0.5*k1[i][j+6*Nx],ferb[i][j]+0.5*k1[i][j+7*Nx],fmr[i][j]+0.5*k1[i][j+8*Nx],ful[i][j]+0.5*k1[i][j+9*Nx],fulb[i][j]+0.5*k1[i][j+10*Nx],ftl[i][j]+0.5*k1[i][j+11*Nx],ftlb[i][j]+0.5*k1[i][j+12*Nx],fdl[i][j]+0.5*k1[i][j+13*Nx],fdlb[i][j]+0.5*k1[i][j+14*Nx],fbl[i][j]+0.5*k1[i][j+15*Nx],fblb[i][j]+0.5*k1[i][j+16*Nx],fur[i][j]+0.5*k1[i][j+17*Nx],furb[i][j]+0.5*k1[i][j+18*Nx],ftr[i][j]+0.5*k1[i][j+19*Nx],ftrb[i][j]+0.5*k1[i][j+20*Nx],fdr[i][j]+0.5*k1[i][j+21*Nx],fdrb[i][j]+0.5*k1[i][j+22*Nx],fbr[i][j]+0.5*k1[i][j+23*Nx],fbrb[i][j]+0.5*k1[i][j+24*Nx],fh[i][j]+0.5*k1[i][j+25*Nx],fwpl[i][j]+0.5*k1[i][j+26*Nx],fwml[i][j]+0.5*k1[i][j+27*Nx],fzl[i][j]+0.5*k1[i][j+28*Nx],fwpp[i][j]+0.5*k1[i][j+29*Nx],fwpm[i][j]+0.5*k1[i][j+30*Nx],fwmp[i][j]+0.5*k1[i][j+31*Nx],fwmm[i][j]+0.5*k1[i][j+32*Nx],fzp[i][j]+0.5*k1[i][j+33*Nx],fzm[i][j]+0.5*k1[i][j+34*Nx],fphp[i][j]+0.5*k1[i][j+35*Nx],fphm[i][j]+0.5*k1[i][j+36*Nx],fzphp[i][j]+0.5*k1[i][j+37*Nx],fzphm[i][j]+0.5*k1[i][j+38*Nx],fgp[i][j]+0.5*k1[i][j+39*Nx],fgm[i][j]+0.5*k1[i][j+40*Nx],fhzl[i][j]+0.5*k1[i][j+41*Nx],ftal[i][j]+0.5*k1[i][j+42*Nx],ftalb[i][j]+0.5*k1[i][j+43*Nx],ftar[i][j]+0.5*k1[i][j+44*Nx],ftarb[i][j]+0.5*k1[i][j+45*Nx],fnuta[i][j]+0.5*k1[i][j+46*Nx],fnutab[i][j]+0.5*k1[i][j+47*Nx],fcl[i][j]+0.5*k1[i][j+48*Nx],fclb[i][j]+0.5*k1[i][j+49*Nx],fcr[i][j]+0.5*k1[i][j+50*Nx],fcrb[i][j]+0.5*k1[i][j+51*Nx],fsl[i][j]+0.5*k1[i][j+52*Nx],fslb[i][j]+0.5*k1[i][j+53*Nx],k2,newterm(grid,deltagrid,j,Nx),talpha);
			k2 = RKjuc(i,j,Nx,t0+(i+0.5)*dt,dt,grid,deltagrid,fel[i][j]+0.5*k1[i][j],felb[i][j]+0.5*k1[i][j+Nx],fml[i][j]+0.5*k1[i][j+2*Nx],fnue[i][j]+0.5*k1[i][j+3*Nx],fnueb[i][j]+0.5*k1[i][j+4*Nx],fnum[i][j]+0.5*k1[i][j+5*Nx],fer[i][j]+0.5*k1[i][j+6*Nx],ferb[i][j]+0.5*k1[i][j+7*Nx],fmr[i][j]+0.5*k1[i][j+8*Nx],ful[i][j]+0.5*k1[i][j+9*Nx],fulb[i][j]+0.5*k1[i][j+10*Nx],ftl[i][j]+0.5*k1[i][j+11*Nx],ftlb[i][j]+0.5*k1[i][j+12*Nx],fdl[i][j]+0.5*k1[i][j+13*Nx],fdlb[i][j]+0.5*k1[i][j+14*Nx],fbl[i][j]+0.5*k1[i][j+15*Nx],fblb[i][j]+0.5*k1[i][j+16*Nx],fur[i][j]+0.5*k1[i][j+17*Nx],furb[i][j]+0.5*k1[i][j+18*Nx],ftr[i][j]+0.5*k1[i][j+19*Nx],ftrb[i][j]+0.5*k1[i][j+20*Nx],fdr[i][j]+0.5*k1[i][j+21*Nx],fdrb[i][j]+0.5*k1[i][j+22*Nx],fbr[i][j]+0.5*k1[i][j+23*Nx],fbrb[i][j]+0.5*k1[i][j+24*Nx],fh[i][j]+0.5*k1[i][j+25*Nx],fwpl[i][j]+0.5*k1[i][j+26*Nx],fwml[i][j]+0.5*k1[i][j+27*Nx],fzl[i][j]+0.5*k1[i][j+28*Nx],fwpp[i][j]+0.5*k1[i][j+29*Nx],fwpm[i][j]+0.5*k1[i][j+30*Nx],fwmp[i][j]+0.5*k1[i][j+31*Nx],fwmm[i][j]+0.5*k1[i][j+32*Nx],fzp[i][j]+0.5*k1[i][j+33*Nx],fzm[i][j]+0.5*k1[i][j+34*Nx],fphp[i][j]+0.5*k1[i][j+35*Nx],fphm[i][j]+0.5*k1[i][j+36*Nx],fzphp[i][j]+0.5*k1[i][j+37*Nx],fzphm[i][j]+0.5*k1[i][j+38*Nx],fgp[i][j]+0.5*k1[i][j+39*Nx],fgm[i][j]+0.5*k1[i][j+40*Nx],fhzl[i][j]+0.5*k1[i][j+41*Nx],ftal[i][j]+0.5*k1[i][j+42*Nx],ftalb[i][j]+0.5*k1[i][j+43*Nx],ftar[i][j]+0.5*k1[i][j+44*Nx],ftarb[i][j]+0.5*k1[i][j+45*Nx],fnuta[i][j]+0.5*k1[i][j+46*Nx],fnutab[i][j]+0.5*k1[i][j+47*Nx],fcl[i][j]+0.5*k1[i][j+48*Nx],fclb[i][j]+0.5*k1[i][j+49*Nx],fcr[i][j]+0.5*k1[i][j+50*Nx],fcrb[i][j]+0.5*k1[i][j+51*Nx],fsl[i][j]+0.5*k1[i][j+52*Nx],fslb[i][j]+0.5*k1[i][j+53*Nx],k2,newterm(grid,deltagrid,j,Nx),talpha);
			for(k=j+1;k<Nx;k++){
				k2 = RK(i,j,k,Nx,t0+(i+0.5)*dt,dt,grid,deltagrid,fel[i][k]+0.5*k1[i][k],felb[i][k]+0.5*k1[i][k+Nx],fml[i][k]+0.5*k1[i][k+2*Nx],fnue[i][k]+0.5*k1[i][k+3*Nx],fnueb[i][k]+0.5*k1[i][k+4*Nx],fnum[i][k]+0.5*k1[i][k+5*Nx],fer[i][k]+0.5*k1[i][k+6*Nx],ferb[i][k]+0.5*k1[i][k+7*Nx],fmr[i][k]+0.5*k1[i][k+8*Nx],ful[i][k]+0.5*k1[i][k+9*Nx],fulb[i][k]+0.5*k1[i][k+10*Nx],ftl[i][k]+0.5*k1[i][k+11*Nx],ftlb[i][k]+0.5*k1[i][k+12*Nx],fdl[i][k]+0.5*k1[i][k+13*Nx],fdlb[i][k]+0.5*k1[i][k+14*Nx],fbl[i][k]+0.5*k1[i][k+15*Nx],fblb[i][k]+0.5*k1[i][k+16*Nx],fur[i][k]+0.5*k1[i][k+17*Nx],furb[i][k]+0.5*k1[i][k+18*Nx],ftr[i][k]+0.5*k1[i][k+19*Nx],ftrb[i][k]+0.5*k1[i][k+20*Nx],fdr[i][k]+0.5*k1[i][k+21*Nx],fdrb[i][k]+0.5*k1[i][k+22*Nx],fbr[i][k]+0.5*k1[i][k+23*Nx],fbrb[i][k]+0.5*k1[i][k+24*Nx],fh[i][k]+0.5*k1[i][k+25*Nx],fwpl[i][k]+0.5*k1[i][k+26*Nx],fwml[i][k]+0.5*k1[i][k+27*Nx],fzl[i][k]+0.5*k1[i][k+28*Nx],fwpp[i][k]+0.5*k1[i][k+29*Nx],fwpm[i][k]+0.5*k1[i][k+30*Nx],fwmp[i][k]+0.5*k1[i][k+31*Nx],fwmm[i][k]+0.5*k1[i][k+32*Nx],fzp[i][k]+0.5*k1[i][k+33*Nx],fzm[i][k]+0.5*k1[i][k+34*Nx],fphp[i][k]+0.5*k1[i][k+35*Nx],fphm[i][k]+0.5*k1[i][k+36*Nx],fzphp[i][k]+0.5*k1[i][k+37*Nx],fzphm[i][k]+0.5*k1[i][k+38*Nx],fgp[i][k]+0.5*k1[i][k+39*Nx],fgm[i][k]+0.5*k1[i][k+40*Nx],fhzl[i][k]+0.5*k1[i][k+41*Nx],ftal[i][k]+0.5*k1[i][k+42*Nx],ftalb[i][k]+0.5*k1[i][k+43*Nx],ftar[i][k]+0.5*k1[i][k+44*Nx],ftarb[i][k]+0.5*k1[i][k+45*Nx],fnuta[i][k]+0.5*k1[i][k+46*Nx],fnutab[i][k]+0.5*k1[i][k+47*Nx],fcl[i][k]+0.5*k1[i][k+48*Nx],fclb[i][k]+0.5*k1[i][k+49*Nx],fcr[i][k]+0.5*k1[i][k+50*Nx],fcrb[i][k]+0.5*k1[i][k+51*Nx],fsl[i][k]+0.5*k1[i][k+52*Nx],fslb[i][k]+0.5*k1[i][k+53*Nx],k2,talpha);
			}
		}
		for(j=0;j<Nx-1;j++){
			k3 = RKj(i,j,Nx,t0+(i+0.5)*dt,dt,grid,deltagrid,fel[i][j]+0.5*k2[i][j],felb[i][j]+0.5*k2[i][j+Nx],fml[i][j]+0.5*k2[i][j+2*Nx],fnue[i][j]+0.5*k2[i][j+3*Nx],fnueb[i][j]+0.5*k2[i][j+4*Nx],fnum[i][j]+0.5*k2[i][j+5*Nx],fer[i][j]+0.5*k2[i][j+6*Nx],ferb[i][j]+0.5*k2[i][j+7*Nx],fmr[i][j]+0.5*k2[i][j+8*Nx],ful[i][j]+0.5*k2[i][j+9*Nx],fulb[i][j]+0.5*k2[i][j+10*Nx],ftl[i][j]+0.5*k2[i][j+11*Nx],ftlb[i][j]+0.5*k2[i][j+12*Nx],fdl[i][j]+0.5*k2[i][j+13*Nx],fdlb[i][j]+0.5*k2[i][j+14*Nx],fbl[i][j]+0.5*k2[i][j+15*Nx],fblb[i][j]+0.5*k2[i][j+16*Nx],fur[i][j]+0.5*k2[i][j+17*Nx],furb[i][j]+0.5*k2[i][j+18*Nx],ftr[i][j]+0.5*k2[i][j+19*Nx],ftrb[i][j]+0.5*k2[i][j+20*Nx],fdr[i][j]+0.5*k2[i][j+21*Nx],fdrb[i][j]+0.5*k2[i][j+22*Nx],fbr[i][j]+0.5*k2[i][j+23*Nx],fbrb[i][j]+0.5*k2[i][j+24*Nx],fh[i][j]+0.5*k2[i][j+25*Nx],fwpl[i][j]+0.5*k2[i][j+26*Nx],fwml[i][j]+0.5*k2[i][j+27*Nx],fzl[i][j]+0.5*k2[i][j+28*Nx],fwpp[i][j]+0.5*k2[i][j+29*Nx],fwpm[i][j]+0.5*k2[i][j+30*Nx],fwmp[i][j]+0.5*k2[i][j+31*Nx],fwmm[i][j]+0.5*k2[i][j+32*Nx],fzp[i][j]+0.5*k2[i][j+33*Nx],fzm[i][j]+0.5*k2[i][j+34*Nx],fphp[i][j]+0.5*k2[i][j+35*Nx],fphm[i][j]+0.5*k2[i][j+36*Nx],fzphp[i][j]+0.5*k2[i][j+37*Nx],fzphm[i][j]+0.5*k2[i][j+38*Nx],fgp[i][j]+0.5*k2[i][j+39*Nx],fgm[i][j]+0.5*k2[i][j+40*Nx],fhzl[i][j]+0.5*k2[i][j+41*Nx],ftal[i][j]+0.5*k2[i][j+42*Nx],ftalb[i][j]+0.5*k2[i][j+43*Nx],ftar[i][j]+0.5*k2[i][j+44*Nx],ftarb[i][j]+0.5*k2[i][j+45*Nx],fnuta[i][j]+0.5*k2[i][j+46*Nx],fnutab[i][j]+0.5*k2[i][j+47*Nx],fcl[i][j]+0.5*k2[i][j+48*Nx],fclb[i][j]+0.5*k2[i][j+49*Nx],fcr[i][j]+0.5*k2[i][j+50*Nx],fcrb[i][j]+0.5*k2[i][j+51*Nx],fsl[i][j]+0.5*k2[i][j+52*Nx],fslb[i][j]+0.5*k2[i][j+53*Nx],k3,newterm(grid,deltagrid,j,Nx),talpha);
			k3 = RKjuc(i,j,Nx,t0+(i+0.5)*dt,dt,grid,deltagrid,fel[i][j]+0.5*k2[i][j],felb[i][j]+0.5*k2[i][j+Nx],fml[i][j]+0.5*k2[i][j+2*Nx],fnue[i][j]+0.5*k2[i][j+3*Nx],fnueb[i][j]+0.5*k2[i][j+4*Nx],fnum[i][j]+0.5*k2[i][j+5*Nx],fer[i][j]+0.5*k2[i][j+6*Nx],ferb[i][j]+0.5*k2[i][j+7*Nx],fmr[i][j]+0.5*k2[i][j+8*Nx],ful[i][j]+0.5*k2[i][j+9*Nx],fulb[i][j]+0.5*k2[i][j+10*Nx],ftl[i][j]+0.5*k2[i][j+11*Nx],ftlb[i][j]+0.5*k2[i][j+12*Nx],fdl[i][j]+0.5*k2[i][j+13*Nx],fdlb[i][j]+0.5*k2[i][j+14*Nx],fbl[i][j]+0.5*k2[i][j+15*Nx],fblb[i][j]+0.5*k2[i][j+16*Nx],fur[i][j]+0.5*k2[i][j+17*Nx],furb[i][j]+0.5*k2[i][j+18*Nx],ftr[i][j]+0.5*k2[i][j+19*Nx],ftrb[i][j]+0.5*k2[i][j+20*Nx],fdr[i][j]+0.5*k2[i][j+21*Nx],fdrb[i][j]+0.5*k2[i][j+22*Nx],fbr[i][j]+0.5*k2[i][j+23*Nx],fbrb[i][j]+0.5*k2[i][j+24*Nx],fh[i][j]+0.5*k2[i][j+25*Nx],fwpl[i][j]+0.5*k2[i][j+26*Nx],fwml[i][j]+0.5*k2[i][j+27*Nx],fzl[i][j]+0.5*k2[i][j+28*Nx],fwpp[i][j]+0.5*k2[i][j+29*Nx],fwpm[i][j]+0.5*k2[i][j+30*Nx],fwmp[i][j]+0.5*k2[i][j+31*Nx],fwmm[i][j]+0.5*k2[i][j+32*Nx],fzp[i][j]+0.5*k2[i][j+33*Nx],fzm[i][j]+0.5*k2[i][j+34*Nx],fphp[i][j]+0.5*k2[i][j+35*Nx],fphm[i][j]+0.5*k2[i][j+36*Nx],fzphp[i][j]+0.5*k2[i][j+37*Nx],fzphm[i][j]+0.5*k2[i][j+38*Nx],fgp[i][j]+0.5*k2[i][j+39*Nx],fgm[i][j]+0.5*k2[i][j+40*Nx],fhzl[i][j]+0.5*k2[i][j+41*Nx],ftal[i][j]+0.5*k2[i][j+42*Nx],ftalb[i][j]+0.5*k2[i][j+43*Nx],ftar[i][j]+0.5*k2[i][j+44*Nx],ftarb[i][j]+0.5*k2[i][j+45*Nx],fnuta[i][j]+0.5*k2[i][j+46*Nx],fnutab[i][j]+0.5*k2[i][j+47*Nx],fcl[i][j]+0.5*k2[i][j+48*Nx],fclb[i][j]+0.5*k2[i][j+49*Nx],fcr[i][j]+0.5*k2[i][j+50*Nx],fcrb[i][j]+0.5*k2[i][j+51*Nx],fsl[i][j]+0.5*k2[i][j+52*Nx],fslb[i][j]+0.5*k2[i][j+53*Nx],k3,newterm(grid,deltagrid,j,Nx),talpha);
			for(k=j+1;k<Nx;k++){
				k3 = RK(i,j,k,Nx,t0+(i+0.5)*dt,dt,grid,deltagrid,fel[i][k]+0.5*k2[i][k],felb[i][k]+0.5*k2[i][k+Nx],fml[i][k]+0.5*k2[i][k+2*Nx],fnue[i][k]+0.5*k2[i][k+3*Nx],fnueb[i][k]+0.5*k2[i][k+4*Nx],fnum[i][k]+0.5*k2[i][k+5*Nx],fer[i][k]+0.5*k2[i][k+6*Nx],ferb[i][k]+0.5*k2[i][k+7*Nx],fmr[i][k]+0.5*k2[i][k+8*Nx],ful[i][k]+0.5*k2[i][k+9*Nx],fulb[i][k]+0.5*k2[i][k+10*Nx],ftl[i][k]+0.5*k2[i][k+11*Nx],ftlb[i][k]+0.5*k2[i][k+12*Nx],fdl[i][k]+0.5*k2[i][k+13*Nx],fdlb[i][k]+0.5*k2[i][k+14*Nx],fbl[i][k]+0.5*k2[i][k+15*Nx],fblb[i][k]+0.5*k2[i][k+16*Nx],fur[i][k]+0.5*k2[i][k+17*Nx],furb[i][k]+0.5*k2[i][k+18*Nx],ftr[i][k]+0.5*k2[i][k+19*Nx],ftrb[i][k]+0.5*k2[i][k+20*Nx],fdr[i][k]+0.5*k2[i][k+21*Nx],fdrb[i][k]+0.5*k2[i][k+22*Nx],fbr[i][k]+0.5*k2[i][k+23*Nx],fbrb[i][k]+0.5*k2[i][k+24*Nx],fh[i][k]+0.5*k2[i][k+25*Nx],fwpl[i][k]+0.5*k2[i][k+26*Nx],fwml[i][k]+0.5*k2[i][k+27*Nx],fzl[i][k]+0.5*k2[i][k+28*Nx],fwpp[i][k]+0.5*k2[i][k+29*Nx],fwpm[i][k]+0.5*k2[i][k+30*Nx],fwmp[i][k]+0.5*k2[i][k+31*Nx],fwmm[i][k]+0.5*k2[i][k+32*Nx],fzp[i][k]+0.5*k2[i][k+33*Nx],fzm[i][k]+0.5*k2[i][k+34*Nx],fphp[i][k]+0.5*k2[i][k+35*Nx],fphm[i][k]+0.5*k2[i][k+36*Nx],fzphp[i][k]+0.5*k2[i][k+37*Nx],fzphm[i][k]+0.5*k2[i][k+38*Nx],fgp[i][k]+0.5*k2[i][k+39*Nx],fgm[i][k]+0.5*k2[i][k+40*Nx],fhzl[i][k]+0.5*k2[i][k+41*Nx],ftal[i][k]+0.5*k2[i][k+42*Nx],ftalb[i][k]+0.5*k2[i][k+43*Nx],ftar[i][k]+0.5*k2[i][k+44*Nx],ftarb[i][k]+0.5*k2[i][k+45*Nx],fnuta[i][k]+0.5*k2[i][k+46*Nx],fnutab[i][k]+0.5*k2[i][k+47*Nx],fcl[i][k]+0.5*k2[i][k+48*Nx],fclb[i][k]+0.5*k2[i][k+49*Nx],fcr[i][k]+0.5*k2[i][k+50*Nx],fcrb[i][k]+0.5*k2[i][k+51*Nx],fsl[i][k]+0.5*k2[i][k+52*Nx],fslb[i][k]+0.5*k2[i][k+53*Nx],k3,talpha);
			}
		}
		for(j=0;j<Nx-1;j++){
			k4 = RKj(i,j,Nx,t0+(i+1)*dt,dt,grid,deltagrid,fel[i][j]+k3[i][j],felb[i][j]+k3[i][j+Nx],fml[i][j]+k3[i][j+2*Nx],fnue[i][j]+k3[i][j+3*Nx],fnueb[i][j]+k3[i][j+4*Nx],fnum[i][j]+k3[i][j+5*Nx],fer[i][j]+k3[i][j+6*Nx],ferb[i][j]+k3[i][j+7*Nx],fmr[i][j]+k3[i][j+8*Nx],ful[i][j]+k3[i][j+9*Nx],fulb[i][j]+k3[i][j+10*Nx],ftl[i][j]+k3[i][j+11*Nx],ftlb[i][j]+k3[i][j+12*Nx],fdl[i][j]+k3[i][j+13*Nx],fdlb[i][j]+k3[i][j+14*Nx],fbl[i][j]+k3[i][j+15*Nx],fblb[i][j]+k3[i][j+16*Nx],fur[i][j]+k3[i][j+17*Nx],furb[i][j]+k3[i][j+18*Nx],ftr[i][j]+k3[i][j+19*Nx],ftrb[i][j]+k3[i][j+20*Nx],fdr[i][j]+k3[i][j+21*Nx],fdrb[i][j]+k3[i][j+22*Nx],fbr[i][j]+k3[i][j+23*Nx],fbrb[i][j]+k3[i][j+24*Nx],fh[i][j]+k3[i][j+25*Nx],fwpl[i][j]+k3[i][j+26*Nx],fwml[i][j]+k3[i][j+27*Nx],fzl[i][j]+k3[i][j+28*Nx],fwpp[i][j]+k3[i][j+29*Nx],fwpm[i][j]+k3[i][j+30*Nx],fwmp[i][j]+k3[i][j+31*Nx],fwmm[i][j]+k3[i][j+32*Nx],fzp[i][j]+k3[i][j+33*Nx],fzm[i][j]+k3[i][j+34*Nx],fphp[i][j]+k3[i][j+35*Nx],fphm[i][j]+k3[i][j+36*Nx],fzphp[i][j]+k3[i][j+37*Nx],fzphm[i][j]+k3[i][j+38*Nx],fgp[i][j]+k3[i][j+39*Nx],fgm[i][j]+k3[i][j+40*Nx],fhzl[i][j]+k3[i][j+41*Nx],ftal[i][j]+k3[i][j+42*Nx],ftalb[i][j]+k3[i][j+43*Nx],ftar[i][j]+k3[i][j+44*Nx],ftarb[i][j]+k3[i][j+45*Nx],fnuta[i][j]+k3[i][j+46*Nx],fnutab[i][j]+k3[i][j+47*Nx],fcl[i][j]+k3[i][j+48*Nx],fclb[i][j]+k3[i][j+49*Nx],fcr[i][j]+k3[i][j+50*Nx],fcrb[i][j]+k3[i][j+51*Nx],fsl[i][j]+k3[i][j+52*Nx],fslb[i][j]+k3[i][j+53*Nx],k4,newterm(grid,deltagrid,j,Nx),talpha);
			k4 = RKjuc(i,j,Nx,t0+(i+1)*dt,dt,grid,deltagrid,fel[i][j]+k3[i][j],felb[i][j]+k3[i][j+Nx],fml[i][j]+k3[i][j+2*Nx],fnue[i][j]+k3[i][j+3*Nx],fnueb[i][j]+k3[i][j+4*Nx],fnum[i][j]+k3[i][j+5*Nx],fer[i][j]+k3[i][j+6*Nx],ferb[i][j]+k3[i][j+7*Nx],fmr[i][j]+k3[i][j+8*Nx],ful[i][j]+k3[i][j+9*Nx],fulb[i][j]+k3[i][j+10*Nx],ftl[i][j]+k3[i][j+11*Nx],ftlb[i][j]+k3[i][j+12*Nx],fdl[i][j]+k3[i][j+13*Nx],fdlb[i][j]+k3[i][j+14*Nx],fbl[i][j]+k3[i][j+15*Nx],fblb[i][j]+k3[i][j+16*Nx],fur[i][j]+k3[i][j+17*Nx],furb[i][j]+k3[i][j+18*Nx],ftr[i][j]+k3[i][j+19*Nx],ftrb[i][j]+k3[i][j+20*Nx],fdr[i][j]+k3[i][j+21*Nx],fdrb[i][j]+k3[i][j+22*Nx],fbr[i][j]+k3[i][j+23*Nx],fbrb[i][j]+k3[i][j+24*Nx],fh[i][j]+k3[i][j+25*Nx],fwpl[i][j]+k3[i][j+26*Nx],fwml[i][j]+k3[i][j+27*Nx],fzl[i][j]+k3[i][j+28*Nx],fwpp[i][j]+k3[i][j+29*Nx],fwpm[i][j]+k3[i][j+30*Nx],fwmp[i][j]+k3[i][j+31*Nx],fwmm[i][j]+k3[i][j+32*Nx],fzp[i][j]+k3[i][j+33*Nx],fzm[i][j]+k3[i][j+34*Nx],fphp[i][j]+k3[i][j+35*Nx],fphm[i][j]+k3[i][j+36*Nx],fzphp[i][j]+k3[i][j+37*Nx],fzphm[i][j]+k3[i][j+38*Nx],fgp[i][j]+k3[i][j+39*Nx],fgm[i][j]+k3[i][j+40*Nx],fhzl[i][j]+k3[i][j+41*Nx],ftal[i][j]+k3[i][j+42*Nx],ftalb[i][j]+k3[i][j+43*Nx],ftar[i][j]+k3[i][j+44*Nx],ftarb[i][j]+k3[i][j+45*Nx],fnuta[i][j]+k3[i][j+46*Nx],fnutab[i][j]+k3[i][j+47*Nx],fcl[i][j]+k3[i][j+48*Nx],fclb[i][j]+k3[i][j+49*Nx],fcr[i][j]+k3[i][j+50*Nx],fcrb[i][j]+k3[i][j+51*Nx],fsl[i][j]+k3[i][j+52*Nx],fslb[i][j]+k3[i][j+53*Nx],k4,newterm(grid,deltagrid,j,Nx),talpha);
			for(k=j+1;k<Nx;k++){
				k4 = RK(i,j,k,Nx,t0+(i+1)*dt,dt,grid,deltagrid,fel[i][k]+k3[i][k],felb[i][k]+k3[i][k+Nx],fml[i][k]+k3[i][k+2*Nx],fnue[i][k]+k3[i][k+3*Nx],fnueb[i][k]+k3[i][k+4*Nx],fnum[i][k]+k3[i][k+5*Nx],fer[i][k]+k3[i][k+6*Nx],ferb[i][k]+k3[i][k+7*Nx],fmr[i][k]+k3[i][k+8*Nx],ful[i][k]+k3[i][k+9*Nx],fulb[i][k]+k3[i][k+10*Nx],ftl[i][k]+k3[i][k+11*Nx],ftlb[i][k]+k3[i][k+12*Nx],fdl[i][k]+k3[i][k+13*Nx],fdlb[i][k]+k3[i][k+14*Nx],fbl[i][k]+k3[i][k+15*Nx],fblb[i][k]+k3[i][k+16*Nx],fur[i][k]+k3[i][k+17*Nx],furb[i][k]+k3[i][k+18*Nx],ftr[i][k]+k3[i][k+19*Nx],ftrb[i][k]+k3[i][k+20*Nx],fdr[i][k]+k3[i][k+21*Nx],fdrb[i][k]+k3[i][k+22*Nx],fbr[i][k]+k3[i][k+23*Nx],fbrb[i][k]+k3[i][k+24*Nx],fh[i][k]+k3[i][k+25*Nx],fwpl[i][k]+k3[i][k+26*Nx],fwml[i][k]+k3[i][k+27*Nx],fzl[i][k]+k3[i][k+28*Nx],fwpp[i][k]+k3[i][k+29*Nx],fwpm[i][k]+k3[i][k+30*Nx],fwmp[i][k]+k3[i][k+31*Nx],fwmm[i][k]+k3[i][k+32*Nx],fzp[i][k]+k3[i][k+33*Nx],fzm[i][k]+k3[i][k+34*Nx],fphp[i][k]+k3[i][k+35*Nx],fphm[i][k]+k3[i][k+36*Nx],fzphp[i][k]+k3[i][k+37*Nx],fzphm[i][k]+k3[i][k+38*Nx],fgp[i][k]+k3[i][k+39*Nx],fgm[i][k]+k3[i][k+40*Nx],fhzl[i][k]+k3[i][k+41*Nx],ftal[i][k]+k3[i][k+42*Nx],ftalb[i][k]+k3[i][k+43*Nx],ftar[i][k]+k3[i][k+44*Nx],ftarb[i][k]+k3[i][k+45*Nx],fnuta[i][k]+k3[i][k+46*Nx],fnutab[i][k]+k3[i][k+47*Nx],fcl[i][k]+k3[i][k+48*Nx],fclb[i][k]+k3[i][k+49*Nx],fcr[i][k]+k3[i][k+50*Nx],fcrb[i][k]+k3[i][k+51*Nx],fsl[i][k]+k3[i][k+52*Nx],fslb[i][k]+k3[i][k+53*Nx],k4,talpha);
			}
		}
		for(j=0;j<Nx-1;j++){
			fel[i+1][j] = fel[i][j] + k1[i][j]/6 + k2[i][j]/3 + k3[i][j]/3 + k4[i][j]/6;
			felb[i+1][j] = felb[i][j] + k1[i][j+Nx]/6 + k2[i][j+Nx]/3 + k3[i][j+Nx]/3 + k4[i][j+Nx]/6;
			fml[i+1][j] = fml[i][j] + k1[i][j+2*Nx]/6 + k2[i][j+2*Nx]/3 + k3[i][j+2*Nx]/3 + k4[i][j+2*Nx]/6;
			fnue[i+1][j] = fnue[i][j] + k1[i][j+3*Nx]/6 + k2[i][j+3*Nx]/3 + k3[i][j+3*Nx]/3 + k4[i][j+3*Nx]/6;
			fnueb[i+1][j] = fnueb[i][j] + k1[i][j+4*Nx]/6 + k2[i][j+4*Nx]/3 + k3[i][j+4*Nx]/3 + k4[i][j+4*Nx]/6;
			fnum[i+1][j] = fnum[i][j] + k1[i][j+5*Nx]/6 + k2[i][j+5*Nx]/3 + k3[i][j+5*Nx]/3 + k4[i][j+5*Nx]/6;
			fer[i+1][j] = fer[i][j] + k1[i][j+6*Nx]/6 + k2[i][j+6*Nx]/3 + k3[i][j+6*Nx]/3 + k4[i][j+6*Nx]/6;
			ferb[i+1][j] = ferb[i][j] + k1[i][j+7*Nx]/6 + k2[i][j+7*Nx]/3 + k3[i][j+7*Nx]/3 + k4[i][j+7*Nx]/6;
			fmr[i+1][j] = fmr[i][j] + k1[i][j+8*Nx]/6 + k2[i][j+8*Nx]/3 + k3[i][j+8*Nx]/3 + k4[i][j+8*Nx]/6;
			ful[i+1][j] = ful[i][j] + k1[i][j+9*Nx]/6 + k2[i][j+9*Nx]/3 + k3[i][j+9*Nx]/3 + k4[i][j+9*Nx]/6;
			fulb[i+1][j] = fulb[i][j] + k1[i][j+10*Nx]/6 + k2[i][j+10*Nx]/3 + k3[i][j+10*Nx]/3 + k4[i][j+10*Nx]/6;
			ftl[i+1][j] = ftl[i][j] + k1[i][j+11*Nx]/6 + k2[i][j+11*Nx]/3 + k3[i][j+11*Nx]/3 + k4[i][j+11*Nx]/6;
			ftlb[i+1][j] = ftlb[i][j] + k1[i][j+12*Nx]/6 + k2[i][j+12*Nx]/3 + k3[i][j+12*Nx]/3 + k4[i][j+12*Nx]/6;
			fdl[i+1][j] = fdl[i][j] + k1[i][j+13*Nx]/6 + k2[i][j+13*Nx]/3 + k3[i][j+13*Nx]/3 + k4[i][j+13*Nx]/6;
			fdlb[i+1][j] = fdlb[i][j] + k1[i][j+14*Nx]/6 + k2[i][j+14*Nx]/3 + k3[i][j+14*Nx]/3 + k4[i][j+14*Nx]/6;
			fbl[i+1][j] = fbl[i][j] + k1[i][j+15*Nx]/6 + k2[i][j+15*Nx]/3 + k3[i][j+15*Nx]/3 + k4[i][j+15*Nx]/6;
			fblb[i+1][j] = fblb[i][j] + k1[i][j+16*Nx]/6 + k2[i][j+16*Nx]/3 + k3[i][j+16*Nx]/3 + k4[i][j+16*Nx]/6;
			fur[i+1][j] = fur[i][j] + k1[i][j+17*Nx]/6 + k2[i][j+17*Nx]/3 + k3[i][j+17*Nx]/3 + k4[i][j+17*Nx]/6;
			furb[i+1][j] = furb[i][j] + k1[i][j+18*Nx]/6 + k2[i][j+18*Nx]/3 + k3[i][j+18*Nx]/3 + k4[i][j+18*Nx]/6;
			ftr[i+1][j] = ftr[i][j] + k1[i][j+19*Nx]/6 + k2[i][j+19*Nx]/3 + k3[i][j+19*Nx]/3 + k4[i][j+19*Nx]/6;
			ftrb[i+1][j] = ftrb[i][j] + k1[i][j+20*Nx]/6 + k2[i][j+20*Nx]/3 + k3[i][j+20*Nx]/3 + k4[i][j+20*Nx]/6;
			fdr[i+1][j] = fdr[i][j] + k1[i][j+21*Nx]/6 + k2[i][j+21*Nx]/3 + k3[i][j+21*Nx]/3 + k4[i][j+21*Nx]/6;
			fdrb[i+1][j] = fdrb[i][j] + k1[i][j+22*Nx]/6 + k2[i][j+22*Nx]/3 + k3[i][j+22*Nx]/3 + k4[i][j+22*Nx]/6;
			fbr[i+1][j] = fbr[i][j] + k1[i][j+23*Nx]/6 + k2[i][j+23*Nx]/3 + k3[i][j+23*Nx]/3 + k4[i][j+23*Nx]/6;
			fbrb[i+1][j] = fbrb[i][j] + k1[i][j+24*Nx]/6 + k2[i][j+24*Nx]/3 + k3[i][j+24*Nx]/3 + k4[i][j+24*Nx]/6;
			fh[i+1][j] = fh[i][j] + k1[i][j+25*Nx]/6 + k2[i][j+25*Nx]/3 + k3[i][j+25*Nx]/3 + k4[i][j+25*Nx]/6;
			fwpl[i+1][j] = fwpl[i][j] + k1[i][j+26*Nx]/6 + k2[i][j+26*Nx]/3 + k3[i][j+26*Nx]/3 + k4[i][j+26*Nx]/6;
			fwml[i+1][j] = fwml[i][j] + k1[i][j+27*Nx]/6 + k2[i][j+27*Nx]/3 + k3[i][j+27*Nx]/3 + k4[i][j+27*Nx]/6;
			fzl[i+1][j] = fzl[i][j] + k1[i][j+28*Nx]/6 + k2[i][j+28*Nx]/3 + k3[i][j+28*Nx]/3 + k4[i][j+28*Nx]/6;
			fwpp[i+1][j] = fwpp[i][j] + k1[i][j+29*Nx]/6 + k2[i][j+29*Nx]/3 + k3[i][j+29*Nx]/3 + k4[i][j+29*Nx]/6;
			fwpm[i+1][j] = fwpm[i][j] + k1[i][j+30*Nx]/6 + k2[i][j+30*Nx]/3 + k3[i][j+30*Nx]/3 + k4[i][j+30*Nx]/6;
			fwmp[i+1][j] = fwmp[i][j] + k1[i][j+31*Nx]/6 + k2[i][j+31*Nx]/3 + k3[i][j+31*Nx]/3 + k4[i][j+31*Nx]/6;
			fwmm[i+1][j] = fwmm[i][j] + k1[i][j+32*Nx]/6 + k2[i][j+32*Nx]/3 + k3[i][j+32*Nx]/3 + k4[i][j+32*Nx]/6;
			fzp[i+1][j] = fzp[i][j] + k1[i][j+33*Nx]/6 + k2[i][j+33*Nx]/3 + k3[i][j+33*Nx]/3 + k4[i][j+33*Nx]/6;
			fzm[i+1][j] = fzm[i][j] + k1[i][j+34*Nx]/6 + k2[i][j+34*Nx]/3 + k3[i][j+34*Nx]/3 + k4[i][j+34*Nx]/6;
			fphp[i+1][j] = fphp[i][j] + k1[i][j+35*Nx]/6 + k2[i][j+35*Nx]/3 + k3[i][j+35*Nx]/3 + k4[i][j+35*Nx]/6;
			fphm[i+1][j] = fphm[i][j] + k1[i][j+36*Nx]/6 + k2[i][j+36*Nx]/3 + k3[i][j+36*Nx]/3 + k4[i][j+36*Nx]/6;
			fzphp[i+1][j] = fzphp[i][j] + k1[i][j+37*Nx]/6 + k2[i][j+37*Nx]/3 + k3[i][j+37*Nx]/3 + k4[i][j+37*Nx]/6;
			fzphm[i+1][j] = fzphm[i][j] + k1[i][j+38*Nx]/6 + k2[i][j+38*Nx]/3 + k3[i][j+38*Nx]/3 + k4[i][j+38*Nx]/6;
			fgp[i+1][j] = fgp[i][j] + k1[i][j+39*Nx]/6 + k2[i][j+39*Nx]/3 + k3[i][j+39*Nx]/3 + k4[i][j+39*Nx]/6;
			fgm[i+1][j] = fgm[i][j] + k1[i][j+40*Nx]/6 + k2[i][j+40*Nx]/3 + k3[i][j+40*Nx]/3 + k4[i][j+40*Nx]/6;
			fhzl[i+1][j] = fhzl[i][j] + k1[i][j+41*Nx]/6 + k2[i][j+41*Nx]/3 + k3[i][j+41*Nx]/3 + k4[i][j+41*Nx]/6;
			ftal[i+1][j] = ftal[i][j] + k1[i][j+42*Nx]/6 + k2[i][j+42*Nx]/3 + k3[i][j+42*Nx]/3 + k4[i][j+42*Nx]/6;
			ftalb[i+1][j] = ftalb[i][j] + k1[i][j+43*Nx]/6 + k2[i][j+43*Nx]/3 + k3[i][j+43*Nx]/3 + k4[i][j+43*Nx]/6;
			ftar[i+1][j] = ftar[i][j] + k1[i][j+44*Nx]/6 + k2[i][j+44*Nx]/3 + k3[i][j+44*Nx]/3 + k4[i][j+44*Nx]/6;
			ftarb[i+1][j] = ftarb[i][j] + k1[i][j+45*Nx]/6 + k2[i][j+45*Nx]/3 + k3[i][j+45*Nx]/3 + k4[i][j+45*Nx]/6;
			fnuta[i+1][j] = fnuta[i][j] + k1[i][j+46*Nx]/6 + k2[i][j+46*Nx]/3 + k3[i][j+46*Nx]/3 + k4[i][j+46*Nx]/6;
			fnutab[i+1][j] = fnutab[i][j] + k1[i][j+47*Nx]/6 + k2[i][j+47*Nx]/3 + k3[i][j+47*Nx]/3 + k4[i][j+47*Nx]/6;
			fcl[i+1][j] = fcl[i][j] + k1[i][j+48*Nx]/6 + k2[i][j+48*Nx]/3 + k3[i][j+48*Nx]/3 + k4[i][j+48*Nx]/6;
			fclb[i+1][j] = fclb[i][j] + k1[i][j+49*Nx]/6 + k2[i][j+49*Nx]/3 + k3[i][j+49*Nx]/3 + k4[i][j+49*Nx]/6;
			fcr[i+1][j] = fcr[i][j] + k1[i][j+50*Nx]/6 + k2[i][j+50*Nx]/3 + k3[i][j+50*Nx]/3 + k4[i][j+50*Nx]/6;
			fcrb[i+1][j] = fcrb[i][j] + k1[i][j+51*Nx]/6 + k2[i][j+51*Nx]/3 + k3[i][j+51*Nx]/3 + k4[i][j+51*Nx]/6;
			fsl[i+1][j] = fsl[i][j] + k1[i][j+52*Nx]/6 + k2[i][j+52*Nx]/3 + k3[i][j+52*Nx]/3 + k4[i][j+52*Nx]/6;
			fslb[i+1][j] = fslb[i][j] + k1[i][j+53*Nx]/6 + k2[i][j+53*Nx]/3 + k3[i][j+53*Nx]/3 + k4[i][j+53*Nx]/6;
		}
		fml[i+1][Nx-1] = Lt(grid,deltagrid,Nx,i+1,fel,felb,fml,ftal,ftalb,fnue,fnueb,fnum,fnuta,fnutab,fer,ferb,fmr,ftar,ftarb,ful,fulb,fcl,fclb,ftl,ftlb,fdl,fdlb,fsl,fslb,fbl,fblb,fur,furb,fcr,fcrb,ftr,ftrb,fdr,fdrb,fbr,fbrb,fh,fwpl,fwml,fzl,fwpp,fwpm,fwmp,fwmm,fzp,fzm,fphp,fphm,fgp,fgm)/deltagrid[Nx-1];
		fmr[i+1][Nx-1] = fml[i+1][Nx-1];
	}
	printf("SM evolution completed.\n");
	
	
	/* ----------------------- Exporting the PDFs in LHAPDF6-inspired format ----------------------- */
	printf("Saving the PDFs in LHAPDF6 format.\n");
	int write; /*This variable is used to call the functions that write on the output files*/
	int itSteps = 6; /*To compress the output file, we save one every "itSteps" in t, default is 6*/
	
	fprintf(lhapdf, "PdfType: central\n");
	fprintf(lhapdf, "Format: lhagrid1\n");
	fprintf(lhapdf, "---\n");
	for(i=0;i<lr;i++){
		fprintf(lhapdf, "%.10e ", grid[rgrid[i]]); /* we use 10 significant digits to save the x-grid, this is needed to describe well the last few values, which are very close to 1 */
	}
	fprintf(lhapdf, "\n");
	for(i=iw;i<Nt0-1;i+=itSteps){
		fprintf(lhapdf, "%e	", qgrid[i]);
	}
	for(i=0;i<Nt;i+=itSteps){
		fprintf(lhapdf, "%e	", qgrid[Nt0+i-1]);
	}
	fprintf(lhapdf, "\n");
	fprintf(lhapdfbar, "PdfType: central\n");
	fprintf(lhapdfbar, "Format: lhagrid1\n");
	fprintf(lhapdfbar, "---\n");
	for(i=0;i<lr;i++){
		fprintf(lhapdfbar, "%.10e	", grid[rgrid[i]]); /* we use 10 significant digits to save the x-grid, this is needed to describe well the last few values, which are very close to 1 */
	}
	fprintf(lhapdfbar, "\n");
	for(i=iw;i<Nt0-1;i+=itSteps){
		fprintf(lhapdfbar, "%e	", qgrid[i]);
	}
	for(i=0;i<Nt;i+=itSteps){
		fprintf(lhapdfbar, "%e	", qgrid[Nt0+i-1]);
	}
	fprintf(lhapdfbar, "\n");
	if(FS==5){
		fprintf(lhapdf, "eL	eR	nue	muL	muR	numu	taL	taR	nuta	eLb	eRb	nueb	muLb	muRb	numub	taLb	taRb	nutab	dL	dR	uL	uR	sL	sR	cL	cR	bL	bR	dLb	dRb	uLb	uRb	sLb	sRb	cLb	cRb	bLb	bRb	gp	gm	gap	gam	Zp	Zm	ZL	Zgap	Zgam	Wpp	Wpm	WpL	Wmp	Wmm	WmL	h	hZL\n");
		fprintf(lhapdf, "11	11	12	13	13	14	15	15	16	-11	-11	-12	-13	-13	-14	-15	-15	-16	1	1	2	2	3	3	4	4	5	5	-1	-1	-2	-2	-3	-3	-4	-4	-5	-5	21	21	22	22	23	23	23	2223	2223	24	24	24	-24	-24	-24	25	2523\n");
		fprintf(lhapdf, "-	+	-	-	+	-	-	+	-	+	-	+	+	-	+	+	-	+	-	+	-	+	-	+	-	+	-	+	+	-	+	-	+	-	+	-	+	-	+	-	+	-	+	-	0	+	-	+	-	0	+	-	0	0	0\n");
		fprintf(lhapdfbar, "eL	eR	nue	muL	muR	numu	taL	taR	nuta	eLb	eRb	nueb	muLb	muRb	numub	taLb	taRb	nutab	dL	dR	uL	uR	sL	sR	cL	cR	bL	bR	dLb	dRb	uLb	uRb	sLb	sRb	cLb	cRb	bLb	bRb	gp	gm	gap	gam	Zp	Zm	ZL	Zgap	Zgam	Wpp	Wpm	WpL	Wmp	Wmm	WmL	h	hZL\n");
		fprintf(lhapdfbar, "11	11	12	13	13	14	15	15	16	-11	-11	-12	-13	-13	-14	-15	-15	-16	1	1	2	2	3	3	4	4	5	5	-1	-1	-2	-2	-3	-3	-4	-4	-5	-5	21	21	22	22	23	23	23	2223	2223	24	24	24	-24	-24	-24	25	2523\n");
		fprintf(lhapdfbar, "-	+	-	-	+	-	-	+	-	+	-	+	+	-	+	+	-	+	-	+	-	+	-	+	-	+	-	+	+	-	+	-	+	-	+	-	+	-	+	-	+	-	+	-	0	+	-	+	-	0	+	-	0	0	0\n");
		
		if(lepton==0){
			printf("Electron 5FS: ");
			for(i=0;i<lr;i++){
				for(j=iw;j<Nt0-1;j+=itSteps){
					write = e5FS0(rgrid[i],j,lhapdf,grid,fe,fm,fu,fd,fb,fph,fg,fc,ftau);
				}
				for(j=0;j<Nt;j+=itSteps){
					write = e5FS(rgrid[i],j,lhapdf,grid,fel,felb,fml,fnue,fnueb,fnum,fer,ferb,fmr,ful,fulb,ftl,ftlb,fdl,fdlb,fbl,fblb,fur,furb,ftr,ftrb,fdr,fdrb,fbr,fbrb,fh,fwpl,fwml,fzl,fwpp,fwpm,fwmp,fwmm,fzp,fzm,fphp,fphm,fzphp,fzphm,fgp,fgm,fhzl,ftal,ftalb,ftar,ftarb,fnuta,fnutab,fcl,fclb,fcr,fcrb,fsl,fslb);
				}
			}
			printf("done!\n");
			printf("Positron 5FS: ");
			for(i=0;i<lr;i++){
				for(j=iw;j<Nt0-1;j+=itSteps){
					write = eb5FS0(rgrid[i],j,lhapdfbar,grid,fe,fm,fu,fd,fb,fph,fg,fc,ftau);
				}
				for(j=0;j<Nt;j+=itSteps){
					write = eb5FS(rgrid[i],j,lhapdfbar,grid,fel,felb,fml,fnue,fnueb,fnum,fer,ferb,fmr,ful,fulb,ftl,ftlb,fdl,fdlb,fbl,fblb,fur,furb,ftr,ftrb,fdr,fdrb,fbr,fbrb,fh,fwpl,fwml,fzl,fwpp,fwpm,fwmp,fwmm,fzp,fzm,fphp,fphm,fzphp,fzphm,fgp,fgm,fhzl,ftal,ftalb,ftar,ftarb,fnuta,fnutab,fcl,fclb,fcr,fcrb,fsl,fslb);
				}
			}
			printf("done!\n");
		}
		else{
			printf("Muon 5FS: ");
			for(i=0;i<lr;i++){
				for(j=iw;j<Nt0-1;j+=itSteps){
					write = mu5FS0(rgrid[i],j,lhapdf,grid,fe,fm,fu,fd,fb,fph,fg,fc,ftau);
				}
				for(j=0;j<Nt;j+=itSteps){
					write = mu5FS(rgrid[i],j,lhapdf,grid,fel,felb,fml,fnue,fnueb,fnum,fer,ferb,fmr,ful,fulb,ftl,ftlb,fdl,fdlb,fbl,fblb,fur,furb,ftr,ftrb,fdr,fdrb,fbr,fbrb,fh,fwpl,fwml,fzl,fwpp,fwpm,fwmp,fwmm,fzp,fzm,fphp,fphm,fzphp,fzphm,fgp,fgm,fhzl,ftal,ftalb,ftar,ftarb,fnuta,fnutab,fcl,fclb,fcr,fcrb,fsl,fslb);
				}
			}
			printf("done!\n");
			printf("Antimuon 5FS: ");
			for(i=0;i<lr;i++){
				for(j=iw;j<Nt0-1;j+=itSteps){
					write = mub5FS0(rgrid[i],j,lhapdfbar,grid,fe,fm,fu,fd,fb,fph,fg,fc,ftau);
				}
				for(j=0;j<Nt;j+=itSteps){
					write = mub5FS(rgrid[i],j,lhapdfbar,grid,fel,felb,fml,fnue,fnueb,fnum,fer,ferb,fmr,ful,fulb,ftl,ftlb,fdl,fdlb,fbl,fblb,fur,furb,ftr,ftrb,fdr,fdrb,fbr,fbrb,fh,fwpl,fwml,fzl,fwpp,fwpm,fwmp,fwmm,fzp,fzm,fphp,fphm,fzphp,fzphm,fgp,fgm,fhzl,ftal,ftalb,ftar,ftarb,fnuta,fnutab,fcl,fclb,fcr,fcrb,fsl,fslb);
				}
			}
			printf("done!\n");
		}
	}
	else{
		fprintf(lhapdf, "eL	eR	nue	muL	muR	numu	taL	taR	nuta	eLb	eRb	nueb	muLb	muRb	numub	taLb	taRb	nutab	dL	dR	uL	uR	sL	sR	cL	cR	bL	bR	tL	tR	dLb	dRb	uLb	uRb	sLb	sRb	cLb	cRb	bLb	bRb	tLb	tRb	gp	gm	gap	gam	Zp	Zm	ZL	Zgap	Zgam	Wpp	Wpm	WpL	Wmp	Wmm	WmL	h	hZL\n");
		fprintf(lhapdf, "11	11	12	13	13	14	15	15	16	-11	-11	-12	-13	-13	-14	-15	-15	-16	1	1	2	2	3	3	4	4	5	5	6	6	-1	-1	-2	-2	-3	-3	-4	-4	-5	-5		-6	-6	21	21	22	22	23	23	23	2223	2223	24	24	24	-24	-24	-24	25	2523\n");
		fprintf(lhapdf, "-	+	-	-	+	-	-	+	-	+	-	+	+	-	+	+	-	+	-	+	-	+	-	+	-	+	-	+	-	+	+	-	+	-	+	-	+	-	+	-	+	-	+	-	+	-	+	-	0	+	-	+	-	0	+	-	0	0	0\n");
		fprintf(lhapdfbar, "eL	eR	nue	muL	muR	numu	taL	taR	nuta	eLb	eRb	nueb	muLb	muRb	numub	taLb	taRb	nutab	dL	dR	uL	uR	sL	sR	cL	cR	bL	bR	tL	tR	dLb	dRb	uLb	uRb	sLb	sRb	cLb	cRb	bLb	bRb	tLb	tRb	gp	gm	gap	gam	Zp	Zm	ZL	Zgap	Zgam	Wpp	Wpm	WpL	Wmp	Wmm	WmL	h	hZL\n");
		fprintf(lhapdfbar, "11	11	12	13	13	14	15	15	16	-11	-11	-12	-13	-13	-14	-15	-15	-16	1	1	2	2	3	3	4	4	5	5	6	6	-1	-1	-2	-2	-3	-3	-4	-4	-5	-5		-6	-6	21	21	22	22	23	23	23	2223	2223	24	24	24	-24	-24	-24	25	2523\n");
		fprintf(lhapdfbar, "-	+	-	-	+	-	-	+	-	+	-	+	+	-	+	+	-	+	-	+	-	+	-	+	-	+	-	+	-	+	+	-	+	-	+	-	+	-	+	-	+	-	+	-	+	-	+	-	0	+	-	+	-	0	+	-	0	0	0\n");
		
		if(lepton==0){
			printf("Electron 6FS: ");
			for(i=0;i<lr;i++){
				for(j=iw;j<Nt0-1;j+=itSteps){
					write = e6FS0(rgrid[i],j,lhapdf,grid,fe,fm,fu,fd,fb,fph,fg,fc,ftau);
				}
				for(j=0;j<Nt;j+=itSteps){
					write = e6FS(rgrid[i],j,lhapdf,grid,fel,felb,fml,fnue,fnueb,fnum,fer,ferb,fmr,ful,fulb,ftl,ftlb,fdl,fdlb,fbl,fblb,fur,furb,ftr,ftrb,fdr,fdrb,fbr,fbrb,fh,fwpl,fwml,fzl,fwpp,fwpm,fwmp,fwmm,fzp,fzm,fphp,fphm,fzphp,fzphm,fgp,fgm,fhzl,ftal,ftalb,ftar,ftarb,fnuta,fnutab,fcl,fclb,fcr,fcrb,fsl,fslb);
				}
			}
			printf("done!\n");
			printf("Positron 6FS: ");
			for(i=0;i<lr;i++){
				for(j=iw;j<Nt0-1;j+=itSteps){
					write = eb6FS0(rgrid[i],j,lhapdfbar,grid,fe,fm,fu,fd,fb,fph,fg,fc,ftau);
				}
				for(j=0;j<Nt;j+=itSteps){
					write = eb6FS(rgrid[i],j,lhapdfbar,grid,fel,felb,fml,fnue,fnueb,fnum,fer,ferb,fmr,ful,fulb,ftl,ftlb,fdl,fdlb,fbl,fblb,fur,furb,ftr,ftrb,fdr,fdrb,fbr,fbrb,fh,fwpl,fwml,fzl,fwpp,fwpm,fwmp,fwmm,fzp,fzm,fphp,fphm,fzphp,fzphm,fgp,fgm,fhzl,ftal,ftalb,ftar,ftarb,fnuta,fnutab,fcl,fclb,fcr,fcrb,fsl,fslb);
				}
			}
			printf("done!\n");
		}
		else{
			printf("Muon 6FS: ");
			for(i=0;i<lr;i++){
				for(j=iw;j<Nt0-1;j+=itSteps){
					write = mu6FS0(rgrid[i],j,lhapdf,grid,fe,fm,fu,fd,fb,fph,fg,fc,ftau);
				}
				for(j=0;j<Nt;j+=itSteps){
					write = mu6FS(rgrid[i],j,lhapdf,grid,fel,felb,fml,fnue,fnueb,fnum,fer,ferb,fmr,ful,fulb,ftl,ftlb,fdl,fdlb,fbl,fblb,fur,furb,ftr,ftrb,fdr,fdrb,fbr,fbrb,fh,fwpl,fwml,fzl,fwpp,fwpm,fwmp,fwmm,fzp,fzm,fphp,fphm,fzphp,fzphm,fgp,fgm,fhzl,ftal,ftalb,ftar,ftarb,fnuta,fnutab,fcl,fclb,fcr,fcrb,fsl,fslb);
				}
			}
			printf("done!\n");
			printf("Antimuon 6FS: ");
			for(i=0;i<lr;i++){
				for(j=iw;j<Nt0-1;j+=itSteps){
					write = mub6FS0(rgrid[i],j,lhapdfbar,grid,fe,fm,fu,fd,fb,fph,fg,fc,ftau);
				}
				for(j=0;j<Nt;j+=itSteps){
					write = mub6FS(rgrid[i],j,lhapdfbar,grid,fel,felb,fml,fnue,fnueb,fnum,fer,ferb,fmr,ful,fulb,ftl,ftlb,fdl,fdlb,fbl,fblb,fur,furb,ftr,ftrb,fdr,fdrb,fbr,fbrb,fh,fwpl,fwml,fzl,fwpp,fwpm,fwmp,fwmm,fzp,fzm,fphp,fphm,fzphp,fzphm,fgp,fgm,fhzl,ftal,ftalb,ftar,ftarb,fnuta,fnutab,fcl,fclb,fcr,fcrb,fsl,fslb);
				}
			}
			printf("done!\n");
		}
	}
	
	/* ----------------------- Freeing memory allocation ----------------------- */
	freemat(fe,Nt0);
	freemat(fm,Nt0);
	freemat(ftau,Nt0);
	freemat(fu,Nt0);
	freemat(fc,Nt0);
	freemat(fd,Nt0);
	freemat(fb,Nt0);
	freemat(fph,Nt0);
	freemat(fg,Nt0);
	freemat(fel,Nt);
	freemat(felb,Nt);
	freemat(fer,Nt);
	freemat(ferb,Nt);
	freemat(fml,Nt);
	freemat(fmr,Nt);
	freemat(ftal,Nt);
	freemat(ftalb,Nt);
	freemat(ftar,Nt);
	freemat(ftarb,Nt);
	freemat(fnue,Nt);
	freemat(fnueb,Nt);
	freemat(fnum,Nt);
	freemat(fnuta,Nt);
	freemat(fnutab,Nt);
	freemat(ful,Nt);
	freemat(fulb,Nt);
	freemat(fur,Nt);
	freemat(furb,Nt);
	freemat(fcl,Nt);
	freemat(fclb,Nt);
	freemat(fcr,Nt);
	freemat(fcrb,Nt);
	freemat(fdl,Nt);
	freemat(fdlb,Nt);
	freemat(fdr,Nt);
	freemat(fdrb,Nt);
	freemat(fsl,Nt);
	freemat(fslb,Nt);
	freemat(ftl,Nt);
	freemat(ftlb,Nt);
	freemat(ftr,Nt);
	freemat(ftrb,Nt);
	freemat(fbl,Nt);
	freemat(fblb,Nt);
	freemat(fbr,Nt);
	freemat(fbrb,Nt);
	freemat(fh,Nt);
	freemat(fwpl,Nt);
	freemat(fwml,Nt);
	freemat(fzl,Nt);
	freemat(fwpp,Nt);
	freemat(fwpm,Nt);
	freemat(fwmp,Nt);
	freemat(fwmm,Nt);
	freemat(fzp,Nt);
	freemat(fzm,Nt);
	freemat(fphp,Nt);
	freemat(fphm,Nt);
	freemat(fzphp,Nt);
	freemat(fzphm,Nt);
	freemat(fgp,Nt);
	freemat(fgm,Nt);
	freemat(fhzl,Nt);
	freemat(k01,Nt0-1);
	freemat(k02,Nt0-1);
	freemat(k03,Nt0-1);
	freemat(k04,Nt0-1);
	freemat(k1,Nt-1);
	freemat(k2,Nt-1);
	freemat(k3,Nt-1);
	freemat(k4,Nt-1);
	freemat(talpha,Na);
	free(grid);
	free(deltagrid);
	free(rgrid);
	
	/* ----------------------- The End ----------------------- */
	printf("Calculation completed. The End.\n");
	return 0;	
}


