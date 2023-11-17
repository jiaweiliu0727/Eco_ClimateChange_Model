double DDEVF(void *Paramstuff,gsl_rng *RandNumsPass,size_t dim,int pop,int maxy_t, double hatch, int q)
{
// DDEVF sets up model and calls 0DE_SOLVER, returns Params->sim_results
//printf("DDEVF: population:%d\n",pop);

STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
		//declarations of stuff for plotting
double PLOT=0.0;

double Fprob;
double p1,p2,p3;	        // simulation results (p_i is the probability of being in each of the three classes)
double Cpar, Cprob, Cprob2, Cprob3;					//CK//
double Opar, Oprob, Oprob2, Oprob3;					//CK//
double r1, r2, r3, c1, c2;	            //CK// simulation results (r1 is resting spore density, c1 is condidia density)
r1 = 0.0; r2=0.0; c1= 0.0; c2=0.0;

double cover_C = Params->PARS[17];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
double cover_R = Params->PARS[20];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
//double cover = Params->PARS[20];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
double open_C = Params->PARS[24];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
double open_R = Params->PARS[25];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals

int FlagDay=0;
int MAXT3=maxy_t;


int DIM = Params->PARS[9]+3;   //number of ODE equations
int m = Params->PARS[9];

double t=h;		double t_next=h;	double t_0=h;	int i;				// time loop and index
double epsilon = pow(10,-6);
double y_ode[DIM];
double rand_nuR[MAXT3];
double rand_nuF[MAXT3];

double ave_R = Params->PARS[50+pop];   //average R(0) for all the sites
double specific_muF = Params->PARS[6];   //general intercept for MAX TEMP decay function for conidia
double specific_nuF = Params->PARS[3];   //Site-specific infection rate for conidia
//double R_end = Params->PARS[16];   //CK//  Flat time for end of resting spore blooming
double rain_P = Params->PARS[21];  //fit param used to scale accumulating rain.
double rain_P2 = Params->PARS[26];  //fit param used to scale accumulating rain.
double rain_P3 = Params->PARS[29];  //fit param used to scale accumulating rain.
double RH_P = Params->PARS[22];   //CK//  Parameter for RH data
double temp_P = Params->PARS[23];   //CK//  Parameter for temperature data
double fourth_size=Params->PARS[28];	//CK// degree day when the bugs reach 4th instar


double DDstart	= Params->PARS[27];  //CK// param used for starting date
double DDstop	= Params->PARS[19];  //CK// param used for stopping date

Params->size_C = 1.0;
double C_end=Params->PARS[16];	  //CK// fit param that turns off new conidia production once a specific size has been reached

double temp_now;  //CIK// used simplify decay functions
double total_rainfall=0.0;  //used to sum up rainfall
int rain_day;  //used to sum up rainfall
int beta;		//used to dictate how many days to go back when accumulating rain
int theta;    //used to determine lag period before calculating accumulated rainfall
//beta=Params->PARS[40+pop];
beta=Params->PARS[7];
//beta=1;
theta=1;
//theta=Params->PARS[13];
//printf("Beta: %d Theta: %d\n", beta, theta);		getc(stdin);


double DDtemp_now;  //CK// used simplify decay functions

double nuF2;  //Transmission rate of conidia
double nuR2;  //Transmission rate of resting spores
double FIO_Cc;
double FIO_Cr;
double FIO_Oc;
double FIO_Or;

double DD10=0;    //accumulated degree days about 10 degrees C

//double specific_Rend = Params->PARS[40+pop];   //Site-specific decay rate for conidia

//printf("specific_nuF: %e\n", specific_muF);		getc(stdin);

// ----------------------------------- Generate Random Numbers -------------------------------------------- //
for (i=0;i<=MAXT3;i++)	{
    //JL: The stochasticity to change the transmission rates of conidia and resting spores vary every day.
	rand_nuR[i]=gsl_ran_gaussian(RandNumsPass,Params->PARS[11]); //2nd entry is stdev
	rand_nuF[i]=gsl_ran_gaussian(RandNumsPass,Params->PARS[12]);
	//rand_nuF[i]=gsl_cdf_gaussian_Pinv(RandNumsPass[i],Params->PARS[12]);
	//printf("i(%d)\t var1=%f\t var2=%f\t rand parts: nuR=%e\t nuF=%e\n",i,Params->PARS[11],Params->PARS[12],rand_nuR[i],rand_nuF[i]);
}//getc(stdin);
// ------------------------------------- Initial Conditions ---------------------------------------------- //
Params->INITS[0] = Params->PARS[30+pop];			// initS
Params->INITS[1] = 0.0;						// initV  SET TO 0 FOR FUNGUS ONLY MODEL!!
//Params->INITS[3] = Params->PARS[40+pop];								// initR
Params->INITS[3] = ave_R;												// initR  //CK// changed to use average R(0), not site-specific

double initS = Params->INITS[0];			// initS
double initV = Params->INITS[1];			// initV
double initR = Params->INITS[3];			// initR

// ----------------------------------------- Fixed Parameters ---------------------------------------------- //
int gstepsV		= (int) Params->PARS[8];	int gstepsF	= (int) Params->PARS[9];
double ratio	= Params->PARS[10];
double neo_v	= 7.0;			// latent period of neonates (days) FUNGUS ONLY MODEL!

double R_end;   //CK//  The day when resting spore germination ends
double R_start;   //CK//  The day when resting spore germination starts


Params->PARS[0]=1.0;


// ------------------------------------- initialize model parameters --------------------------------------- //
int FlagWeek;	int FlagV=0;	int FlagR=0;	int FlagR_end=0;		// keep track of end of week

double ConiBefore=0.0;    //CK// Thing to store conidia the day before the feral collection
double RestBefore=0.0;    //CK// Thing to store resting spores the day before the feral collection

int day = 0; int week = 0;							// keeps track of day and week number

int line_ticker=0;   //CK// Ticker used to associate t in function with numbered days.
int line_ticker2;
int test_day=0;	//CK// used to find the line in the weather data that corresponds to the starting day of collections
//int num_day =  Params->DATA[pop][0][4];  //CK// Starting day number

int num_day =  hatch;  //CK// Starting day number

//printf("starting day number: %d\n", num_day);		getc(stdin);

line_ticker = num_day;

line_ticker=line_ticker-1;
line_ticker2=line_ticker-1;

//printf("corresponding line in WDATA: %d\n", line_ticker);		getc(stdin);

double S,V,F,R;	double IV=0, IF=0;
double E_V[gstepsV+1]; double E_F[gstepsF+1];
int num_weeks=MAXT3/7;

// -----------------------------------//CK// calculating ending blooming times //CK//--------------------------------------- //

DD10=0.0;		R_start = 0.0;
test_day = hatch;

while(DD10 <= DDstart){

	DDtemp_now = Params->CCDATA[q][test_day][3]-10.0;  //CK// begin calculation of accumulated Degree Days
	if(DDtemp_now<0.0){DDtemp_now=0.0;}
	DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time
	R_start++;
	test_day++;
//printf("temp now: %f DD10: %f R_start: %f test_day: %d\n", DDtemp_now, DD10, R_start, test_day);		//getc(stdin);
}
//printf("temp now: %f DD10: %f R_start: %f test_day: %d\n", DDtemp_now, DD10, R_start, test_day);		//getc(stdin);


//R_start=0.0;

DD10=0.0;		R_end = 0.0;
test_day = hatch;

while(DD10 <= DDstop){

	DDtemp_now = Params->CCDATA[q][test_day][3]-10.0;  //CK// begin calculation of accumulated Degree Days
	if(DDtemp_now<0.0){DDtemp_now=0.0;}
	DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time
	R_end++;
	test_day++;
//printf("temp now: %f DD10: %f R_end: %f test_day: %d\n", DDtemp_now, DD10, R_end, test_day);		//getc(stdin);
}
//printf("temp now: %f DD10: %f R_end: %f test_day: %d\n", DDtemp_now, DD10, R_end, test_day);		//getc(stdin);

DD10=0.0;

//printf("Checking blooming params: pop=%d\t Rstart= %f\t lat=%f\t r_time= %f\t Rend= %f\n", pop,Params->DAY_F[pop],lats[pop-1], r_time, R_end);		getc(stdin);



//-----------------------------------//CK// Infections on day 0!!!! -------------------------------//
if(R_start<1.0){  //JL: Usually R_start>1.0, which means the starting days of collection are in advance of resting spore blooming.


	total_rainfall = 0.0;

	//printf("stochasticity: %f\t initR: %f\n", rand_nuF[0], initR);

	for (rain_day= (line_ticker - beta - theta - 1);rain_day <= line_ticker - theta -1;rain_day++){
		total_rainfall = Params->CCDATA[q][rain_day][0]+total_rainfall;
	//printf("line: %d\t accumulating rain: %lf\n", rain_day, total_rainfall);		getc(stdin);
	}

	nuR2=(rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall)) - rain_P/(rain_P2+1))*exp(rand_nuR[0]);
	Params->nuR = (DD10/fourth_size)*nuR2;

	r2=Params->nuR*initR;

	total_rainfall = 0.0;


	for (rain_day= (line_ticker2 - beta - theta - 1);rain_day <= line_ticker2 - theta -1;rain_day++){
		total_rainfall = Params->CCDATA[q][rain_day][0]+total_rainfall;
	//printf("line: %d\t accumulating rain: %lf\n", rain_day, total_rainfall);		getc(stdin);
	}

	nuR2=(rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall)) - rain_P/(rain_P2+1))*exp(rand_nuR[0]);
	Params->nuR = (DD10/fourth_size)*nuR2;


	c1 = 0.0;	c2 = 0.0;
	r1 = Params->nuR*initR;

	Cpar = cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2);	//covered cages
	Opar = open_C*((c1+c2)/2) + open_R*((r1+r2)/2);	//uncovered cages

	FIO_Cc = cover_C*((c1+c2)/2);
	FIO_Cr =  cover_R*((r1+r2)/2);
	FIO_Oc = open_C*((c1+c2)/2);
	FIO_Or = open_R*((r1+r2)/2);

	Cprob = exp(-Cpar);	Oprob = exp(-Opar);

//		printf("%d\t %d\t %f\t %f\t %f\t %f\t\n",pop,day,prob,1.0-prob, prob2, prob3); // getc(stdin);

	Cprob2 = (( cover_C*((c1+c2)/2))/(cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2)))*((1 - Cprob));  //OPEN fraction of F infection due to conidia
	Cprob3 = ((cover_R*((r1+r2)/2))/(cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2)))*((1 - Cprob));  //fraction of F infection due to resting spores

	//Cprob2 = 0.0;  //fraction of F infection due to conidia   COVERED CAGES!!!
	//Cprob3 = 1-Cprob;  //fraction of F infection due to resting spores

	Oprob2 = ((open_C*((c1+c2)/2))/(open_C*((c1+c2)/2) + open_R*((r1+r2)/2)))*((1 - Oprob));  //OPEN fraction of F infection due to conidia
	Oprob3 = ((open_R*((r1+r2)/2))/(open_C*((c1+c2)/2) + open_R*((r1+r2)/2)))*((1 - Oprob));  //fraction of F infection due to resting spores

	if(Cprob == 1.0){Cprob2=0.0; Cprob3=0.0;}
	if(Oprob == 1.0){Oprob2=0.0; Oprob3=0.0;}


	printf("%d\t %d\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",pop,0,0.0,0.0,0.0, Cprob,1.0-Cprob,Cprob2,Cprob3, Oprob,1-Oprob,Oprob2,Oprob3, FIO_Cc, FIO_Cr, FIO_Oc, FIO_Or); //getc(stdin);

	r2= r1;
}
else{
	//sim_results[0][5]=0.0; sim_results[0][6]=0.0;
	if(PLOT==1.0){printf("%d\t %d\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",pop,0,0.0,0.0,0.0, 1.0,0.0,0.0,0.0, 1.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0); //getc(stdin);
	}
}

// ---------------------------- Calculated Parameters  and Population Sizes ------------------------------- //
double Vstart = ratio*initV;						// viral cadavers after infected neonates die

//printf("Starting Virus Conditions: initV=%f\t ratio=%f\t Vstart=%f\n", initV, ratio, Vstart);	getc(stdin);

//double r_germ = Params->DAY_F[pop]-r_time;		// calculated day of resting spore germination
double r_germ = R_start;		//CK// nixing r_time because I made all germ dates start at beginning of the collections
if (r_germ<0)	r_germ=0;
if		(initS>Vstart)	{	S = initS-Vstart;	}	// S(0) (number of healthy neonates)
else					{	S = 0;				}
V=0.0; F=0.0; R=0.0;
//for (i=1;i<=gstepsV;i++)		{	E_V[i]=0;			}
for (i=1;i<=gstepsF;i++)		{	E_F[i]=0;			}
//double timing[6]={neo_v,r_germ,R_end,MAXT3};
double timing[6]={r_germ,R_end,MAXT3};
// ----------------------------------------- initialize results ------------------------------------------- //
double P[4] = {0,0,0,0};				// model results (0,healthy,viral kill,fungal kill)
double Vkill=0, Fkill=0;				// total individuals infected (used in calculating fraction infected)
double VirFI[100], FungFI[100];			// fraction of infected individals (by week number)
VirFI[0] = 0.0; VirFI[1] = 0.0; FungFI[0] = 0.0; FungFI[1] = 0.0;		// fraction infected first two weeks

// -------------------- MAIN LOOP!! (calculate populations as time is increased) -------------------------- //
//printf("pop:%d\t t_0=%f\t neo_v=%f\t r_germ=%f R_end=%f week_end=%d maxt=%f\n",pop,t_0,timing[0],timing[1],timing[2],7*(week+1),timing[3]);
//printf("begin DDEVF main loop\n");		getc(stdin);
while (t_0<MAXT3+h)	{
	//if (week>=num_weeks)	break;
	FlagWeek=0;
    //printf("t_0=%f\t day=%d\t week=%d\t num_weeks=%d max_time=%d\n",t_0,day,week,num_weeks,Params->MAXT[pop]);//getc(stdin);
	//printf("neo=%f\t rgerm=%f\t rend=%f\t day+1=%d\n",timing[0],timing[1],timing[2],day+1);	getc(stdin);
	// ------------------------------- Find Stoppage Event ----------------------------------------------- //
	// --------------------- end of day -------------------- //
	if (day+1<timing[0] && day+1<timing[1])	{

		FlagDay=1;

		t_next=day+1;
		day++;
		num_day++;
		line_ticker++;

		//printf("end of day: t=%f\n",t_next);

		if (day%7==0)	{
		  //printf("end of week!!: t=%f\n",t_next);
			week++;
			FlagWeek=1;
		}
		//getc(stdin);
	}

	// --------------------- resting spores bloom -------------------- //
	else if (timing[0]<timing[1] && timing[0]<day+1)	{
		FlagR=1;
		t_next=timing[0];
		timing[0]=999.9;
		//printf("Resting spores bloom: t=%f\n",t_next);		//getc(stdin);
	}
	// --------------------- resting spores done -------------------- //
	else if (timing[1]<timing[0] && timing[1]<day+1)	{
		FlagR_end=1;
		t_next=timing[1];
		timing[0]=999.9;
		timing[1]=999.9;
		//printf("Resting spores done: t=%f\n",t_next);		//getc(stdin);
	}

	// --------------------- resting spores bloom and end of day -------------------- //
	else if (abs(day+1-timing[0])<epsilon)	{

		FlagDay=1;

		FlagR=1;
		t_next=day+1;
		timing[0]=999.9;
		day++;
		num_day++;
		line_ticker++;

		if (day%7==0)	{
			week++;
			FlagWeek=1;
		}
		//printf("end of day and Resting Spores bloom: t=%f\n",t_next);	//getc(stdin);
	}
	// --------------------- resting spores done and end of day -------------------- //
	else if (abs(day+1-timing[1])<epsilon)	{

		FlagDay=1;

		FlagR_end=1;
		t_next=day+1;
		timing[0]=999.9;
		timing[1]=999.9;
		day++;
		num_day++;
		line_ticker++;

		if (day%7==0)	{
			week++;
			FlagWeek=1;
		}
		//printf("end of day and Resting Spores done: t=%f\n",t_next);	//getc(stdin);
	}
	else {
		printf("ERROR: NO EVENT IS NEXT IN TIME???\n");
		getc(stdin);
	}
	// -------------------------- integrate until next stoppage event ---------------------------------- //
	while (t<t_next)	{
		y_ode[0]=S;	y_ode[m+1]=F;	Params->POPS[3]=R;
		//printf("pop:%d t=%f\t S=%4.3e\t V=%4.3e\t F=%4.3e\t R=%4.3e\n",pop,t,S,V,F,R);

		for (i=1;i<=gstepsF;i++)	{
			//y_ode[1+i]=E_V[i];
			y_ode[i]=E_F[i];
		}
		y_ode[m+2]=Fkill;
		//printf("DDEVF:\n parm 2=%f parm 3=%f parm 4=%f parm 5=%f parm 6=%f\n",Params->PARS[2],Params->PARS[3],Params->PARS[4],Params->PARS[5],Params->PARS[6]);
		//getc(stdin);

		//Params->nuF = specific_nuF;
		//Params->nuR = initR;
		Params->nuV = Params->PARS[2];
		//Params->muF = specific_muF;

		DDtemp_now = Params->CCDATA[q][line_ticker - 1][3]-10.0;  //CK// begin calculation of accumulated Degree Days
		if(DDtemp_now<0.0){DDtemp_now=0.0;}
		DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time

		if(DD10>=C_end){Params->size_C=0.0;}	//CK// stops new conidia production once a fit size has been reached.


		//Params->nuF = specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6]);
		//Params->nuF = specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[(int)t]);
		Params->nuF = (DD10/fourth_size)*specific_nuF*exp(RH_P*Params->CCDATA[q][line_ticker - 1][1]) * exp(rand_nuF[(int)t]);
		nuF2 = specific_nuF*exp(RH_P*Params->CCDATA[q][line_ticker - 1][1]) * exp(rand_nuF[(int)t]);

		temp_now = Params->CCDATA[q][line_ticker - 1][2];  //CK// putting max temp in smaller object

		//Params->muF = specific_muF;	//CK// Conidia Decay Response #2.2  BEST SO FAR!!
		//Params->muF = specific_muF*temp_now;	//CK// Conidia Decay Response #2.2  BEST SO FAR!!
		Params->muF = specific_muF*exp(temp_P*temp_now);	//CK// Conidia Decay Response #2.2  BEST SO FAR!!

		total_rainfall = 0.0;
		rain_day= line_ticker - beta - 1;
//printf("Start rain calc! day: %d\n", line_ticker-1);		getc(stdin);
		for (rain_day= (line_ticker - beta - theta - 1);rain_day < line_ticker - theta -1;rain_day++)	{
			total_rainfall = Params->CCDATA[q][rain_day][0]+total_rainfall;
//printf("line: %d\t accumulating rain: %lf\n", rain_day, total_rainfall);		getc(stdin);
		}

		nuR2=(rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall)) - rain_P/(rain_P2+1))*exp(rand_nuR[(int)t]);
		Params->nuR = (DD10/fourth_size)*nuR2;

		if(Params->POPS[3] == 0.0){Params->nuR = 0.0; nuR2 = 0.0;}


		t=ODE_Solver(t,t_next,Params,y_ode); //Numerical integration from t to t_next, y_ode holds group densities
		//Update the values for each class//
		S=y_ode[0];	F=y_ode[m+1]; //set group densities equal to ode output
		IF=0;

//printf("pop: %d day: %d conidia: %f\n", pop, day, y_ode[7]*Params->nuF);



	if ((day+1)%7==0)	{
		ConiBefore=y_ode[m+1]*nuF2;  //CK// Saving the conidia 24 hours before feral collection
		RestBefore = Params->POPS[3]*nuR2;
		//ConiBefore=y_ode[m+1]*Params->nuF;  //CK// Saving the conidia 24 hours before feral collection
		//RestBefore = Params->POPS[3]*Params->nuR;
		//FlagConidia=2;
	}
		//getc(stdin);

		for (i=1;i<=gstepsF;i++)	{
			//E_V[i]=y_ode[1+i];
			E_F[i]=y_ode[i];
			//IV += E_V[i];
			IF += E_F[i];
			//printf("E_V(%d)=%e\t",i,E_V[i]);	//printf("IV=%f\n",IV);
		}//printf("\n");
		if(IF>initS)		IF=initS;
		Fkill=y_ode[m+2];			// killed populations
	}
	if (FlagV==1)			{	//printf("VIRAL INFECTED NEONATES DIE: t=%f\n",t_next);		//getc(stdin);
		V=Vstart;
		FlagV=2;
	}
	if (FlagR==1)			{	//printf("Resting spores bloom: t=%f\n",t_next);			//	getc(stdin);
		R=initR;
		FlagR=2;
	}

	else if (FlagR_end==1)	{	//printf("Resting spores done: t=%f\n",t_next);				// getc(stdin);
		R=0;
		FlagR_end=2;
	}

//Add daily plotting output here
	if (FlagDay==1)	{
         //printf("%d\t %d\t %d\t %d\t %d\t %e\t %e\n",spot,year,repeat,pop,day,y_ode[0], y_ode[m+1]);

		if ((day+1)%7==0)	{
			ConiBefore=y_ode[m+1]*nuF2;  //CK// Saving the conidia 24 hours before feral collection
			RestBefore = Params->POPS[3]*nuR2;
			//FlagConidia=2;
		}
	if(PLOT==1.0){
		if ((S+IF)==0)	{
			Fprob=0;
		}
		else					{
			//VirFI[week]  = IV/(S+IV+IF);		// fraction infected at the end of each week
			Fprob = IF/(S+IF);
		}
		p1 = y_ode[m+1]*Params->nuF;
		p2 = Params->POPS[3]*Params->nuR;
		//p1 = 1-(Fprob);
		//p2 = 0;  			//CK//	formerly virus infected
		p3 = Fprob;

//		printf("%d\t %d\t %f\t %f\t %f\n",pop,i,p1,p2,p3); //getc(stdin);

		c1 = y_ode[m+1]*nuF2;
		r1 = Params->POPS[3]*nuR2;

		Cpar = cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2);	//covered cages
		Opar = open_C*((c1+c2)/2) + open_R*((r1+r2)/2);	//uncovered cages

		FIO_Cc = cover_C*((c1+c2)/2);
		FIO_Cr =  cover_R*((r1+r2)/2);
		FIO_Oc = open_C*((c1+c2)/2);
		FIO_Or = open_R*((r1+r2)/2);

		Cprob = exp(-Cpar);	Oprob = exp(-Opar);

//		printf("%d\t %d\t %f\t %f\t %f\t %f\t\n",pop,day,prob,1.0-prob, prob2, prob3); // getc(stdin);

		Cprob2 = (( cover_C*((c1+c2)/2))/(cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2)))*((1 - Cprob));  //OPEN fraction of F infection due to conidia
		Cprob3 = ((cover_R*((r1+r2)/2))/(cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2)))*((1 - Cprob));  //fraction of F infection due to resting spores

		//Cprob2 = 0.0;  //fraction of F infection due to conidia   COVERED CAGES!!!
		//Cprob3 = 1-Cprob;  //fraction of F infection due to resting spores

		Oprob2 = ((open_C*((c1+c2)/2))/(open_C*((c1+c2)/2) + open_R*((r1+r2)/2)))*((1 - Oprob));  //OPEN fraction of F infection due to conidia
		Oprob3 = ((open_R*((r1+r2)/2))/(open_C*((c1+c2)/2) + open_R*((r1+r2)/2)))*((1 - Oprob));  //fraction of F infection due to resting spores

		if(Cprob == 1.0){Cprob2=0.0; Cprob3=0.0;}
		if(Oprob == 1.0){Oprob2=0.0; Oprob3=0.0;}

		//printf("%d\t %d\t %f\t %f\t %f\t %f\n",pop, day, c1, c2, r1, r2);

		printf("%d\t %d\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",pop,day,p1,p2,p3, Cprob,1.0-Cprob,Cprob2,Cprob3, Oprob,1-Oprob,Oprob2,Oprob3, FIO_Cc, FIO_Cr, FIO_Oc, FIO_Or); //getc(stdin);

		c2=c1;	r2=r1;   //make today's C and R yesterday's C and R
	}
		FlagDay=0;
	}
//Done with daily plotting output here


	// ---------------------------- end of the week updates ---------------------------------- //

	t_0=t_next;
}

Params->survivors += y_ode[0]; //JL: Update the cumulative number of hosts survived after the epizootic. Used to calculate the average fraction infected at this location.
Params->survive = y_ode[0]; //JL: Update the number of hosts survived after the epizootic at the certain location at the certain year, for this specific realization. Used in recording the exact values in each realization.

//printf("Resting spores day -1: %f\n", sim_results[0][5]);

//printf("end DDEVF\n");	getc(stdin);
return 0;
}
