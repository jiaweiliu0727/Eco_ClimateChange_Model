#include "headCC_Multi_avevar.h"

//JL: The general structure of code is similar to the code structure in Kyle et al. 2020.
//See Supplemental data at https://www.journals.uchicago.edu/doi/suppl/10.1086/707138 for the previous version of code in model fitting and parameterization.
int main(int argc, char *argv[])
{
//int test = atoi(argv[1]);                  // test is 99 for checking, 0 otherwise

//int test = 99;
int test = 66; //CK// Second test mode.  Like the full program but less runs and less MISER calls

int View = 0;  //CK// turn to 1 to print line search progress.  turn to 0 to run for real.

STRUCTURE Params;
int pro = 1;//atoi(argv[1]);						// pro and argv[1] are the inputs (argv[i] is the i^th input)
//printf("Profile Parameter is %d\n",pro);	fflush(stdout);
// ------------------------------------- Adustable accuracy vs. speed ------------------------------------------------ //
int num_runs	 = 20;
double parm_inc, host_inc, initR_inc;	//int inc_gamma_box= 1;

//if (pro==1)	{	parm_inc=200.0;		host_inc=100.0;	initR_inc=100.0;	}
if (pro==1)	{	parm_inc=50.0;		host_inc=10.0;	initR_inc=15.0;	}
//if (pro==1)	{	parm_inc=34.0;		host_inc=18.0;	initR_inc=20.0;	}
else		{	parm_inc=15.0;		host_inc=10.0;	initR_inc=10.0;	}

//if (test==66)	{	printf("for checking CK MODE2!!!\n");        num_runs=5;} //CK// New test mode!

if (test==66)	{      num_runs=5;	parm_inc=12.0;	host_inc=6.0;	initR_inc=8.0;}

//JL: Best parameters of the full model, from Kyle et al. 2020
double bestparams[30]={1.0, 1.0, 0.000000e+00, 0.000241071699421562, 0, 0, 0.00962435749864498, 10, 5, 50, 0, 0.37161719994828, 0.913699399999732, 2.2223804999527, 0.945435549999967, 0, 525.015699999847, 8.32036899999904, 0.119701349994476, 267.034499999981, 7.88482749903281, 3.80285399989692, 0.070488499999861, 0.233982799999915, 7.05116449999956, 6.38002749970359, 3.54725448752468, 100.157149999888, 291.2745, 0.166585199947054};

//Initial host densities at each location (larvae/m^2)
double initialS[9] = {25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0};

//Initial resting spore densities (resting spores/m^2), taken from Kyle et al. 2020. The 1st to 3rd values represent low to high fungus densities.
double initialR[9] = {0.000691100449204947, 0.00571302599998285, 0.018277349998676,0.000691100449204947, 0.00571302599998285, 0.018277349998676,0.000691100449204947, 0.00571302599998285, 0.018277349998676}; //sites:6,2,1,6,2,1,6,2,1



// ------------------------------------------------------------------------------------------------------------------ //
int i=0; int j;int ii; int jj; int k; int q = 0; int qq=0;
//int run;	            int changer;	    double index, tot_index;

int num_adj_pars=29;			// number of adjustable parameters


int pop;
//double log_pop;
//Params.th_id=0;
// -------------------------------------------- MISTER STUFF --------------------------------------------------------- //
inputdata(&Params);				// gets Params->CCDATA[0-9][0-19799][0-3] from inputdata.h

int calls=50;					// number of stochastic simulations for each parameter and IC set
if (test==99)	calls=50;
if (test==66)	calls=2; //CK// second test mode
//if (test==66)	calls=5; //CK// second test mode
size_t dim;
//size_t dim2;
// --------------------------------------- Name for Output Files ----------------------------------------------------- //
char strFileName[129];					// from filenames.h; name for the output file of mean of fraction infected at each location
GetString(pro,0,strFileName,128);		fflush(stdout);		//getc(stdin);
FILE *fp_results;

char strFileNameV[129]; //Name for the output file of variance of fraction infected at each location
GetStringV(pro,0,strFileNameV,128);       fflush(stdout);		//getc(stdin);
FILE *fp_results_v;
// ---------------------------------------- Random Number Stuff ------------------------------------------------------ //
gsl_rng *r_seed;
r_seed=random_setup();
//printf("Random Seed: %f\n", r_seed); //getc(stdin);
// -------------------- parameter high/low values and increments and fixed parameter values ------------------------- //
global_fixed_parms(&Params);  // gets Params.PARS[i] for fixed parameters from bounds.h
parm_range_inc(&Params,parm_inc,host_inc,initR_inc,num_adj_pars); // gets Params.parm_set,low,high,R_END from bounds.h
// ------------------------------------ Declare Likelidhood Quanitites ----------------------------------------------- //
double total_lhood;						// sum of pop_best_lhood over all patches


// -----------------------Declaring things for CC simulation------------------------------------ //

double DD10=0.0;		double S_start = 0.0;	double S_end = 0.0;
int test_day;
//int line_ticker;
double DDtemp_now;
double hatch =317.0; //Degree day needed for larvae hatch, from Russo et al. 1993
double Hlim1 = 3.0;
double Hlim2 = 40.0-Hlim1;
double pupate = 586.5;  //Degree day needed for pupation, from Carter et al. 1992
double Plim1=7.65;
double Plim2=41.0-Plim1;
double MAXT3;
int limit;


int reps;	//number of stochastic simulations to do per
double INFECTED[10]={0,0,0,0,0,0,0,0,0,0};
double sumINFECTED=0;
double aveINFECTED=0;
double varINFECTED=0;

// ------------------------------ OUTER MAIN LOOP (used for profile lhood ------------------------------------------ //
//for (outer_parm=Params.parm_low[pro]; outer_parm<=Params.parm_high[pro]; outer_parm+=Params.parm_step[pro])	{ //pro=20 for S(0)
//printf("profile=%d\t low=%f\t high=%f\t step=%f\n",pro,Params.parm_low[pro],Params.parm_high[pro],Params.parm_step[pro]);
//if      (pro==1) printf("max lhood!\n");
//else if (pro==2||pro==3||pro==5||pro==6||pro==10)	    {
//	Params.PARS[pro] = pow(10,outer_parm);
//}
//else if (pro==4||pro==15||pro==16||pro==17||pro==18)   {
//	Params.PARS[pro] = outer_parm;
//}
//else {	printf("bad profile input\n");	getc(stdin);	}

//while (1==1)	{

//---------------------Write over the initial params with known fit params --------------------------//

for (k=0;k<=num_adj_pars;k++)	{
	Params.PARS[k] = bestparams[k];
	//printf("%e\n",Params.PARS[k]);
}


reps = 100;  //was 500

for(k=1; k<2; k++){  //k=0: low fungus density; k=1: intermediate fungus density; k=2: high fungus density

	pop=k;

	Params.PARS[30+pop]=initialS[pop];
	Params.PARS[50+pop] = initialR[pop];


    FILE *fp;
    char name[100];
    sprintf(name,"max_lhood/full_mf_s25_hist_reps.txt");  //File name for the output of all values of fraction infected across all realizations
    fp=fopen(name,"a+");

    for(i=0;i<19800;i++){       //JL: Whole gypsy moth range 132*150 grids
        //for (i=2550;i<2558;i++){
        spot=i+1;

        sumINFECTED=0;
		aveINFECTED=0;
		varINFECTED=0;

		for(q = 0; q < 10; q++){  	//begin loop to go through 10 years of CC data FOR EACH GRID CELL!!!  QUICK AND DIRTY AVERAGES ACROSS YEARS
        //printf("Location %d, year %d\n",i,q);

		    year=q;

			DD10=0.0;
			test_day = i*365;  //The starting line of location i in year q
			S_start = test_day;
			limit = (i+1)*365;  //The ending line of location i in year q

			while(DD10 <= hatch & test_day<limit){  //Getting the starting day of epizootic (larvae hatch)

				DDtemp_now = Params.CCDATA[q][test_day][3]-Hlim1;  //CK// begin calculation of accumulated Degree Days
				if(DDtemp_now<0.0){DDtemp_now=0.0;}
				if(DDtemp_now> Hlim2){DDtemp_now=Hlim2;}
				DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time
				S_start++;
				test_day++;
			//printf("temp now: %lf DD10: %f S_start: %f test_day: %d\n", Params.CCDATA[q][test_day][3], DD10, S_start, test_day);		//getc(stdin);
			}
			//printf("S_start in that year: %lf\n", S_start-i*365);
			//printf("S_start: %lf\n", S_start);
			//getc(stdin);
			//printf("S_start: %lf\n", S_start-i*365);		//getc(stdin);


			DD10=0.0;
			test_day = S_start;
			S_end = test_day;

			while(DD10 <= pupate & test_day<limit){  //Getting the ending day of epizootic (pupation)

				DDtemp_now = Params.CCDATA[q][test_day][3]-Plim1;  //CK// begin calculation of accumulated Degree Days
				if(DDtemp_now<0.0){DDtemp_now=0.0;}
				if(DDtemp_now> Plim2){DDtemp_now=Plim2;}
				DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time
				S_end++;
				test_day++;
			//printf("temp now: %f DD10: %f S_end: %f test_day: %d\n", Params.CCDATA[q][test_day-1][3], DD10, S_end, test_day);		getc(stdin);
			}
			//if (S_end-i*365>300){
            //    S_end=300+i*365;
			//}
			//getc(stdin);
		//printf("test 2 here\n");//getc(stdin);

			MAXT3 = S_end - S_start;	//number of days the bugs are active

			//printf("S_end in that year: %lf\n", S_end-i*365);		//getc(stdin);
			//printf("MAXT3: %lf S_start: %lf\t S_end: %lf\n", MAXT3, S_start-i*365, S_end-i*365);		//getc(stdin);
			//printf("MAXT3: %lf S_start: %lf\t S_end: %lf\n", MAXT3, S_start, S_end);


			Params.survivors = 0.0;
			Params.survive=0.0;
			Params.total = Params.PARS[30+pop];

			dim = 2*MAXT3;     //CK//  holds the random variables for each day.  Need 2 per day

			//printf("MAXT3: %lf S_start: %lf\t S_end: %lf\n", MAXT3, S_start-i*365, S_end-i*365);		//getc(stdin);

			for(j=0; j<reps; j++){
                repeat=j;
				DDEVF(&Params,r_seed,dim,pop,MAXT3, S_start, q);   //NEED TO ADD q TO THE DDEVF FUNCTION.  UPDATE ALL WEATHER READING LINES
				//if (pop==1){
				fprintf(fp,"%e\t",1-Params.survive/Params.PARS[30+pop]);
				//}
			}


			Params.total = reps*Params.PARS[30+pop];

			INFECTED[q] = 1.0 - (Params.survivors/Params.total);  //Average fraction infected across realizations
			//printf("Fraction of infected in year %d: %e,\n",q,INFECTED[q]);
			sumINFECTED += INFECTED[q];
			//INFECTED += 1.0 - (Params.survivors/Params.total);
            //printf("Accumulative INFECTED: %lf\t q: %d\t i: %d\n",INFECTED,q,i);
			//}
		}		//END OF 10 YEAR DATA LOOP
		fprintf(fp,"\n");
        aveINFECTED=sumINFECTED/10;  //Average fraction infected across years

        for (qq=0;qq<10;qq++){
            varINFECTED += (INFECTED[qq]-aveINFECTED)*(INFECTED[qq]-aveINFECTED);  //Variance of fraction infected across years
		}

		varINFECTED=varINFECTED/10;
		//printf("Variance: %e\n",varINFECTED);

		Params.CCoutV[i] = varINFECTED; //record variance of infection for each grid cell here

		Params.CCout[i] = aveINFECTED;	//record average infection for each grid cell here

		//printf("INFECTED: %lf\t q: %d\t i: %d\n",Params.CCout[i],q,i);//getc(stdin);

		//printf("survivors= %lf\t total= %lf\t INFECTED=%lf\n",Params.survivors, Params.total, INFECTED);//getc(stdin);
		//printf("%lf\n",INFECTED);//getc(stdin);


		//Params.CCout[i] = INFECTED/10.0;	//calculate average infection for each grid cell here
	//}
    }
    fclose(fp);
	// ------------------------------------------ output results to file  --------------------------------------- //
	fp_results = fopen(strFileName,"a");
	output_file(&Params,fp_results); // prints to output file (filenames.h)
	fclose(fp_results);

	fp_results_v = fopen(strFileNameV,"a");
	output_file_v(&Params,fp_results_v); // prints to output file (filenames.h)
	fclose(fp_results_v);

	// ------------------------------------------------------------------------------------------------------- //

//	printf("test 3 here");getc(stdin);
//	}
}


free_i3tensor(Params.DATA,0,DATA_SETS,0,MAX_WEEKS,0,3);
free_i3tensor(Params.EXPDATA,0,DATA_SETS,0,MAX_WEEKS,0,3);
free_d3tensor(Params.CCDATA,0,10,0,MAX_WEEKS2,0,4);
//free(Params.Rain);
//free(Params.MaxT);
//free(Params.MinRH);
//printf("DONE!!!\n");

return 0;
}
