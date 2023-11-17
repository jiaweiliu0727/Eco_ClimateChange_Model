#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_errno.h>	// GSL_SUCCESS ...
#include <gsl/gsl_odeiv.h>	// ODE solver

#include <gsl/gsl_types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_min.h>

#include "nrutil.c"

char *strFileNameDate;

#define DATA_SETS 1		        // number of data sets
#define NUM_PARS 100		        // number of parameters to be passed from main to hood


const double h = 0.01;		        // time step
int year=0;
int repeat=0;
int spot=0;

//CK Structure for experimental data!!!

struct dataset{  //  building the structure, to be declared later
	int Date;
	int Tree;
	int Covered;
	double total;
	double fungus;
};

typedef struct
{
	double PARS[NUM_PARS];

	double nuV;
	double nuF;
	double nuR;
	double muF;
	double size_C;

	double R_END[DATA_SETS+1];
	double INITS[DATA_SETS+1];
	double POPS[4];
	double EV[100];
	double EF[100];
	double DAY_F[DATA_SETS+1];	    // day where fungal infection first happens
	int MAXT[DATA_SETS+1];		    // number of days in data set (different for each data set)
	int MAXT2[DATA_SETS+1];		    // number of days in EXPERIMENTAL data set (different for each EXPERIMENTAL data set)
	int MAXT3[DATA_SETS+1];		    // number of days in WEATHER data set (different for each WEATHER data set)
	int ***DATA;				    // 3-dimensional array that holds all the data
	int ***EXPDATA;				    // 3-dimensional array that holds all the EXPERIMENTAL data
	double ***CCDATA;				// 3-dimensional array that holds all the WEATHER data

    double CCout[19800];  //JL: Whole gypsy moth range 132*150 grids
    double CCoutV[19800];
	double survivors;
	double survive;
	double total;

	double test_data[1000][36];

	double AcceptedVect[NUM_PARS];
	double LoopVect[NUM_PARS];

	double parm_low[NUM_PARS];
	double parm_high[NUM_PARS];
	double parm_step[NUM_PARS];

	double MLE[NUM_PARS];
	double MLE_host[DATA_SETS+1];
	double MLE_initR[DATA_SETS+1];

	double best_initS[DATA_SETS+1];
	double best_initR[DATA_SETS+1];

	int parm_inc;
	double sim_results[55][4];		// 1st entry larger than the number of weeks in any data set
	int th_id;
	int pop;
}STRUCTURE;


#include "inputdataCC2_multi2.h"
#include "random_setup.h"
#include "filenamesCC_avevar.h"
#include "bounds_stochWEATHER_window10_AveR_logistic2.h"
#include "prob_dists.h"
#include "fast_odes2.h"
#include "DDEVF_GCC_R+T+RH_multi.h"
