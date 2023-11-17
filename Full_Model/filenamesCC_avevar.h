
void GetString(int j,int i,char* strOut, unsigned int strSize)
{
// ---------------------------------------- Name for Output Files ------------------------------------------- //
char *Prefix=" ";
char *max_lhood_prefix="2023_full_mf_s25_average_future_rep100";
char *mcmc_prefix="big8june_";
char *profile_prefix="intrinsic_";

char *FileName	= " ";
char *Path		= " ";
char *Name		= " ";
char *Type		= ".dat";						// type of output file
char profile[33];
// ---------------------------------------- create time suffixes -------------------------------------------- //
struct tm *ptr;									// struct defined in time.h;  used to hold the time and date
time_t lt;										// type defined in time.h; for storing the calendar time
lt = time(NULL);
ptr = localtime(&lt);							// converts calendar time to local time
char Date[30];									// string to hold date and time
strftime(Date, 30, "%m_%d_%H_%M_%S", ptr);		// adds to Date month, day, year, hour, minute, second from ptr
// ---------------------------------- create paths and allocate memory ------------------------------------- //
if	(j==0)	{										// MCMC
	Path="/home/jiaweiliu/WholeGM_CCModel/mcmc_results/";
	Prefix=mcmc_prefix;
	if		(i==1)	{	Name = "parms_";											}
	else if	(i==2)	{	Name = "pc_";												}
	else if (i==3)	{	Name = "acc_";												}
	else if (i==4)	{	Name = "lastpars_";											}
	else			{	printf("bad i value in filenames!!!\n");	getc(stdin);	}
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name)+strlen(Type)+strlen(Date)+1), sizeof(char));
}
else if (j==1)	{									// Line Search MLE
	Prefix=max_lhood_prefix;
	if		(i==0)	{	Path="/home/jiaweiliu/WholeGM_CCModel/max_lhood/";			Name="max_lhood";			}
	else if (i==4)	{	Path="/home/jiaweiliu/WholeGM_CCModel/max_lhood/";			Name="L_";				}
	else			{	printf("bad i value in filenames!!!\n");		getc(stdin);	}
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name) + strlen(Type)+1), sizeof(char));
}
else if (j>1)	{									// Profile Likelihoods
	Path="/home/jiaweiliu/WholeGM_CCModel/profile_data/";
	Prefix=profile_prefix;
	sprintf(profile,"%d",j);
	Name="pro_";
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name) + strlen(Type)+strlen(profile) +1), sizeof(char));
}
// ---------------------------------- add strings onto the file name ------------------------------------------ //
strcat(FileName,Path);
strcat(FileName,Prefix);
strcat(FileName,Name);

if		(j==0||(j==1 && i==4))	strcat(FileName,Date);		// add Date to output file name		(MCMC & max lhood)
else if (j>1)					strcat(FileName,profile);	// add profile number to file name	(profile likelihood)

strcat(FileName,Type);

if (j==0)			strncpy(strOut,FileName,strSize);
else if (j>=1)		strncpy(strOut,FileName,strSize);

return;
}


void GetStringV(int j,int i,char* strOut, unsigned int strSize)
{
// ---------------------------------------- Name for Output Files ------------------------------------------- //
char *Prefix=" ";
char *max_lhood_prefix="2023_full_mf_s25_variance_future_rep100";
char *mcmc_prefix="big8june_";
char *profile_prefix="intrinsic_";

char *FileName	= " ";
char *Path		= " ";
char *Name		= " ";
char *Type		= ".dat";						// type of output file
char profile[33];
// ---------------------------------------- create time suffixes -------------------------------------------- //
struct tm *ptr;									// struct defined in time.h;  used to hold the time and date
time_t lt;										// type defined in time.h; for storing the calendar time
lt = time(NULL);
ptr = localtime(&lt);							// converts calendar time to local time
char Date[30];									// string to hold date and time
strftime(Date, 30, "%m_%d_%H_%M_%S", ptr);		// adds to Date month, day, year, hour, minute, second from ptr
// ---------------------------------- create paths and allocate memory ------------------------------------- //
if	(j==0)	{										// MCMC
	Path="/home/jiaweiliu/WholeGM_CCModel/mcmc_results/";
	Prefix=mcmc_prefix;
	if		(i==1)	{	Name = "parms_";											}
	else if	(i==2)	{	Name = "pc_";												}
	else if (i==3)	{	Name = "acc_";												}
	else if (i==4)	{	Name = "lastpars_";											}
	else			{	printf("bad i value in filenames!!!\n");	getc(stdin);	}
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name)+strlen(Type)+strlen(Date)+1), sizeof(char));
}
else if (j==1)	{									// Line Search MLE
	Prefix=max_lhood_prefix;
	if		(i==0)	{	Path="/home/jiaweiliu/WholeGM_CCModel/max_lhood/";			Name="max_lhood";			}
	else if (i==4)	{	Path="/home/jiaweiliu/WholeGM_CCModel/max_lhood/";			Name="L_";				}
	else			{	printf("bad i value in filenames!!!\n");		getc(stdin);	}
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name) + strlen(Type)+1), sizeof(char));
}
else if (j>1)	{									// Profile Likelihoods
	Path="/home/jiaweiliu/WholeGM_CCModel/profile_data/";
	Prefix=profile_prefix;
	sprintf(profile,"%d",j);
	Name="pro_";
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name) + strlen(Type)+strlen(profile) +1), sizeof(char));
}
// ---------------------------------- add strings onto the file name ------------------------------------------ //
strcat(FileName,Path);
strcat(FileName,Prefix);
strcat(FileName,Name);

if		(j==0||(j==1 && i==4))	strcat(FileName,Date);		// add Date to output file name		(MCMC & max lhood)
else if (j>1)					strcat(FileName,profile);	// add profile number to file name	(profile likelihood)

strcat(FileName,Type);

if (j==0)			strncpy(strOut,FileName,strSize);
else if (j>=1)		strncpy(strOut,FileName,strSize);

return;
}

void output_file(void *Paramstuff,FILE *fp_results) //Output file of average
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
int i,pop;

//for (i=0;i<1008;i++)	{
for(i=0;i<19800;i++){       //JL: Whole gypsy moth range 132*150 grids
//for (i=2550;i<2558;i++){
	fprintf(fp_results,"%e\t",Params->CCout[i]);
	//printf("%e\n",Params->CCout[i]);
}
fprintf(fp_results,"\n");


}
void output_file_v(void *Paramstuff,FILE *fp_results_v) //Output file of variance
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
int i,pop;

//for (i=0;i<1008;i++)	{
for(i=0;i<19800;i++){       //JL: Whole gypsy moth range 132*150 grids
//for (i=2550;i<2558;i++){
	fprintf(fp_results_v,"%e\t",Params->CCoutV[i]);
	//printf("%e\n",Params->CCout[i]);
}
fprintf(fp_results_v,"\n");


}




