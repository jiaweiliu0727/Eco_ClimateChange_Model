int inputdata(void *Paramstuff)
{
#define MAX_WEEKS 100	// larger than the number of weeks in any data set
#define MAX_WEEKS2 8000000	// larger than the number of lines in the climate change data set


STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
// loads the data into a matrix and finds the number of weeks in the data set
int Sdata[MAX_WEEKS]; int Vdata[MAX_WEEKS]; int Fdata[MAX_WEEKS]; int Ddata[MAX_WEEKS]; int Cdata[MAX_WEEKS]; int D2data[MAX_WEEKS];
int weeks;	int i; int j; int q;
int num_weeks[DATA_SETS+1];			// used for output to file
int num_weeks2[DATA_SETS+1];			//CK// used for output to file
int num_weeks3[DATA_SETS+1];			//CK// used for output to file
int total_days=0;					// the number of days summed over all data sets (for MISER)
Params->DATA = i3tensor(0,20,0,MAX_WEEKS,0,5);
Params->EXPDATA = i3tensor(0,20,0,MAX_WEEKS,0,5);
Params->CCDATA = d3tensor(0,10,0,MAX_WEEKS2,0,4);

double rain;
double maxT;
double aveT;
double minRH;

int FlagF;

char *file;
char *file_name="data";
char *file_name2="expdata";  //CK// name for inputing the experimental data
//////////////////////////////////////////////////////////////////////////////////////////////
char *file_name3="highRCP_BC_ccENDall";  //CK// SWITCH BETWEEN HISTORIC DATA AND FUTURE DATA HERE! EITHER USE "ccHIST" OR "ccFUTURE"
//////////////////////////////////////////////////////////////////////////////////////////////
char *file_type=".txt";

char *code;
char *code_name="ftp";

char numbs[5];
/*------------------------------- Data Sets in Kyle et al. 2020---------------------------------*/
for (j=1;j<=DATA_SETS;j++)	{
	weeks=0;	i=0;	FlagF=0;
	FILE *ftp_data;

//printf("just before numbs...\n");
	sprintf(numbs,"%d",j);

	file = (char*)calloc((strlen(file_name)+strlen(file_type)+strlen(numbs)+1),sizeof(char));
	code = (char*)calloc((strlen(code_name)+strlen(numbs)+1),sizeof(char));

	strcat(file,file_name);
	strcat(file,numbs);
	strcat(file,file_type);

	strcat(code,code_name);
	strcat(code,numbs);

//printf("just before fopen...\n");

	ftp_data=fopen(file,"r");
	if (ftp_data==0)	{printf("file %d open error \n",j);		getc(stdin);	}
	//else				{printf("open data %d success \n",j);	fflush(stdout);	}

	while (fscanf(ftp_data,"%d %d %d %d %d\n",&Sdata[i],&Vdata[i],&Fdata[i],&Ddata[i],&D2data[i])!= EOF)			{
		Params->DATA[j][i][0]=Sdata[i]; Params->DATA[j][i][1]=Vdata[i]; Params->DATA[j][i][2]=Fdata[i]; Params->DATA[j][i][3]=Ddata[i]; Params->DATA[j][i][4]=D2data[i];
		//printf("FERALS: i=%d\t pop:%d\t healthy:%d\t viral:%d\t fungal:%d\t week:%d\t week2:%d\n",i,j,Params->DATA[j][i][0],Params->DATA[j][i][1],Params->DATA[j][i][2], Params->DATA[j][i][3], Params->DATA[j][i][4]);
		//if ((Fdata[i]>0) && (FlagF<1))	{
		//	FlagF=2;
		//	Params->DAY_F[j]=7*i;
		//}

		weeks++; i++;
	}
	fclose(ftp_data);

	Params->DAY_F[j]=0.0;  //CK// making resting spores start blooming on first day

	Params->MAXT[j]=7*(weeks-1);				// number of days
	total_days += Params->MAXT[j];
	//printf("data set %d has %d days\t total=%d\n",j,Params->MAXT[j],total_days);getc(stdin);
	num_weeks[j]=i;

	//if (FlagF==0)	Params->DAY_F[j]=Params->MAXT[j];
}
//getc(stdin);

//printf("just before Experimental Data...\n");getc(stdin);

/*-------------------------------Experimental Data Sets in Kyle et al. 2020---------------------------------*/
for (j=1;j<=DATA_SETS;j++)	{
	weeks=0;	i=0;	FlagF=0;
	FILE *ftp_data;

//printf("just before numbs...\n");

	sprintf(numbs,"%d",j);

	file = (char*)calloc((strlen(file_name2)+strlen(file_type)+strlen(numbs)+1),sizeof(char));
	code = (char*)calloc((strlen(code_name)+strlen(numbs)+1),sizeof(char));

	strcat(file,file_name2);
	strcat(file,numbs);
	strcat(file,file_type);

	strcat(code,code_name);
	strcat(code,numbs);

//printf("just before fopen...\n");getc(stdin);

	ftp_data=fopen(file,"r");
	if (ftp_data==0)	{printf("file %d open error \n",j);		getc(stdin);	}
	//else				{printf("open data %d success \n",j);	fflush(stdout);	}

	while (fscanf(ftp_data,"%d %d %d %d %d \n",&Sdata[i],&Fdata[i],&Ddata[i],&Cdata[i],&D2data[i])!= EOF)			{
		Params->EXPDATA[j][i][0]=Sdata[i]; Params->EXPDATA[j][i][1]=Fdata[i]; Params->EXPDATA[j][i][2]=Ddata[i]; Params->EXPDATA[j][i][3]=Cdata[i]; Params->EXPDATA[j][i][4]=D2data[i];
		//printf("EXPERIMENTALS: i=%d\t wk_number:%d\t healthy:%d\t fungal:%d\t week:%d\t covered:%d\t date2:%d\n",i,weeks,Params->EXPDATA[j][i][0],Params->EXPDATA[j][i][1],Params->EXPDATA[j][i][2],Params->EXPDATA[j][i][3],Params->EXPDATA[j][i][4]);
		//if ((Fdata[i]>0) && (FlagF<1))	{
		//	FlagF=2;
		//	Params->DAY_F[j]=7*i;
		//}
		weeks++; i++;
	}
	fclose(ftp_data);

	Params->MAXT2[j]=weeks-1;				// number of days
	total_days += Params->MAXT2[j];
	//printf("TRUE data set %d has %d days\n",j,Params->MAXT2[j]);getc(stdin);
	num_weeks2[j]=i;

	//if (FlagF==0)	Params->DAY_F[j]=Params->MAXT2[j];
}

//getc(stdin);

//printf("just before Weater Data...\n");

/*-------------------------------Projected Weather Data Sets ---------------------------------*/
//* Generated by Jiali Wang and Rao Kotamarthi. Available upon request since files can be quite large.*//
for (j=1;j<=DATA_SETS;j++)	{
	weeks=0;	i=0;	FlagF=0; 	q=0;
	FILE *ftp_data;

	sprintf(numbs,"%d",j);

	file = (char*)calloc((strlen(file_name3)+strlen(file_type)+strlen(numbs)+1),sizeof(char));
	code = (char*)calloc((strlen(code_name)+strlen(numbs)+1),sizeof(char));

	strcat(file,file_name3);
	strcat(file,numbs);
	strcat(file,file_type);

	strcat(code,code_name);
	strcat(code,numbs);
//printf("within weather data sets \n");

	ftp_data=fopen(file,"r");
	//ftp_data = fopen("weather1.txt","r");
	if (ftp_data==0)	{printf("file %d open error \n",j);		getc(stdin);	}
	//else				{printf("open data %d success \n",j);	fflush(stdout);	}

	//printf("just before while statement...\n");
	while (fscanf(ftp_data,"%lf %lf %lf %lf \n",&rain,&minRH,&maxT,&aveT)!= EOF){
//	while (fscanf(ftp_data,"%lf %lf %lf \n",&Params.WDATA[i][0],&Params.WDATA[i][1],&Params.WDATA[i][2])!= EOF)			{
		//printf("inside while statement...\n");getc(stdin);

		if(i==7227000){ weeks=0;	i=0;	FlagF=0; q++;}         //JL: Whole gypsy moth range; 132*150 grids *365 days=7227000 lines for each year.
        //JL: q indicates year (1-10); i has information of location id and day (location 1 day 1-365; location 2 day 1-365 ...); the last dimension is for different weather variables.
		Params->CCDATA[q][i][0]=rain; Params->CCDATA[q][i][1]=minRH; Params->CCDATA[q][i][2]=maxT; Params->CCDATA[q][i][3]=aveT;

		//printf("WEATHER: rain:%lf\t  minRH:%lf\t maxT:%lf\t aveT:%lf\n",Params->Rain[i],Params->MinRH[i],Params->MaxT[i],Params->AveT[i]);getc(stdin);
		//printf("WEATHER: q=%d\t i=%d\t rain:%lf\t  minRH:%lf\t maxT:%lf\t aveT:%lf\n",q,i,Params->CCDATA[q][i][0],Params->CCDATA[q][i][1],Params->CCDATA[q][i][2], Params->CCDATA[q][i][3]); //getc(stdin);
		//if ((Fdata[i]>0) && (FlagF<1))	{
		//	FlagF=2;
		//	Params->DAY_F[j]=7*i;
		//}
		weeks++; i++;
		//printf("weeks: %d i:%d \n",weeks,i);

	}
	fclose(ftp_data);

	//Params->MAXT3[j]=weeks-1;				// number of days
	//total_days += Params->MAXT3[j];    //don't think weather data should count for miser...
	//printf("data set %d has %d days\t total=%d\n",j,Params->MAXT2[j],total_days);getc(stdin);
	num_weeks3[j]=i;

	//getc(stdin);

	//if (FlagF==0)	Params->DAY_F[j]=Params->MAXT2[j];
}

	//printf("ok?\n");
	//getc(stdin);
	//exit(1);

/*----------------------------end of data sets------------------------------*/
FILE *fp_weeks;
fp_weeks=fopen("weeks.dat","w");

for (j=1;j<=DATA_SETS;j++)	{
	fprintf(fp_weeks,"%d\t",num_weeks[j]);
}
fclose(fp_weeks);

return 0;
}
