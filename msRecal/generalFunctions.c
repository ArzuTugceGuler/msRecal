#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "msRecal.h"

void showHelp(){
    
    printf("Usage: msRecal.exe [parameters]\n" ); 
    printf("Compulsory Parameters:\n"); 
    printf("-p<string>\tpepXML file location. Multiple locations can be specified, all must be preceeded by -p.\n"); 
    printf("-m<string>\tmzXML file location.\n");
    printf("\n");
    printf("Optional Parameters:\n");
    printf("-o<string>\tmzXML output file location. May not be equal to input location.\n"); 
    printf("-i<string>\tInstrument type that will be calibrated.\n");
    printf("-a<string>\tMass analyzer type of the instrument that will be calibrated.\n");
    printf("-e<float>\tMaximum mass measurement error.\n"); 
    printf("-t<string>\tScore name. Name of the score in the pepXML file, which value will be considered in comparison to bounds.\n");
    printf("-L<float>\tLower peptide score bound.\n"); 
    printf("-U<float>\tUpper peptide score bound.\n"); 
    printf("-s<int,int>\tmzXML scan window. The first number is the beginscan, the last one the endscan.\n");
    printf("-r<float,float>\tRetention time window. The first number is the -delta on a peptide with respect to a scan spectrum, the second number the +delta.\n");
    printf("-c\tCrop flag. By specifying this, spectra that cannot be recalibrated will be emptied.\n");
    printf("-C<int>\tMinimum number of calibrants needed for spectrum recalibration.\n"); 
    printf("-b\tBackground intensity. Only peaks with an intensity higher than this will be used in recalibration.\n"); 
    
    fflush(stdout);

}

msrecal_params* readParameters(int argc, char *argv[])
{
	char temp[100];
	char* p;
       
	int i;			
	
	msrecal_params* params = (msrecal_params*) malloc(sizeof(msrecal_params));
        printf("\n>>Parameters being initialized..");
	initParameters(params);
	
        for(i=1; i<argc; i++) {	  
		if (argv[i][0] == '-') {	
			
			p = &argv[strlen(argv[i])>2? i: i+1][strlen(argv[i])>2? 2: 0];

			
			if (argv[i][1] == 'p') {
				params->pepxml_file = strclone(p);
			}
						
			
			else if (argv[i][1] == 'm') {
				params->mzxml_file = strclone(p);
			}

			else if (argv[i][1]=='o') {
				params->output_mzXML_file = strclone(p);
			}

			else if (argv[i][1] == 'i') {
				params->instrument = strclone(p);
			}

			else if (argv[i][1] == 'a') {
				params->mass_analyzer = strclone(p);
			}

			else if (argv[i][1] == 'e') {
				params->mmme = atof(p);
			}

			else if (argv[i][1] == 't') {
				params->score_name = strclone(p);
			}

			else if (argv[i][1] == 'L'){
				params->min_score_threshold = atof(p);
			}

			else if (argv[i][1] == 'U'){
				params->max_score_threshold = atof(p);
			}
			
			else if (argv[i][1]=='s') {
				strcpy(temp, p); 
				p = strtok(temp,","); 
				params->ms_start_scan = atoi(p); 
				p = strtok('\0',","); 
				params->ms_end_scan = atoi(p);
			}

			else if (argv[i][1]=='r') {
				strcpy(temp, p); 
				p = strtok(temp,","); 
				params->lower_rel_bnd_rt = atof(p); 
				p = strtok('\0',","); 
				params->upper_rel_bnd_rt = atof(p);
			}

			else if (argv[i][1] == 'f') {
				params->recal_offset = atoi(p);
			}

			else if (argv[i][1] == 'c') {
				params->crop = 1;
			}
			
			else if (argv[i][1] == 'C') {
				params->min_cal = atoi(p);
			}
			
			else if (argv[i][1] == 'b') {
				params->bg = atof(p);
			}

			else if (argv[i][1]=='h') {
				showHelp();
				exit(-1);
			}
  		}// if               
	}// for
        
        printf("\npepXML file: %s", params->pepxml_file);
        printf("\nmzXML file: %s", params->mzxml_file);
        printf("\nInstrument type: %s", params->instrument);
        printf("\nMass analyzer type: %s", params->mass_analyzer);
        printf("\nOutput: %s", params->output_mzXML_file);
        printf("\nMaximum allowed mme: %g ppm", params->mmme);
        printf("\nScore: %s", params->score_name);
        printf("\nScore interval: [%g , %g]", params->min_score_threshold, params->max_score_threshold);
        printf("\nMS scan interval: [%i , %i]", params->ms_start_scan, params->ms_end_scan);
        printf("\nCrop flag: %i", params->crop);
        printf("\nMinimum # of calibrants: %i", params->min_cal);
        printf("\nBackground intensity: %g", params->bg);
        printf("\nRetention time window: [-%g , +%g]", params->lower_rel_bnd_rt, params->upper_rel_bnd_rt);
        printf("\nRetention time offset: %g", params->recal_offset);
        
	// Checking integrity of the parameter file
	/*if (params->pepxml_file == NULL  || params->mzxml_file == NULL)
		return NULL;	*/

	return params;

}// msrecal_params readParameters(int argc, char *argv[])


//Initialization routine for parameters struct 
void initParameters(msrecal_params* params)
{
	int i;

	params->mzxml_file = NULL;
    params->pepxml_file = NULL;
	params->output_mzXML_file = NULL;
	params->score_name = DEFAULT_SCORE;

	params->ms_start_scan = DEFAULT_START_SCAN;
	params->ms_end_scan = INT_MAX;		

	params->crop = DEFAULT_MODE;
	params->min_cal = DEFAULT_MIN_CALIBRANTS;
	params->bg = DEFAULT_MIN_INTENSITY;

	params->lower_rel_bnd_rt = DEFAULT_RT_LOWER;
	params->upper_rel_bnd_rt = DEFAULT_RT_UPPER;
	params->recal_offset = DEFAULT_RECAL_OFFSET;
	
	params->mmme = DEFAULT_MMME;
	params->min_score_threshold = DBL_MIN;
	params->max_score_threshold = DBL_MAX;

	params->instrument = NULL;
	params->mass_analyzer = NULL;
        
        printf("\ndone.");
	
} //void initParameters(parameters* params)
 

