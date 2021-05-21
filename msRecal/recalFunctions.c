#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "msRecal.h"
#include "string.h"
#include "PeptideProphetDelegate.h"

static int candidate_list[MAX_ROWS][MAX_CALIBRANTS];
static int scan_window;
static int* scan_cal_index;

static int seqlen;
static int fragment[5];
static calibrant calibrant_list[MAX_CALIBRANTS];

static double mimass[5] = 
{
	1.0078250321,
	12,
	14.0030740052,
	15.9949146221,
	31.97207069
};

static double cyclosiloxanes[5] = 
{
	593.157605,
	667.176396,
	741.195187,
	815.213979,
	889.232770
};

static unsigned char aa[20][5] = 
{
	{5,3,1,1,0}, /* amino acid compositions */
	{12,6,4,1,0},
	{6,4,2,2,0},
	{5,4,1,3,0},
	{5,3,1,1,1}, 
	{7,5,1,3,0},
	{8,5,2,2,0},
	{3,2,1,1,0},
	{7,6,3,1,0},
	{11,6,1,1,0},
	{11,6,1,1,0},
	{12,6,2,1,0},
	{9,5,1,1,1},
	{9,9,1,1,0},
	{7,5,1,1,0},
	{5,3,1,2,0},
	{7,4,1,2,0},
	{10,11,2,1,0},
	{9,9,1,2,0},
	{9,5,1,1,0}
};

static int n_calibrants;
static double mz, Ca, Cb;

  
/* Function that determines if the peptide should be processed or not */
int process_peptide(search_hit sh, msrecal_params* params)
{
	int process = 0, j;
	search_score ss;
	//double* score = NULL;
	double* score = malloc(sizeof(double));
	void* hookstruct;	

	for (j=0; j<strlen(sh.peptide); j++) {
	    if (strchr(ANTI_ACIDS, sh.peptide[j])) {
                return 0;
            }
	}

	for (j=0; j<sh.search_score_count; j++) {
            ss = sh.search_score_array[j];
            if (strstr(ss.name, params->score_name)== NULL) {
                continue;
            }
            else{

                if (ss.value >= params->min_score_threshold && ss.value <= params->max_score_threshold) {
                    process = 1;
                }
                return process;
            }
            
	}
        
	/* None of the regular score measures applied, now looking for hooked ones from peptide prophet */
	if (!process) {
        for (j=0; j<sh.analysis_result_count; j++) {
        	hookstruct = (void*) sh.analysis_result_array[j].hook;

            if (strstr(sh.analysis_result_array[j].analysis, "peptideprophet") && hookstruct){
            	score = (double*) peptide_prophet_result_property(params->score_name, hookstruct);
            	//printf("Score: %f\n", *score); fflush(stdout);
            }

            if (score && *score >= params->min_score_threshold && *score <= params->max_score_threshold) {
            	//printf("OK!\n");
            	process = 1;

            }
        }
	}

	return process;

}

peptideset* build_peptide_set(pmsms_pipeline_analysis pepfile, msrecal_params* params, int* pepnum)
{    
	int i, j;
	int peptide_count = 0;
	peptideset *peptide, *retval;
	spectrum_query sq;
	search_hit sh;

	// Counting the total number of peptides
	*pepnum = 0;
	retval = NULL;
	
    peptide_count = 0;
	for (i=0; i<pepfile->run_summary_count; i++) {
            peptide_count += pepfile->run_summary_array[i].spectrum_query_count;
	}
        
    printf("\n%i MS/MS run(s) with %i spectra.\n",pepfile->run_summary_count , peptide_count); fflush(stdout);

	peptide = (peptideset*) malloc(sizeof(peptideset)*peptide_count);

	
	for (i=0; i<pepfile->run_summary_count; i++) {
            for (j=0; j<pepfile->run_summary_array[i].spectrum_query_count; j++) {
            	sq = pepfile->run_summary_array[i].spectrum_query_array[j];	/* ith search hit */
            	sh = sq.search_result_array[0].search_hit_array[0];
            	//Check if the score satisfies the boundaries given as cl arg
            	//printf("Spectrum query: %s\n", sh.peptide); fflush(stdout);
            	if (!process_peptide(sh, params)){
            		//printf("Not selected\n"); fflush(stdout);
                    continue;
                }
         
		/* Found valid peptide */
		peptide[*pepnum].sequence = strclone(sh.peptide);
		peptide[*pepnum].retention = sq.retention_time_sec;
		*pepnum += 1;
	    }
	}

	if (!retval)
	    retval = (peptideset*) malloc(sizeof(peptideset) * (*pepnum));
	else
	    retval = (peptideset*) realloc(retval, sizeof(peptideset) * (*pepnum));

	for (i=0; i<(*pepnum); i++) {
	    retval[i] = peptide[i];
	}/* for */
	free(peptide);

	return retval;
}


/* Function that builds the collection of internal calibrants */
void build_internal_calibrants(pmzxml_file mzXML_file, peptideset* peptide_set, int set_size, msrecal_params* params)
{
	int i, j, k, unique;
	double min_rt, max_rt;
	scan_attributes satt;

    scan_window = (params->ms_end_scan - params->ms_start_scan)+1;
    scan_cal_index = (int*) malloc(sizeof(int) * scan_window);

	for (i=0; i<scan_window; i++)
		scan_cal_index[i]=0;

	for(i=0; i<set_size; i++) {

		printf("Count: %i\n", i); fflush(stdout);

        min_rt = (peptide_set[i].retention + params->recal_offset) - params->lower_rel_bnd_rt;

        max_rt = (peptide_set[i].retention + params->recal_offset) + params->upper_rel_bnd_rt;

        for(j=params->ms_start_scan; j<=params->ms_end_scan; j++) {

        	//printf("\tScan: %i\n", j); fflush(stdout);
            satt = get_scan_attributes(mzXML_file, j);
            //printf("Count: %i\n", i); fflush(stdout);
            if (satt.retentionTime < min_rt){
            	continue;
            }
            else if (satt.retentionTime > max_rt){
                break;
            }

            unique = 1;
            // Check if the same sequence isn't already in the list for this spectrum
            for (k=0; k<scan_cal_index[j-(params->ms_start_scan)]; k++) {
            	if (strcmp(peptide_set[candidate_list[j-1][k]].sequence, peptide_set[i].sequence) == 0) {
            		unique = 0;
            		break;
            	}
            }

            if (unique) {
            	candidate_list[j-1][scan_cal_index[j-(params->ms_start_scan)]] = i; /* copy peptide sequence to internal calibrant candidate lists for nearby MS spectra... [row][calibrant]*/
                scan_cal_index[j-(params->ms_start_scan)]++; //  number of calibrants for one ms scan
            }

        }
        //printf("\tNumber: %i\n", scan_cal_index[j-(params->ms_start_scan)]); fflush(stdout);

	}

        for(j=params->ms_start_scan; j<=params->ms_end_scan; j++) {
            printf("\nMS scan: %i\t%f\n", j, get_scan_attributes(mzXML_file, j).retentionTime); fflush(stdout);
            for(k=0; k<scan_cal_index[j-(params->ms_start_scan)]; k++) {
                printf("%i\t%s\t%f\n", k+1, peptide_set[candidate_list[j-1][k]].sequence, peptide_set[candidate_list[j-1][k]].retention); fflush(stdout);
            }
        }


}

/* Monisotopic mass calculator */
double calc_mass(char* sequence,int charge)
{
	int i,a,b;
	double mass;

	seqlen=strlen(sequence);


	for(i=0;i<5;i++)
            fragment[i]=0;

	for(i=0; i<seqlen; i++){ // generate molecular formula for fragment
            a=20-strlen(strchr(AMINO_ACIDS,sequence[i]));

            for(b=0;b<5;b++) {
		fragment[b]=fragment[b]+aa[a][b];
            }
	}

	fragment[0]=fragment[0]+2;
	fragment[3]++;  // add H2O

	// calculate and return mass based on molecular formula
	mass=0.0;
	for(b=0;b<5;b++)
		mass += fragment[b]*mimass[b];

	return (mass+charge*HPLUS_MASS)/charge;

}

int sort_type_comp_inv_int(const void *i, const void *j){
	if (((*(calibrant*)j).intensity - (*(calibrant*)i).intensity) > 0)
		return 1;
	else if (((*(calibrant*)j).intensity - (*(calibrant*)i).intensity) < 0)
		return - 1;
	return 0;

}

void makeCalibrantList(int scan, pscan_peaks mzpeaks, peptideset* peptide_set, msrecal_params* params){

        int i, j, k, z;
        n_calibrants = 0;
        // Find internal calibrants for the peaks of the spectrum
        for (i=0; i<mzpeaks->count; i++) {
            if (mzpeaks->intensities[i] < params->bg){
                continue;
            }

            for(z=1; z<=4; z++) {
		for(j=0; j<scan_cal_index[scan-(params->ms_start_scan)]; j++) {
                    mz = calc_mass(peptide_set[candidate_list[scan-1][j]].sequence, z);
                    if(fabs((mz - mzpeaks->mzs[i])/mz)<=params->mmme) {
                        //printf("\nTheoretical mz: %f", mz); fflush(stdout);
                        //printf("\nMeasured mz: %f", mzpeaks->mzs[i]); fflush(stdout);
                        //printf("\nMeasured intensity: %f\n", mzpeaks->intensities[i]); fflush(stdout);
			calibrant_list[n_calibrants].mz = mz;
			calibrant_list[n_calibrants].peak = mzpeaks->mzs[i];
			calibrant_list[n_calibrants].intensity = mzpeaks->intensities[i];
			n_calibrants++;
                    }
                }

                // add cyclosiloxane peaks as potential internal calibrants in row i
		if(z==1) { // all cyclosiloxanes are singly charged
                    for(j=0; j<5; j++) {
			if(fabs((cyclosiloxanes[j] - mzpeaks->mzs[i])/ cyclosiloxanes[j]) <= params->mmme/1000000) {
                            calibrant_list[n_calibrants].mz = cyclosiloxanes[j];
                            calibrant_list[n_calibrants].peak = mzpeaks->mzs[i];
                            calibrant_list[n_calibrants].intensity = mzpeaks->intensities[i];
                            n_calibrants++;
			}
                    }
		}
            }
	}
        qsort(calibrant_list, n_calibrants, sizeof(calibrant), sort_type_comp_inv_int);
}

