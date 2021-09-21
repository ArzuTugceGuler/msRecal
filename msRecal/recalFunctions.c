#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "msRecal.h"
#include "string.h"
#include "PeptideProphetDelegate.h"

static int candidate_list[MAX_ROWS][MAX_CALIBRANTS];
static int scan_window;
static int* scan_cal_index;

static int seqlen;
static int fragment[5];
static calibrant calibrant_list[MAX_CALIBRANTS];

static int calibrated_scans[MAX_ROWS];
//For each MS1, Ca, Cb, Cc (if applicable)
static double calibration_coefficients[MAX_ROWS][2];

static int calib_mode;

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
static int calibration_counter = 0;
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
            	//Add try-catch for pepXML files with empty <search_result>
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

		//printf("Count: %i\n", i); fflush(stdout);

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
            	candidate_list[j-1][scan_cal_index[j-(params->ms_start_scan)]] = i; // copy peptide sequence to internal calibrant candidate lists for nearby MS spectra... [ms_scan_no-1][ms_scan_calibrant_no-1]
                scan_cal_index[j-(params->ms_start_scan)]++; //  number of candidate calibrants for each ms scan within the RT window
            }
        }


	}

	//Print candidate peptides for each MS scan
    /*for(j=params->ms_start_scan; j<=params->ms_end_scan; j++) {
    	printf("\nMS scan: %i\t retention time: %f\n", j, get_scan_attributes(mzXML_file, j).retentionTime);
        for(k=0; k<scan_cal_index[j-(params->ms_start_scan)]; k++) {
        	printf("%i\t%s\t%f\n", k+1, peptide_set[candidate_list[j-1][k]].sequence, peptide_set[candidate_list[j-1][k]].retention); fflush(stdout);
        }
    }*/
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
	int i, j, z, cyclo;
    n_calibrants = 0;
    cyclo = 0;

    // Find internal calibrants for the peaks of the spectrum
    for (i=0; i<mzpeaks->count; i++) {
    	//do not take the intensities below the background threshold into account
    	if (mzpeaks->intensities[i] < params->bg){
    		//printf("\nPeak intensity:%f", mzpeaks->intensities[i]); fflush(stdout);
    		continue;
        }
    	//charges from 1 to 4
    	for(z=1; z<=4; z++) {
    		//loop through the all the candidate calibrants of this scan selected based on rt
    		for(j=0; j<scan_cal_index[scan-(params->ms_start_scan)]; j++) {
    			//calculate the m/z of the candidate calibrant
    			mz = calc_mass(peptide_set[candidate_list[scan-1][j]].sequence, z);
    			//printf("\nmz calculated:%f measured:%f", mz, mzpeaks->mzs[i]); fflush(stdout);
    			//Check if each potential calibrant fits the mmme given
                if(fabs( (mz - mzpeaks->mzs[i]) / mz ) <= params->mmme / 1000000 ) {
            		/*
            		printf("\nTheoretical mz: %f", mz); fflush(stdout);
            		printf("\nMeasured mz: %f", mzpeaks->mzs[i]); fflush(stdout);
            		printf("\nMeasured intensity: %f\n", mzpeaks->intensities[i]); fflush(stdout);
            		*/
            		calibrant_list[n_calibrants].mz = mz;
            		calibrant_list[n_calibrants].peak = mzpeaks->mzs[i];
            		calibrant_list[n_calibrants].intensity = mzpeaks->intensities[i];
            		n_calibrants++;
            	}
            }
            // add cyclosiloxane peaks as potential internal calibrants in row i
            if(z==1) { // all cyclosiloxanes are singly charged
            	for(j=0; j<5; j++) {
            		if(fabs( (cyclosiloxanes[j] - mzpeaks->mzs[i]) / cyclosiloxanes[j]) <= params->mmme / 1000000 ) {
            			calibrant_list[n_calibrants].mz = cyclosiloxanes[j];
            			calibrant_list[n_calibrants].peak = mzpeaks->mzs[i];
            			calibrant_list[n_calibrants].intensity = mzpeaks->intensities[i];
            			n_calibrants++;
            			cyclo++;
            		}
            	}
            }
    	}
    }

    // --------------------------------------------------------PRINT
    printf("\tcyclosiloxanes\t%i", cyclo); fflush(stdout);
    printf("\tinit_calibrants\t%i", n_calibrants); fflush(stdout);

    qsort(calibrant_list, n_calibrants, sizeof(calibrant), sort_type_comp_inv_int);
}

// Calibration function CAL2 Inverted for FTICR
int calib_f_FTICR(const gsl_vector *x, void *params, gsl_vector *f)
{
	double *y = ((struct data *)params)->y;
	double *mz = ((struct data *)params)->mz2;
	double a = gsl_vector_get (x, 0);
	double b = gsl_vector_get (x, 1);
	double M;
    size_t i;

	for (i=0;i<n_calibrants;i++) {
		//Model m = a/(f-b) (CAL2 inverted)
		M = a/(y[i]-b);
		gsl_vector_set (f, i, (M-mz[i]));
	}// for

    return GSL_SUCCESS;

}

// DF calibrator for FTICR
int calib_df_FTICR(const gsl_vector *x, void *params, gsl_matrix *J)
{
	double *y = ((struct data *)params)->y;
	double a = gsl_vector_get (x, 0);
	double b = gsl_vector_get (x, 1);
	size_t i;

	for (i=0;i<n_calibrants;i++) {
		gsl_matrix_set (J,i,0, 1/(y[i]-b) );
		gsl_matrix_set (J,i,1, a/((y[i]-b)*(y[i]-b)) );
	}// for

	return GSL_SUCCESS;

}

// FDF Calibrator for FTICR
int calib_fdf_FTICR(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
	calib_f_FTICR (x,params,f);
	calib_df_FTICR (x,params,J);

	return GSL_SUCCESS;
}

// Calibration function for Orbitrap FTMS
int calib_f_FTMS(const gsl_vector *x, void *params, gsl_vector *f)
{
	double *y = ((struct data *)params)->y;
	double *mz = ((struct data *)params)->mz2;
	double a = gsl_vector_get (x, 0);
	double M;
    size_t i;

	for (i=0;i<n_calibrants;i++) {
		//Model m = A/(f*f)
		M = a/(y[i]*y[i]);
		gsl_vector_set (f, i, (M-mz[i]));
	}

    return GSL_SUCCESS;
}

// DF calibrator Orbitrap FTMS
int calib_df_FTMS(const gsl_vector *x, void *params, gsl_matrix *J)
{
	double *y = ((struct data *)params)->y;
	//double a = gsl_vector_get (x, 0);
	size_t i;

	for (i=0;i<n_calibrants;i++) {
		gsl_matrix_set (J,i,0, 1/(y[i]*y[i]) );
	}

	return GSL_SUCCESS;

}

// FDF Calibrator for FTMS
int calib_fdf_FTMS(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
	calib_f_FTMS (x,params,f);
	calib_df_FTMS (x,params,J);

	return GSL_SUCCESS;
}

// Calibration function inverted for TOF
int calib_f_TOF(const gsl_vector *x, void *params, gsl_vector *f)
{
	double *y = ((struct data *)params)->y;
	double *mz = ((struct data *)params)->mz2;
	double a = gsl_vector_get (x, 0);
	double b = gsl_vector_get (x, 1);
	double M;
    size_t i;

	for (i=0;i<n_calibrants;i++) {
		//Model m = ((t-b)/a)^2
		M = ( (y[i]-b) / a ) * ( (y[i]-b) / a );
		gsl_vector_set (f, i, (M-mz[i]));
	}// for

    return GSL_SUCCESS;

}

// DF calibrator for TOF
int calib_df_TOF(const gsl_vector *x, void *params, gsl_matrix *J)
{
	double *y = ((struct data *)params)->y;
	double a = gsl_vector_get (x, 0);
	double b = gsl_vector_get (x, 1);
	size_t i;

	for (i=0;i<n_calibrants;i++) {
		gsl_matrix_set (J,i,0, -2 * (y[i]-b) * (y[i]-b) / (a*a*a) );
		gsl_matrix_set (J,i,1, 2 * (b-y[i]) / (a*a) );
	}// for

	return GSL_SUCCESS;

}

// FDF Calibrator for TOF
int calib_fdf_TOF(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
	calib_f_TOF (x,params,f);
	calib_df_TOF (x,params,J);

	return GSL_SUCCESS;
}


double mz_recal(double peak)
{
	//FTICR calibration
	if(calib_mode == 1){
		return Ca/((1/peak)-Cb);
	}
	//Orbitrap calibration
	else if(calib_mode == 2){
		return Ca*peak;
	}
	//TOF calibration
	else
		return ( (sqrt(peak)-Cb) / Ca ) * ( (sqrt(peak)-Cb) / Ca );
}/* double mz_recal(double peak) */


/* compare the integers */
int sort_type_comp_inv_err(const void *i, const void *j)
{
	calibrant *ip, *jp;
	double erri, errj;

	ip = (calibrant*)i;
	jp = (calibrant*)j;

	errj = fabs((jp->mz- mz_recal(jp->peak))/jp->mz);
	erri = fabs((ip->mz - mz_recal(ip->peak))/ip->mz);

	if ((errj - erri) > 0)
		return 1;
	else if ((errj - erri) < 0)
		return -1;
	return 0;

}// int comp(const void *i, const void *j)

//int recalibratePeaks(msrecal_params* params){
int recalibratePeaks(msrecal_params* params){
	int status, SATISFIED, j;
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	double chi;

    size_t iter=0;
    //number of free parameters in calibration function
    //const size_t pp = 1;
    size_t pp;
	double y[MAX_CALIBRANTS];
	double mz2[MAX_CALIBRANTS];
	struct data d={y,mz2};
	double x_init[2]={1.0,0.0}; // start here, close to minimum if reasonably calibrated beforehand

	calib_mode = params->match;

	//set the number of free parameters based on the calibration function
	//FTICR and TOF calibration
	if(calib_mode == 1 || calib_mode == 3){
		pp=2;
	}
	//Orbitrap calibration
	else if(calib_mode == 2){
		pp=1;
	}

    gsl_multifit_function_fdf func;
    gsl_vector_view x=gsl_vector_view_array(x_init,pp);

	//Functions set based on the mass analyzer
	//FTICR calibration
	if(calib_mode == 1){
	    func.f = &calib_f_FTICR;
	    func.df = &calib_df_FTICR;
	    func.fdf = &calib_fdf_FTICR;
	}
	//Orbitrap calibration
	else if(calib_mode == 2){
	    func.f = &calib_f_FTMS;
	    func.df = &calib_df_FTMS;
		func.fdf = &calib_fdf_FTMS;
	}
	//TOF calibration
	else if(calib_mode == 3){
		func.f = &calib_f_TOF;
		func.df = &calib_df_TOF;
		func.fdf = &calib_fdf_TOF;
	}

    SATISFIED=0;
    while (n_calibrants >=params->min_cal && !SATISFIED) {
    	// least-squares fit first using all peaks, than removing those that don't fit
    	for (j=0;j<n_calibrants;j++) {
    		//Dummy units based on the mass analyzer
    		//FTICR calibration
    		if(calib_mode == 1){
    			d.y[j] = 1 / calibrant_list[j].peak;
    		}
    		//Orbitrap calibration
    		else if(calib_mode == 2){
    			d.y[j] = 1 / sqrt(calibrant_list[j].peak);
    		}
    		//TOF calibration
    		else if(calib_mode == 3){
    			d.y[j] = sqrt(calibrant_list[j].peak);
    		}
        	d.mz2[j] = calibrant_list[j].mz;
        }// for

        iter=0;
        T = gsl_multifit_fdfsolver_lmder;
        s = gsl_multifit_fdfsolver_alloc (T, n_calibrants, pp);
        func.n = n_calibrants;
        func.p = pp;
        func.params = &d;
        gsl_multifit_fdfsolver_set(s,&func,&x.vector);

        do {
        	iter++;
        	status = gsl_multifit_fdfsolver_iterate (s);

        	if (status)
        		break;
        	status=gsl_multifit_test_delta (s->dx, s->x, 1e-9, 1e-9);
        } while (status==GSL_CONTINUE && iter<500);

		//Dummy units based on the mass analyzer
		//FTICR calibration
		if(calib_mode == 1){
	        Ca = gsl_vector_get(s->x,0);
	        Cb = gsl_vector_get(s->x,1);
		}
		//Orbitrap calibration
		else if(calib_mode == 2){
			Ca = gsl_vector_get(s->x,0);
		}
		//TOF calibration
		else if(calib_mode == 3){
			Ca = gsl_vector_get(s->x,0);
			Cb = gsl_vector_get(s->x,1);
		}
        chi = gsl_blas_dnrm2(s->f);
        gsl_multifit_fdfsolver_free(s);

        // OK, that was one internal recalibration, now lets check if all calibrants are < INTERNAL_CALIBRATION_TARGET, if not, throw these out
        // and recalibrate (as long as we have at least three peaks)
        qsort(calibrant_list, n_calibrants, sizeof(calibrant), sort_type_comp_inv_err);



        for(j=n_calibrants-1; j>=0; j--)
        	if (fabs((calibrant_list[j].mz-mz_recal(calibrant_list[j].mz))/calibrant_list[j].mz)<INTERNAL_CALIBRATION_TARGET)
        		break;
        if (j==n_calibrants-1){
        	SATISFIED=1;
        }// all calibrants < INTERNAL_CALIBRATION_TARGET (e.g. 2.5 ppm)
        n_calibrants=j+1; // remove calibrants that doesn't fit CAL2 better than e.g. 2 ppm
    }

    if(SATISFIED){
    	printf("\tcalibrated\t1"); fflush(stdout);
    	printf("\tCa: %f\t", Ca); fflush(stdout);
    	if(calib_mode == 1 || calib_mode == 3){
    		printf("\tCb: %f\t", Cb); fflush(stdout);
    	}
    }

	return SATISFIED;
}

void applyCalibration(int scan, pscan_peaks mzpeaks)
{
	int j;
	/*
	printf("\nFinal calibration for scan %i:\n", scan);
	for(j=0;j<n_calibrants;j++) {
		printf("\t%f %f %f %.4f\n", calibrant_list[j].peak, calibrant_list[j].mz, mz_recal(calibrant_list[j].peak), 1e6*(mz_recal(calibrant_list[j].peak)-calibrant_list[j].mz)/calibrant_list[j].mz); fflush(stdout);
	}// for
	*/

	for (j=0; j<mzpeaks->count; j++) {
		//FTICR calibration
		if(calib_mode == 1){
			mzpeaks->mzs[j] = Ca/((1/mzpeaks->mzs[j])-Cb);
		}
		//Orbitrap calibration
		else if(calib_mode == 2){
			mzpeaks->mzs[j] = Ca * mzpeaks->mzs[j];
		}
		//TOF calibration
		else if(calib_mode == 3){
			mzpeaks->mzs[j] = ( (sqrt(mzpeaks->mzs[j])-Cb) / Ca ) * ( (sqrt(mzpeaks->mzs[j])-Cb) / Ca );
		}

	}// for

	// --------------------------------------------------------PRINT
	//printf("\tcalibrated\t1"); fflush(stdout);
	//printf("\tCa\t%.43f", Ca); fflush(stdout);

	calibrated_scans[calibration_counter] = scan;
	//FTICR, Orbitrap, TOF all have a Ca coefficient
	calibration_coefficients[calibration_counter][1] = Ca;
	//FTICR and TOF have a second calibrant
	if(calib_mode == 1 || calib_mode == 3 )
		calibration_coefficients[calibration_counter][2] = Cb;
	calibration_counter++;

}

double returnCoefficientCa(int scan){

	for(int i = 0; i <= MAX_ROWS; i++){
		if(calibrated_scans[i] == scan){
			return calibration_coefficients[i][1];
		}
	}
	return 0;
}

double returnCoefficientCb(int scan){

	for(int i = 0; i <= MAX_ROWS; i++){
		if(calibrated_scans[i] == scan){
			return calibration_coefficients[i][2];
		}
	}
	return 0;
}

