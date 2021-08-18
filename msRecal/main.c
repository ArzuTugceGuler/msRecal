/* 
 * File:   main.c
 * Author: Arzu Tugce Guler
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "msRecal.h"
#include "pepXMLReader.h"
#include "PeptideProphetDelegate.h"
#include "mzXMLStructures.h"
#include "mzXMLWriter.h"



int main(int argc, char *argv[]) {

    msrecal_params* params;
    pdelegate_list dlgl;
    pdelegate_type dlg;
    pmsms_pipeline_analysis pepXML_file;
    pmzxml_file mzXML_file;  
    pscan mzscan;
    scan_peaks mzpeaks;
    peptideset* peptide_set;


    int pepnum;
    long scan;
    
    int SATISFIED;

    double calibration_coefficient;
    double calibrated_precursor_mz;

    params = readParameters(argc, argv);
    
    dlgl = make_delegate_list();
    dlg = make_peptide_prophet_result_delegate();
    add_delegate_to_list(dlgl, dlg);

    printf("\n\n>>Reading pepXML file %s...", params->pepxml_file);
    fflush(stdout);

    pepXML_file = read_pepxml_file(params->pepxml_file, 0, 0, dlgl);

    if (pepXML_file == NULL) {
    	printf("\nerror opening pepXML file\n"); fflush(stdout);
    	return -1;
    }// if

    if (!params){
    	showHelp();
        return -1;
    }
    printf("\ndone.\n" ); fflush(stdout);

	// Filtering out peptides and making them unique
    printf("\n>>Filtering peptides from pepXML... "); fflush(stdout);
    peptide_set = build_peptide_set(pepXML_file, params, &pepnum);
    printf("done.\n%i peptides that are potential calibrants.\n", pepnum); fflush(stdout);	

    //prints the potential peptide sequences and retention times
    /*
    for(int i=0; i<pepnum; i++){
    	printf("\nPep no: %i.\n", i+1); fflush(stdout);
    	printf("\t%s\n", peptide_set[i].sequence); fflush(stdout);
    	printf("\t%f\n", peptide_set[i].retention); fflush(stdout);
    }
    */

    printf("\n>>Reading mzXML file %s...",params->mzxml_file); fflush(stdout);
    //mzXML_file = read_mzxml_file_spectrum(params->mzxml_file, 0, 0, params->ms_start_scan, params->ms_end_scan);
    mzXML_file = read_mzxml_file_spectrum(params->mzxml_file, 0, scan_all_flag, params->ms_start_scan, params->ms_end_scan);

    if (mzXML_file == NULL) {
    	printf("error opening mzXML file\n"); fflush(stdout);
    	return -1;
    }

    if (mzXML_file->scan_num < params->ms_end_scan)
    	params->ms_end_scan = mzXML_file->scan_num;
    printf("\ndone.\n%i scans within the scan window: [%i, %i].\n", params->ms_end_scan - params->ms_start_scan + 1 , params->ms_start_scan , params->ms_end_scan); fflush(stdout);


    printf("\n>>Making internal calibrant candidate list for MS scans based on retention time...\n"); fflush(stdout);
    build_internal_calibrants(mzXML_file, peptide_set, pepnum, params);
    printf("done.\n"); fflush(stdout);

    printf("\n>>Recalibrating mzXML data..."); fflush(stdout);

    int calibrated [params->ms_end_scan - params->ms_start_scan + 1][2];
    int counter = 0;
    int all_calibrated = 0;
    int ms1_calibrated = 0;
    int ms2_calibrated = 0;

    for(int i = 0; i <= params->ms_end_scan - params->ms_start_scan; i++){
    	calibrated[i][0] = -1;
    	calibrated[i][1] = -1;
    }

    // iterate each ms scan for calibration
    for(scan = params->ms_start_scan; scan <= params->ms_end_scan; scan++) {
    	//printf("\nscan\t%li", scan); fflush(stdout);
    	mzscan = get_scan(mzXML_file, scan, 0);

        if(mzscan->attributes.msLvl != 1 ) {
        	if(mzscan->attributes.msLvl == 2){
        		// --------------------------------------------------------PRINT
        		printf("\nscan\t%li\tms_level\t2", scan); fflush(stdout);
        		printf("\tprecursor_scan\t%i\tprecursor_mz\t%12f", mzXML_file ->scan_array[scan-1]->precursor_array[0].precursorScanNum, mzXML_file ->scan_array[scan-1]->precursor_array[0].value); fflush(stdout);
        		printf("\tcyclosiloxanes\t-1\tinit_calibrants\t-1"); fflush(stdout);

        		for(int i = 0; i <= params->ms_end_scan - params->ms_start_scan; i++){
        			if(calibrated[i][0] == mzXML_file ->scan_array[scan-1]->precursor_array[0].precursorScanNum){
        				if(calibrated[i][1] == 1){
        					// --------------------------------------------------------PRINT
        					printf("\tcalibrated\t1"); fflush(stdout);
        					calibrated[counter][0] = scan;
        					calibrated[counter][1] = 1;
        					calibration_coefficient = returnCoefficient(mzXML_file ->scan_array[scan-1]->precursor_array[0].precursorScanNum);
        					printf("\tCa\t%.43f", calibration_coefficient); fflush(stdout);
        					calibrated_precursor_mz = calibration_coefficient * mzXML_file ->scan_array[scan-1]->precursor_array[0].value;
        					printf("\tprecursor_before\t%.43f\tprecursor_after\t%.43f", mzXML_file ->scan_array[scan-1]->precursor_array[0].value, calibrated_precursor_mz); fflush(stdout);
        					//update the precursor mass of the ms2 scan
        					mzXML_file ->scan_array[scan-1]->precursor_array[0].value = calibrated_precursor_mz;
        					ms2_calibrated++;
        					all_calibrated++;
        				}
        				else{
        					// --------------------------------------------------------PRINT
        					printf("\tcalibrated\t0"); fflush(stdout);
        					calibrated[counter][0] = scan;
        					calibrated[counter][1] = 0;
        				}
        			counter++;
        			break;
        			}
        		}
        	}
            continue;
        }

        // If crop param is 1, empty invalid scans
        if(!mzscan->peaks || mzscan->peaks->count < 0) {
        	//printf("\nInvalid scan!\n"); fflush(stdout);
        	if (params->crop)
        		empty_scan(mzXML_file, scan);
        	continue;
        }

        //Loading and filtering the peaks
        mzpeaks = load_scan_peaks(mzXML_file, scan);

        if (mzpeaks.count == 0) {
        	//printf("\nScan: %li empty", scan); fflush(stdout);
            unload_scan_peaks(mzXML_file, scan);
            continue;
        }

        // Make the list of calibrants for this scan and sort them in descending order of intensity
        // --------------------------------------------------------PRINT
        printf("\nscan\t%li\tms_level\t1", scan); fflush(stdout);
        printf("\tprecursor_scan\t-1\tprecursor_mz\t-1"); fflush(stdout);
        makeCalibrantList(scan, &mzpeaks, peptide_set, params);
        // Recalibrate peaks
        SATISFIED = recalibratePeaks(params);

        if (SATISFIED) {
            applyCalibration(scan, &mzpeaks);
            update_scan_peaks(mzXML_file, scan, mzpeaks.count, 32, mzpeaks.mzs, mzpeaks.intensities);
            calibrated[counter][0] = scan;
            calibrated[counter][1] = 1;
            all_calibrated++;
            ms1_calibrated++;
		}
		else {
			// --------------------------------------------------------PRINT
			printf("\tcalibrated\t0"); fflush(stdout);
			calibrated[counter][0] = scan;
			calibrated[counter][1] = 0;
            unload_scan_peaks(mzXML_file, scan);
            if (params->crop)
                empty_scan(mzXML_file, scan);						
		}

        counter ++;

    }

	printf("\n\ndone.\n\n>>Writing recalibrated mzXML file to destination %s...", params->output_mzXML_file); fflush(stdout);
    write_mzxml_file(mzXML_file, params->output_mzXML_file);
    printf("\ndone.\n"); fflush(stdout);

    printf("\nTotal number of scans checked for calibration: %i\n", counter); fflush(stdout);
    printf("\nTotal number of calibrated ms1 scans: %i\n", ms1_calibrated ); fflush(stdout);
    printf("\nTotal number of calibrated ms2 precursors: %i\n", ms2_calibrated ); fflush(stdout);
    printf("\nTotal number of calibrated scans: %i\n", all_calibrated ); fflush(stdout);

    return 0;	
}

