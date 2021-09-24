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

    double calibration_coefficient_Ca;
    double calibration_coefficient_Cb;

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
    mzXML_file = read_mzxml_file_spectrum(params->mzxml_file, 2, scan_all_flag, params->ms_start_scan, params->ms_end_scan);

    if (mzXML_file == NULL) {
    	printf("error opening mzXML file\n"); fflush(stdout);
    	return -1;
    }

    if (mzXML_file->scan_num < params->ms_end_scan)
    	params->ms_end_scan = mzXML_file->scan_num;
    printf("\ndone.\n%i scans within the scan window: [%i, %i].\n", params->ms_end_scan - params->ms_start_scan + 1 , params->ms_start_scan , params->ms_end_scan); fflush(stdout);

    printf("\n>>Extracting instrument information from metadata..."); fflush(stdout);
    printf("\n%s: %s", mzXML_file->msinstrument_array[0].model_category, mzXML_file->msinstrument_array[0].model_value); fflush(stdout);
    printf("\n%s: %s\n", mzXML_file->msinstrument_array[0].massanalyzer_category, mzXML_file->msinstrument_array[0].massanalyzer_value); fflush(stdout);

    if(!params->mass_analyzer){
    	params -> mass_analyzer = mzXML_file->msinstrument_array[0].massanalyzer_value;
    }

    if(!params->instrument){
    	params -> instrument = mzXML_file->msinstrument_array[0].model_value;
    }

    printf("\n>>Choosing the most suitable calibration function for the data..."); fflush(stdout);
    processInstrument(params);

	if (!params->mass_analyzer && !params->instrument) {
		printf("\nNo valid mass analyzer type found!\n"); fflush(stdout);
	    return -1;
	}

    printf("\n>>Making internal calibrant candidate list for MS scans based on retention time...\n"); fflush(stdout);
    build_internal_calibrants(mzXML_file, peptide_set, pepnum, params);
    printf("done.\n"); fflush(stdout);

    printf("\n>>Recalibrating mzXML data..."); fflush(stdout);

    //int calibrated [params->ms_end_scan - params->ms_start_scan + 1][2];
    int counter = 0;
    int all_calibrated = 0;
    int ms1_calibrated = 0;
    int ms2_calibrated = 0;
    int previous_ms1 = 0;

    /*for(int i = 0; i <= params->ms_end_scan - params->ms_start_scan; i++){
    	calibrated[i][0] = -1;
    	calibrated[i][1] = -1;
    }*/

    // iterate each scan for calibration
    for(scan = params->ms_start_scan; scan <= params->ms_end_scan; scan++) {
    	//printf("\nscan\t%li", scan); fflush(stdout);
    	//load scan
    	mzscan = get_scan(mzXML_file, scan, 0);

    	//Calibrating ms2 precursors, it will only happen if the ms1 scan where the precursor is coming from is already calibrated
    	//Make sure not ms=1 level
        if(mzscan->attributes.msLvl != 1 ) {
        	//ms2 level
        	if(mzscan->attributes.msLvl == 2){
        		// --------------------------------------------------------PRINT
        		printf("\nscan\t%li\tms_level\t2", scan); fflush(stdout);
        		//printf("\tprecursor_scan\t%i\tprecursor_mz\t%12f", mzXML_file ->scan_array[scan-1]->precursor_array[0].precursorScanNum, mzXML_file ->scan_array[scan-1]->precursor_array[0].value); fflush(stdout);
        		//ms2 scans do not have their own calibrants, they get it from the respective ms1 scans
        		//only precursor masses are recalibrated with the respective ms1 scan coefficients
        		printf("\tcyclosiloxanes\t-1\tinit_calibrants\t-1"); fflush(stdout);

        		//This check is omited because some files do not have precursor scan info, so we assume the precursor is selected from the previous ms1 scan
        		//for(int i = 0; i <= params->ms_end_scan - params->ms_start_scan; i++){
        			//find the matching ms1 scan where the precursor is selected
        			//if(calibrated[i][0] == mzXML_file ->scan_array[scan-1]->precursor_array[0].precursorScanNum){
        				//check if that ms1 scan is calibrated (we will use it to calibrate the ms2 precursor mass)
        				if(previous_ms1 == 1){
        					// --------------------------------------------------------PRINT
        					printf("\tcalibrated\t1"); fflush(stdout);
        					//calibrated[counter][0] = scan;
        					//calibrated[counter][1] = 1;
        					//FTICR calibration
        					if(params->match == 1){
        						//calibration_coefficient_Ca = returnCoefficientCa(mzXML_file ->scan_array[scan-1]->precursor_array[0].precursorScanNum);
        						//calibration_coefficient_Cb = returnCoefficientCb(mzXML_file ->scan_array[scan-1]->precursor_array[0].precursorScanNum);
           						printf("\tCa\t%f", calibration_coefficient_Ca); fflush(stdout);
        						printf("\tCb\t%f", calibration_coefficient_Cb); fflush(stdout);
        						calibrated_precursor_mz = calibration_coefficient_Ca / ((1/mzXML_file ->scan_array[scan-1]->precursor_array[0].value)-calibration_coefficient_Cb);
        						printf("\tprecursor_before\t%f\tprecursor_after\t%f", mzXML_file ->scan_array[scan-1]->precursor_array[0].value, calibrated_precursor_mz); fflush(stdout);
        					}
        					//Orbitrap calibration
        					else if(params->match == 2){
        						//calibration_coefficient_Ca = returnCoefficientCa(mzXML_file ->scan_array[scan-1]->precursor_array[0].precursorScanNum);
        						printf("\tCa\t%f", calibration_coefficient_Ca); fflush(stdout);
        						calibrated_precursor_mz = calibration_coefficient_Ca * mzXML_file ->scan_array[scan-1]->precursor_array[0].value;
        						printf("\tprecursor_before\t%f\tprecursor_after\t%f", mzXML_file ->scan_array[scan-1]->precursor_array[0].value, calibrated_precursor_mz); fflush(stdout);
        					}
        					//TOF calibration
        					else if(params->match == 3){
        						//calibration_coefficient_Ca = returnCoefficientCa(mzXML_file ->scan_array[scan-1]->precursor_array[0].precursorScanNum);
        						//calibration_coefficient_Cb = returnCoefficientCb(mzXML_file ->scan_array[scan-1]->precursor_array[0].precursorScanNum);
           						printf("\tCa\t%f", calibration_coefficient_Ca); fflush(stdout);
        						printf("\tCb\t%f", calibration_coefficient_Cb); fflush(stdout);
        						calibrated_precursor_mz = ( (sqrt(mzXML_file ->scan_array[scan-1]->precursor_array[0].value) - calibration_coefficient_Cb) / calibration_coefficient_Ca ) * ( (sqrt(mzXML_file ->scan_array[scan-1]->precursor_array[0].value) - calibration_coefficient_Cb) / calibration_coefficient_Ca );
        						printf("\tprecursor_before\t%f\tprecursor_after\t%f", mzXML_file ->scan_array[scan-1]->precursor_array[0].value, calibrated_precursor_mz); fflush(stdout);
        					}
        					//update the precursor mass of the ms2 scan after applying the calibration function with the coefficients of the ms1 scan
        					mzXML_file ->scan_array[scan-1]->precursor_array[0].value = calibrated_precursor_mz;
        					ms2_calibrated++;
        					all_calibrated++;
        				}
        				//if the ms1 scan where the precursor is selected is not calibrated then do not calibrate ms2 precursor
        				else{
        					// --------------------------------------------------------PRINT
        					printf("\tcalibrated\t0"); fflush(stdout);
        					//calibrated[counter][0] = scan;
        					//calibrated[counter][1] = 0;
        				}
        			counter++;
        			//break;
        			//}// if precursor scan number match
        		//} for loop for all the previous scans
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
        //Loads the m/z-intensity pair of an ms1 spectrum
        mzpeaks = load_scan_peaks(mzXML_file, scan);

        if (mzpeaks.count == 0) {
        	//printf("\nScan: %li empty", scan); fflush(stdout);
            unload_scan_peaks(mzXML_file, scan);
            continue;
        }

        // Make the list of calibrants for this scan and sort them in descending order of intensity
        // --------------------------------------------------------PRINT
        printf("\nscan\t%li\tms_level\t1", scan); fflush(stdout);
        //Make the calibrant list for that ms1 scan
        makeCalibrantList(scan, &mzpeaks, peptide_set, params);
        // Recalibrate ms1 peaks
        SATISFIED = recalibratePeaks(params);

        if (SATISFIED) {
            applyCalibration(scan, &mzpeaks);
            update_scan_peaks(mzXML_file, scan, mzpeaks.count, 32, mzpeaks.mzs, mzpeaks.intensities);
            //calibrated[counter][0] = scan;
            //calibrated[counter][1] = 1;
            all_calibrated++;
            ms1_calibrated++;
            previous_ms1 = 1;
            calibration_coefficient_Ca = returnCoefficientCa(scan);
			if(params->match == 1 || params->match == 3){
				calibration_coefficient_Cb = returnCoefficientCb(scan);
			}
		}
		else {
			// --------------------------------------------------------PRINT
			//printf("\nTEST"); fflush(stdout);
			printf("\tcalibrated\t0"); fflush(stdout);
			//printf("\nNOT CALIBRATED"); fflush(stdout);
			//calibrated[counter][0] = scan;
			//calibrated[counter][1] = 0;
			previous_ms1 = 0;
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

