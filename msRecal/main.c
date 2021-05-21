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
    printf("\n>>Filtering peptides... "); fflush(stdout);
    peptide_set = build_peptide_set(pepXML_file, params, &pepnum);
    printf("done.\n%i peptides that are potential calibrants.\n", pepnum); fflush(stdout);	

    /*for(int i=0; i<pepnum; i++){
    	printf("\nPep no: %i.\n", i); fflush(stdout);
    	printf("\t%s\n", peptide_set[i].sequence); fflush(stdout);
    	printf("\t%f\n", peptide_set[i].retention); fflush(stdout);
    }*/

    printf("\n>>Reading mzXML file %s...",params->mzxml_file); fflush(stdout);
    mzXML_file = read_mzxml_file_spectrum(params->mzxml_file, 0, 0, params->ms_start_scan, params->ms_end_scan);
    if (mzXML_file == NULL) {
    	printf("error opening mzXML file\n"); fflush(stdout);
    	return -1;
    }
    if (mzXML_file->scan_num < params->ms_end_scan)
    	params->ms_end_scan = mzXML_file->scan_num;
    printf("\ndone.\n%i MS scans to be calibrated within the scan window: [%i, %i].\n", params->ms_end_scan - params->ms_start_scan + 1 , params->ms_start_scan , params->ms_end_scan); fflush(stdout);    

    /*
    printf("\n>>Making internal calibrant candidate list...\n"); fflush(stdout);
    build_internal_calibrants(mzXML_file, peptide_set, pepnum, params);

    printf("\ndone.\n"); fflush(stdout);

    /*
    printf("\n>>Recalibrating mzXML data...\n"); fflush(stdout);

    for(scan = params->ms_start_scan; scan <= params->ms_end_scan; scan++) {
    //for(scan = params->ms_start_scan; scan <= params->ms_start_scan +1; scan++) {
        mzscan = get_scan(mzXML_file, scan, 0);
        if(mzscan->attributes.msLvl != 1 ) 
            continue;
        // If crop param is 1, empty invalid scans
	if(!mzscan->peaks || mzscan->peaks->count < 0) {
            if (params->crop)
		empty_scan(mzXML_file, scan);
	    continue;
	}
        // Loading and filtering the peaks
        mzpeaks = load_scan_peaks(mzXML_file, scan);
        //printf("\nLoaded scan: %i\t Peak count: %i\n", scan, mzpeaks.count); fflush(stdout);
	if (mzpeaks.count == 0) {
            unload_scan_peaks(mzXML_file, scan);
            continue;
	}
        // Make the list of calibrants for this spectrum  and sort them in descending order of intensity  
        makeCalibrantList(scan, &mzpeaks, peptide_set, params);
        
        // Recalibrate peaks
	SATISFIED = recalibratePeaks(params);
        
        if (SATISFIED) {
            applyCalibration(scan, &mzpeaks);	
            update_scan_peaks(mzXML_file, scan, mzpeaks.count, 32, mzpeaks.mzs, mzpeaks.intensities);
            printf("-------------------------------------------------------------\n"); fflush(stdout);
	}// if
	else {
            unload_scan_peaks(mzXML_file, scan);
            if (params->crop)
                empty_scan(mzXML_file, scan);						
	}

    } 
    
    printf("done.\n\n>>Writing recalibrated mzXML file to destination %s...", params->output_mzXML_file); fflush(stdout);
    write_mzxml_file(mzXML_file, params->output_mzXML_file);
    printf("\ndone.\n"); fflush(stdout);
*/
    return 0;	
}

