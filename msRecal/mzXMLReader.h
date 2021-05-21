/* 
 * File:   mzXMLReader.h
 * Author: Arzu Tugce Guler
 *
 * Created on 15 April 2016, 16:34
 */

#ifndef MZXMLREADER_H
#define	MZXMLREADER_H

#include <stdio.h>
#include "mzXMLStructures.h"
#include "mzXMLReaderDefs.h"

/* Main interface functions */
/* Parses mzxml file fully */
pmzxml_file read_mzxml_file(char* file_name, file_flags fflags, scan_config_flags sflags);

/* Parses subset of scans */
pmzxml_file read_mzxml_file_spectrum(char* file_name, file_flags fflags, scan_config_flags sflags, int begin_scan, int end_scan);

/* Unloads an mzXML file from memory */
void unload_mzxml_file(pmzxml_file mzXML_file);

/**************************************** File specific readers ***************************************************/
/* Gets the parentfile structure */
pparentfile get_parentfile_structure(pmzxml_file file, int* parentfile_count);

/* Gets the msinstrument structure */
pmsinstrument get_msinstrument_structure(pmzxml_file file, int* instrument_count);

/* Gets the dataprocessing structure */
pdataprocessing get_dataprocessing_structure(pmzxml_file file, int* dataprocessing_count);

/* Gets the separation structure */
pseparation get_separation_structure(pmzxml_file file);

/* Gets the spotting structure */
pspotting get_spotting_structure(pmzxml_file file);

/* Function that retrieves the whole scan */
pscan get_scan(pmzxml_file file, int scan_nr, scan_config_flags sflags);


/**************************************** Scan specific readers ***************************************************/
/* Function that retrieves a specific scan's attributes */
scan_attributes get_scan_attributes(pmzxml_file file, int scan_nr);

/* Function that retrieves a specific scan's origin structure */
pscan_origin get_scan_origin(pmzxml_file file, int scan_nr, int* origin_count);

/* Function that retrieves a specific scan's precursor array */
pscan_precursor	get_scan_precursor(pmzxml_file file, int scan_nr, int* precursor_count);

/* Function that retrieves a specific scan's maldi structure */
pscan_maldi get_scan_maldi(pmzxml_file file, int scan_nr);

/* Function that retrieves peaks from an mzxml */
scan_peaks load_scan_peaks(pmzxml_file file, int scan_nr);

/* Function that unloads the peaks in a scan */
void unload_scan_peaks(pmzxml_file file, int scan_nr);

#endif	/* MZXMLREADER_H */

