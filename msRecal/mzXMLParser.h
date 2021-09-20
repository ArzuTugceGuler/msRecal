/* 
 * File:   mzXMLParser.h
 * Author: Arzu Tugce Guler
 *
 * Created on 15 April 2016, 16:36
 */

#ifndef MZXMLPARSER_H
#define	MZXMLPARSER_H

#include <stdio.h>
#include "mzXMLStructures.h"
#include "mzXMLReaderDefs.h"

/*SIZES DOUBLED!*/
#define READ_TAIL_BUFF_SIZE		16384		/* 256 bytes tail size assumption */
#define SHORT_HEADER_BUFF_SIZE	64768	/* 4KB buffer size for short headers */
#define INDEX_REV_BUFF_SIZE		480		/* 24 bytes reverse side size assumption */
#define INDEX_ELT_BUFF_SIZE		1200		/* 60 bytes to find last index element */
#define READ_BUFF_SIZE			131072	/* 32KB reading buffer as standard */


/****************************************** Hidden parsing functions ***********************************************/
/* Parses the comlete mzXML file */;
void parse_mzxml_file(pmzxml_file file, FILE* finput, file_flags fflags, scan_config_flags sflags, int begin_scan, int end_scan);

/* Parses the header of the msrun section. also records the offsets of the parentfile and dataprocessing structures */
void parse_msrun_header(pmzxml_file file, FILE* finput, file_flags fflags);

/* Parses the entire parentfile structure */
void parse_parentfile_structure(pmzxml_file file, char* beginptr, FILE* finput);

/* Parses the entire msinstrument structure */
void parse_msinstrument_structure(pmzxml_file file, char* beginptr, FILE* finput);

/* Parses a single dataprocessing element */
msinstrument parse_msinstrument_element(char* beginptr, FILE* finput);

/* Parses the entire dataprocessing structure. */
void parse_dataprocessing_structure(pmzxml_file file, char* beginptr, FILE* finput);

/* Parses a single dataprocessing element */
dataprocessing parse_dataprocessing_element(char* beginptr, FILE* finput);

/* Parses the processingoperation structure */
void parse_processingoperation_structure(pdataprocessing dataproc, char* beginptr, FILE* finput);

/* Parses the entire separation structure */
void parse_separation_structure(pmzxml_file file, char* beginptr, FILE* finput);

/* Parses the entire spotting structure */
void parse_spotting_structure(pmzxml_file file, char* beginptr, FILE* finput);

/* Function that parses scan headers */
void parse_scan_header(pmzxml_file file, int scan_number, FILE* finput, scan_config_flags sflags);

/* Function that parses scan attributes */
void parse_scan_attributes(pscan scan, char* content);

/* Function that scans the origin structure */
void parse_scanorigin_structure(pscan scan, char* contents, FILE* finput);

/* Filters the info of a single scanorigin entry */
scan_origin parse_origin(char* content);

/* Function that scans the precursor structure */
void parse_precursor_structure(pscan scan, char* beginptr, FILE* finput);

/* Filters the info of a single precursor entry */
scan_precursor parse_precursor(char* beginptr, FILE* finput);

/* Filter the optional maldi part */
void parse_maldi_structure(pscan scan, char* contents, FILE* finput);

/* Function that filters peak attributes, and records the offset of the peaks in a file */
void parse_scan_peaks(pscan scan, char* beginptr, FILE* finput);

/* Function that filters peak attributes */
void parse_scan_peaks_header(pscan scan, char* beginptr);

/* Filter the optional maldi part */
void parse_namevalue_structure(pscan scan, char* beginptr, FILE* finput);

/* Parses the end of the scan segment */
void parse_scan_end(pmzxml_file mzxml_file, FILE* finput);

/* Parses search score sequence */
void parse_index_sequence(pmzxml_file mzxml_file, FILE* finput);

/* Handles file tail of the mzxml file */
void parse_file_tail(pmzxml_file mzxml_file, FILE* finput);

#endif	/* MZXMLPARSER_H */

