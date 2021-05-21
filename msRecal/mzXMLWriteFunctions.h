/* 
 * File:   mzXMLWriteFunctions.h
 * Author: Arzu Tugce Guler
 *
 * Created on 15 April 2016, 16:42
 */

#ifndef MZXMLWRITEFUNCTIONS_H
#define	MZXMLWRITEFUNCTIONS_H

#include <stdio.h>
#include "mzXMLStructures.h"

#define WRITE_BUFF_SIZE 129472	/* 128KB reading buffer as standard */
#ifdef WIN32
  #define endline "\r\n"
#else
  #define endline "\n"
#endif

/* Writes the entire file. Copies all residual info from the file recorded in mzxml_file */
void print_mzxml(pmzxml_file mzxml_file_hndl, char* output_file);

/* Writes all parentfile structures */
void print_parentfile_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin);

/* Writes all msinstrument structures */
void print_msinstrument_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin);

/* Writes a single msinstrument element */
void print_msinstrument_element(msinstrument instr, FILE* fout);

/* Writes all dataprocessing structures */
void print_dataprocessing_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin);

/* Writes a single dataprocessing element */ 
void print_dataprocessing_element(dataprocessing datapr, FILE* fout);

/* Writes all separation structures */
void print_separation_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin);

/* Writes all spotting structures */
void print_spotting_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin);

/* Writes all scan structures */
void print_scan_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin);

/* Prints a single scan element */
void print_scan_element(pscan mzscan, int scan_no, FILE* fout, FILE* fin);

/* Writes all scan origin structures */
void print_scan_origin_structure(pscan mzscan, FILE* fout, FILE* fin);

/* Writes all precursormz structures */
void print_precursormz_structure(pscan mzscan, FILE* fout, FILE* fin);

/* Writes all maldi structures */
void print_maldi_structure(pscan mzscan, FILE* fout, FILE* fin);

/* Writes all peaks structures */
void print_peaks_structure(pscan mzscan, FILE* fout, FILE* fin);

/* Writes all maldi structures */
void print_namevalue_structure(pscan mzscan, FILE* fout, FILE* fin);

/* Writes the sha1 structure */
void print_sha1(FILE** foutput, char* output_file);	

/* Writes the scan index structures */
void print_indices(pmzxml_file mzxml_file_hndl, FILE* fout);	

/* Writes the index offset structure */
void print_index_offset(pmzxml_file mzxml_file_hndl, FILE* fout);	

#endif	/* MZXMLWRITEFUNCTIONS_H */

