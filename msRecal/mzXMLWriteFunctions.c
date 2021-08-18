#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "mzXMLDefs.h"
#include "SHA.h"
#include "Base64Lib.h"
#include "mzXMLParser.h"
#include "mzXMLWriteFunctions.h"
#include "XMLFunctions.h"

static char write_buffer[WRITE_BUFF_SIZE];
static long* index_offset = NULL;
static long index_array_offset = 0;
static long tot_written_chars = 0;

/* Function that converts network order to host order */
static unsigned int ntohl(int nLongNumber)
{
   short word = 0x4321;
   if((*(char *)& word) != 0x21)
     return nLongNumber;
   else {
		return (((nLongNumber&0x000000FF)<<24)+((nLongNumber&0x0000FF00)<<8)+ ((nLongNumber&0x00FF0000)>>8)+((nLongNumber&0xFF000000)>>24));
   }
}

/* Calculates the length of the encoded string */
static int calc_encoding_length(pscan_peaks sp)
{
	int len;
	//if (sp->compressedLen >= 0){
	//CONDITION CHANGED!
	/*We do not work with compressed files, if the file is not compressed then
	 * this value is 0, so it has to be calculated. If it is compressed the
	 * meaning of this value has to be checked*/
	if (sp->compressedLen > 0)
		return sp->compressedLen;
	else {
		len = ((sp->count) * 2) * (sp->precision/8);
		
		return (((2 + len - ((len + 2) % 3)) * 4) / 3);
	}/* else */

}/* int calc_encoding_length(pscan_peaks sp) */


static char* read_xml_segment(FILE* fin, FILE* fout, long offset, long length)
{
	int read_chars = 0;
	char *beginptr, *endptr;

	fseek(fin, offset, SEEK_SET);
	read_chars = fread(write_buffer, sizeof(char), length, fin);
	write_buffer[read_chars] = '\0';

	beginptr = write_buffer;
	endptr = write_buffer + (read_chars);
	while (isspace(beginptr[0]))
		beginptr++;
	while (isspace(endptr[0]) || iscntrl(endptr[0]))
		endptr--;
	endptr[1] = '\0';

	return beginptr;

}/* char* read_xml_segment(FILE* fin, FILE* fout, long offset, long length) */

/* Copies an XML segment from input to output */
static void copy_xml_segment(FILE* fin, FILE* fout, long offset, long length)
{
	char* beginptr;
	
	beginptr = read_xml_segment(fin, fout, offset, length);
		
	tot_written_chars += fprintf(fout, "%s", beginptr);
	tot_written_chars += fprintf(fout, endline);

}/* char* copy_xml_segment(FILE* fin, FILE* fout, long offset, long length) */


/* Writes the entire file. Copies all residual info from the file recorded in mzxml_file */
void print_mzxml(pmzxml_file mzxml_file_hndl, char* output_file)
{
	FILE *fout, *fin; 
	int i;
	int newnum;
	
	fout = fopen(output_file, "wb");
	if (!fout)
		return;

	fin = fopen(mzxml_file_hndl->file_name, "rb");
	if (!fin)
		return;

	/* Print the same header that is present in the source file  */
	copy_xml_segment(fin, fout, 0, mzxml_file_hndl->msrun_offset);	

	/* Print the msrun header and the attributes */
	/* updating scan number attribute */
	newnum = mzxml_file_hndl->scan_num;
	for (i=0; i<mzxml_file_hndl->scan_num; i++) {
		if (mzxml_file_hndl->scan_array[i] && mzxml_file_hndl->scan_array[i]->deprecated)
		  newnum--;
	}// for

	tot_written_chars += fprintf(fout, "%s %s=\"%i\"", MZXML_MSRUN_OTAG, MZXML_MSRUN_ATTRIB_SCANCOUNT, newnum);
	if (mzxml_file_hndl->start_time >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_MSRUN_ATTRIB_STARTTIME, seconds_to_xml_duration(mzxml_file_hndl->start_time));
	if (mzxml_file_hndl->end_time >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_MSRUN_ATTRIB_ENDTIME, seconds_to_xml_duration(mzxml_file_hndl->end_time));
	tot_written_chars += fprintf(fout, "%s", XML_TAG_ETAG);
	tot_written_chars += fprintf(fout, endline);

	/* Printing the individual elements */
	print_parentfile_structure(mzxml_file_hndl, fout, fin);
	print_msinstrument_structure(mzxml_file_hndl, fout, fin);
	print_dataprocessing_structure(mzxml_file_hndl, fout, fin);
	print_separation_structure(mzxml_file_hndl, fout, fin);
	print_spotting_structure(mzxml_file_hndl, fout, fin);
	print_scan_structure(mzxml_file_hndl, fout, fin);

	/* Done with the input file */
	fclose(fin);

	/* End of msrun, writing closing tags */
	tot_written_chars += fprintf(fout, "%s%s", MZXML_MSRUN_CTAG, XML_TAG_ETAG);
	tot_written_chars += fprintf(fout, endline);

	/* Since this is the section where the indices come now, we save this as the offset and write the indices */
	index_array_offset = tot_written_chars;
	print_indices(mzxml_file_hndl, fout);

	/* Writing of the index */
	print_index_offset(mzxml_file_hndl, fout);

	/* Finish with the sa1 for the entire file */
	print_sha1(&fout, output_file); 

	/* Write the endtag for the mzxml file */
	tot_written_chars += fprintf(fout, "%s%s", MZXML_CTAG, XML_TAG_ETAG);
	tot_written_chars += fprintf(fout, endline);
	fclose(fout);
	
}/* void print_mzxml(pmzxml_file mzxml_file_hndl, char* output_file) */


/* Writes all parentfile structures */
void print_parentfile_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin)
{
	int i;

	/* Loaded in memory */
	if (mzxml_file_hndl->parentfile_count > 0) {
		for (i=0; i<mzxml_file_hndl->parentfile_count; i++) {
			tot_written_chars += fprintf(fout, "%s %s=\"%s\" %s=\"%s\" %s=\"%s\"%s%s", MZXML_PARENTFILE_OTAG,
				MZXML_PARENTFILE_ATTRIB_FILENAME, mzxml_file_hndl->parentfile_array[i].file_name,
				MZXML_PARENTFILE_ATTRIB_FILETYPE, mzxml_file_hndl->parentfile_array[i].file_type,
				MZXML_PARENTFILE_ATTRIB_FILESHA1, mzxml_file_hndl->parentfile_array[i].sha1,
				MZXML_PARENTFILE_CTAG, XML_TAG_ETAG);			
			tot_written_chars += fprintf(fout, endline);
		}/* for */
	}/* if */
	/* Present in the original file */
	else if (mzxml_file_hndl->parentfile_offset >= 0){
		copy_xml_segment(fin, fout, mzxml_file_hndl->parentfile_offset, mzxml_file_hndl->parentfile_length);
	}/* else */

}/* void print_parentfile_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin) */


/* Writes all msinstrument structures */
void print_msinstrument_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin)
{
	int i;

	/* Loaded in memory */
	if (mzxml_file_hndl->msinstrument_count > 0) {
		for (i=0; i<mzxml_file_hndl->msinstrument_count; i++) {
			print_msinstrument_element(mzxml_file_hndl->msinstrument_array[i], fout);		
		}/* for */
	}/* if */
	/* Present in the original file */
	else if (mzxml_file_hndl->msinstrument_offset >= 0){
		copy_xml_segment(fin, fout, mzxml_file_hndl->msinstrument_offset, mzxml_file_hndl->msinstrument_length);
	}/* else */

}/* void print_msinstrument_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin) */


/* Writes a single msinstrument element TO BE IMPLEMENTED */
void print_msinstrument_element(msinstrument instr, FILE* fout)
{


}


/* Writes all dataprocessing structures */
void print_dataprocessing_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin)
{
	int i;

	/* Loaded in memory */
	if (mzxml_file_hndl->dataprocessing_count > 0) {
		for (i=0; i<mzxml_file_hndl->dataprocessing_count; i++) {
			print_dataprocessing_element(mzxml_file_hndl->dataprocessing_array[i], fout);	
		}/* for */
	}/* if */
	/* Present in the original file */
	else if (mzxml_file_hndl->dataprocessing_offset >= 0){
		copy_xml_segment(fin, fout, mzxml_file_hndl->dataprocessing_offset, mzxml_file_hndl->dataprocessing_length);
	}/* else */

}/* void print_dataprocessing_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin) */


/* Writes a single dataprocessing element */
void print_dataprocessing_element(dataprocessing datapr, FILE* fout)
{
	int i;

	/* Writing the opening tag and required attributes */
	tot_written_chars += fprintf(fout, "%s", MZXML_DATAPROCESSING_OTAG);			

	/* Writing optional attributes */
	if (datapr.intensity_cutoff >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_DATAPROCESSING_ATTRIB_INTENSITYCUTOFF, datapr.intensity_cutoff);	
	if (datapr.centroided)
		tot_written_chars += fprintf(fout, " %s=\"1\"", MZXML_DATAPROCESSING_ATTRIB_CENTROIDED);	
	if (datapr.deisotoped)
		tot_written_chars += fprintf(fout, " %s=\"1\"", MZXML_DATAPROCESSING_ATTRIB_DEISOTOPED);	
	if (datapr.deconvoluted)
		tot_written_chars += fprintf(fout, " %s=\"1\"", MZXML_DATAPROCESSING_ATTRIB_CHARGEDECONVOLUTED);	
	if (datapr.spotintegration)
		tot_written_chars += fprintf(fout, " %s=\"1\"", MZXML_DATAPROCESSING_ATTRIB_SPOTINTEGRATION);	
	tot_written_chars += fprintf(fout, "%s", XML_TAG_ETAG);	
	tot_written_chars += fprintf(fout, endline);

	/* Writing elements */
	/* Writing Software structure: opening tag and required attributes */
	tot_written_chars += fprintf(fout, "%s %s=\"%s\" %s=\"%s\" %s=\"%s\"", MZXML_SOFTWARE_OTAG, 
		MZXML_SOFTWARE_ATTRIB_TYPE, datapr.sw_type,
		MZXML_SOFTWARE_ATTRIB_NAME, datapr.sw_name,
		MZXML_SOFTWARE_ATTRIB_VERSION, datapr.sw_version);

	/* Writing Software structure: optional attributes */
	if (datapr.sw_completion_time > 0)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_SOFTWARE_ATTRIB_COMPLETIONTIME, seconds_to_xml_duration(datapr.sw_completion_time));	
	tot_written_chars += fprintf(fout, "%s", XML_TAG_CTAG);
	tot_written_chars += fprintf(fout, endline);

	/* Printing processing operation array */
	for (i=0; i<datapr.proc_op_count; i++) {
		tot_written_chars += fprintf(fout, "%s", MZXML_PROCESSINGOPERATION_OTAG);
		if (datapr.proc_op_array[i].name)
			tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_PROCESSINGOPERATION_ATTRIB_NAME, datapr.proc_op_array[i].name);	
		if (datapr.proc_op_array[i].value)
			tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_PROCESSINGOPERATION_ATTRIB_VALUE, datapr.proc_op_array[i].value);
		if (datapr.proc_op_array[i].type)
			tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_PROCESSINGOPERATION_ATTRIB_TYPE, datapr.proc_op_array[i].type);
		tot_written_chars += fprintf(fout, "%s", XML_TAG_CTAG);		
		tot_written_chars += fprintf(fout, endline);
		if (datapr.proc_op_array[i].comment) {
			tot_written_chars += fprintf(fout, "%s%s%s%s%s", MZXML_COMMENT_OTAG, XML_TAG_ETAG, 
			datapr.proc_op_array[i].comment,
			MZXML_COMMENT_CTAG, XML_TAG_ETAG);	
			tot_written_chars += fprintf(fout, endline);
		}// if
	}/* for */

	/* Writing endtag */
	tot_written_chars += fprintf(fout, "%s%s", MZXML_DATAPROCESSING_CTAG, XML_TAG_ETAG);	
	tot_written_chars += fprintf(fout, endline);

}/* void print_dataprocessing_element(dataprocessing datapr, FILE* fout) */


/* Writes all separation structures */
void print_separation_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin)
{
	/* Loaded in memory */
	if (mzxml_file_hndl->separation != NULL) {
		/* TO BE IMPLEMENTED */		
	}/* if */
	/* Present in the original file */
	else if (mzxml_file_hndl->separation_offset >= 0){
		copy_xml_segment(fin, fout, mzxml_file_hndl->dataprocessing_offset, mzxml_file_hndl->dataprocessing_length);
	}/* else */

}/* void print_separation_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin) */ 

/* Writes all spotting structures */
void print_spotting_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin)
{
	/* Loaded in memory */
	if (mzxml_file_hndl->spotting != NULL) {
		/* TO BE IMPLEMENTED */		
	}/* if */
	/* Present in the original file */
	else if (mzxml_file_hndl->spotting_offset >= 0){
		copy_xml_segment(fin, fout, mzxml_file_hndl->spotting_offset, mzxml_file_hndl->spotting_length);
	}/* else */

}/* void print_spotting_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin) */


static void print_scan_closure(int prev_scan_lvl, int curr_scan_lvl, FILE* fout)
{
	int lvl_diff, i;

	if (prev_scan_lvl > 0) {
		lvl_diff = curr_scan_lvl - prev_scan_lvl;
		
		if (lvl_diff <= 0) {
			for (i=0; i<(-lvl_diff)+1; i++) {
				tot_written_chars += fprintf(fout, "%s%s", MZXML_SCAN_CTAG, XML_TAG_ETAG);	
				tot_written_chars += fprintf(fout, endline);
			}// for
		}// if
	}// if
}// static void print_scan_closure(int prev_scan_lvl, int curr_scan_lvl, FILE* fout)

static char* copy_whole_scan(pmzxml_file mzxml_file_hndl, int scan_index, FILE* fout, FILE* fin)
{
	char* beginptr, *endptr;
	int readlen;
	int endtag = 1;

	/* read the entire scan block */
	if (scan_index < mzxml_file_hndl->scan_num - 1)
		readlen = mzxml_file_hndl->index_array[scan_index+1] - mzxml_file_hndl->index_array[scan_index];
	else
		readlen = mzxml_file_hndl->scan_end_offset - mzxml_file_hndl->index_array[scan_index];
	beginptr = read_xml_segment(fin, fout, mzxml_file_hndl->index_array[scan_index], readlen);

	endptr = &(beginptr[strlen(beginptr)-1]); 
	while (endtag && (endptr > beginptr)) {
		if (endptr[0] == MZXML_ELEMENT_OTAG[0]) {
			if (strncmp(endptr, MZXML_SCAN_CTAG, strlen(MZXML_SCAN_CTAG))==0) {
				endptr--;
				while (isspace(endptr[0]) || iscntrl(endptr[0]))
					endptr--;
				endptr[1] = '\0';
			}// if
			else {
				endtag = 0;
			}
		}// if
		else {
			endptr--;
		}// else	
	}// while

	return beginptr;
	
}// static void print_scan_closure(int prev_scan_lvl, int curr_scan_lvl, FILE* fout)

/* Writes all scan structures */
void print_scan_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin)
{
	int i;
	int prevscan_lvl = -1, thisscan_lvl;	
	char* scantext, *mslvl;

	/* initializing the offset array */
	index_offset = (long*) malloc(sizeof(long) * mzxml_file_hndl->scan_num);

	for (i=0; i<mzxml_file_hndl->scan_num; i++) {
		if (!mzxml_file_hndl->scan_array[i] || mzxml_file_hndl->scan_array[i]->deprecated) {
			scantext = copy_whole_scan(mzxml_file_hndl, i, fout, fin);
			mslvl = get_xml_attribute_value(scantext, MZXML_SCAN_ATTRIB_MSLEVEL);
			thisscan_lvl = atoi(mslvl);
			print_scan_closure(prevscan_lvl, thisscan_lvl, fout);
			index_offset[i] = tot_written_chars;
			tot_written_chars += fprintf(fout, "%s", scantext);
			tot_written_chars += fprintf(fout, endline);
		}
		else {
			thisscan_lvl = mzxml_file_hndl->scan_array[i]->attributes.msLvl;
			print_scan_closure(prevscan_lvl, thisscan_lvl, fout);
			index_offset[i] = tot_written_chars;
			print_scan_element(mzxml_file_hndl->scan_array[i], i+1, fout, fin);
		}// else
		prevscan_lvl = thisscan_lvl;
	}/* for */
	print_scan_closure(prevscan_lvl, 1, fout);

}/* void print_scan_structure(pmzxml_file mzxml_file_hndl, FILE* fout, FILE* fin) */


/* Prints a single scan element */
void print_scan_element(pscan mzscan, int scan_no, FILE* fout, FILE* fin)
{
	/* Writing the opening tag and required attributes */
	tot_written_chars += fprintf(fout, "%s %s=\"%i\"", MZXML_SCAN_OTAG, MZXML_SCAN_ATTRIB_NUM, scan_no);	

	/* Writing optional attributes */
	if (mzscan->attributes.msLvl >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%i\"", MZXML_SCAN_ATTRIB_MSLEVEL, mzscan->attributes.msLvl);	
	if (mzscan->attributes.peakscount >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%i\"", MZXML_SCAN_ATTRIB_PEAKSCOUNT, mzscan->attributes.peakscount);	
	if (mzscan->attributes.polarity)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_SCAN_ATTRIB_POLARITY, mzscan->attributes.polarity);	
	if (mzscan->attributes.scanType)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_SCAN_ATTRIB_SCANTYPE, mzscan->attributes.scanType);	
	if (mzscan->attributes.filterLine)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_SCAN_ATTRIB_FILTERLINE, mzscan->attributes.filterLine);	
	if (mzscan->attributes.centroided)
		tot_written_chars += fprintf(fout, " %s=\"1\"", MZXML_SCAN_ATTRIB_CENTROIDED);	
	if (mzscan->attributes.deisotoped)
		tot_written_chars += fprintf(fout, " %s=\"1\"", MZXML_SCAN_ATTRIB_DEISOTOPED);	
	if (mzscan->attributes.chargeDeconvoluted)
		tot_written_chars += fprintf(fout, " %s=\"1\"", MZXML_SCAN_ATTRIB_CHARGEDECONVOLUTED);	
	if (mzscan->attributes.retentionTime >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_SCAN_ATTRIB_RETENTIONTIME, seconds_to_xml_duration(mzscan->attributes.retentionTime));
	if (mzscan->attributes.ionisationEnergy >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_SCAN_ATTRIB_IONISATIONENERGY, mzscan->attributes.ionisationEnergy);
	if (mzscan->attributes.collisionEnergy >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_SCAN_ATTRIB_COLLISIONENERGY, mzscan->attributes.collisionEnergy);
	if (mzscan->attributes.collisonGasPressure >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_SCAN_ATTRIB_CIDGASPRESSURE, mzscan->attributes.collisonGasPressure);
	if (mzscan->attributes.startMz >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_SCAN_ATTRIB_STARTMZ, mzscan->attributes.startMz);
	if (mzscan->attributes.endMz >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_SCAN_ATTRIB_ENDMZ, mzscan->attributes.endMz);
	if (mzscan->attributes.lowMz >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_SCAN_ATTRIB_LOWMZ, mzscan->attributes.lowMz);
	if (mzscan->attributes.highMz >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_SCAN_ATTRIB_HIGHMZ, mzscan->attributes.highMz);
	if (mzscan->attributes.basePeakMz >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_SCAN_ATTRIB_BASEPEAKMZ, mzscan->attributes.basePeakMz);
	if (mzscan->attributes.basePeakIntensity >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_SCAN_ATTRIB_BASEPEAKINTENSITY, mzscan->attributes.basePeakIntensity);
	if (mzscan->attributes.totIonCurrent >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_SCAN_ATTRIB_TOTIONCURRENT, mzscan->attributes.totIonCurrent);
	if (mzscan->attributes.msInstrumentID >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%i\"", MZXML_SCAN_ATTRIB_MSINSTRUMENTID, mzscan->attributes.msInstrumentID);
	tot_written_chars += fprintf(fout, "%s", XML_TAG_ETAG);
	tot_written_chars += fprintf(fout, endline);

	/* Writing elements */	
	print_scan_origin_structure(mzscan, fout, fin);
	print_precursormz_structure(mzscan, fout, fin);
	print_maldi_structure(mzscan, fout, fin);
	print_peaks_structure(mzscan, fout, fin);
	print_namevalue_structure(mzscan, fout, fin);

}/* void print_scan_element(pscan mzscan, FILE* fout, FILE* fin) */


/* Writes all scan origin structures */
void print_scan_origin_structure(pscan mzscan, FILE* fout, FILE* fin)
{
	int i;

	/* Loaded in memory */
	if (mzscan->origin_count > 0) {
		for (i=0; i<mzscan->origin_count; i++) {
			tot_written_chars += fprintf(fout, "%s %s=\"%s\" %s=\"%i\"%s", MZXML_SCANORIGIN_OTAG,
				MZXML_SCANORIGIN_ATTRIB_PARENTFILEID, mzscan->origin_array[i].parentFileID,
				MZXML_SCANORIGIN_ATTRIB_NUM, mzscan->origin_array[i].num,
				XML_TAG_CTAG);			
			tot_written_chars += fprintf(fout, endline);
		}/* for */
	}/* if */
	/* Present in the original file */
	else if (mzscan->origin_offset >= 0){
		copy_xml_segment(fin, fout, mzscan->origin_offset, mzscan->origin_length);
	}/* else */

}/* void print_scan_origin_structure(pscan mzscan, FILE* fout, FILE* fin) */

/* Prints a single pecursor segment */
void print_precursormz_element(scan_precursor prec, FILE* fout)
{
	tot_written_chars += fprintf(fout, "%s", MZXML_PRECURSORMZ_OTAG);

	/* Writing optional attributes */
	if (prec.precursorScanNum >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%i\"", MZXML_PRECURSORMZ_ATTRIB_PRECURSORSCANNUM, prec.precursorScanNum);	
	if (prec.precursorIntensity >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.16f\"", MZXML_PRECURSORMZ_ATTRIB_PRECURSORINTENSITY, prec.precursorIntensity);	
	if (prec.precursorCharge >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%i\"", MZXML_PRECURSORMZ_ATTRIB_PRECURSORCHARGE, prec.precursorCharge);	
	if (prec.possibleCharges)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_PRECURSORMZ_ATTRIB_POSSIBLECHARGES, prec.possibleCharges);	
	if (prec.windowWideness >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%.8f\"", MZXML_PRECURSORMZ_ATTRIB_WINDOWWIDENESS, prec.windowWideness);	
	if (prec.activationMethods)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_PRECURSORMZ_ATTRIB_ACTIVATIONMETHOD, prec.activationMethods);	
	tot_written_chars += fprintf(fout, "%s%.8f%s%s", XML_TAG_ETAG, prec.value, MZXML_PRECURSORMZ_CTAG, XML_TAG_ETAG);	
	tot_written_chars += fprintf(fout, endline);

}/* void print_precursormz_element(scan_precursor prec, FILE* fout) */


/* Writes all precursormz structures */
void print_precursormz_structure(pscan mzscan, FILE* fout, FILE* fin)
{
	int i; 

	if (mzscan->precursor_count > 0) {
		for (i=0; i<mzscan->precursor_count; i++) {
			print_precursormz_element(mzscan->precursor_array[i], fout);
		}/* for */
	}/* if */
	/* Present in the precursoral file */
	else if (mzscan->precursor_offset >= 0) {
		copy_xml_segment(fin, fout, mzscan->precursor_offset, mzscan->precursor_length);
	}/* else */

}/* void print_precursormz_structure(pscan mzscan, FILE* fout, FILE* fin) */


/* Writes all maldi structures */
void print_maldi_structure(pscan mzscan, FILE* fout, FILE* fin)
{
	if (mzscan->maldi) {
		/* NOT IMPLEMENTED YET */
	}/* if */
	else if (mzscan->maldi_offset >= 0) {
		copy_xml_segment(fin, fout, mzscan->maldi_offset, mzscan->maldi_length);
	}/* else */

}/* void print_maldi_structure(pscan mzscan, FILE* fout, FILE* fin) */


/* Prints the encoded peaks that are present in the peaks */
void print_encoded_peaks(pscan_peaks peaks, FILE* fout)
{
	float* encode_array_f = NULL;
	int* encode_array_i;
	double* encode_array_d = NULL;
	
	char* encoded_array;
	int i, j, len = 0, newlen, mod;

	/* Create and encode merged array */
	if (peaks->precision == 32) {
		len = 2*peaks->count * sizeof(float);
		mod = sizeof(float) / sizeof(int);
		encode_array_f = malloc(len);
		encode_array_i = (int*)encode_array_f;
		for (i=0; i<peaks->count; i++) {		
			encode_array_f[2*i] = (float)(peaks->mzs[i]);
			encode_array_f[(2*i)+1] = (float)(peaks->intensities[i]);
			for (j = 2*i*mod; j<2*(i+1)*mod; j++) {
				encode_array_i[j] = ntohl(encode_array_i[j]);
			}// for
		}/* for */
	}// if
	else if (peaks->precision == 64) {
		len = 2*peaks->count * sizeof(double);
		mod = sizeof(double) / sizeof(int);
		encode_array_d = malloc(len);
		encode_array_i = (int*)encode_array_d;
		for (i=0; i<peaks->count; i++) {		
			encode_array_d[2*i] = peaks->mzs[i];
			encode_array_d[(2*i)+1] = peaks->intensities[i];
			for (j = 2*i*mod; j<2*(i+1)*mod; j++) {
				encode_array_i[j] = ntohl(encode_array_i[j]);
			}// for
		}/* for */
	}// else if

	newlen = calc_encoding_length(peaks);
	encoded_array = malloc(newlen+1 * sizeof(char));
	encode_b64((unsigned char*)encode_array_f, encoded_array, len, newlen);
	tot_written_chars += fprintf(fout, "%s%s%s", encoded_array, MZXML_PEAKS_CTAG, XML_TAG_ETAG);
	tot_written_chars += fprintf(fout, endline);

	if (peaks->precision == 32) {
		free(encode_array_f);
	}// if
	else if (peaks->precision == 64) {
		free(encode_array_d);
	}// else if

	free(encoded_array);

}/* print_encoded_peaks(pscan_peaks peaks, FILE* fout) */


/* Writes all peaks structures */
void print_peaks_structure(pscan mzscan, FILE* fout, FILE* fin)
{
	/* Printing the header */
	tot_written_chars += fprintf(fout, "%s", MZXML_PEAKS_OTAG);

	/* Writing optional attributes */
	if (mzscan->peaks->precision >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%i\"", MZXML_PEAKS_ATTRIB_PRECISION, mzscan->peaks->precision);	
	if (mzscan->peaks->byteOrder)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_PEAKS_ATTRIB_BYTEORDER, mzscan->peaks->byteOrder);	
	if (mzscan->peaks->pairOrder)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_PEAKS_ATTRIB_PAIRORDER, mzscan->peaks->pairOrder);	
	if (mzscan->peaks->contentType)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_PEAKS_ATTRIB_CONTENTTYPE, mzscan->peaks->contentType);	
	if (mzscan->peaks->compressionType)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_PEAKS_ATTRIB_COMPRESSIONTYPE, mzscan->peaks->compressionType);	
	if (mzscan->peaks->compressedLen >= 0)
		tot_written_chars += fprintf(fout, " %s=\"%i\"", MZXML_PEAKS_ATTRIB_COMPRESSEDLEN, mzscan->peaks->compressedLen);	
	tot_written_chars += fprintf(fout, "%s", XML_TAG_ETAG);
	
	// no peaks, just print the endtag
	if (mzscan->attributes.peakscount == 0) {
		tot_written_chars += fprintf(fout, "%s%s", MZXML_PEAKS_CTAG, XML_TAG_ETAG);
		tot_written_chars += fprintf(fout, endline);
	}//if 
	else if (mzscan->peaks->mzs && mzscan->peaks->intensities && mzscan->peaks->precision) {
		print_encoded_peaks(mzscan->peaks, fout);
	}/* if */
	else if (mzscan->peak_content_offset) {
		copy_xml_segment(fin, fout, mzscan->peak_content_offset, calc_encoding_length(mzscan->peaks) + strlen(MZXML_PEAKS_CTAG) + 1);
	}/* else */


}/* void print_peaks_structure(pscan mzscan, FILE* fout, FILE* fin) */


/* Prints a single namevalue element */
void print_namevalue_element(namevalue nmval, FILE* fout)
{
	/* TO BE IMPLEMENTED */

}/* void print_namevalue_element(namevalue nmval, FILE* fout) */


/* Writes all maldi structures */
void print_namevalue_structure(pscan mzscan, FILE* fout, FILE* fin)
{
	int i;

	/* Loaded in memory */
	if (mzscan->namevalue_array_len > 0) {
		for (i=0; i<mzscan->namevalue_array_len; i++) {
			print_namevalue_element(mzscan->namevalue_array[i], fout);
		}/* for */
	}/* if */
	/* Present in the original file */
	else if (mzscan->namevalue_offset >= 0){
		copy_xml_segment(fin, fout, mzscan->namevalue_offset, mzscan->namevalue_length);
	}/* else */

}/* void print_namevalue_structure(pscan mzscan, FILE* fout, FILE* fin) */


/* Writes the sha1 structure */
void print_sha1(FILE** foutput, char* output_file)
{
	char* retval;

	/* Write the opening tag first */
	tot_written_chars += fprintf(*foutput, "%s%s", MZXML_SHA1_OTAG, XML_TAG_ETAG);

	/* Close the file, determine the SHA1 tage and reopen */
	fclose(*foutput);
	retval = get_sha1_string(output_file);

	/* Reopening file and printing the correct value */
	*foutput = fopen(output_file, "ab");
	tot_written_chars += fprintf(*foutput, "%s", retval);	
	free(retval);

	/* Writing endtag */
	tot_written_chars += fprintf(*foutput, "%s%s", MZXML_SHA1_CTAG, XML_TAG_ETAG);
	tot_written_chars += fprintf(*foutput, endline);

}/* void print_sha1(FILE** foutput, char* output_file) */


/* Writes the scan index structures */
void print_indices(pmzxml_file mzxml_file_hndl, FILE* fout)
{
	int i;

	if (mzxml_file_hndl->scan_num == 0)
		return;

	/* Printing opening tag and required attributes */
	tot_written_chars += fprintf(fout, "%s", MZXML_INDEX_OTAG);
	if (mzxml_file_hndl->index_name)
		tot_written_chars += fprintf(fout, " %s=\"%s\"", MZXML_INDEX_ATTRIB_NAME, mzxml_file_hndl->index_name);
	tot_written_chars += fprintf(fout, "%s", XML_TAG_ETAG);
	tot_written_chars += fprintf(fout, endline);

	/* Printing all indices */
	for (i=0; i<mzxml_file_hndl->scan_num; i++) {
		tot_written_chars += fprintf(fout, "%s %s=\"%i\"%s%d%s%s", MZXML_OFFSET_OTAG, 
			MZXML_OFFSET_ATTRIB_ID, i+1, XML_TAG_ETAG,
			index_offset[i], MZXML_OFFSET_CTAG, XML_TAG_ETAG);
			tot_written_chars += fprintf(fout, endline);
	}/* for */

	/* Printing closing tag and freeing up the memory of the offsets */
	tot_written_chars += fprintf(fout, "%s%s", MZXML_INDEX_CTAG, XML_TAG_ETAG);
	tot_written_chars += fprintf(fout, endline);
	free(index_offset);

}/* void print_indices(pmzxml_file mzxml_file_hndl, FILE* fout) */


/* Writes the index offset structure */
void print_index_offset(pmzxml_file mzxml_file_hndl, FILE* fout)
{
	if (mzxml_file_hndl->scan_num == 0) {
		tot_written_chars += fprintf(fout, "%s xsi:nil=\"1\"%s", MZXML_INDEXOFFSET_OTAG, XML_TAG_CTAG);
		tot_written_chars += fprintf(fout, endline);
	}/* if */
	else {
		tot_written_chars += fprintf(fout, "%s%s%d%s%s", MZXML_INDEXOFFSET_OTAG, XML_TAG_ETAG,
			index_array_offset, MZXML_INDEXOFFSET_CTAG, XML_TAG_ETAG);
		tot_written_chars += fprintf(fout, endline);
	}/* else */		
}/* void print_index_offset(FILE* fout) */

