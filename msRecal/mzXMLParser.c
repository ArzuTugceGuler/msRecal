#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "mzXMLDefs.h"
#include "mzXMLParser.h"
#include "XMLFunctions.h"
#include "Base64Lib.h"

static char read_buffer[READ_BUFF_SIZE];
static long  offset = 0;
static long  read_chars = 0;
static char* walkptr = NULL;

static int parentfile_alloc_count = 1;
static int dataprocessing_alloc_count = 1;
static int proc_op_alloc_count = 1;
static int origin_alloc_count = 1;
static int precursor_alloc_count = 1;
static int namevalue_alloc_count = 1;
static int scan_index_begin = -1;


/* Calculates the length of the encoded string */
//same function also in mzXMLReader.c
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


/* Parses the comlete mzXML file */;
void parse_mzxml_file(pmzxml_file file, FILE* finput, file_flags fflags, scan_config_flags sflags, int begin_scan, int end_scan)
{
	int topscan, i;

	parse_file_tail(file, finput);
	parse_index_sequence(file, finput);
	parse_scan_end(file, finput);;
	parse_msrun_header(file, finput, fflags);

	file->scan_array = (pscan*) malloc(file->scan_num * sizeof(pscan));
	memset(file->scan_array, 0, file->scan_num * sizeof(pscan));

	topscan = (end_scan < file->scan_num? end_scan: file->scan_num);

	for (i=begin_scan; i<=topscan; i++) {
		file->scan_array[i-1] = (pscan) malloc(sizeof(scan));
		parse_scan_header(file, i, finput, sflags);
	}/* for */

}/* void parse_mzxml_file(pmzxml_file file, FILE* finput, file_flags fflags, scan_config_flags sflags, int begin_scan, int end_scan) */


/* Sets all arrays and attribute values to undefined (-1) or 0(NULL) */
void set_msrun_header_defaults(pmzxml_file file)
{
	/* Setting the defaults for the optional attributes and elements */
	file->msrun_offset = -1;
	file->start_time = -1;
	file->end_time = -1;

	file->parentfile_offset = -1;
	file->parentfile_length = 0;
	file->parentfile_array = NULL;
	file->parentfile_count = 0;

	file->msinstrument_offset = -1;
	file->msinstrument_length = 0;
	file->msinstrument_array = NULL;
	file->msinstrument_count = 0;

	file->dataprocessing_offset = -1;
	file->dataprocessing_length = 0;
	file->dataprocessing_array = NULL;
	file->dataprocessing_count = 0;

	file->separation_offset = -1;
	file->separation_length = 0;
	file->separation = NULL;

	file->spotting_offset = -1;
	file->spotting_length = 0;
	file->spotting = NULL;

}/* void set_msrun_header_defaults(pmzxml_file file) */


/* Calculates the lengths of the segments */
void set_msrun_lengths(pmzxml_file file, long endpos)
{
	long offset_array[6];
	long length_array[5];
	int i, j;

	offset_array[0] = file->parentfile_offset;
	offset_array[1] = file->msinstrument_offset;
	offset_array[2] = file->dataprocessing_offset;
	offset_array[3] = file->separation_offset;
	offset_array[4] = file->spotting_offset;
	offset_array[5] = endpos;

	for (i=0; i<5; i++) {
		if (offset_array[i] == -1)
			length_array[i] = 0;
		else {
			j = i+1;
			while (offset_array[j] == -1 && j < 6)
				j += 1;
			length_array[i] = (offset_array[j] - offset_array[i])-1;
		}/* else */
	}/* for */

	file->parentfile_length = length_array[0];
	file->msinstrument_length = length_array[1];
	file->dataprocessing_length = length_array[2];
	file->separation_length = length_array[3];
	file->spotting_length = length_array[4];	

}/* void set_msrun_lengths(pmzxml_file file) */


/* Parses the header of the msrun section. also records the offsets of the parentfile and dataprocessing structures */
void parse_msrun_header(pmzxml_file file, FILE* finput, file_flags fflags)
{	
	char *tag, *tmpbuffer;
	int test;

	rewind(finput);
	offset = 0;
	read_chars = fread(read_buffer, sizeof(char), SHORT_HEADER_BUFF_SIZE-1, finput);
	read_buffer[read_chars] = '\0';	
	walkptr = read_buffer;

	/* Setting the defaults of the struct */
	set_msrun_header_defaults(file);

	tag = get_xml_tag(read_buffer, &walkptr, finput, SHORT_HEADER_BUFF_SIZE, &offset);
	test = offset + ((int)tag - (int)read_buffer);
	while (strstr(tag, MZXML_SCAN_OTAG) == NULL) {
		/* Find msrun tag, so filling attributes */
		if (strstr(tag, MZXML_MSRUN_OTAG)) {	
			file->msrun_offset = offset + ((int)tag - (int)read_buffer);
			tmpbuffer = get_xml_attribute_value(tag, MZXML_MSRUN_ATTRIB_STARTTIME);
			if (tmpbuffer) {
				file->start_time = xml_duration_to_seconds(tmpbuffer);
				free(tmpbuffer);
			}/* if */
			tmpbuffer = get_xml_attribute_value(tag, MZXML_MSRUN_ATTRIB_ENDTIME);
			if (tmpbuffer) {
				file->end_time = xml_duration_to_seconds(tmpbuffer);
				free(tmpbuffer);
			}/* if */
		}/* if */
		/* Parsing the parentfile structure */
		else if (strstr(tag, MZXML_PARENTFILE_OTAG)) {
			if (file->parentfile_offset < 0)
				file->parentfile_offset = offset + ((int)tag - (int)read_buffer);			
			if (fflags & parentfile_flag) 
				parse_parentfile_structure(file, tag, finput);			
		}/* else if */
		else if (strstr(tag, MZXML_MSINSTRUMENT_OTAG)) {
			if (file->msinstrument_offset < 0)
				file->msinstrument_offset = offset + ((int)tag - (int)read_buffer);
			if (fflags & instrument_flag) 
				parse_msinstrument_structure(file, tag, finput);
		}/* else if */
		else if (strstr(tag, MZXML_DATAPROCESSING_OTAG)) {
			if (file->dataprocessing_offset < 0)
				file->dataprocessing_offset = offset + ((int)tag - (int)read_buffer);
			if (fflags & dataprocessing_flag) 
				parse_dataprocessing_structure(file, tag, finput);
		}/* else if */
		else if (strstr(tag, MZXML_SEPARATION_OTAG)) {
			if (file->separation_offset < 0)
				file->separation_offset = offset + ((int)tag - (int)read_buffer);		
			if (fflags & separation_flag) 
				parse_separation_structure(file, tag, finput);
		}/* else if */
		else if (strstr(tag, MZXML_SPOTTING_OTAG)) {
			if (file->spotting_offset < 0)
				file->spotting_offset = offset + ((int)tag - (int)read_buffer);
			if (fflags & spotting_flag) 
				parse_spotting_structure(file, tag, finput);
		}/* else if */
		tag = get_xml_tag(read_buffer, &walkptr, finput, SHORT_HEADER_BUFF_SIZE, &offset);
	}/* while */

	/* Calculating lengths of segments */
	set_msrun_lengths(file, offset + ((int)tag - (int)read_buffer));

	/* Adjusting the size of all created arrays to reflect their actual length */
	if (file->parentfile_count > 0)
		file->parentfile_array = (pparentfile) realloc(file->parentfile_array, file->parentfile_count * sizeof(parentfile));
	if (file->msinstrument_count > 0)
		file->msinstrument_array = (pmsinstrument) realloc(file->msinstrument_array, file->msinstrument_count * sizeof(msinstrument));
	if (file->dataprocessing_count > 0)
		file->dataprocessing_array = (pdataprocessing) realloc(file->dataprocessing_array, file->dataprocessing_count * sizeof(dataprocessing));

}/* void parse_msrun_header(pmzxml_file file, FILE* finput, file_flags fflags) */


/* Parses the parentfile structure */
void parse_parentfile_structure(pmzxml_file file, char* beginptr, FILE* finput)
{
	if (!(file->parentfile_array)) {
		file->parentfile_array = malloc(sizeof(parentfile));
		file->parentfile_count = 0;	
	}/* if */

	/* Retrieving the parentfile attributes */
	file->parentfile_array[file->parentfile_count].file_name = get_xml_attribute_value(beginptr, MZXML_PARENTFILE_ATTRIB_FILENAME);
	file->parentfile_array[file->parentfile_count].file_type = get_xml_attribute_value(beginptr, MZXML_PARENTFILE_ATTRIB_FILETYPE);
	file->parentfile_array[file->parentfile_count].sha1 = get_xml_attribute_value(beginptr, MZXML_PARENTFILE_ATTRIB_FILESHA1);
				
	/* Updating counters and making more room for future structures */
	file->parentfile_count += 1;
	if (file->parentfile_count == parentfile_alloc_count) {
		parentfile_alloc_count *= 2;
		file->parentfile_array = realloc(file->parentfile_array, parentfile_alloc_count * sizeof(parentfile));
	}/* if */		

}/* pparentfile get_parentfile_structure(pmzxml_file file, int* parentfile_count) */


/* Parses the entire msinstrument structure TO BE IMPLEMENTED */
void parse_msinstrument_structure(pmzxml_file file, char* beginptr, FILE* finput)
{

}/* pmsinstrument get_msinstrument_structure(pmzxml_file file, int* msinstrument_count) */


/* Parses the dataprocessing structure */
void parse_dataprocessing_structure(pmzxml_file file, char* beginptr, FILE* finput)
{
	if (!(file->dataprocessing_array)) {
		file->dataprocessing_array = malloc(sizeof(dataprocessing));
		file->dataprocessing_count = 0;	
	}/* if */

	/* Retrieving the dataprocessing attributes */
	file->dataprocessing_array[file->dataprocessing_count] = parse_dataprocessing_element(beginptr, finput);
				
	/* Updating counters and making more room for future structures */
	file->dataprocessing_count += 1;
	if (file->dataprocessing_count == dataprocessing_alloc_count) {
		dataprocessing_alloc_count *= 2;
		file->dataprocessing_array = realloc(file->dataprocessing_array, dataprocessing_alloc_count * sizeof(dataprocessing));
	}/* if */		

}/* void parse_dataprocessing_structure(pmzxml_file file, char* contents) */


/* Parses a single dataprocessing element */
dataprocessing parse_dataprocessing_element(char* beginptr, FILE* finput)
{
	dataprocessing dp;
	char* tmpbuffer, *tag;

	/* Setting optional values to 0 or NULL */
	dp.intensity_cutoff = -1;
	dp.centroided = 0;
	dp.deisotoped = 0;
	dp.deconvoluted = 0;
	dp.spotintegration = 0;

	dp.proc_op_array = NULL;
	dp.proc_op_count = 0;

	/* Filter intensity cutoff attribute */
	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_DATAPROCESSING_ATTRIB_INTENSITYCUTOFF);
	if (tmpbuffer) {
		dp.intensity_cutoff = strtod(tmpbuffer, NULL);
		free(tmpbuffer);
	}/* if */

	/* Filter centroided attribute */
	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_DATAPROCESSING_ATTRIB_CENTROIDED);
	if (tmpbuffer) {
		dp.centroided = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */			

	/* Filter deisotoped attribute */
	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_DATAPROCESSING_ATTRIB_DEISOTOPED);
	if (tmpbuffer) {
		dp.deisotoped = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */

	/* Filter charge deconvoluted attribute */
	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_DATAPROCESSING_ATTRIB_CHARGEDECONVOLUTED);
	if (tmpbuffer) {
		dp.deconvoluted = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */

	/* Filter spot integration attribute */
	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_DATAPROCESSING_ATTRIB_SPOTINTEGRATION);
	if (tmpbuffer) {
		dp.spotintegration = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */

	tag = get_xml_tag(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
	while (strstr(tag, MZXML_DATAPROCESSING_CTAG) == NULL) {
		if (strstr(tag, MZXML_SOFTWARE_OTAG)) {
			dp.sw_type = get_xml_attribute_value(tag, MZXML_SOFTWARE_ATTRIB_TYPE);
			dp.sw_name = get_xml_attribute_value(tag, MZXML_SOFTWARE_ATTRIB_NAME);
			dp.sw_version = get_xml_attribute_value(tag, MZXML_SOFTWARE_ATTRIB_VERSION);
			tmpbuffer = get_xml_attribute_value(tag, MZXML_SOFTWARE_ATTRIB_COMPLETIONTIME);
			if (tmpbuffer) {
				dp.sw_completion_time = xml_duration_to_seconds(tmpbuffer);
				free(tmpbuffer);
			}/* if */
			else
				dp.sw_completion_time = -1;
		}/* if */
		else if (strstr(tag, MZXML_PROCESSINGOPERATION_OTAG)) {
			parse_processingoperation_structure(&dp, tag, finput);
		}/* else if */
		else if (strstr(tag, MZXML_COMMENT_OTAG)) {
			dp.proc_op_array[dp.proc_op_count].comment = get_xml_tag_value(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
		}/* else if */
		tag = get_xml_tag(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
	}/* while */
	
	return dp;	

}/* dataprocessing parse_dataprocessing(char* content) */


/* Parses the processingoperation structure */
void parse_processingoperation_structure(pdataprocessing dataproc, char* beginptr, FILE* finput)
{
	if (!(dataproc->proc_op_array)) {
		dataproc->proc_op_array = malloc(sizeof(processingoperation));
		dataproc->proc_op_count = 0;	
	}/* if */

	/* Retrieving the processingoperation attributes */
	dataproc->proc_op_array[dataproc->proc_op_count].name = get_xml_attribute_value(beginptr, MZXML_PROCESSINGOPERATION_ATTRIB_NAME);
	dataproc->proc_op_array[dataproc->proc_op_count].type = get_xml_attribute_value(beginptr, MZXML_PROCESSINGOPERATION_ATTRIB_TYPE);
	dataproc->proc_op_array[dataproc->proc_op_count].value = get_xml_attribute_value(beginptr, MZXML_PROCESSINGOPERATION_ATTRIB_VALUE);
				
	/* Updating counters and making more room for future structures */
	dataproc->proc_op_count += 1;
	if (dataproc->proc_op_count == proc_op_alloc_count) {
		proc_op_alloc_count *= 2;
		dataproc->proc_op_array = realloc(dataproc->proc_op_array, proc_op_alloc_count * sizeof(processingoperation));
	}/* if */		

}/* void parse_dataprocessing_structure(pmzxml_file file, char* contents) */


/* Parses the entire separation structure TO BE IMPLEMENTED */
void parse_separation_structure(pmzxml_file file, char* beginptr, FILE* finput)
{


}/* void parse_separation_structure(pmzxml_file file, char* contents) */


/* Parses the entire spotting structure TO BE IMPLEMENTED */
void parse_spotting_structure(pmzxml_file file, char* beginptr, FILE* finput)
{
	

}/* void parse_spotting_structure(pmzxml_file file, char* contents) */


/* Sets all arrays and attribute values to undefined (-1) or 0(NULL) */
void set_scan_header_defaults(pscan curr_scan)
{
	/* Put default values for elements and attributes */	
	curr_scan->deprecated = 0;
	curr_scan->origin_offset = -1;
	curr_scan->origin_length = 0;
	curr_scan->origin_array = NULL;
	curr_scan->origin_count = 0;

	curr_scan->precursor_offset = -1;
	curr_scan->precursor_length = 0;
	curr_scan->precursor_array = NULL;
	curr_scan->precursor_count = 0;

	curr_scan->maldi_offset = -1;
	curr_scan->maldi_length = 0;
	curr_scan->maldi = NULL;

	curr_scan->peak_offset = -1;
	curr_scan->peak_content_offset = -1;	
	curr_scan->peaks = NULL;

	curr_scan->namevalue_offset = -1;
	curr_scan->namevalue_length = 0;
	curr_scan->namevalue_array = NULL;
	curr_scan->namevalue_array_len = 0;

}/* void set_scan_header_defaults(pscan curr_scan) */


/* Calculates the lengths of the segments */
void set_scan_lengths(pscan curr_scan, long endpos)
{
	long offset_array[6];
	long length_array[5];
	int i, j;

	offset_array[0] = curr_scan->origin_offset;
	offset_array[1] = curr_scan->precursor_offset;
	offset_array[2] = curr_scan->maldi_offset;
	offset_array[3] = curr_scan->peak_offset;
	offset_array[4] = curr_scan->namevalue_offset;
	offset_array[5] = endpos;

	for (i=0; i<5; i++) {
		if (offset_array[i] == -1)
			length_array[i] = 0;
		else {
			j = i+1;
			while (offset_array[j] == -1 && j < 6)
				j += 1;
			length_array[i] = (offset_array[j] - offset_array[i])-1;
		}/* else */
	}/* for */

	curr_scan->origin_length = length_array[0];
	curr_scan->precursor_length = length_array[1];
	curr_scan->maldi_length = length_array[2];
	curr_scan->namevalue_length = length_array[4];	

}/* void set_scan_lengths(pscan curr_scan, long endpos) */


/* Function that initialises scan by determining offsets of all the scan parts */
void parse_scan_header(pmzxml_file file, int scan_number, FILE* finput, scan_config_flags sflags)
{
	char *tag;	
	pscan curr_scan;
	long eopeaks_offset;

	offset = file->index_array[scan_number-1];
	fseek(finput, offset, SEEK_SET);
	read_chars = fread(read_buffer, sizeof(char), READ_BUFF_SIZE-1, finput);
	read_buffer[read_chars] = '\0';
	walkptr = read_buffer;

	if (!file->scan_array[scan_number-1])
		file->scan_array[scan_number-1] = malloc(sizeof(scan));
	curr_scan = file->scan_array[scan_number-1];

	/* Put default values for elements and attributes */	
	set_scan_header_defaults(curr_scan);

	/* Scan until done */
	tag = get_xml_tag(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
	while (strstr(tag, MZXML_PEAKS_OTAG) == NULL) {
		if (strstr(tag, MZXML_SCAN_OTAG)) {
			curr_scan->attributes.number = scan_number;
			parse_scan_attributes(curr_scan, tag);
		}/* if */
		else if (strstr(tag, MZXML_SCANORIGIN_OTAG)) {
			if (curr_scan->origin_offset < 0)
				curr_scan->origin_offset = offset + ((int)tag - (int)read_buffer);
			if (sflags & scan_origin_flag) 
				parse_scanorigin_structure(curr_scan, tag, finput);
		}/* else if */
		else if (strstr(tag, MZXML_PRECURSORMZ_OTAG)) {
			if (curr_scan->precursor_offset < 0)
				curr_scan->precursor_offset = offset + ((int)tag - (int)read_buffer);			
			if (sflags & scan_precursor_flag) 
				parse_precursor_structure(curr_scan, tag, finput);
		}/* else if */
		else if (strstr(tag, MZXML_MALDI_OTAG)) {
			if (curr_scan->maldi_offset < 0) 
				curr_scan->maldi_offset = offset + ((int)tag - (int)read_buffer);								
			if (sflags & scan_maldi_flag)
				parse_maldi_structure(curr_scan, tag, finput);
		}/* else if */		
		tag = get_xml_tag(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
	}/* while */

	/* Doing some bookkeeping for offsets and lengths */
	curr_scan->peak_offset = offset + ((int)tag - (int)read_buffer);		
	curr_scan->peak_content_offset = offset + ((int)walkptr - (int)read_buffer);	

	/* Parsing the scan header */
	parse_scan_peaks_header(curr_scan, tag);

	/* Move on to the part where the potential namevalues are */
	eopeaks_offset = curr_scan->peak_content_offset + calc_encoding_length(curr_scan->peaks) + strlen(MZXML_PEAKS_CTAG) + 1;	
	if (eopeaks_offset > offset + READ_BUFF_SIZE) {
		offset = eopeaks_offset;
		fseek(finput, offset, SEEK_SET);
		read_chars = fread(read_buffer, sizeof(char), READ_BUFF_SIZE-1, finput);
		read_buffer[read_chars] = '\0';
		walkptr = read_buffer;
	}/* if */
	else {
		walkptr = &(read_buffer[eopeaks_offset-offset]);
	}/* else */

	/* Read until we find the end */
	tag = get_xml_tag(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
	while (strstr(tag, MZXML_SCAN_CTAG) == NULL && strstr(tag, MZXML_SCAN_OTAG) == NULL) {
		if (strstr(tag, MZXML_NAMEVALUE_OTAG)) {
			if (curr_scan->namevalue_offset < 0) 
				curr_scan->namevalue_offset = offset + ((int)tag - (int)read_buffer);				
			parse_namevalue_structure(curr_scan, tag, finput);
		}/* if */	
		else if (strstr(tag, MZXML_COMMENT_OTAG)) {			
			curr_scan->namevalue_array[curr_scan->namevalue_array_len].comment = get_xml_tag_value(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
		}/* else if */	
		tag = get_xml_tag(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
	}/* while */

	/* Calculating the lengths of all segments */
	set_scan_lengths(curr_scan, offset + ((int)tag - (int)read_buffer));

	/* Doing some bookkeeping on the lengths of the arrarys created */
	if (curr_scan->origin_count > 0)
		curr_scan->origin_array = realloc(curr_scan->origin_array, curr_scan->origin_count * sizeof(scan_origin));
	if (curr_scan->precursor_count > 0)
		curr_scan->precursor_array = realloc(curr_scan->precursor_array, curr_scan->precursor_count * sizeof(scan_precursor));
	if (curr_scan->namevalue_array_len > 0)
		curr_scan->namevalue_array = realloc(curr_scan->namevalue_array, curr_scan->namevalue_array_len * sizeof(namevalue));

}/* void init_scan(pmzxml_file file, int scan_number) */


/* Function that parses scan attributes */
void parse_scan_attributes(pscan scan, char* content)
{
	char *tmpbuffer;	// compare/work buffer, and pointer for interm. processing	

	/* Setting msLevel */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_MSLEVEL);	
	if (!tmpbuffer)
		return;
	scan->attributes.msLvl = atoi(tmpbuffer);
	free(tmpbuffer);
	
	/* Setting peakcount */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_PEAKSCOUNT);
	if (!tmpbuffer)
		return;
	scan->attributes.peakscount = atoi(tmpbuffer);
	free(tmpbuffer);
	
	/* Setting polarity */
	scan->attributes.polarity = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_POLARITY);						
	
	/* Setting scan type */
	scan->attributes.scanType = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_SCANTYPE);			

	/* Setting filterline */
	scan->attributes.filterLine = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_FILTERLINE);	

	/* Setting centroided indication */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_CENTROIDED);	
	if (tmpbuffer) {			
		scan->attributes.centroided = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.centroided = 0;

	/* Setting isotope indication */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_DEISOTOPED);	
	if (tmpbuffer) {			
		scan->attributes.deisotoped = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.deisotoped = 0;

	/* Setting deconvolution indication */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_CHARGEDECONVOLUTED);	
	if (tmpbuffer) {			
		scan->attributes.chargeDeconvoluted = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.chargeDeconvoluted = 0;

	/* Setting retention time */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_RETENTIONTIME);	
	if (tmpbuffer) {			
		scan->attributes.retentionTime = xml_duration_to_seconds(tmpbuffer);
		free(tmpbuffer);	
	}/* if */
	else 
		scan->attributes.retentionTime = -1;

	/* Setting ionisation energy */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_IONISATIONENERGY);	
	if (tmpbuffer) {			
		scan->attributes.ionisationEnergy = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.ionisationEnergy = -1;

	/* Setting collision energy */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_COLLISIONENERGY);	
	if (tmpbuffer) {			
		scan->attributes.collisionEnergy = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else 
		scan->attributes.collisionEnergy = -1;

	/* Setting collision gas indicator */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_COLLISIONGAS);	
	if (tmpbuffer) {			
		scan->attributes.collisonGas = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.collisonGas = 0;

	/* Setting collision gas pressure */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_COLLISIONGASPRESSURE);	
	if (!tmpbuffer)
		tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_CIDGASPRESSURE);	
	if (tmpbuffer)	{
		scan->attributes.collisonGasPressure = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.collisonGasPressure = -1;

	/* Setting startMz */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_STARTMZ);	
	if (tmpbuffer) {			
		scan->attributes.startMz = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.startMz = -1;

	/* Setting endMz */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_ENDMZ);	
	if (tmpbuffer) {			
		scan->attributes.endMz = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.endMz = -1;

	/* Setting lowMz */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_LOWMZ);	
	if (tmpbuffer) {			
		scan->attributes.lowMz = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.lowMz = -1;

	/* Setting highMz */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_HIGHMZ);	
	if (tmpbuffer) {			
		scan->attributes.highMz = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.highMz = -1;

	/* Setting basePeak M/z */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_BASEPEAKMZ);	
	if (tmpbuffer) {			
		scan->attributes.basePeakMz = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.basePeakMz = -1;

	/* Setting basePeak Intensity */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_BASEPEAKINTENSITY);	
	if (tmpbuffer) {			
		scan->attributes.basePeakIntensity = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.basePeakIntensity = -1;
	
	/* Setting basePeak M/z */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_TOTIONCURRENT);	
	if (tmpbuffer) {			
		scan->attributes.totIonCurrent = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.totIonCurrent = -1;

	/* Setting scan ID */
	tmpbuffer = get_xml_attribute_value(content, MZXML_SCAN_ATTRIB_MSINSTRUMENTID);	
	if (tmpbuffer) {			
		scan->attributes.msInstrumentID = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->attributes.msInstrumentID = -1;

}/* void parse_scan_attributes(pscan scan, char* header_content) */


/* Parses the scan origin structure */
void parse_scanorigin_structure(pscan scan, char* beginptr, FILE* finput)
{
	char* tmpbuffer;

	if (!(scan->origin_array)) {
		scan->origin_array = malloc(sizeof(scan_origin));
		scan->origin_count = 0;	
	}/* if */
	
	/* Retrieving the scanorigin attributes */
	scan->origin_array[scan->origin_count].num = -1;	
	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_SCANORIGIN_ATTRIB_NUM);
	if (tmpbuffer) {
		scan->origin_array[scan->origin_count].num = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	scan->origin_array[scan->origin_count].parentFileID = get_xml_attribute_value(beginptr, MZXML_SCANORIGIN_ATTRIB_PARENTFILEID);	
				
	/* Updating counters and making more room for future structures */
	scan->origin_count += 1;
	if (scan->origin_count == origin_alloc_count) {
		origin_alloc_count *= 2;
		scan->origin_array = realloc(scan->origin_array, origin_alloc_count * sizeof(scan_origin));
	}/* if */		

}/* void parse_scanorigin_structure(pscan scan, char* beginptr, FILE* finput) */


/* Function that scans the precursor structure */
void parse_precursor_structure(pscan scan, char* beginptr, FILE* finput)
{	
	if (!(scan->precursor_array)) {
		scan->precursor_array = malloc(sizeof(scan_precursor));
		scan->precursor_count = 0;	
	}/* if */
	
	/* Retrieving the precursor attributes */
	scan->precursor_array[scan->precursor_count] = parse_precursor(beginptr, finput);
				
	/* Updating counters and making more room for future structures */
	scan->precursor_count += 1;
	if (scan->precursor_count == precursor_alloc_count) {
		precursor_alloc_count *= 2;
		scan->precursor_array = realloc(scan->precursor_array, precursor_alloc_count * sizeof(scan_precursor));
	}/* if */

}/* void parse_precursor_structure(pscan scan, char* content) */


/* Filters the info of a single precursor entry */
scan_precursor parse_precursor(char* beginptr, FILE* finput)
{
	char *tmpbuffer;
	scan_precursor prec;

	/* Setting precursor scan number */
	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_PRECURSORMZ_ATTRIB_PRECURSORSCANNUM);	
	prec.precursorScanNum = -1;
	if (tmpbuffer) {			
		prec.precursorScanNum = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */

	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_PRECURSORMZ_ATTRIB_PRECURSORINTENSITY);	
	prec.precursorIntensity = -1;
	if (tmpbuffer) {			
		prec.precursorIntensity = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */	
		
	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_PRECURSORMZ_ATTRIB_PRECURSORCHARGE);	
	prec.precursorCharge = -1;
	if (tmpbuffer) {			
		prec.precursorCharge = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */

	/* Filters possible charges */
	prec.possibleCharges = get_xml_attribute_value(beginptr, MZXML_PRECURSORMZ_ATTRIB_POSSIBLECHARGES);	

	/* Filters window wideness */
	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_PRECURSORMZ_ATTRIB_WINDOWWIDENESS);	
	prec.windowWideness = -1;
	if (tmpbuffer) {			
		prec.windowWideness = atof(tmpbuffer);
		free(tmpbuffer);
	}/* if */			

	/* Filters activation method */
	prec.activationMethods = get_xml_attribute_value(beginptr, MZXML_PRECURSORMZ_ATTRIB_ACTIVATIONMETHOD);	
	
	tmpbuffer = get_xml_tag_value(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
	prec.value = -1;
	if (tmpbuffer)
		prec.value = strtod(tmpbuffer, NULL);

	return prec;

}/* pprecursor parse_precursor(char* content) */


/* Filter the optional maldi part */
void parse_maldi_structure(pscan scan, char* contents, FILE* finput)
{

}/* void parse_maldi(pscan scan, char* content) */


/* Function that filters peak attributes and the peak intensities */
void parse_scan_peaks(pscan scan, char* beginptr, FILE* finput)
{
	float* flcastpointer;	
	int* intcastpointer;
	char* tmpbuffer;
	int	i, len, newlen;
	
	scan->peaks->mzs = malloc(scan->peaks->count * sizeof(double));
	scan->peaks->intensities = malloc(scan->peaks->count * sizeof(double));
	len = calc_encoding_length(scan->peaks);
	//printf("\nlen: %i", len); fflush(stdout);
	//PROBLEM PROBABLY HERE!
	//function call: parse_scan_peaks(file->scan_array[scan_nr-1], contents, 0);
	//beginptr -> contents
	tmpbuffer = decode_b64(beginptr, len+1, &newlen);
	//printf("\ntmpbuffer: %s", tmpbuffer); fflush(stdout);
	if (scan->peaks->precision == 32) {
		flcastpointer = (float*)tmpbuffer;
		intcastpointer = (int*)tmpbuffer;

		for (i=0; i<scan->peaks->count; i++) {
			//printf("\nPeak no: %i", i); fflush(stdout);
			intcastpointer[2*i] = ntohl((unsigned int)intcastpointer[2*i]);
			intcastpointer[(2*i)+1] = ntohl((unsigned int)intcastpointer[(2*i)+1]);
			scan->peaks->mzs[i] = (double)(flcastpointer[(2*i)]);
			scan->peaks->intensities[i] = (double)(flcastpointer[(2*i)+1]);
			//print the m/z-intensity pairs for each scan
			//printf("\nmz: %f", scan->peaks->mzs[i]); fflush(stdout);
			//printf("\nintensity: %f\n", scan->peaks->intensities[i]); fflush(stdout);

		}// for
		free(tmpbuffer);
	}// if

}/* void parse_peaks_header(pscan scan, char* content, long offset) */


/* Function that filters peak header attributes */
void parse_scan_peaks_header(pscan scan, char *beginptr)
{
	char* tmpbuffer;

	if (!scan->peaks)
		scan->peaks = malloc(sizeof(scan_peaks));

	/* Setting the count */
	scan->peaks->mzs = NULL;
	scan->peaks->intensities = NULL;
	scan->peaks->count = scan->attributes.peakscount;
	scan->peaks->compressionType = get_xml_attribute_value(beginptr, MZXML_PEAKS_ATTRIB_COMPRESSIONTYPE);
	scan->peaks->contentType = get_xml_attribute_value(beginptr, MZXML_PEAKS_ATTRIB_CONTENTTYPE);
	scan->peaks->byteOrder = get_xml_attribute_value(beginptr, MZXML_PEAKS_ATTRIB_BYTEORDER);
	scan->peaks->pairOrder = get_xml_attribute_value(beginptr, MZXML_PEAKS_ATTRIB_PAIRORDER);

	/* Setting peak precision */
	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_PEAKS_ATTRIB_PRECISION);	
	if (tmpbuffer) {			
		scan->peaks->precision = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->peaks->precision = -1;

	/* Setting compression length */
	tmpbuffer = get_xml_attribute_value(beginptr, MZXML_PEAKS_ATTRIB_COMPRESSEDLEN);	
	if (tmpbuffer) {			
		scan->peaks->compressedLen = atoi(tmpbuffer);
		free(tmpbuffer);
	}/* if */
	else
		scan->peaks->compressedLen = -1;

}/* void parse_scan_peaks_header(pscan scan, char *beginptr) */


/* Filter the optional maldi part */
void parse_namevalue_structure(pscan scan, char* beginptr, FILE* finput)
{
	if (!(scan->namevalue_array)) {
		scan->namevalue_array = malloc(sizeof(namevalue));
		scan->namevalue_array_len = 0;	
	}/* if */

	/* Retrieving the processingoperation attributes */
	scan->namevalue_array[scan->namevalue_array_len].name = get_xml_attribute_value(beginptr, MZXML_NAMEVALUE_ATTRIB_NAME);
	scan->namevalue_array[scan->namevalue_array_len].type = get_xml_attribute_value(beginptr, MZXML_NAMEVALUE_ATTRIB_TYPE);
	scan->namevalue_array[scan->namevalue_array_len].value = get_xml_attribute_value(beginptr, MZXML_NAMEVALUE_ATTRIB_VALUE);
				
	/* Updating counters and making more room for future structures */
	scan->namevalue_array_len += 1;
	if (scan->namevalue_array_len == namevalue_alloc_count) {
		namevalue_alloc_count *= 2;
		scan->namevalue_array = realloc(scan->namevalue_array, namevalue_alloc_count * sizeof(namevalue));
	}/* if */			

}/* void parse_namevalue_structure(pscan scan, char* begin, FILE* finput) */


/* Parses the end of the scan segment */
void parse_scan_end(pmzxml_file mzxml_file, FILE* finput)
{	
	char* tmpbuffer, *prevbuffer=NULL;

	/* Positioning */
	fseek(finput, scan_index_begin-SHORT_HEADER_BUFF_SIZE, SEEK_SET);
	read_chars = fread(read_buffer, sizeof(char), SHORT_HEADER_BUFF_SIZE-1, finput);
	read_buffer[read_chars] = '\0';
	
	/* Finding the last scanclosure tag */
	tmpbuffer = strstr(read_buffer, MZXML_SCAN_CTAG);
	while (tmpbuffer) {
		prevbuffer = tmpbuffer;
		tmpbuffer = strstr(tmpbuffer+1, MZXML_SCAN_CTAG);
	}/* while */

	mzxml_file->scan_end_offset =  scan_index_begin-SHORT_HEADER_BUFF_SIZE + ((int)prevbuffer - (int)read_buffer) + strlen(MZXML_SCAN_CTAG) + 1;

}/* void parse_scan_end(pmzxml_file mzxml_file, FILE* finput) */


/* Parses search score sequence */
void parse_index_sequence(pmzxml_file mzxml_file, FILE* finput)
{	 
	int alloc_count = 100;	
	char* tmpbuffer, *tag;

	/* Positioning */
	fseek(finput, scan_index_begin, SEEK_SET);
	read_chars = fread(read_buffer, sizeof(char), READ_BUFF_SIZE-1, finput);
	read_buffer[read_chars] = '\0';
	walkptr = read_buffer;

	mzxml_file->index_array = malloc(alloc_count * sizeof(long));
	mzxml_file->scan_num = 0;

	tag = get_xml_tag(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
	while (strstr(tag, MZXML_CTAG) == NULL) {
		//printf("\n%s", tag); fflush(stdout);
		/* finding index tag */
		if (strstr(tag, MZXML_INDEX_OTAG)) {
			mzxml_file->index_name =  get_xml_attribute_value(tag, MZXML_INDEX_ATTRIB_NAME);
			//printf("\n***%s", mzxml_file->index_name); fflush(stdout);
		}/* if */
		else if (strstr(tag, MZXML_OFFSET_OTAG)) {
			tmpbuffer = get_xml_tag_value(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);			
			mzxml_file->index_array[mzxml_file->scan_num] = atol(tmpbuffer);	
			//printf("\n%i", mzxml_file->index_array[mzxml_file->scan_num]); fflush(stdout);
			mzxml_file->scan_num += 1;
			free(tmpbuffer);
			if (mzxml_file->scan_num == alloc_count) {
				alloc_count *= 2;
				mzxml_file->index_array = realloc(mzxml_file->index_array, alloc_count * sizeof(long));
			}/* if */
		}/* else if */						
		tag = get_xml_tag(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
	}/* while */
	//printf("\n%i", mzxml_file->scan_num); fflush(stdout);
	mzxml_file->index_array = realloc(mzxml_file->index_array, mzxml_file->scan_num * sizeof(long));
	
}/* void parse_index_sequence(pmzxml_file mzxml_file, FILE* finput) */


/* Handles file tail of the mzxml file */
void parse_file_tail(pmzxml_file mzxml_file, FILE* finput)
{	
	char *tmpbuffer, *tag;

	/* Set defaults of optional items */
	mzxml_file->sha1 = NULL;
	mzxml_file->scan_end_offset = 0;

	/* Positioning */
	fseek(finput, -READ_TAIL_BUFF_SIZE+1, SEEK_END);
	read_chars = fread(read_buffer, sizeof(char), READ_TAIL_BUFF_SIZE-1, finput);
	read_buffer[read_chars] = '\0';
	walkptr = read_buffer;

	tag = get_xml_tag(read_buffer, &walkptr, finput, READ_TAIL_BUFF_SIZE, &offset);
	while (strstr(tag, MZXML_CTAG) == NULL) {		
		/* Retrieving SHA1 content (if any) */
		if (strstr(tag, MZXML_SHA1_OTAG))
			mzxml_file->sha1 = get_xml_tag_value(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
		/* Retrieving index structure offset */
		else if (strstr(tag, MZXML_INDEXOFFSET_OTAG)) {
			tmpbuffer = get_xml_tag_value(read_buffer, &walkptr, finput, READ_BUFF_SIZE, &offset);
			if (tmpbuffer) {
				scan_index_begin = atol(tmpbuffer);			
				free(tmpbuffer);
			}/* if */
			else 
				return;	
		}/* else */		
		tag = get_xml_tag(read_buffer, &walkptr, finput, READ_TAIL_BUFF_SIZE, &offset);
	}/* while */
			
}/* void parse_file_tail(pmzxml_file mzxml_file, FILE* finput) */
