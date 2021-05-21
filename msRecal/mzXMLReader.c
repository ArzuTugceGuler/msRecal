#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "StringFunctions.h"
#include "mzXMLDefs.h"
#include "mzXMLReader.h"
#include "mzXMLParser.h"
#include "XMLFunctions.h"


/* Calculates the length of the encoded string */
static int calc_encoding_length(pscan_peaks sp)
{
	int len;
	if (sp->compressedLen >= 0)
		return sp->compressedLen;
	else {
		len = ((sp->count) * 2) * (sp->precision/8);
		
		return (((2 + len - ((len + 2) % 3)) * 4) / 3);
	}/* else */

}/* int calc_encoding_length(pscan_peaks sp) */


/* Only parses mzxml scans */
pmzxml_file read_mzxml_file(char* file_name, file_flags fflags, scan_config_flags sflags)
{
	return read_mzxml_file_spectrum(file_name, fflags, sflags, 1, INT_MAX);

}/* pmzxml_file read_mzxml_file(char* file_name, file_flags fflags, scan_config_flags sflags) */


/* Parses subset of scans */
pmzxml_file read_mzxml_file_spectrum(char* file_name, file_flags fflags, scan_config_flags sflags, int begin_scan, int end_scan)
{
	FILE* input;

	pmzxml_file retval = (pmzxml_file) malloc(sizeof(mzxml_file));
	retval->file_name = file_name;
	input = fopen(file_name, "rb");

	parse_mzxml_file(retval, input, fflags, sflags, begin_scan, end_scan);
	
	fclose(input);
	return retval;

}/* pmzxml_file read_mzxml_file_spectrum(char* file_name, file_flags fflags, scan_config_flags sflags, int begin_scan, int end_scan) */


void unload_namevalue_struct(namevalue nmv)
{
	free(nmv.name);
	free(nmv.value);
	free(nmv.type);
	free(nmv.comment);

}// void unload_namevalue_struct(namevalue nmv)

/* Function that unloads the peaks in a scan */
void unload_scanpeaks_struct(pscan_peaks pks)
{
	if (!pks)
		return;

	free (pks->mzs);
	free (pks->intensities);
	free (pks->byteOrder);
	free (pks->pairOrder);
	free (pks->contentType);
	free (pks->compressionType);

}// void unload_scan_peaks(pscan_peaks pks)

/* Unloads the scan maldi structure */
void unload_scanmaldi_struct(pscan_maldi mal)
{
	if (mal)
		free(mal->contents);
	
}// void unload_scanorigin_struct(scan_origin org)


/* Unloads the scan origin structure */
void unload_scanprecursor_struct(scan_precursor prec)
{
	free(prec.possibleCharges);
	free(prec.activationMethods);

}// void unload_scanorigin_struct(scan_origin org)


/* Unloads the scan origin structure */
void unload_scanorigin_struct(scan_origin org)
{
	free(org.parentFileID);

}// void unload_scanorigin_struct(scan_origin org)

/* Unloads the scan structure */
void unload_scan_struct(pscan scan_ptr)
{
	int i;

	if (!scan_ptr)
		return;

	// Attributes first
	free(scan_ptr->attributes.polarity);
	free(scan_ptr->attributes.scanType);
	free(scan_ptr->attributes.filterLine);

		// Scanorigin part
	for (i=0; i<scan_ptr->origin_count; i++)
		unload_scanorigin_struct((scan_ptr->origin_array[i]));
	free(scan_ptr->origin_array);

	// Scanorigin part
	for (i=0; i<scan_ptr->precursor_count; i++)
		unload_scanprecursor_struct((scan_ptr->precursor_array[i]));
	free(scan_ptr->precursor_array);

	// Maldi part	
	unload_scanmaldi_struct(scan_ptr->maldi);
	free(scan_ptr->maldi);

	// Scan peaks part
	unload_scanpeaks_struct(scan_ptr->peaks);
	free(scan_ptr->peaks);

	// Namevalue part
	for (i=0; i<scan_ptr->namevalue_array_len; i++)
		unload_namevalue_struct((scan_ptr->namevalue_array[i]));
	free(scan_ptr->namevalue_array);

}// void unload_scan_struct(pscan scan_ptr)

/* Unloads the spotting structure */
void unload_spotting_struct(pspotting spo)
{
	if (spo) 
		free(spo->content);

}// void unload_spotting_struct(spotting spo)

/* Unloads the separation structure */
void unload_separation_struct(pseparation sep)
{
	if (sep)
		free(sep->content);

}// void unload_separation_struct(separation sep)

/* Unloads the dataprocessing structure */
void unload_dataprocessing_struct(dataprocessing dp)
{
	int i;

	free(dp.sw_name);
	free(dp.sw_type);
	free(dp.sw_version);

	for (i=0; i<dp.proc_op_count; i++) {
		free(dp.proc_op_array[i].name);
		free(dp.proc_op_array[i].value);
		free(dp.proc_op_array[i].type);
		free(dp.proc_op_array[i].comment);
	}// for

}// void unload_dataprocessing_struct(dataprocessing dp)

/* Unloads the msinstrument structure */
void unload_msinstrument_struct(msinstrument msi)
{
	free(msi.content);

}// void unload_msinstrument_struct(msinstrument msi)


/* Unloads the parentfile structure */
void unload_parentfile_struct(parentfile pf)
{
	free(pf.file_name);
	free(pf.file_type);
	free(pf.sha1);

}// void unload_parentfile_struct(pparentfile pf)


/* Unloads an mzXML file from memory */
void unload_mzxml_file(pmzxml_file mzXML_file)
{
	int i;

	for (i=0; i<mzXML_file->parentfile_count; i++)
		unload_parentfile_struct((mzXML_file->parentfile_array[i]));
	free(mzXML_file->parentfile_array);

	for (i=0; i<mzXML_file->msinstrument_count; i++)
		unload_msinstrument_struct((mzXML_file->msinstrument_array[i]));
	free(mzXML_file->msinstrument_array);

	for (i=0; i<mzXML_file->dataprocessing_count; i++)
		unload_dataprocessing_struct((mzXML_file->dataprocessing_array[i]));
	free(mzXML_file->dataprocessing_array);

	unload_separation_struct((mzXML_file->separation));
	free(mzXML_file->separation);

	unload_spotting_struct((mzXML_file->spotting));
	free(mzXML_file->spotting);

	for (i=0; i<mzXML_file->scan_num; i++) {
		unload_scan_struct(mzXML_file->scan_array[i]);
		free(mzXML_file->scan_array[i]);
	}// for
	
	free(mzXML_file->index_array);
	free(mzXML_file->index_name);
	free(mzXML_file->sha1);
	free(mzXML_file);

}// void unload_mzxml_file(pmzxml_file mzXML_file)


/* Gets the parentfile structure */
pparentfile get_parentfile_structure(pmzxml_file file, int* parentfile_count)
{
	FILE* finput;
	int read;
	long offset;
	char* readbuffer, *tag, *walkptr;

	/* In case it's not there */
	if (!file->parentfile_array && file->parentfile_offset >= 0) {
		finput = fopen(file->file_name, "rb");
		fseek(finput, file->parentfile_offset, SEEK_SET);
		readbuffer = malloc(file->parentfile_length + 1);
		read = fread(readbuffer, sizeof(char), file->parentfile_length + 1, finput);
		readbuffer[read] = '\0';
		walkptr = readbuffer;
		fclose(finput);

		tag = get_xml_tag(readbuffer, &walkptr, NULL, file->parentfile_length+1, &offset);
		while (tag) {
			if (strstr(tag, MZXML_PARENTFILE_OTAG)) {
				parse_parentfile_structure(file, tag, NULL);			
			}/* if */
			tag = get_xml_tag(readbuffer, &walkptr, NULL, file->parentfile_length+1, &offset);
		}/* while */	
		file->parentfile_array = realloc(file->parentfile_array, file->parentfile_count * sizeof(parentfile));
	}/* if */

	*parentfile_count = file->parentfile_count;	
	return file->parentfile_array;

}/* pparentfile get_parentfile_structure(pmzxml_file file, int* parentfile_count) */


/* Gets the msinstrument structure */
pmsinstrument get_msinstrument_structure(pmzxml_file file, int* instrument_count)
{
	FILE* finput;
	int read;
	long offset;
	char* readbuffer, *tag, *walkptr;

	/* In case it's not there */
	if (!file->msinstrument_array && file->msinstrument_offset >= 0) {
		finput = fopen(file->file_name, "rb");
		fseek(finput, file->msinstrument_offset, SEEK_SET);
		readbuffer = malloc((file->dataprocessing_offset - file->msinstrument_offset) + 1);
		read = fread(readbuffer, sizeof(char), (file->dataprocessing_offset - file->msinstrument_offset) + 1, finput);
		readbuffer[read] = '\0';
		walkptr = readbuffer;
		fclose(finput);

		tag = get_xml_tag(readbuffer, &walkptr, NULL, file->dataprocessing_offset - file->msinstrument_offset, &offset);
		while (tag) {
			if (strstr(tag, MZXML_MSINSTRUMENT_OTAG)) {
				parse_msinstrument_structure(file, tag, NULL);	
			}/* if */
			else
				break;			
			tag = get_xml_tag(readbuffer, &walkptr, NULL, file->dataprocessing_offset - file->msinstrument_offset, &offset);
		}/* while */	
		file->msinstrument_array = realloc(file->msinstrument_array, file->msinstrument_count * sizeof(msinstrument));
	}/* if */

	*instrument_count = file->msinstrument_count;	
	return file->msinstrument_array;

}/* pparentfile get_parentfile_structure(pmzxml_file file, int* parentfile_count) */


/* Gets the dataprocessing structure */
pdataprocessing get_dataprocessing_structure(pmzxml_file file, int* dataprocessing_count)
{
	FILE* finput;
	int read;
	long offset;
	char* readbuffer, *tag, *walkptr;

	/* In case it's not there */
	if (!file->dataprocessing_array && file->dataprocessing_offset >= 0) {
		finput = fopen(file->file_name, "rb");
		fseek(finput, file->dataprocessing_offset, SEEK_SET);
		readbuffer = malloc((file->index_array[0] - file->dataprocessing_offset) + 1);
		read = fread(readbuffer, sizeof(char), (file->index_array[0] - file->dataprocessing_offset) + 1, finput);
		readbuffer[read] = '\0';
		walkptr = readbuffer;
		fclose(finput);

		tag = get_xml_tag(readbuffer, &walkptr, NULL, file->index_array[0] - file->dataprocessing_offset, &offset);
		while (tag) {
			if (strstr(tag, MZXML_DATAPROCESSING_OTAG)) {
				parse_dataprocessing_structure(file, tag, NULL);	
			}/* if */
			else
				break;			
			tag = get_xml_tag(readbuffer, &walkptr, NULL, file->index_array[0] - file->dataprocessing_offset, &offset);
		}/* while */	
		file->dataprocessing_array = realloc(file->dataprocessing_array, file->dataprocessing_count * sizeof(dataprocessing));
	}/* if */

	*dataprocessing_count = file->dataprocessing_count;	
	return file->dataprocessing_array;

}/* pdataprocessing get_dataprocessing_structure(pmzxml_file file,  dataprocessing_array, int* dataprocessing_count) */


/* Gets the separation structure */
pseparation get_separation_structure(pmzxml_file file)
{
	return NULL;

}/* pseparation get_separation_structure(pmzxml_file file) */


/* Gets the spotting structure */
pspotting get_spotting_structure(pmzxml_file file)
{
	return NULL;

}/* pspotting get_spotting_structure(pmzxml_file file) */


/* Function that retrieves the whole scan */
pscan get_scan(pmzxml_file file, int scan_nr, scan_config_flags sflags)
{
	FILE* input;
	int count;	
	
	if (!file->scan_array[scan_nr-1]) {		
		file->scan_array[scan_nr-1] = malloc(sizeof(scan));		
		input = fopen(file->file_name, "rb");
		parse_scan_header(file, scan_nr, input, sflags);	
		fclose(input);
	}/* if */
	else {
		if (sflags & scan_origin_flag) 
			get_scan_origin(file, scan_nr, &count);
		if (sflags & scan_precursor_flag) 
			get_scan_precursor(file, scan_nr, &count);
		if (sflags & scan_maldi_flag)
			get_scan_maldi(file, scan_nr);
	}/* else */

	return file->scan_array[scan_nr-1];

}/* pscan get_scan(pmzxml_file file, int scan_nr) */


/* Function that retrieves a specific scan's attributes */
scan_attributes get_scan_attributes(pmzxml_file file, int scan_nr)
{
	FILE* input;

	if (!file->scan_array[scan_nr-1]) {
		file->scan_array[scan_nr-1] = malloc(sizeof(scan));		
		input = fopen(file->file_name, "rb");
		parse_scan_header(file, scan_nr, input, 0);	
		fclose(input);
	}/* if */		

	return file->scan_array[scan_nr-1]->attributes;

}/* pscan_attributes get_scan_attributes(pmzxml_file file, int scan_nr) */


/* Function that retrieves a specific scan's origin structure */
pscan_origin get_scan_origin(pmzxml_file file, int scan_nr, int* origin_count)
{
	FILE* input;
	long offset;
	char* readbuff, *tag, *walkptr;
	int read;

	if (!file->scan_array[scan_nr-1]) {
		file->scan_array[scan_nr-1] = malloc(sizeof(scan));		
		input = fopen(file->file_name, "rb");
		parse_scan_header(file, scan_nr, input, scan_origin_flag);	
		fclose(input);
	}/* if */	
	else if (!file->scan_array[scan_nr-1]->origin_array && file->scan_array[scan_nr-1]->origin_offset >= 0) {
		input = fopen(file->file_name, "rb");
		readbuff = malloc(sizeof(char) * ((file->scan_array[scan_nr-1]->peak_offset - file->scan_array[scan_nr-1]->origin_offset)+1));
		fseek(input, file->scan_array[scan_nr-1]->origin_offset, SEEK_SET);
		read = fread(readbuff, sizeof(char), (file->scan_array[scan_nr-1]->peak_offset - file->scan_array[scan_nr-1]->origin_offset)+1, input);
		readbuff[read] = '\0';
		walkptr = readbuff;
		fclose(input);

		tag = get_xml_tag(readbuff, &walkptr, NULL, file->scan_array[scan_nr-1]->peak_offset - file->scan_array[scan_nr-1]->origin_offset, &offset);
		while (tag != NULL) {
			if (strstr(tag, MZXML_SCANORIGIN_OTAG)) {
				parse_scanorigin_structure(file->scan_array[scan_nr-1], tag, NULL);	
			}/* if */
			else
				break;			
			tag = get_xml_tag(readbuff, &walkptr, NULL, file->scan_array[scan_nr-1]->peak_offset - file->scan_array[scan_nr-1]->origin_offset, &offset);
		}/* while */	
		file->scan_array[scan_nr-1]->origin_array = realloc(file->scan_array[scan_nr-1]->origin_array, file->scan_array[scan_nr-1]->origin_count * sizeof(scan_origin));		
		free(readbuff);
	}/* else if */

	*origin_count = file->scan_array[scan_nr-1]->origin_count;
	return file->scan_array[scan_nr-1]->origin_array;

}/* pscan_origin get_scan_origin(pmzxml_file file, int scan_nr) */


/* Function that retrieves a specific scan's precursor array */
pscan_precursor get_scan_precursor(pmzxml_file file, int scan_nr, int* precursor_count)
{
	FILE* input;
	long offset;
	char *readbuff, *tag, *walkptr;
	int read;

	if (!file->scan_array[scan_nr-1]) {
		file->scan_array[scan_nr-1] = malloc(sizeof(scan));		
		input = fopen(file->file_name, "rb");
		parse_scan_header(file, scan_nr, input, scan_precursor_flag);	
		fclose(input);
	}/* if */	
	else if (!file->scan_array[scan_nr-1]->precursor_array && file->scan_array[scan_nr-1]->precursor_offset) {
		input = fopen(file->file_name, "rb");
		readbuff = malloc(sizeof(char) * ((file->scan_array[scan_nr-1]->peak_offset - file->scan_array[scan_nr-1]->precursor_offset) + 1));
		fseek(input, file->scan_array[scan_nr-1]->precursor_offset, SEEK_SET);
		read = fread(readbuff, sizeof(char), (file->scan_array[scan_nr-1]->peak_offset - file->scan_array[scan_nr-1]->precursor_offset), input);
		readbuff[read] = '\0';
		fclose(input);

		tag = get_xml_tag(readbuff, &walkptr, NULL, file->scan_array[scan_nr-1]->peak_offset - file->scan_array[scan_nr-1]->precursor_offset, &offset);
		while (tag != NULL) {
			if (strstr(tag, MZXML_PRECURSORMZ_OTAG)) {
				parse_precursor_structure(file->scan_array[scan_nr-1], tag, NULL);	
			}/* if */
			else
				break;			
			tag = get_xml_tag(readbuff, &walkptr, NULL, file->scan_array[scan_nr-1]->peak_offset - file->scan_array[scan_nr-1]->precursor_offset, &offset);
		}/* while */	
		file->scan_array[scan_nr-1]->precursor_array = realloc(file->scan_array[scan_nr-1]->precursor_array, file->scan_array[scan_nr-1]->precursor_count * sizeof(scan_precursor));				
		free(readbuff);		
	}/* else if */

	*precursor_count = file->scan_array[scan_nr-1]->precursor_count;
	return file->scan_array[scan_nr-1]->precursor_array;

}/* pscan_precursor get_scan_precursor(pmzxml_file file, int scan_nr, int* precursor_count) */


/* Function that retrieves a specific scan's maldi structure TO BE IMPLEMENTED */
pscan_maldi get_scan_maldi(pmzxml_file file, int scan_nr)
{
	return NULL;

}/* pscan_maldi get_scan_maldi(pmzxml_file file, int scan_nr) */

/* Function that retrieves peaks from an mzxml */
scan_peaks load_scan_peaks(pmzxml_file file, int scan_nr)
{
	FILE* input;
	long begin, end;
	char* contents;
	int read;

	input = fopen(file->file_name, "rb");
	if (!file->scan_array[scan_nr-1]) {
		file->scan_array[scan_nr-1] = malloc(sizeof(scan));
		parse_scan_header(file, scan_nr, input, 0);
	}/* if */	

	if (file->scan_array[scan_nr-1]->peaks->mzs)
		free(file->scan_array[scan_nr-1]->peaks->mzs);
	if (file->scan_array[scan_nr-1]->peaks->intensities)
		free(file->scan_array[scan_nr-1]->peaks->intensities);

	begin = file->scan_array[scan_nr-1]->peak_content_offset;	
	end = begin + calc_encoding_length(file->scan_array[scan_nr-1]->peaks);

	fseek(input, begin, SEEK_SET);
	contents = malloc(((end-begin)+1)*sizeof(char));
	read = fread(contents, sizeof(char), (end-begin)*sizeof(char), input);
	contents[read] = '\0';
	fclose(input);
	parse_scan_peaks(file->scan_array[scan_nr-1], contents, 0);
	free(contents);
				
	return (*(file->scan_array[scan_nr-1]->peaks));
}


/* Function that unloads the peaks in a scan */
void unload_scan_peaks(pmzxml_file file, int scan_nr)
{
	if (file->scan_array[scan_nr-1]) {
		free(file->scan_array[scan_nr-1]->peaks->intensities);
		file->scan_array[scan_nr-1]->peaks->intensities = NULL;
		free(file->scan_array[scan_nr-1]->peaks->mzs);
		file->scan_array[scan_nr-1]->peaks->mzs = NULL;
	}// if

}// void unload_scan_peaks(pmzxml_file file, int scan_nr)
