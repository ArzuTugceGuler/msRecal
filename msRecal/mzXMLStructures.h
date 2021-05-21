/* 
 * File:   mzXMLStructures.h
 * Author: Arzu Tugce Guler
 *
 * Created on 15 April 2016, 16:33
 */

#ifndef MZXMLSTRUCTURES_H
#define	MZXMLSTRUCTURES_H

/* Processingoperation name value type */
typedef struct mzxml_namevaluetype_type
{
	char* name;
	char* value;
	char* type;
	char* comment;

} processingoperation, *pprocessingoperation, namevalue, *pnamevalue;


/* Precursor structure, records scan precursors */
typedef struct mzxml_precursor_type
{
	int precursorScanNum;		/* Scan number */
	double precursorIntensity;	/* Ion intensity */
	int precursorCharge;		/* Ion charge state */
	char* possibleCharges;		/* Comma-seperated values of possible estimated precursor charges */
	double windowWideness;		/* Selection window wideness */
	double value;				/* Value of the element */
	char* activationMethods;	/* Activation or fragmentation method: ETD, ECD, CID, etc */

} scan_precursor,*pscan_precursor;


/* Peak structure, used to return peaks that were decoded */
typedef struct mzxml_peak_type
{
	double* mzs;
	double* intensities;
	int     count;
	int		precision;
	char*	byteOrder;
	char*	pairOrder;
	char*	contentType;
	char*	compressionType;
	int		compressedLen;
	
} scan_peaks, *pscan_peaks;


/* Scan origin structure, used to represent scan origin */
typedef struct mzxml_scan_maldi_type
{
	char* contents;

} scan_maldi, *pscan_maldi;


/* Scan origin structure, used to represent scan origin */
typedef struct mzxml_scan_origin_type
{
	char* parentFileID;
	int num;

} scan_origin, *pscan_origin;


/* Attribute structure, used to record attributes of a scan */
typedef struct mzxml_scan_attribute_type
{
	/* Attributes */		
	int number;					/* Scan number */
	int msLvl;					/* MS level: 1= MS, 2= MS/MS, etc */
	int peakscount;				/* number of peaks this scan contains */
	char* polarity;				/* Polarity of the scan */
	char* scanType;				/* Scantype of the current scan */
	char* filterLine;			/* For Thermo instruments, a space to record the "filter line" string for the scan. */
	int centroided;				/* Centroid indication */
	int deisotoped;				/* Deisotoped indication */
	int chargeDeconvoluted;		/* Deconvolution indication */
	double retentionTime;		/* Retention time in seconds */	
	double ionisationEnergy;	/* Ionisation energy for current scan */
	double collisionEnergy;		/* Collision energy for current scan */
	int collisonGas;			/* Collision gas indicator */
	double collisonGasPressure;	/* Collision gas pressure */
	double startMz;				/* Lowest mz reading bound */
	double endMz;				/* Highest mz readind bound */
	double lowMz;				/* Lowest mz actually measured */
	double highMz;				/* Highest mz actually measured */
	double basePeakMz;			/* M/z value for highest intensity */
	double basePeakIntensity;	/* Value of the base peak intensity */
	double totIonCurrent;		/* Total ion current of the scan */
	int msInstrumentID;			/* ID of instrument used to procure the scan data */

} scan_attributes, *pscan_attributes;


/* Scan structure, records individual scans */
typedef struct mzxml_scan_type
{
	short deprecated;			/* indication that scan should be deleted */
	scan_attributes attributes;

	/* scanOrigin part */
	long origin_offset;
	long origin_length;
	pscan_origin origin_array;	
	int	origin_count;					

	/* precursorMz part */	
	long precursor_offset;
	long precursor_length;
	pscan_precursor	precursor_array;	/* array of all precursors */
	int	precursor_count;				/* size of precursor array */

	/* maldi part */
	long maldi_offset;
	long maldi_length;
	pscan_maldi maldi;			

	/* peaks part */
	long peak_offset;					/* offset of the peak */
	long peak_content_offset;			/* offset of the encoded peak content segment*/	
	pscan_peaks	peaks;					/* peaks structures */	
	
	/* namevalue part */
	long namevalue_offset;				/* offset of extra comments */
	long namevalue_length;				/* length in the original file */
	pnamevalue namevalue_array;			/* array of commments */
	int namevalue_array_len;			/* length of array of comments */ 
	
} scan, *pscan;


/* Spotting structure */
typedef struct mzxml_spotting_type
{
	char* content;

} spotting, *pspotting;


/* Separation structure */
typedef struct mzxml_separation_type
{
	char* content;

} separation, *pseparation;


/* Dataprocessing structure, records data processing steps in chronological order */
typedef struct mzxml_dataprocessing_type
{
	/* Attributes */
	double		  intensity_cutoff;	/* minimal intensity value for m/z - intensity pair */
	unsigned char centroided;		/* intean yes - no attribute */
	unsigned char deisotoped;		/* intean yes - no attribute */
	unsigned char deconvoluted;		/* intean yes - no attribute */
	unsigned char spotintegration;	/* intean yes - no attribute specific to LC-MALDI */

	/* Structures */
	char*  sw_type;				/* type of the software operation, either "conversion" or "processing" */
	char*  sw_name;				/* name of the software tool used */
	char*  sw_version;			/* version of the software package used */
	double sw_completion_time;	/* xs:duration completion time entry. optional */

	/* Ignored for now: */
	pprocessingoperation	proc_op_array;	/* array of processing operation structures */
	int						proc_op_count;	/* length of proc_op_array */

} dataprocessing, *pdataprocessing;


/* msinstrument section */
typedef struct mzxml_msinstrument_type
{
	char* content;

} msinstrument, *pmsinstrument;


/* Parentfile structure, records source of data in chronological order */
typedef struct mzxml_parentfile_type
{
	char* file_name;	/* file name of the source file */
	char* file_type;	/* file type of the source file */
	char* sha1;		/* sha1 crc of the source file */

} parentfile, *pparentfile;


/* mzxml file structure */
typedef struct mzxml_file
{
	char*			file_name;				/* source file name/location of the mzxml file */

	/*************************** MS RUN BEGIN **************************/
	double			start_time;
	double			end_time;
	long			msrun_offset;

	/* Structures */	
	long			parentfile_offset;		/* Offset of the section */
	long			parentfile_length;		/* Length of the section */
	pparentfile		parentfile_array;		/* Array of the parentfile structrure */
	int				parentfile_count;		/* Parentfile array length */	

	long			msinstrument_offset;	/* Offset of potential msinstrument section */
	long			msinstrument_length;	/* Length of the section */
	pmsinstrument	msinstrument_array;		/* pointer to the msinstument section, if loaded */
	int				msinstrument_count;		/* Msinstrument array length */

	long			dataprocessing_offset;	/* Offset of the section */
	long			dataprocessing_length;	/* Length of the section */	
	pdataprocessing	dataprocessing_array;	/* Array of the dataprocessing structrure */
	int				dataprocessing_count;	/* Dataprocessing array length */

	long			separation_offset;		/* Offset of the separation section */
	long			separation_length;		/* Length of the separation section */
	pseparation		separation;				/* pointer to the separation section, if loaded */

	long			spotting_offset;		/* Offset of the spotting section */
	long			spotting_length;		/* Length of the spotting section */
	pspotting		spotting;				/* pointer to the spotting section, if loaded */
	
	pscan*			scan_array;				/* Array of the scan structrures, size recorded in scan_count */
	/*************************** MS RUN END **************************/	
	
	long			scan_end_offset;		/* offset at which the scan section ends */
	long*			index_array;			/* array of indices. Arraylength is equal to the number of scans */
	int				scan_num;				/* number of scans */	
	char*			index_name;				/* name of the index */
	char*			sha1;					/* sha1 part */

} mzxml_file, *pmzxml_file;

#endif	/* MZXMLSTRUCTURES_H */

