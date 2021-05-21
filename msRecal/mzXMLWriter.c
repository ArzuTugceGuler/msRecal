#include <float.h>

#include "mzXMLWriter.h"
#include "mzXMLWriteFunctions.h"


/* Writes an mzXML file from a loaded / modified file */
void write_mzxml_file(pmzxml_file mzxml_file_hndl, char* output_file)
{
	print_mzxml(mzxml_file_hndl, output_file);

}/* void write_mzxml_file(pmzxml_file mzxml_file_hndl, char* output_file) */


/* Function that sets scan data, and updates all scan attributes relevant to the data */
void update_scan_peaks(pmzxml_file mzxml_file_hndl, int scan_num, int peak_num, int precision, double* mzs, double* intensities)
{
	int i;
	double low_mz = DBL_MAX, high_mz = -DBL_MAX, base_mz = -DBL_MAX, base_int = -DBL_MAX;
	pscan peaks = mzxml_file_hndl->scan_array[scan_num-1];

	/* Checking scan validity */
	if (!peaks) {
		printf("Error, scan %i does not exist or is not loaded.\n", scan_num);
		return;
	}/* if */
	
	peaks->attributes.peakscount = peak_num;
	peaks->peaks->count = peak_num;
	peaks->peaks->precision = precision;

	/* Looking for the delimiters */
	for (i=0; i<peak_num; i++) {
		/* Scanning for lowest mz */
		if (mzs[i] < low_mz)
			low_mz = mzs[i];
		/* Scanning for highest mz */
		if (mzs[i] > high_mz)
			high_mz = mzs[i];
		/* Scanning for base peak */
		if (intensities[i] > base_int) {
			base_int = intensities[i];
			base_mz = mzs[i];
		}/* if */
	}/* for */

	if (mzs != peaks->peaks->mzs)
		peaks->peaks->mzs = mzs;
	if (intensities != peaks->peaks->intensities)
		peaks->peaks->intensities = intensities;

	peaks->attributes.highMz = high_mz;
	peaks->attributes.lowMz = low_mz;
	peaks->attributes.basePeakMz = base_mz;
	peaks->attributes.basePeakIntensity = base_int;	

}/* void update_scan_peaks(pmzxml_file mzxml_file_hndl, int scan_num, int peak_num, int precision, double* mzs, double* intensities) */


/* Function that empties a scan, resetting all data except scan number, retentiontime, polarity, ms level and peakscont set to 0 */
void empty_scan(pmzxml_file mzxml_file_hndl, int scan_num)
{
	pscan mzscan = mzxml_file_hndl->scan_array[scan_num-1];
	if (!mzscan) {
		printf("Error, scan %i does not exist or is not loaded.\n", scan_num);
		return;
	}/* if */

	/* Modding attributes */		
	mzscan->attributes.peakscount = 0;				
	mzscan->peaks->count = 0;				
	mzscan->attributes.centroided = 0;	
	mzscan->attributes.deisotoped = 0;
	mzscan->attributes.chargeDeconvoluted = 0;				
	mzscan->attributes.ionisationEnergy = -1;
	mzscan->attributes.collisionEnergy = -1;	
	mzscan->attributes.collisonGasPressure = -1;
	mzscan->attributes.startMz = -1;
	mzscan->attributes.endMz = -1;
	mzscan->attributes.lowMz = -1;
	mzscan->attributes.highMz = -1;	
	mzscan->attributes.basePeakMz = -1;
	mzscan->attributes.basePeakIntensity = -1;
	mzscan->attributes.totIonCurrent = -1;
	mzscan->attributes.msInstrumentID = -1;
	
}// void empty_scan(pmzxml_file mzxml_file_hndl, int scan_num)


