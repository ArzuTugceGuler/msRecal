/* 
 * File:   mzXMLDefs.h
 * Author: Arzu Tugce Guler
 *
 * Created on 15 April 2016, 16:38
 */

#ifndef MZXMLDEFS_H
#define	MZXMLDEFS_H

/* Global tags and chars */
#define MZXML_ATTRIB_OTAG "\""
#define MZXML_ATTRIB_CTAG "\""
#define MZXML_ELEMENT_OTAG "<"
#define MZXML_ELEMENT_ETAG ">"
#define MZXML_ELEMENT_CTAG "/>"
#define MZXML_ATTRIB_NIL "xsi:nil=\"1\""


/* mzXML section */
#define MZXML_OTAG "<mzXML"
#define MZXML_CTAG "</mzXML"

/* msRun section */
#define MZXML_MSRUN_OTAG "<msRun"
#define MZXML_MSRUN_CTAG "</msRun"
#define MZXML_MSRUN_ATTRIB_SCANCOUNT "scanCount"
#define MZXML_MSRUN_ATTRIB_STARTTIME "startTime"
#define MZXML_MSRUN_ATTRIB_ENDTIME "endTime"

/* parentFile section */
#define MZXML_PARENTFILE_OTAG "<parentFile"
#define MZXML_PARENTFILE_CTAG "/"
#define MZXML_PARENTFILE_ATTRIB_FILENAME "fileName"
#define MZXML_PARENTFILE_ATTRIB_FILETYPE "fileType"
#define MZXML_PARENTFILE_ATTRIB_FILETYPE_RAWDATA "RAWData"
#define MZXML_PARENTFILE_ATTRIB_FILETYPE_PROCESSEDDATA "processedData"
#define MZXML_PARENTFILE_ATTRIB_FILESHA1 "fileSha1"

/* msInstrument delimiters */
#define MZXML_MSINSTRUMENT_OTAG "<msInstrument"
#define MZXML_MSINSTRUMENT_CTAG "</msInstrument"

/* dataProcessing section */
#define MZXML_DATAPROCESSING_OTAG "<dataProcessing"
#define MZXML_DATAPROCESSING_CTAG "</dataProcessing"
#define MZXML_DATAPROCESSING_ATTRIB_INTENSITYCUTOFF "intensityCutoff"
#define MZXML_DATAPROCESSING_ATTRIB_CENTROIDED "centroided"
#define MZXML_DATAPROCESSING_ATTRIB_DEISOTOPED "deisotoped"
#define MZXML_DATAPROCESSING_ATTRIB_CHARGEDECONVOLUTED "chargeDeconvoluted"
#define MZXML_DATAPROCESSING_ATTRIB_SPOTINTEGRATION "spotIntegration"

/* software section */
#define MZXML_SOFTWARE_OTAG "<software"
#define MZXML_SOFTWARE_ATTRIB_TYPE "type"
#define MZXML_SOFTWARE_ATTRIB_TYPE_ACQUISITION "acquisition"
#define MZXML_SOFTWARE_ATTRIB_TYPE_CONVERSION "conversion"
#define MZXML_SOFTWARE_ATTRIB_TYPE_PROCESSING "processing"
#define MZXML_SOFTWARE_ATTRIB_NAME "name"
#define MZXML_SOFTWARE_ATTRIB_VERSION "version"
#define MZXML_SOFTWARE_ATTRIB_COMPLETIONTIME "completionTime"

/* processingOperation segment */
#define MZXML_PROCESSINGOPERATION_OTAG "<processingOperation"
#define MZXML_PROCESSINGOPERATION_ATTRIB_NAME "name"
#define MZXML_PROCESSINGOPERATION_ATTRIB_VALUE "value"
#define MZXML_PROCESSINGOPERATION_ATTRIB_TYPE "type"

/* separation delimiters */
#define MZXML_SEPARATION_OTAG "<separation"
#define MZXML_SEPARATION_CTAG "</separation"

/* spotting delimiters */
#define MZXML_SPOTTING_OTAG "<spotting"
#define MZXML_SPOTTING_CTAG "</spotting"

/* comment segment */
#define MZXML_COMMENT_OTAG "<comment"
#define MZXML_COMMENT_CTAG "</comment"

/* scan segment */
#define MZXML_SCAN_OTAG "<scan"
#define MZXML_SCAN_CTAG "</scan"
#define MZXML_SCAN_ATTRIB_NUM "num"
#define MZXML_SCAN_ATTRIB_MSLEVEL "msLevel"
#define MZXML_SCAN_ATTRIB_PEAKSCOUNT "peaksCount"
#define MZXML_SCAN_ATTRIB_POLARITY "polarity"
#define MZXML_SCAN_ATTRIB_SCANTYPE "scanType"
#define MZXML_SCAN_ATTRIB_FILTERLINE "filterLine"
#define MZXML_SCAN_ATTRIB_CENTROIDED "centroided"
#define MZXML_SCAN_ATTRIB_DEISOTOPED "deisotoped"
#define MZXML_SCAN_ATTRIB_CHARGEDECONVOLUTED "chargeDeconvoluted"
#define MZXML_SCAN_ATTRIB_RETENTIONTIME "retentionTime"
#define MZXML_SCAN_ATTRIB_IONISATIONENERGY "ionisationEnergy"
#define MZXML_SCAN_ATTRIB_COLLISIONENERGY "collisionEnergy"
#define MZXML_SCAN_ATTRIB_COLLISIONGAS "collisionGas"					/* 2.1 version */
#define MZXML_SCAN_ATTRIB_CIDGASPRESSURE "cidGasPressure"				/* 3.1 version */
#define MZXML_SCAN_ATTRIB_COLLISIONGASPRESSURE "collisionGasPressure"	/* 2.1 version */
#define MZXML_SCAN_ATTRIB_STARTMZ "startMz\""
#define MZXML_SCAN_ATTRIB_ENDMZ "endMz"
#define MZXML_SCAN_ATTRIB_LOWMZ "lowMz"
#define MZXML_SCAN_ATTRIB_HIGHMZ "highMz"
#define MZXML_SCAN_ATTRIB_BASEPEAKMZ "basePeakMz"
#define MZXML_SCAN_ATTRIB_BASEPEAKINTENSITY "basePeakIntensity"
#define MZXML_SCAN_ATTRIB_TOTIONCURRENT "totIonCurrent"
#define MZXML_SCAN_ATTRIB_MSINSTRUMENTID "msInstrumentID"

/* scnaOrigin segment */
#define MZXML_SCANORIGIN_OTAG "<scanOrigin"
#define MZXML_SCANORIGIN_ATTRIB_PARENTFILEID "parentFileID"
#define MZXML_SCANORIGIN_ATTRIB_NUM "num"

/* precursorMz segment */
#define MZXML_PRECURSORMZ_OTAG "<precursorMz"
#define MZXML_PRECURSORMZ_CTAG "</precursorMz"
#define MZXML_PRECURSORMZ_ATTRIB_PRECURSORSCANNUM "precursorScanNum"
#define MZXML_PRECURSORMZ_ATTRIB_PRECURSORINTENSITY "precursorIntensity"
#define MZXML_PRECURSORMZ_ATTRIB_PRECURSORCHARGE "precursorCharge"
#define MZXML_PRECURSORMZ_ATTRIB_POSSIBLECHARGES "possibleCharges"
#define MZXML_PRECURSORMZ_ATTRIB_WINDOWWIDENESS "windowWideness"
#define MZXML_PRECURSORMZ_ATTRIB_ACTIVATIONMETHOD "activationMethod"

/* maldi segment */
#define MZXML_MALDI_OTAG "<maldi"
#define MZXML_MALDI_CTAG "</maldi"

/* peaks segment */
#define MZXML_PEAKS_OTAG "<peaks"
#define MZXML_PEAKS_CTAG "</peaks"
#define MZXML_PEAKS_ATTRIB_PRECISION "precision"
#define MZXML_PEAKS_ATTRIB_BYTEORDER "byteOrder"
#define MZXML_PEAKS_ATTRIB_PAIRORDER "pairOrder"
#define MZXML_PEAKS_ATTRIB_CONTENTTYPE "contentType"
#define MZXML_PEAKS_ATTRIB_COMPRESSIONTYPE "compressionType"
#define MZXML_PEAKS_ATTRIB_COMPRESSEDLEN "compressedLen"

/* nameValue segment */
#define MZXML_NAMEVALUE_OTAG "<nameValue"
#define MZXML_NAMEVALUE_ATTRIB_NAME "name"
#define MZXML_NAMEVALUE_ATTRIB_VALUE "value"
#define MZXML_NAMEVALUE_ATTRIB_TYPE "type"

/* sha1 segment */
#define MZXML_SHA1_OTAG "<sha1"
#define MZXML_SHA1_CTAG "</sha1"

/* index segment */
#define MZXML_INDEX_OTAG "<index"
#define MZXML_INDEX_CTAG "</index"
#define MZXML_INDEX_ATTRIB_NAME "name"

#define MZXML_OFFSET_OTAG "<offset"
#define MZXML_OFFSET_CTAG "</offset"
#define MZXML_OFFSET_ATTRIB_ID "id"

/* indexOffset segment */
#define MZXML_INDEXOFFSET_OTAG "<indexOffset"
#define MZXML_INDEXOFFSET_CTAG "</indexOffset"


#endif	/* MZXMLDEFS_H */

