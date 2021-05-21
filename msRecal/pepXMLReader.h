/* 
 * File:   pepXMLReader.h
 * Author: Arzu Tugce Guler
 *
 * Created on 15 April 2016, 15:17
 */

#ifndef PEPXMLREADER_H
#define	PEPXMLREADER_H

#include <stdio.h>

#include "pepXMLDelegates.h"
#include "pepXMLReaderDefs.h"
#include "pepXMLStructures.h"

/* Main interface functions */
/* Initial reader, read the entire file */
pmsms_pipeline_analysis read_pepxml_file(char* file_name, pipeline_analysis_flags paflags, run_summary_flags rsflags, pdelegate_list dlist);

/* Gets the analysis summary element of the pipeline summary */
panalysis_summary get_analysis_summary(pmsms_pipeline_analysis analysis, int* array_length, pdelegate_list dlist);

/* Gets the dataset derivation element of the pipeline summary */
pdataset_derivation get_dataset_derivation(pmsms_pipeline_analysis analysis);

/* Gets the sample enzyme element from a specific run summary */
psample_enzyme run_summary_sample_enzyme(pmsms_run_summary run_summary);

/* Gets the search_summary element array from a specific run summary */
psearch_summary run_summary_search_summary(pmsms_run_summary run_summary, int* array_length);

/* Gets the analysis timestamp element array from a specific run summary */
panalysis_timestamp run_summary_analysis_timestamp(pmsms_run_summary run_summary, int* array_length);

#endif	/* PEPXMLREADER_H */

