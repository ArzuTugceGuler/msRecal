/* 
 * File:   pepXMLReaderDefs.h
 * Author: Arzu Tugce Guler
 *
 * Created on 15 April 2016, 15:17
 */

#ifndef PEPXMLREADERDEFS_H
#define	PEPXMLREADERDEFS_H

/* Configuration flags for pre-loading the run summary */ 
typedef unsigned char pipeline_analysis_flags;	
#define load_analysis_summary	1
#define load_dataset_derivation	2
#define load_all_analysis		3


/* Configuration flags for pre-loading the run summary */ 
typedef unsigned char run_summary_flags;	
#define load_sample_enzyme		1
#define load_search_summary		2
#define load_analysis_timestamp	4
#define load_all_summary		7


#endif	/* PEPXMLREADERDEFS_H */

