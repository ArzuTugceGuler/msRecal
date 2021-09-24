/* 
 * File:   msRecalDefs.h
 * Author: Arzu Tugce Guler
 *
 * Created on 14 April 2016, 14:13
 */

#ifndef MSRECALDEFS_H
#define	MSRECALDEFS_H

#define HPLUS_MASS 1.00727646688
//#define DEFAULT_MIN_INTENSITY 100000
#define DEFAULT_MIN_INTENSITY 15	// do not use peaks below MIN_INTENSITY for calibration
#define AMINO_ACIDS "ARNDCEQGHILKMFPSTWYV"
#define ANTI_ACIDS "BJOUXZ"
//#define MAX_ROWS 8192
#define MAX_ROWS 450000
//#define MAX_ROWS 480000
#define MAX_CALIBRANTS 1000
#define DEFAULT_MODE 0						// lossless by default
#define DEFAULT_MIN_CALIBRANTS 3			// minimum number of internal calibrants to recalibrate
#define INTERNAL_CALIBRATION_TARGET 3	// 1.5e-6 discard internal calibrants that do not fit better than this */
#define DEFAULT_START_SCAN 1
#define DEFAULT_RT_LOWER 30
#define DEFAULT_RT_UPPER 90
#define DEFAULT_RECAL_OFFSET 0
#define DEFAULT_MMME 20
#define DEFAULT_SCORE "expect"

#endif	/* MSRECALDEFS_H */

