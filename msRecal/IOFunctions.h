/* 
 * File:   IOFunctions.h
 * Author: Arzu Tugce Guler
 *
 * Created on 15 April 2016, 14:53
 */

#ifndef IOFUNCTIONS_H
#define	IOFUNCTIONS_H

/* Updates the contents of a buffer by copying the bufpos residu to the beginning, and adding the rest from file */
int update_pepxml_buffer(char* buffer, char* bufpos, FILE* file_handle, int buffersize);

/* Copies data from the input file to the output file */
void copy_input_output(char* buffer, int buf_len, FILE* inputf, FILE* outputf, int copy_len);

#endif	/* IOFUNCTIONS_H */

