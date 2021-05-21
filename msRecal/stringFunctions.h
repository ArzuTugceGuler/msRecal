/* 
 * File:   stringFunctions.h
 * Author: Arzu Tugce Guler
 *
 * Created on 14 April 2016, 15:19
 */

#ifndef STRINGFUNCTIONS_H
#define	STRINGFUNCTIONS_H

/* Returns the first position of the substring in the string */
int strpos(char* str, char* substr);

/* Returns the first position of the substring in the string, bounded by the window */
int strnpos(char* str, char* substr, int window);

/* Clones a string */
char* strclone(char* src);

#endif	/* STRINGFUNCTIONS_H */

