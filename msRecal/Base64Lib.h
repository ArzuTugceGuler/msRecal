/* 
 * File:   Base64Lib.h
 * Author: Arzu Tugce Guler
 *
 * Created on 15 April 2016, 16:39
 */

#ifndef BASE64LIB_H
#define	BASE64LIB_H

#define MAX_ENC_CHAR_VALUE	64
#define MAX_DEC_CHAR_VALUE	256

/* Encoding function. Extends and overwrites existing buffer, and puts bytes in network order */
char* encode_b64(unsigned char *input_buffer, char *output_buffer, int ibuff_len, int obufflen);

/* Encoding function. Extends and overwrites existing buffer */
char* decode_b64(unsigned char *input_buffer, int buff_len, int *new_len);
//char* decode_b64(char *input_buffer, int buff_len, int *new_len);

#endif	/* BASE64LIB_H */

