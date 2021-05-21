#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Base64Lib.h"

static int encode_char_val[MAX_ENC_CHAR_VALUE] = {
	65, 66, 67, 68, 69, 70, 71, 72,
	73, 74, 75, 76, 77, 78, 79, 80,
	81, 82, 83, 84, 85, 86, 87, 88,
	89, 90, 97, 98, 99, 100, 101, 102,
	103, 104, 105, 106, 107, 108, 109, 110,
	111, 112, 113, 114, 115, 116, 117, 118,
	119, 120, 121, 122, 48, 49, 50, 51,
	52, 53, 54, 55, 56, 57, 43, 47
};

static int decode_char_val[MAX_DEC_CHAR_VALUE] = {
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 62, -1, -1, -1, 63,
	52, 53, 54, 55, 56, 57, 58, 59, 60, 61, -1, -1, -1, -1, -1, -1,
	-1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
	15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, -1, -1, -1, -1, -1,
	-1, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
	41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};

/* Function that tests if the system representation is already network_order */
unsigned char test_network_order()
{
	short word = 0x4321;

	if((*(char *)& word) != 0x21)
		return 1;
	return 0;

}/* unsigned char test_network_order() */


/* Encoding function. Returns pointer to a new buffer, and puts bytes in network order */
char* encode_b64(unsigned char *input_buffer, char *output_buffer, int ibuff_len, int obufflen)
{
	int int_len, rest;
	int i, j, val, *addr;

	addr = (int*)output_buffer;
	rest = ibuff_len%3;
	int_len = ibuff_len - rest;

	i=0; j=0;
	while (i < int_len) {
		/* Reading 3 characters in the int and decoding them into 4 separate numbers */
		val = input_buffer[i] << 16 | input_buffer[i+1] << 8 | input_buffer[i+2];
		i += 3;
		addr[j] = (encode_char_val[val >> 18]) | (encode_char_val[val >> 12 & 0x3f] << 8) | (encode_char_val[val >> 6 & 0x3f] << 16 ) | (encode_char_val[val & 0x3f] << 24);
		j++;
	} /* while */
	if (rest == 1) {
		val = input_buffer[i];
		addr[j] = (encode_char_val[val >> 2]) | (encode_char_val[(val & 0x03) << 4] << 8) | 61 << 16 | 61 << 24;
		j++;
	}/* if */
	else if (rest == 2) {
		val = input_buffer[i] << 8 | input_buffer[i+1];
		i += 2;
		addr[j] = (encode_char_val[val >> 10]) | (encode_char_val[val >> 4 & 0x3f] << 8) | (encode_char_val[(val & 0x0f) << 2] << 16 ) | 61 << 24;
		j++;
	}/* else if */

	output_buffer[obufflen] = '\0';
	return output_buffer;

}/* char* encode_b64_network(char *input_buffer, int *new_len) */


/* Decoding function. returns new buffer with decoded material and puts bytes in host order */
char* decode_b64(unsigned char *input_buffer, int buff_len, int *new_len)
{
	char val1, val2, val3, val4;
	char *retval;
	int i, j;

	*new_len = ((buff_len * 3) / 4)+1;
	retval = malloc(*new_len * sizeof(char));

	i = 0;
	j = 0;

	while (i < buff_len) {
		do {
			val1 = decode_char_val[input_buffer[i]];
			i++;
		} while (i < buff_len && val1 == -1);
		if (val1 == -1)
			break;

		do {
			val2 = decode_char_val[input_buffer[i]];
			i++;
		} while (i < buff_len && val2 == -1);
		if (val2 == -1)
			break;

		retval[j] = ((val1 << 2) | ((val2 & 0x30) >> 4));
		j++;

		do {
			val3 = input_buffer[i];
			i++;
			if (val3 == 61) {
				val3 = -1;
				break;
			}/* if */
			val3 = decode_char_val[(unsigned char)val3];
		} while (i < buff_len && val3 == -1);
		if (val3 == -1)
			break;

		retval[j] = ((val2 & 0x0f) << 4) | ((val3 & 0x3c) >> 2);
		j++;

		do {
			val4 = input_buffer[i];
			i++;
			if (val4 == 61) {
				val4 = -1;
				break;
			}/* if */
			val4 = decode_char_val[(unsigned char)val4];
		} while (i < buff_len && val4 == -1);
		if (val4 == -1)
			break;

		retval[j] = ((val3 & 0x03) << 6) | val4;
		j++;
	}/* while */

	*new_len = j+1;
	retval[j] = '\0';

	return retval;

}

