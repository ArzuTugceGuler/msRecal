/*
 *  sha.cpp
 *
 *  Copyright (C) 1998, 2009
 *  Paul E. Jones <paulej@packetizer.com>
 *  All Rights Reserved
 *
 *****************************************************************************
 *  $Id: sha.c 12 2009-06-22 19:34:25Z paulej $
 *****************************************************************************
 *
 *  Description:
 *      This utility will display the message digest (fingerprint) for
 *      the specified file(s).
 *
 *  Portability Issues:
 *      None.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#include <io.h>
#endif
#include <fcntl.h>
#include "SHA.h"
#include "SHA1.h"


/* Returns a sha1 string of 40 characters long */
char* get_sha1_string(char* file_name)
{
    SHA1Context sha;                /* SHA-1 context                 */
    FILE        *fp;                /* File pointer for reading files*/
    char        c;                  /* Character read from file      */
	char		*result;

	/* Working the sha1 routine */
	result = malloc(41 * sizeof(char));
	fp = fopen(file_name,"rb");
    if (fp) {
        SHA1Reset(&sha);
		c = fgetc(fp);
        while(!feof(fp)) {
            SHA1Input(&sha, &c, 1);
            c = fgetc(fp);
        }/* while */
		fclose(fp);
		sprintf(result, "%08X%08X%08X%08X%08X", sha.Message_Digest[0], sha.Message_Digest[1], sha.Message_Digest[2], sha.Message_Digest[3], sha.Message_Digest[4]);
	}

   return result;
}

