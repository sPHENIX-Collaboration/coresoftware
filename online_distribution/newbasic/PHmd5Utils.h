#ifndef __PHMD5UTILS_H__
#define __PHMD5UTILS_H__

#include "stdio.h"


#ifdef __cplusplus
extern "C" 
{
#endif


int PHmd5Stream(FILE *stream,  unsigned char  *digest,  int * filesize);

int PHmd5File(const char * filename,  unsigned char  *digest, int * filesize);



#ifdef __cplusplus
}  /* end extern "C"  */
#endif


#endif /* __PHMD5UTILS_H__ */
