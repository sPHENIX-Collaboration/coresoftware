#ifndef _phenixTypesIncludeProtection_
#define _phenixTypesIncludeProtection_
/* 
** __DECCXX is defined by the DEC compiler on the alpha.  The definitions
** here deal with the differences in sizes between 32 and 64 bit machines.
*/


#ifdef __DECCXX
typedef unsigned int PHDWORD;
typedef unsigned short SWORD;
typedef unsigned char BYTE;
typedef unsigned int UINT;
#else

/* #if defined(_WIN32_WINNT) */
#ifdef WIN32
typedef unsigned long	PHDWORD;
#else
typedef unsigned int	PHDWORD;
#endif


typedef unsigned short	SWORD;
typedef unsigned char	BYTE;
#ifndef VXWORKS
typedef unsigned int	UINT;
#else
#include <types/vxTypesOld.h>
#endif
#endif /* __DECCXX */


#endif 

     /* end of ifndef _phenixTypesIncludeProtection_ */

