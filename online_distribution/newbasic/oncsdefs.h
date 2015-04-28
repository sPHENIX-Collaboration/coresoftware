// oncsdefs.h
// $Log: oncsdefs.h,v $
// Revision 1.4  2005/12/30 14:17:50  desmond
// ejd remove unused variable nameserver
//
// Revision 1.3  2004/11/15 16:28:02  phnxoncs
// chp: declare initialized char * as const char * (more apropriate and the Sun C compiler stops complaining
//
// Revision 1.2  2002/07/29 15:18:58  pinkenbu
// Make COMMON/inc and newbasic have a consistent oncsdef.h
//
// Revision 1.9  2001/09/10 14:06:32  phoncs
// ejd increase file limit
//
// Revision 1.8  2000/05/31 21:50:06  phoncs
// Changed RESULT to RESULTS. Needed to relieve name clash with /export/software/oncs/R2.6.3/online_distribution/Dcm code. Sorry about that. Adler, May 31st, 2000
//
// Revision 1.7  2000/04/07 15:31:07  phoncs
// ejd increase max objects definition
//
// Revision 1.6  2000/02/01 02:55:08  phoncs
// mlp - check in modifications
//
// Revision 1.5  1999/11/09 15:05:37  phoncs
// ejd add crate partition limit
//
/* general ONCS definitions
 * 
 * created:	Oct 1, 1996
 * 
 * author:	Ed Desmond
 * Module = %M%	Date = %D%	Version = %I%
 *
 * modifications:
 *		none
 *
 */

#ifndef ONCSDEFS_H
#define ONCSDEFS_H

#define MAX_OBJ_REF	2000
#define MAX_EVENTID	1000

// message queue parameters
#define MAX_MSGS 	10
#define MAX_MSG_LEN	1000
#define DEF_MSG_TIMEOUT  400000


// partition parameters
#define MAX_PARTITIONS_CRATE 5
// exit status


#define SUCCESS 0
#define FAILURE 1

//
// CHANGE; by CW Feb 17,98
// typedef RESULT  long
typedef long RESULTS;
typedef long HRESULT;
typedef int MSG_TIMEOUT;

// offsets of device information in constructor arguements
#define DEVICETYPE		0
#define DEVICEEVENTID		1
#define DEVICESUBEVENTID	2

#define MAXDCM_PER_DCB 5


#endif

