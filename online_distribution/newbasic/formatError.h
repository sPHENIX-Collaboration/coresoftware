/* 
** formatError.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:13 $ 
** 
** $Log: formatError.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:13  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.5  1998/12/17 21:57:06  phoncs
** (stephen markacs) better bounds checking in version gets
**
** Revision 1.4  1998/12/11 22:01:18  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  formatError.h
**
**   Definitions of various error-handling macros, prototypes, variables etc
*/

#ifndef _FORMATERROR_
#define _FORMATERROR_

#include "phenixOnline.h"

/*
**  Use C linkage for error routines
*/
#ifdef __cplusplus
extern "C" {
#endif

typedef UINT ERRORVALUE;

#define FORMAT_ERROR_SUBTYPE_FRAME 1
#define FORMAT_ERROR_SUBTYPE_PACKET 2
#define FORMAT_ERROR_SUBTYPE_USER 3

CONSTANT UINT errorTypeFrame = FORMAT_ERROR_SUBTYPE_FRAME;
CONSTANT UINT errorTypePacket = FORMAT_ERROR_SUBTYPE_PACKET;
CONSTANT UINT errorTypeUser = FORMAT_ERROR_SUBTYPE_USER;

/*
**  Some function prototypes
*/
void setFrameError(ERRORVALUE, PHDWORD*, PHDWORD);

void setPacketError(ERRORVALUE, PHDWORD*, PHDWORD);

void setUserError(ERRORVALUE, PHDWORD);

void setFrameSuccess ();

void setPacketSuccess ();

ERRORVALUE formatGetError (UINT*, PHDWORD**, PHDWORD*);

ERRORVALUE formatGetErrorNumber ( );

PHDWORD* formatGetErrorPointer ( );

PHDWORD formatGetErrorAdditionalData ( );

/*
**  Error code definitions
*/
enum formatErrorCodes {
  FORAMT_ERR_SUCCESS = 0,
  FORMAT_ERR_INVALID_HEADER,
  FORMAT_ERR_INVALID_HDRVERSION,
  FORMAT_ERR_INVALID_PACKET_HDRVERSION,
  FORMAT_ERR_INVALID_DATA,
  FORMAT_ERR_HISTORY_OVERFLOW,
  FORMAT_ERR_ERROR_OVERFLOW,
  FORMAT_ERR_BUFFER_OVERFLOW,
  FORMAT_ERR_INVALID_MODIFY,
  FORMAT_ERR_INVALID_FRAMEMARK,
  FORMAT_ERR_LENGTH_OVERFLOW,
  FORMAT_ERR_INVALID_HDRLENGTH,
  FORMAT_ERR_INVALID_APPEND,
  FORMAT_ERR_WRONG_STRUCTURE,
  FORMAT_ERR_HDR_INCONSISTENCY,
  FORMAT_ERR_DATA_INCONSISTENCY,
  FORMAT_ERR_NONEMPTY_PACKET
};

/*
**  Use C linkage for below structures
*/
#ifdef __cplusplus
}
#endif

#endif





