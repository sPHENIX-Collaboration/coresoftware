/* 
** formatError.C
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:13 $ 
** 
** $Log: formatError.C,v $
** Revision 1.1.1.1  2000/07/21 01:51:13  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:02:02  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**	formatError.C
**
**		Contains error-handling routines for "raw format" library
*/

#include "formatError.h"

/*
**  For now the error data will be kept in a statically 
**	  allocated structure. The actual allocation is done here.
*/
#define ERROR_PACKAGE_RAWFMT 1

typedef struct formatError {
	UINT subType ;
	ERRORVALUE errorNumber;

	PHDWORD* frameOrPacket_ptr;
	PHDWORD  additionalData;
} FORMATERROR;

FORMATERROR rawfmtLastError;

/*
**  Set an error status from a "frame" routine
*/
void setFrameError(ERRORVALUE errorNumber, PHDWORD* pointer, PHDWORD additionalData) {
	rawfmtLastError.subType = errorTypeFrame; 
	rawfmtLastError.errorNumber = errorNumber ;
	rawfmtLastError.frameOrPacket_ptr = pointer; 
	rawfmtLastError.additionalData = additionalData;
}

/*
**  Set an error status from a "packet" routine
*/
void setPacketError(ERRORVALUE errorNumber, PHDWORD* pointer, PHDWORD additionalData) {
	rawfmtLastError.subType = errorTypePacket;
	rawfmtLastError.errorNumber = errorNumber;
	rawfmtLastError.frameOrPacket_ptr = pointer;
	rawfmtLastError.additionalData = additionalData;
}

/*
**  Set a "user" error status 
*/
void setUserError(ERRORVALUE errorNumber, PHDWORD additionalData) {
	rawfmtLastError.subType = errorTypeUser;
	rawfmtLastError.errorNumber = errorNumber;
	rawfmtLastError.additionalData = additionalData;
}

/*
** Set success status for a "frame" routine.
*/
void setFrameSuccess () {
	rawfmtLastError.errorNumber = 0;
}

/*
**  Set success for a "packet" routine.
*/
void setPacketSuccess () {
	rawfmtLastError.errorNumber = 0;
}

/*
//  Returns the error code for the "last" recorded error
*/
ERRORVALUE formatGetErrorNumber ( ) {
	return rawfmtLastError.errorNumber;
}

/*
**  Returns the error number and unpacks error data for the last recorded error
*/
ERRORVALUE formatGetError (UINT* subType, PHDWORD** errorPointer, PHDWORD* additionalData) {
	*subType = rawfmtLastError.subType;
	*errorPointer = rawfmtLastError.frameOrPacket_ptr;
	*additionalData = rawfmtLastError.additionalData;

	return rawfmtLastError.errorNumber;
}

/*
**  Returns the pointer to the last failed frame or packet
*/
PHDWORD* formatGetErrorPointer () {
	if (rawfmtLastError.errorNumber == 0) 
		return 0;
	else 
		return rawfmtLastError.frameOrPacket_ptr;
}

/*
** Returns additional data about the last error
*/
PHDWORD formatGetErrorAdditionalData () {
	if (rawfmtLastError.errorNumber == 0) 
		return 0;
	else 
		return rawfmtLastError.additionalData;
}
