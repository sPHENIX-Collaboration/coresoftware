/* 
** Cframe.C
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:10 $ 
** 
** $Log: Cframe.C,v $
** Revision 1.1.1.1  2000/07/21 01:51:10  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.10  1999/05/17 19:47:16  markacs
** (stephen) added findNextError and whole indirection chain
**
** Revision 1.9  1999/04/16 20:12:11  markacs
** (stephen) fixing findFrameErrorStart, findFrameHistoryStart, adjustFrameLength to deal correctly with padding
**
** Revision 1.8  1999/04/13 18:27:38  markacs
** (stephen) put history block access in
**
** Revision 1.7  1998/12/16 15:40:52  markacs
** (stephen markacs) changes to rawFormat code to stop using bool,true,false
**
** Revision 1.6  1998/12/11 22:01:54  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
#include "phenixOnline.h"
#include "frameHdr.h"

#include "Cframe.h"
#include "CframeV1.h"
#include <stdlib.h>

/* Returns a pointer to the start of the frame data block.
** 
** This function does not depend on the version of the frame header since header
** length is in the same location for all versions.
*/
PTR_ret findFrameDataStart ( FRAME_ptr frame_ptr) {
#ifdef STRONG_CHECK
  if (!validFrameHdr(frame_ptr)) return ptrFailure;
#endif
  return (PHDWORD*)frame_ptr + getFrameHdrLength(frame_ptr);
}

/*  Return a pointer to the end of the frame data block/
 */
PTR_ret findFrameDataEnd ( FRAME_ptr frame_ptr) {
  PHDWORD* start_ptr = findFrameDataStart(frame_ptr);
  PHDWORD length = getFrameDataLength(frame_ptr);
  if ((start_ptr == ptrFailure)||(length==valueFailure))
    return ptrFailure;
  else 
    return start_ptr + length - 1;
}


/* 
** checkFrameHdrVersion checks that FrameMark and FrameHdrLength
**  are correct for a frame of the given version and, if so,
**  returns FrameHdrVersion.  Elsewise it returns valueFailure.
*/
VALUE_ret checkFrameHdrVersion ( FRAME_ptr frame_ptr)
{
  PHDWORD frameMark;
  PHDWORD frameHdrLength;
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure) {
    frameMark = getFrameMark(frame_ptr);
    frameHdrLength = getFrameHdrLength(frame_ptr);
    if ( (frameMark == frameMarkV[version]) &&
         (frameHdrLength == frameHdrLengthV[version])) {
      setFrameSuccess();
      return version;
    }
    else {
      setFrameError(FORMAT_ERR_INVALID_FRAMEMARK, frame_ptr, frameMark);
      return valueFailure;
    }
  }
  else return valueFailure;
}

/*  Return a pointer to the last word in the frame.
 */
PTR_ret findFrameEnd ( FRAME_ptr frame_ptr)
{
  if (validFrameHdr(frame_ptr))
  {
    PHDWORD frameLength = getFrameLength(frame_ptr);
    if (frameLength != valueFailure)
    {
      setFrameSuccess();
      return ( (PHDWORD*)frame_ptr + frameLength - 1);
    } 
    else return ptrFailure;
  }
  else return ptrFailure;
}

/*
** checkFrameEndianism returns 1 if the frame header is valid
**  as it is, returns 0 if the frame header would be valid if swapped,
**  and returns valueFailure otherwise.
*/
VALUE_ret checkFrameEndianism(FRAME_ptr frame_ptr)
{
  if (validFrameHdr(frame_ptr)) return 1;
  else 
    { PHDWORD * buffer = (PHDWORD*)malloc(4*maxFrameHdrLength);
      dwordByteSwap(buffer, frame_ptr, maxFrameHdrLength);
      if (validFrameHdr(buffer)) return 0;
      else return valueFailure;
    }
}



/* adjustFrameLength adds addDwords to the frame, with padding 
   if pad is TRUE, and returns the new frame length */
VALUE_ret adjustFrameLength (FRAME_ptr frame_ptr, UINT maxFrameLength, 
                             UINT addDwords     , LOGIC_ret pad)
{
  UINT padDwords;
  UINT currentLength;
  UINT newDwords;
  UINT newLength;

#ifdef STRONGCHECK
  if (!validFrameHdr(frame_ptr)) return valueFailure;
#endif

  if (pad) removeFramePadding(frame_ptr);
  currentLength = getFrameLength(frame_ptr);
			 
  newDwords = addDwords;
  newLength = currentLength + newDwords;

  /*  Perform the padding calculation if requested
   */
  if (pad) {
    UINT modulo = newLength % currentFrameQuantum;
    if (modulo > 0) padDwords = currentFrameQuantum - modulo;
    else padDwords = 0;

    newLength += padDwords;
    newDwords += padDwords;
  }
  else padDwords = 0;

  /*  Check for buffer overflow
   */
  if (newLength > maxFrameLength) {
    setFrameError(FORMAT_ERR_BUFFER_OVERFLOW, frame_ptr, newLength);
    return valueFailure;
  }

  /* 
  ** below:Check to make sure we don't overflow the 32 bit length field.
  ** Note that this calculation may look "funny", but it is necessary
  ** to evaluate it like this to prevent the calculation from
  ** suffering a 32-bit overflow itself.  Also note that 
  ** a length that fills the 32-bit length field corresponds
  ** to a frame taking up over 16 gigabytes of memory.               
  */
  if (maxDwordValue - currentLength < newDwords) {
    setFrameError(FORMAT_ERR_LENGTH_OVERFLOW, frame_ptr, newLength);
    return valueFailure;
  }

  setFrameLength(frame_ptr, newLength);
  if (padDwords != 0)
  {
    PHDWORD* pad_ptr = frame_ptr + newLength - padDwords;
    dwordClear (pad_ptr, padDwords);  /* clear padding words */
    setFramePadding(frame_ptr, padDwords);
  }
  return newLength;
}


/* removeFramePadding strips padding off a frame and
   returns the length of the frame                      */

VALUE_ret removeFramePadding (FRAME_ptr frame_ptr)
{
  UINT length;
  UINT padDwords = getFramePadding(frame_ptr);
  length = getFrameLength(frame_ptr);
  length-=padDwords;
  setFrameLength(frame_ptr, length);
  setFramePadding(frame_ptr, 0);
  return length;
}


/* Here begins the code for the version-DEPENDENT functions */
/* (only if not in DCM mode)                                */
/* ======================================================== */

#ifndef DCM

typedef LOGIC_ret CHECKFUNCTION ( FRAME_ptr);
typedef VALUE_ret ACCESSFUNCTION ( FRAME_ptr);
typedef PTR_ret   PTRACCESSFUNCTION ( FRAME_ptr);
typedef VALUE_ret MODIFYFUNCTION (FRAME_ptr, UINT);
typedef LOGIC_ret LOGICALMODIFYFUNCTION (FRAME_ptr, UINT);
typedef VALUE_ret VARRAYCOPYFUNCTION(FRAME_ptr, PHDWORD*, UINT);
typedef LOGIC_ret LARRAYCOPYFUNCTION(FRAME_ptr, PHDWORD*, UINT);
typedef VALUE_ret MAKEHDRFUNCTION(FRAME_ptr, UINT, UINT, UINT, UINT);
typedef VALUE_ret PARAMACCESSFUNCTION (FRAME_ptr, UINT);
typedef PTR_ret   FINDNEXTFUNCTION (FRAME_ptr, FRAME_ptr);

typedef ACCESSFUNCTION* ACCESSFUNCTIONPTR_arr[NUM_FRAME_VERSIONS]; 
typedef PTRACCESSFUNCTION* PTRACCESSFUNCTIONPTR_arr[NUM_FRAME_VERSIONS];
typedef CHECKFUNCTION*  CHECKFUNCTIONPTR_arr[NUM_FRAME_VERSIONS];
typedef MODIFYFUNCTION* MODIFYFUNCTIONPTR_arr[NUM_FRAME_VERSIONS]; 
typedef LOGICALMODIFYFUNCTION* LOGICALMODIFYFUNCTIONPTR_arr[NUM_FRAME_VERSIONS];
typedef VARRAYCOPYFUNCTION* VARRAYCOPYFUNCTIONPTR_arr[NUM_FRAME_VERSIONS];
typedef LARRAYCOPYFUNCTION* LARRAYCOPYFUNCTIONPTR_arr[NUM_FRAME_VERSIONS];
typedef MAKEHDRFUNCTION* MAKEHDRFUNCTIONPTR_arr[NUM_FRAME_VERSIONS];
typedef PARAMACCESSFUNCTION* PARAMACCESSFUNCTIONPTR_arr[NUM_FRAME_VERSIONS];
typedef FINDNEXTFUNCTION* FINDNEXTFUNCTIONPTR_arr[NUM_FRAME_VERSIONS];

/*
** IMPORTANT: The zeroth element in the indirection arrays is defined in the
**            lines below to be a null pointer.  It is thus very important
**            that anywhere the array is deferenced, it is checked that
**            it is not the zeroth element that is dereferenced.  If it
**            is, a null pointer will be dereferenced and s segmentation
**            fault or some other ugliness will occur.
*/
ACCESSFUNCTIONPTR_arr CONSTANT getFrameSourceIdV = {0, &getFrameSourceIdV1};
ACCESSFUNCTIONPTR_arr CONSTANT getFrameTypeV = {0, &getFrameTypeV1};
ACCESSFUNCTIONPTR_arr CONSTANT getFrameStatusV = {0, &getFrameStatusV1};
ACCESSFUNCTIONPTR_arr CONSTANT getFrameDataTypeV = {0, &getFrameDataTypeV1};
ACCESSFUNCTIONPTR_arr CONSTANT getFramePaddingV = {0, &getFramePaddingV1};
ACCESSFUNCTIONPTR_arr CONSTANT getFrameHistoryLengthV = {0, &getFrameHistoryLengthV1};
ACCESSFUNCTIONPTR_arr CONSTANT getFrameErrorLengthV = {0, &getFrameErrorLengthV1};
ACCESSFUNCTIONPTR_arr CONSTANT getFrameDataLengthV = {0, &getFrameDataLengthV1};
ACCESSFUNCTIONPTR_arr CONSTANT getFrameAlignLengthV = {0, &getFrameAlignLengthV1};

PTRACCESSFUNCTIONPTR_arr CONSTANT findFrameErrorStartV = {0, &findFrameErrorStartV1};
PTRACCESSFUNCTIONPTR_arr CONSTANT findFrameAlignBlockV = {0, &findFrameAlignBlockV1};
PTRACCESSFUNCTIONPTR_arr CONSTANT findFrameHistoryStartV = {0, &findFrameHistoryStartV1};

CHECKFUNCTIONPTR_arr  CONSTANT validFrameHdrV = {0, &validFrameHdrV1};
CHECKFUNCTIONPTR_arr  CONSTANT emptyFrameV = {0, &emptyFrameV1};

MODIFYFUNCTIONPTR_arr CONSTANT orFrameStatusV = {0, &orFrameStatusV1};
MODIFYFUNCTIONPTR_arr CONSTANT setFramePaddingV = {0, &setFramePaddingV1};
MODIFYFUNCTIONPTR_arr CONSTANT adjustFrameDataLengthV = {0, &adjustFrameDataLengthV1};
MODIFYFUNCTIONPTR_arr CONSTANT adjustFrameHistoryLengthV = {0, &adjustFrameHistoryLengthV1};
MODIFYFUNCTIONPTR_arr CONSTANT adjustFrameErrorLengthV = {0, &adjustFrameErrorLengthV1};

LOGICALMODIFYFUNCTIONPTR_arr CONSTANT setDataTypeV = {0, &setDataTypeV1};
LOGICALMODIFYFUNCTIONPTR_arr CONSTANT setFrameTypeV = {0, &setFrameTypeV1};
LOGICALMODIFYFUNCTIONPTR_arr CONSTANT setSourceIdV = {0, &setSourceIdV1};
LOGICALMODIFYFUNCTIONPTR_arr CONSTANT setFrameHistoryLengthV = {0, &setFrameHistoryLengthV1};
LOGICALMODIFYFUNCTIONPTR_arr CONSTANT setFrameErrorLengthV = {0, &setFrameErrorLengthV1};
LOGICALMODIFYFUNCTIONPTR_arr CONSTANT setFrameAlignLengthV = {0, &setFrameAlignLengthV1}; 
LOGICALMODIFYFUNCTIONPTR_arr CONSTANT setFrameStatusV = {0, &setFrameStatusV1};

VARRAYCOPYFUNCTIONPTR_arr CONSTANT getAlignBlockV = {0, &getAlignBlockV1};
LARRAYCOPYFUNCTIONPTR_arr CONSTANT setAlignBlockV = {0, &setAlignBlockV1};

MAKEHDRFUNCTIONPTR_arr CONSTANT makeFrameHdrV = {0, &makeFrameHdrV1};

PARAMACCESSFUNCTIONPTR_arr CONSTANT getHistoryEntryV = {0, &getHistoryEntryV1};
PARAMACCESSFUNCTIONPTR_arr CONSTANT getHistoryStageV = {0, &getHistoryStageV1};
PARAMACCESSFUNCTIONPTR_arr CONSTANT getHistorySourceIndexV = {0, &getHistorySourceIndexV1};
PARAMACCESSFUNCTIONPTR_arr CONSTANT getHistoryStatusV = {0, &getHistoryStatusV1};
FINDNEXTFUNCTIONPTR_arr CONSTANT findNextErrorV = {0, &findNextErrorV1};

VALUE_ret getFrameStatus ( FRAME_ptr frame_ptr) 
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getFrameStatusV[version])(frame_ptr);
  else 
    return valueFailure;
}

VALUE_ret getFrameType ( FRAME_ptr frame_ptr) 
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getFrameTypeV[version])(frame_ptr);
  else 
    return valueFailure;
}


VALUE_ret getFrameSourceId ( FRAME_ptr frame_ptr) 
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getFrameSourceIdV[version])(frame_ptr);
  else 
    return valueFailure;
}


VALUE_ret getFrameDataType ( FRAME_ptr frame_ptr) 
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getFrameDataTypeV[version])(frame_ptr);
  else 
    return valueFailure;
}


VALUE_ret getFrameDataLength ( FRAME_ptr frame_ptr) 
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getFrameDataLengthV[version])(frame_ptr);
  else 
    return valueFailure;
}


VALUE_ret getFrameHistoryLength ( FRAME_ptr frame_ptr)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getFrameHistoryLengthV[version])(frame_ptr);
  else 
    return valueFailure;
}


VALUE_ret getFrameErrorLength ( FRAME_ptr frame_ptr)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getFrameErrorLengthV[version])(frame_ptr);
  else 
    return valueFailure;
}


VALUE_ret getFrameAlignLength ( FRAME_ptr frame_ptr)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getFrameAlignLengthV[version])(frame_ptr);
  else 
    return valueFailure;
}


PTR_ret findFrameAlignBlock ( FRAME_ptr frame_ptr) 
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*findFrameAlignBlockV[version])(frame_ptr);
  else 
    return ptrFailure;
}


PTR_ret findFrameHistoryStart ( FRAME_ptr frame_ptr)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*findFrameHistoryStartV[version])(frame_ptr);
  else 
    return ptrFailure;
}


PTR_ret findFrameErrorStart ( FRAME_ptr frame_ptr)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*findFrameErrorStartV[version])(frame_ptr);
  else
    return ptrFailure;
}


VALUE_ret orFrameStatus (FRAME_ptr frame_ptr, UINT statusBits)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure) 
    return (*orFrameStatusV[version])(frame_ptr, statusBits);
  else 
    return valueFailure;
}


VALUE_ret getFramePadding ( FRAME_ptr frame_ptr)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getFramePaddingV[version])(frame_ptr);
  else 
    return valueFailure;
}


/*  Construct a header for a frame using the "current" header version
 */
VALUE_ret makeFrameHdr (FRAME_ptr frame_ptr,  UINT maxFrameLen,  UINT dataType, 
			UINT frameType,  UINT sourceId)
{
  return (*makeFrameHdrV[currentFrameHdrVersion])(frame_ptr, maxFrameLen, dataType,
						  frameType, sourceId);
}

VALUE_ret adjustFrameDataLength (FRAME_ptr frame_ptr, UINT newDwords)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*adjustFrameDataLengthV[version])(frame_ptr, newDwords);
  else
    return valueFailure;
}


VALUE_ret adjustFrameHistoryLength (FRAME_ptr frame_ptr, UINT newDwords)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*adjustFrameHistoryLengthV[version])(frame_ptr, newDwords);
  else
    return valueFailure;
}

VALUE_ret adjustFrameErrorLength (FRAME_ptr frame_ptr, UINT newDwords)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*adjustFrameErrorLengthV[version])(frame_ptr, newDwords);
  else
    return valueFailure;
}

VALUE_ret setFramePadding (FRAME_ptr frame_ptr, UINT padDwords)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*setFramePaddingV[version])(frame_ptr, padDwords);
  else
    return valueFailure;
}

LOGIC_ret setDataType (FRAME_ptr frame_ptr, UINT dataType)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*setDataTypeV[version])(frame_ptr, dataType);
  else 
    return FALSE;
}

LOGIC_ret setFrameType (FRAME_ptr frame_ptr, UINT frameType)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*setFrameTypeV[version])(frame_ptr, frameType);
  else
    return FALSE;
}

LOGIC_ret setSourceId (FRAME_ptr frame_ptr, UINT sourceId)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*setSourceIdV[version])(frame_ptr, sourceId);
  else
    return FALSE;
}

LOGIC_ret setFrameHistoryLength (FRAME_ptr frame_ptr, UINT historyLength)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*setFrameHistoryLengthV[version])(frame_ptr, historyLength);
  else
    return FALSE;
}

LOGIC_ret setFrameErrorLength (FRAME_ptr frame_ptr, UINT errorLength)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*setFrameErrorLengthV[version])(frame_ptr, errorLength);
  else
    return FALSE;
}

LOGIC_ret setFrameAlignLength (FRAME_ptr frame_ptr, UINT alignLength)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*setFrameAlignLengthV[version])(frame_ptr, alignLength);
  else 
    return FALSE;
}

LOGIC_ret setFrameStatus (FRAME_ptr frame_ptr, UINT status)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*setFrameStatusV[version])(frame_ptr, status);
  else 
    return FALSE;
}

VALUE_ret getAlignBlock (FRAME_ptr frame_ptr, PHDWORD* alignDestination, UINT maxNumDwords)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getAlignBlockV[version])(frame_ptr, alignDestination, maxNumDwords);
  else
    return valueFailure;
}

LOGIC_ret setAlignBlock (FRAME_ptr frame_ptr,PHDWORD* alignSource, UINT numDwords)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*setAlignBlockV[version])(frame_ptr, alignSource, numDwords);
  else
    return FALSE;
}

VALUE_ret getHistoryEntry (FRAME_ptr frame_ptr, UINT index)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getHistoryEntryV[version])(frame_ptr, index);
  else
    return valueFailure;
}

VALUE_ret getHistoryStage (FRAME_ptr frame_ptr, UINT index)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getHistoryStageV[version])(frame_ptr, index);
  else
    return valueFailure;
}

VALUE_ret getHistorySourceIndex (FRAME_ptr frame_ptr, UINT index)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getHistorySourceIndexV[version])(frame_ptr, index);
  else
    return valueFailure;
}

VALUE_ret getHistoryStatus (FRAME_ptr frame_ptr, UINT index)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*getHistoryStatusV[version])(frame_ptr, index);
  else
    return valueFailure;
}

PTR_ret findNextError (FRAME_ptr frame_ptr, FRAME_ptr thisError)
{
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (*findNextErrorV[version])(frame_ptr, thisError);
  else
    return ptrFailure;
}

#endif /* end of ifndef DCM */

/* ====================================================== */
/* Here ends the code for the version-DEPENDENT functions */












