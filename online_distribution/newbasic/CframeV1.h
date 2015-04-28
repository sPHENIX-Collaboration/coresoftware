/* 
** CframeV1.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:10 $ 
** 
** $Log: CframeV1.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:10  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.9  1999/05/17 19:47:06  markacs
** (stephen) added findNextError and whole indirection chain
**
** Revision 1.8  1999/04/16 20:12:10  markacs
** (stephen) fixing findFrameErrorStart, findFrameHistoryStart, adjustFrameLength to deal correctly with padding
**
** Revision 1.7  1999/04/13 18:27:37  markacs
** (stephen) put history block access in
**
** Revision 1.6  1998/12/16 15:40:51  markacs
** (stephen markacs) changes to rawFormat code to stop using bool,true,false
**
** Revision 1.5  1998/12/11 22:01:12  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
** CframeV1.h
**
** Now contains combination of old CframeV1.h and CframeV1.C.
** This change was made to allow inlining
** 
** The layout of the file is as follows:
**
** (1) definition of version-dependent field parameters
** (2) includes of phenixOnline.h, etc.
** (3) prototypes of functions defined later in this file
** (4) include of Cframe.h
** (5) defintion of V1 functions
**
**
** For additional comments, see Cframe.h
**
*/

#ifndef _CFRAMEV1_
#define _CFRAMEV1_

#ifdef __cplusplus
extern "C" {
#endif

#include "phenixOnline.h"
#include "formatError.h"
#include "framePublic.h"
#include "frameV1Public.h"
#include "frameHdrV1.h"
#include "Cframe.h"

/*
**  Function prototype definitions
**  ============================================================
*/

VALUE_ret makeFrameHdrV1 (PHDWORD*, UINT, UINT, UINT, UINT);
VALUE_ret getFrameDataLengthV1 (FRAME_ptr);
VALUE_ret adjustFrameHistoryLengthV1 (FRAME_ptr, UINT);
VALUE_ret adjustFrameErrorLengthV1 (FRAME_ptr, UINT);
VALUE_ret orFrameStatusV1 (FRAME_ptr, UINT);

INLINE_P LOGIC_ret validFrameHdrV1 (FRAME_ptr);
INLINE_P LOGIC_ret emptyFrameV1 (FRAME_ptr);

INLINE_P VALUE_ret getFrameErrorLengthV1 (FRAME_ptr);
INLINE_P VALUE_ret getFrameHistoryLengthV1 (FRAME_ptr);
INLINE_P VALUE_ret getFrameAlignLengthV1 (FRAME_ptr);
INLINE_P VALUE_ret getFramePaddingV1 (FRAME_ptr);

INLINE_P VALUE_ret getFrameSourceIdV1 (FRAME_ptr);
INLINE_P VALUE_ret getFrameTypeV1 (FRAME_ptr);
INLINE_P VALUE_ret getFrameDataTypeV1 (FRAME_ptr);
INLINE_P VALUE_ret getFrameStatusV1 (FRAME_ptr);

INLINE_P VALUE_ret adjustFrameDataLengthV1 (FRAME_ptr, UINT);
INLINE_P VALUE_ret setFramePaddingV1 (FRAME_ptr, UINT);

INLINE_P LOGIC_ret setDataTypeV1(FRAME_ptr, UINT);
INLINE_P LOGIC_ret setFrameTypeV1(FRAME_ptr, UINT);
INLINE_P LOGIC_ret setSourceIdV1(FRAME_ptr, UINT);
INLINE_P LOGIC_ret setFrameHistoryLengthV1(FRAME_ptr, UINT);
INLINE_P LOGIC_ret setFrameErrorLengthV1(FRAME_ptr, UINT);
INLINE_P LOGIC_ret setFrameAlignLengthV1(FRAME_ptr, UINT);
INLINE_P LOGIC_ret setFrameStatusV1(FRAME_ptr, UINT);

INLINE_P PTR_ret findFrameAlignBlockV1 (FRAME_ptr);
INLINE_P PTR_ret findFrameErrorStartV1 (FRAME_ptr);
INLINE_P PTR_ret findFrameHistoryStartV1 (FRAME_ptr);

INLINE_P VALUE_ret getAlignBlockV1 (FRAME_ptr, PHDWORD*, UINT);
INLINE_P LOGIC_ret setAlignBlockV1 (FRAME_ptr, PHDWORD*, UINT);

INLINE_P VALUE_ret getHistoryEntryV1 (FRAME_ptr, UINT);
INLINE_P VALUE_ret getHistoryStageV1 (FRAME_ptr, UINT);
INLINE_P VALUE_ret getHistorySourceIndexV1 (FRAME_ptr, UINT);
INLINE_P VALUE_ret getHistoryStatusV1 (FRAME_ptr, UINT);

/* note: the following three are not referenced by indirection
   from outside */
INLINE_P VALUE_ret getStageFromHistoryEntryV1 (PHDWORD);  
INLINE_P VALUE_ret getSourceIndexFromHistoryEntryV1 (PHDWORD);
INLINE_P VALUE_ret getStatusFromHistoryEntryV1 (PHDWORD);

INLINE_P PTR_ret findNextErrorV1 (FRAME_ptr, FRAME_ptr);

/*
** Inlined routines
** ============================================================
*/

/* Check to see if this can really be a valid V1 header 
*/
INLINE_D LOGIC_ret validFrameHdrV1 (FRAME_ptr frame_ptr)
{
  if (getFrameHdrLength(frame_ptr) == V1_FRAMEHDR_LENGTH) {
    setFrameSuccess();
    return TRUE;
  }
  else {
    setFrameError (FORMAT_ERR_INVALID_HEADER, frame_ptr, 0);
    return FALSE;
  }
}


/* Check to see if this is an empty frame.
*/
INLINE_D LOGIC_ret emptyFrameV1 (FRAME_ptr frame_ptr)
{
  return getFrameHdrLength(frame_ptr) == getFrameLength(frame_ptr);
}


INLINE_D VALUE_ret getFramePaddingV1 (FRAME_ptr frame_ptr)
{
  return getBitsMACRO(frame_ptr,               
		      PADDING_OFFSET_OF_DWORD,
                      PADDING_OFFSET_IN_DWORD, 
		      PADDING_MASK            );
}


INLINE_D VALUE_ret getFrameSourceIdV1 (FRAME_ptr frame_ptr)
{
  return getBitsMACRO(frame_ptr,                 
		      SOURCE_ID_OFFSET_OF_DWORD,
                      SOURCE_ID_OFFSET_IN_DWORD, 
		      SOURCE_ID_MASK            );
}


INLINE_D VALUE_ret getFrameDataTypeV1 (FRAME_ptr frame_ptr)
{
  return getBitsMACRO(frame_ptr,                 
		      DATA_TYPE_OFFSET_OF_DWORD,
                      DATA_TYPE_OFFSET_IN_DWORD, 
		      DATA_TYPE_MASK            );
}


INLINE_D VALUE_ret getFrameTypeV1 (FRAME_ptr frame_ptr)
{
  return getBitsMACRO(frame_ptr,                  
		      FRAME_TYPE_OFFSET_OF_DWORD,
                      FRAME_TYPE_OFFSET_IN_DWORD, 
		      FRAME_TYPE_MASK           );
}


INLINE_D VALUE_ret getFrameStatusV1 (FRAME_ptr frame_ptr)
{
  return getBitsMACRO(frame_ptr,              
		      STATUS_OFFSET_OF_DWORD,
                      STATUS_OFFSET_IN_DWORD, 
		      STATUS_MASK            );
}


INLINE_D VALUE_ret getFrameErrorLengthV1 (FRAME_ptr frame_ptr) 
{
  return getBitsMACRO(frame_ptr,                    
		      ERROR_LENGTH_OFFSET_OF_DWORD,
                      ERROR_LENGTH_OFFSET_IN_DWORD, 
		      ERROR_LENGTH_MASK           );
}


INLINE_D VALUE_ret getFrameHistoryLengthV1 (FRAME_ptr frame_ptr) 
{
  return getBitsMACRO(frame_ptr,                       
		      HISTORY_LENGTH_OFFSET_OF_DWORD,
                      HISTORY_LENGTH_OFFSET_IN_DWORD,  
		      HISTORY_LENGTH_MASK           );
}


INLINE_D VALUE_ret getFrameAlignLengthV1 (FRAME_ptr frame_ptr) 
{
  return getBitsMACRO(frame_ptr,                      
		      ALIGN_LENGTH_OFFSET_OF_DWORD,
                      ALIGN_LENGTH_OFFSET_IN_DWORD,   
		      ALIGN_LENGTH_MASK            );
}


/* Return pointer to start of error block in frame 
*/
INLINE_D PTR_ret findFrameErrorStartV1 (FRAME_ptr frame_ptr)
{
  return frame_ptr + getFrameLength(frame_ptr) 
                   - getFrameErrorLengthV1(frame_ptr)
                   - getFramePadding(frame_ptr); 
}


/* Return pointer to start of history block in frame 
*/
INLINE_D PTR_ret findFrameHistoryStartV1 (FRAME_ptr frame_ptr)
{
  return frame_ptr + getFrameLength(frame_ptr) 
                   - getFrameErrorLengthV1(frame_ptr)
                   - getFrameHistoryLengthV1(frame_ptr)
                   - getFramePadding(frame_ptr);
}


/* Return pointer to start of alignment block */
INLINE_D PTR_ret findFrameAlignBlockV1 (FRAME_ptr frame_ptr)
{
  return frame_ptr + ALIGN_OFFSET_OF_DWORD;
}

INLINE_D VALUE_ret setFramePaddingV1 (FRAME_ptr frame_ptr, UINT numPadDwords)
{
  setBitsMACRO(frame_ptr,
               PADDING_OFFSET_OF_DWORD,
               PADDING_OFFSET_IN_DWORD,
               PADDING_MASK,
               numPadDwords            );
  return TRUE;
}


/*
**  Adjust the size of the data block in the frame header.
**
**    Note: this routine is only intended to change an explicit
**          size that appears in a header. For V1 frames no such
**          field is present so we simply return the current data
**          block size.
*/

INLINE_D UINT adjustFrameDataLengthV1 (FRAME_ptr frame_ptr, UINT newDwords)
{
  UINT currentDataLength = getFrameDataLengthV1 (frame_ptr);

  if (currentDataLength == valueFailure)
    return valueFailure;
  else
    return currentDataLength + newDwords; 
}

INLINE_D LOGIC_ret setDataTypeV1(FRAME_ptr frame_ptr, UINT dataType)
{
#ifndef DCM_CHECK
 if (dataType > maxByteValue) return FALSE;
#endif
 setBitsMACRO(frame_ptr,
              DATA_TYPE_OFFSET_OF_DWORD,
	      DATA_TYPE_OFFSET_IN_DWORD,
              DATA_TYPE_MASK,
              dataType                  );
 return TRUE;
}


INLINE_D LOGIC_ret setFrameTypeV1(FRAME_ptr frame_ptr, UINT frameType)
{
#ifndef DCM_CHECK
  if (frameType > maxByteValue) return FALSE;
#endif
  setBitsMACRO(frame_ptr,
               FRAME_TYPE_OFFSET_OF_DWORD,
               FRAME_TYPE_OFFSET_IN_DWORD,
               FRAME_TYPE_MASK,
               frameType                  );
  return TRUE;
}


INLINE_D LOGIC_ret setSourceIdV1(FRAME_ptr frame_ptr, UINT sourceId)
{
#ifndef DCM_CHECK
  if (sourceId > maxSwordValue) return FALSE;
#endif
  setBitsMACRO(frame_ptr,
               SOURCE_ID_OFFSET_OF_DWORD,
               SOURCE_ID_OFFSET_IN_DWORD,
               SOURCE_ID_MASK,
               sourceId                   );
  return TRUE;
}


INLINE_D LOGIC_ret setFrameHistoryLengthV1(FRAME_ptr frame_ptr, UINT historyLength)
{
#ifndef DCM_CHECK
  if (historyLength > maxSwordValue) return FALSE;
#endif
  setBitsMACRO(frame_ptr,
               HISTORY_LENGTH_OFFSET_OF_DWORD,
               HISTORY_LENGTH_OFFSET_IN_DWORD,
               HISTORY_LENGTH_MASK,
               historyLength                   );
  return TRUE;
}


INLINE_D LOGIC_ret setFrameErrorLengthV1(FRAME_ptr frame_ptr, UINT errorLength)
{
#ifndef DCM_CHECK
  if (errorLength > maxSwordValue) return FALSE;
#endif
  setBitsMACRO(frame_ptr,
               ERROR_LENGTH_OFFSET_OF_DWORD,
               ERROR_LENGTH_OFFSET_IN_DWORD,
               ERROR_LENGTH_MASK,
               errorLength                    );
  return TRUE;
}


INLINE_D LOGIC_ret setFrameAlignLengthV1(FRAME_ptr frame_ptr, UINT alignLength)
{
#ifndef DCM_CHECK
  if (alignLength > maxByteValue) return FALSE;
#endif
  setBitsMACRO(frame_ptr,
               ALIGN_LENGTH_OFFSET_OF_DWORD,
	       ALIGN_LENGTH_OFFSET_IN_DWORD,
               ALIGN_LENGTH_MASK,
               alignLength                   );
  return TRUE;
}

INLINE_D LOGIC_ret setFrameStatusV1(FRAME_ptr frame_ptr, UINT status)
{
#ifndef DCM_CHECK
  if (status > maxSwordValue) return FALSE;
#endif
  setBitsMACRO(frame_ptr,
	       STATUS_OFFSET_OF_DWORD,
	       STATUS_OFFSET_IN_DWORD,
	       STATUS_MASK,
	       status                  );
  return TRUE;
}


INLINE_D VALUE_ret getAlignBlockV1 (FRAME_ptr frame_ptr, PHDWORD* alignDestination, 
			   UINT maxNumDwords)
{
  UINT alignLength = getFrameAlignLength(frame_ptr);
  if (maxNumDwords < alignLength) return valueFailure;
  dwordCopy(alignDestination, (frame_ptr+ALIGN_OFFSET_OF_DWORD), alignLength);
  return alignLength;
}


INLINE_D LOGIC_ret setAlignBlockV1 (FRAME_ptr frame_ptr, PHDWORD* alignSource, 
			   UINT numDwords)
{
  UINT alignLength = getFrameAlignLength(frame_ptr);
  if (numDwords != alignLength) return FALSE;
  dwordCopy((frame_ptr+ALIGN_OFFSET_OF_DWORD), alignSource, alignLength);
  return TRUE;
}

/* get entry in history block by entry index, where indices count
   from THE BOTTOM */

INLINE_D VALUE_ret getHistoryEntryV1 (FRAME_ptr frame_ptr, UINT index)
{
  UINT historyLength = getFrameHistoryLength(frame_ptr);
  if (index>=historyLength) return valueFailure;
  return *(findFrameHistoryStart(frame_ptr)+(historyLength-index-1));
}

INLINE_D VALUE_ret getStageFromHistoryEntryV1(PHDWORD historyEntry)
{
  return (historyEntry&0xf0000000)>>28;
}

INLINE_D VALUE_ret getSourceIndexFromHistoryEntryV1(PHDWORD historyEntry)
{
  return (historyEntry&0x0fff0000)>>16;
}

INLINE_D VALUE_ret getStatusFromHistoryEntryV1(PHDWORD historyEntry)
{
  return historyEntry&0x0000ffff;
}

INLINE_D VALUE_ret getHistoryStageV1(FRAME_ptr frame_ptr, UINT index)
{
  return getStageFromHistoryEntryV1(getHistoryEntryV1(frame_ptr,index));
}

INLINE_D VALUE_ret getHistorySourceIndexV1(FRAME_ptr frame_ptr, UINT index)
{
  return getSourceIndexFromHistoryEntryV1(getHistoryEntryV1(frame_ptr,index));
}

INLINE_D VALUE_ret getHistoryStatusV1(FRAME_ptr frame_ptr, UINT index)
{
  return getStatusFromHistoryEntryV1(getHistoryEntryV1(frame_ptr,index));
}

INLINE_D PTR_ret findNextErrorV1(FRAME_ptr frame_ptr, FRAME_ptr thisError)
{
  if (thisError < findFrameErrorStart(frame_ptr) || 
      thisError > findFrameEnd(frame_ptr) ) return ptrFailure;
  PHDWORD errorType = *thisError;
  FRAME_ptr temp_ptr;
  switch (errorType)
  {
  case 1111573574: temp_ptr = thisError+3; break;             // "BADF" : bad frame
  case 1296651078: temp_ptr = thisError+2; break;              // "MISF" : missing frame
  case 1111573584: temp_ptr = thisError+4; break;               // "BADP" : bad packet
  case 1296651088: temp_ptr = thisError+3; break;                // "MISP" : missing packet
  case 1430537542: temp_ptr = thisError+(*(thisError+1)); break; // "UDEF" : undefined error
          default: return ptrFailure;
  }
  if (temp_ptr > findFrameEnd(frame_ptr)) return ptrFailure;
  return temp_ptr;
}

#ifdef __cplusplus
}
#endif

#endif /* end of ifndef _CFRAMEV1_ */












