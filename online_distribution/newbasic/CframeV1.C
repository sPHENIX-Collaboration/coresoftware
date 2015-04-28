/* 
** CframeV1.C
** 
** Author: $Author: phnxbld $  
**   Date: $Date: 2009/09/19 14:34:27 $ 
** 
** $Log: CframeV1.C,v $
** Revision 1.2  2009/09/19 14:34:27  phnxbld
** fix compiler warning
**
** Revision 1.1.1.1  2000/07/21 01:51:10  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:58  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
#include "CframeV1.h"

/*
**  makeFrameHdrV1 makes a Version 1 frame header.  It takes the address pointed to by
**    the frame_ptr passed to it, clears the right amount of space for the header
**    (sets all bytes to 0), sets the frameMark, hdrVersion, hdrLength, and
**    frameLength fields to the appropriate values for an empty Version 1 header,  
**    sets the dataType, frameType and sourceId fields to the values passes to it,
**    and copies the alignment block from the location pointed to by the 
**    alignment pointer passed to it.
*/
VALUE_ret makeFrameHdrV1 (PHDWORD* frame_ptr, UINT maxFrameLen, UINT dataType, 
			  UINT frameType, UINT sourceId) 
{
  if (maxFrameLen < currentFrameHdrLength) return valueFailure;

  dwordClear (frame_ptr, currentFrameHdrLength);

  setFrameMark(frame_ptr,currentFrameMark);
  setFrameHdrVersion(frame_ptr,currentFrameHdrVersion);
  setFrameHdrLength(frame_ptr,currentFrameHdrLength); 

  setDataType(frame_ptr, dataType);
  setFrameType(frame_ptr, frameType);
  setSourceId(frame_ptr, sourceId);
  setFrameLength(frame_ptr, V1_FRAMEHDR_LENGTH);
  setFrameAlignLength(frame_ptr, V1_ALIGN_LENGTH);

  return 0; 
}

VALUE_ret getFrameDataLengthV1 (FRAME_ptr frame_ptr) 
{
  PHDWORD dataLength =   getFrameLength(frame_ptr) 
                     - getFrameHdrLength(frame_ptr)
                     - getFrameErrorLengthV1(frame_ptr)
                     - getFrameHistoryLengthV1(frame_ptr)
                     - getFramePaddingV1(frame_ptr)       ;

  if (dataLength > getFrameLength(frame_ptr) )
  {
    setFrameError(FORMAT_ERR_INVALID_HEADER, frame_ptr, 0);
    return valueFailure;
  }
  else 
  {
    setFrameSuccess();
    return dataLength;
  }
}

/* set bits in the frame status word */
VALUE_ret orFrameStatusV1 (FRAME_ptr frame_ptr, UINT statusBits)
{
  UINT status;
  /*
  **  Make sure intended bits fit within field length
  */
  if ((statusBits & ((1<<STATUS_NUM_BITS)-1)) == statusBits)
    {
      status = getBitsMACRO(frame_ptr,              STATUS_OFFSET_OF_DWORD,
                            STATUS_OFFSET_IN_DWORD, STATUS_MASK        );
      status|=statusBits;
      setBitsMACRO(frame_ptr, 
                   STATUS_OFFSET_OF_DWORD,
                   STATUS_OFFSET_IN_DWORD,
                   STATUS_MASK,
                   status                 );
      return status;
    }
  else 
    {
      setUserError (FORMAT_ERR_INVALID_DATA, statusBits);
      return valueFailure;
    }
}



VALUE_ret adjustFrameHistoryLengthV1 (FRAME_ptr frame_ptr, UINT newDwords)
{
  UINT newLength = getFrameHistoryLengthV1(frame_ptr) + newDwords;

  /*
  ** Make sure it fits in the current field size
  */
  if (newLength > maxSwordValue) 
  {
    setFrameError(FORMAT_ERR_HISTORY_OVERFLOW, frame_ptr, newDwords);
    return valueFailure;
  }
  else	
  {
    setBitsMACRO(frame_ptr, 
                 HISTORY_LENGTH_OFFSET_OF_DWORD,
                 HISTORY_LENGTH_OFFSET_IN_DWORD,
                 HISTORY_LENGTH_MASK,
                 newLength                      );
    setFrameSuccess();    
    return newLength;
  }
}



VALUE_ret adjustFrameErrorLengthV1 (FRAME_ptr frame_ptr, UINT newDwords)
{
  UINT newLength = getFrameErrorLengthV1(frame_ptr) + newDwords;

  /*
  ** Make sure it fits in the current field size
  */
  if (newLength > maxSwordValue) 
  {
    setFrameError(FORMAT_ERR_ERROR_OVERFLOW, frame_ptr, newDwords);
    return valueFailure;
  }
  else	
  {
    setBitsMACRO(frame_ptr, 
                 ERROR_LENGTH_OFFSET_OF_DWORD,
                 ERROR_LENGTH_OFFSET_IN_DWORD,
                 ERROR_LENGTH_MASK,
                 newLength                      );
    setFrameSuccess();    
    return newLength;
  }
}



