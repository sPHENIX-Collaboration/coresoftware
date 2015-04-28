/* 
** frameRoutines.C
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:13 $ 
** 
** $Log: frameRoutines.C,v $
** Revision 1.1.1.1  2000/07/21 01:51:13  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.4  1999/10/07 19:53:35  steinber
** changed storeFrameHistory to not use padding.  breaks EvB code otherwise
**
** Revision 1.3  1998/12/11 22:02:05  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
** makeFrameHdr
**
**		Routine to make a new frame header in a buffer pointed 
**		to by "newFramePtr". The header is created with "empty"
**		data, history, and error blocks.
*/

#include "phenixOnline.h"
#include "frameRoutines.h"
#include "framePublic.h"
#include "Cframe.h"

/*
** storeFrameData
**
**		Routine to store data in new frame. 
**
*/
VALUE_ret storeFrameData (PHDWORD* frame_ptr, UINT maxFrameLen, 
			  PHDWORD* frameData, UINT dataDwords)
{
  PHDWORD* output_ptr;
  UINT finalLength;

  /*
  **  Make sure we're pointing to a real frame and 
  **	that the header version is the current one
  **  (we don't write old frame headers)
  */
  if (currentFrameHdr(frame_ptr))
    return valueFailure;

  /*
  **	Only store data in empty frames
  */
  if (!emptyFrame (frame_ptr)) return valueFailure;
		
  /*
  **  Strip off any padding in the frame
  */
  removeFramePadding (frame_ptr);

  /*
  **  Find out where to write data
  */
  output_ptr = findFrameDataStart(frame_ptr); 
  if (output_ptr == ptrFailure)
    return valueFailure;

  /*
  **  Now extend the frame to hold the data
  */
  finalLength = extendFrameData (frame_ptr, maxFrameLen, dataDwords);
  if (finalLength == valueFailure)
    return valueFailure;

  /*
  **	Now transfer the data into the frame.
  */
  dwordCopy (output_ptr, frameData, dataDwords);

  return finalLength;
}

/*
**  Store a history block in a frame that doesn't already have one
*/
VALUE_ret storeFrameHistory (PHDWORD* frame_ptr, UINT maxFrameLen, 
			     PHDWORD* frameHistory, UINT historyDwords)
{
  PHDWORD* history_ptr;
  UINT finalLength;

  /*
  **  Make sure we're pointing to a real frame and 
  **	that the header version is the current one
  **  (we don't write old frame headers)
  */
  if (!validFrameHdr (frame_ptr) || !currentFrameHdr (frame_ptr))
    return valueFailure;

  /*
  **  Only store history in empty history block
  */
  if (getFrameHistoryLength (frame_ptr) != 0) return valueFailure;
		
  /*
  **  Strip off any frame padding
  */
  removeFramePadding(frame_ptr);

  /*
  **  Get pointer to history block
  */
  history_ptr = findFrameHistoryStart (frame_ptr);
  if (history_ptr == ptrFailure)
    return valueFailure;

  /*
  **  Extend the frame to hold the history data
  */
  finalLength = extendFrameHistoryNopad (frame_ptr, maxFrameLen, historyDwords);
  if ( finalLength == valueFailure)
    return valueFailure;

  /*
  **  Now copy in history words
  */
  dwordCopy (history_ptr, frameHistory, historyDwords);

  return finalLength;
}

/*
**  Extend the length of the data block in a frame
*/
VALUE_ret extendFrameData (FRAME_ptr frame_ptr, UINT maxFrameLength, UINT dataDwords)
{
  if (adjustFrameDataLength (frame_ptr, dataDwords) != valueFailure) {
    /*
    **  We only need to extend the length of the frame.
    */
    return adjustFrameLength (frame_ptr, maxFrameLength, dataDwords, TRUE);
  }
  else return valueFailure;
}

/*
**  Extend the length of the data block in a frame
*/
VALUE_ret extendFrameDataNopad (FRAME_ptr frame_ptr, UINT maxFrameLength, UINT dataDwords)
{
  if (adjustFrameDataLength (frame_ptr, dataDwords) != valueFailure) {
    /*
    **  We only need to extend the length of the frame.
    */
    return adjustFrameLength (frame_ptr, maxFrameLength, dataDwords, FALSE);
  }
  else return valueFailure;
}


/*
**  Extend the length of the history block in a frame
*/
VALUE_ret extendFrameHistory (FRAME_ptr frame_ptr, UINT maxFrameLength, UINT historyDwords)
{
  if (adjustFrameHistoryLength (frame_ptr, historyDwords) != valueFailure) {
    /*
    **  We only need to extend the length of the frame.
    */
    return adjustFrameLength (frame_ptr, maxFrameLength, historyDwords, TRUE);
  }
  else return valueFailure;
}

/*
**  Extend the length of the history block in a frame
*/
VALUE_ret extendFrameHistoryNopad (FRAME_ptr frame_ptr, UINT maxFrameLength, UINT historyDwords)
{
  if (adjustFrameHistoryLength (frame_ptr, historyDwords) != valueFailure) {
    /*
    **  We only need to extend the length of the frame.
    */
    return adjustFrameLength (frame_ptr, maxFrameLength, historyDwords, FALSE);
  }
  else return valueFailure;
}










