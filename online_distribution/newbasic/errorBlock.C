/* 
** errorBlock.C
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:12 $ 
** 
** $Log: errorBlock.C,v $
** Revision 1.1.1.1  2000/07/21 01:51:12  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:02:01  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
#include "errorBlock.h"

VALUE_ret calcNumErrorsV1 (UINT numErrorDwords)
{	
  if (numErrorDwords % errorEntryV1Length ) {
    return valueFailure;
  }
  else return numErrorDwords/errorEntryV1Length;
}

void endianSwapErrorV1 (ERRORENTRYV1_ptr outError, ERRORENTRYV1_ptr inError)
{
  /*
  **  For now do nothing
  */
  *outError = *inError;
}
