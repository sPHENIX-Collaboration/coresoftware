/* 
** errorBlock.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:12 $ 
** 
** $Log: errorBlock.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:12  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:17  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  formatError.h
**
**		Definitions of various error-handling macros, prototypes, variables etc
*/

#ifndef _ERRORBLOCK_
#define _ERRORBLOCK_

#include "phenixOnline.h"

/*
**  Use C linkage
*/
#ifdef __cplusplus
extern "C" {
#endif


  struct errorEntryV1 {
    BYTE	severity;
    BYTE	deviceType;
    SWORD	deviceId;
    
    SWORD	errorCode;
    SWORD	detectCode;
    
    PHDWORD	addData[2];
  };
  
  typedef struct errorEntryV1 ERRORENTRYV1, *ERRORENTRYV1_ptr;
  
#define ERROR_ENTRY_V1_LENGTH sizeof(ERRORENTRYV1)/4
  CONSTANT UINT errorEntryV1Length = ERROR_ENTRY_V1_LENGTH;
  
  enum daqErrorCodes {
    InvalidFrameHeader = 1,
    InvalidSourceId = 2,
    InvalidFrameType = 3
  };
  
  VALUE_ret calcNumErrorsV1 (UINT);
  void endianSwapErrorV1 (ERRORENTRYV1_ptr, ERRORENTRYV1_ptr);
  
#ifdef __cplusplus
}
#endif

#endif
