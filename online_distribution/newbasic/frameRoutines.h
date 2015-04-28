/* 
** frameRoutines.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:14 $ 
** 
** $Log: frameRoutines.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:14  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:43  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  Prototype definitions for "public" routines in frameRoutines.C
*/

#ifndef _FRAME_ROUTINES_
#define _FRAME_ROUTINES_

#include "phenixOnline.h"
#include "framePublic.h"

/*
**  Use C linkage for below structures
*/
#ifdef __cplusplus
extern "C" {
#endif

/*
**  Routines that perform more general and sophistcated operations on frames
**    than the member-function-like routines defined in frames.h
*/

  VALUE_ret makeFrameHdr (PHDWORD*, UINT, UINT, UINT, UINT);
  
  VALUE_ret storeFrameData (FRAME_ptr, UINT, PHDWORD*, UINT);
  VALUE_ret storeFrameHistory (FRAME_ptr, UINT, PHDWORD*, UINT);
  
  VALUE_ret extendFrameData (FRAME_ptr, UINT, UINT);
  VALUE_ret extendFrameDataNopad (FRAME_ptr, UINT, UINT);
  VALUE_ret extendFrameHistory (FRAME_ptr, UINT, UINT);
  VALUE_ret extendFrameHistoryNopad (FRAME_ptr, UINT, UINT);
  
#ifdef __cplusplus
}
#endif

#endif
	
