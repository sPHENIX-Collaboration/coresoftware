/* 
** frameV1Public.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:14 $ 
** 
** $Log: frameV1Public.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:14  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:44  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  Public definitions for V1 frames.
**
*/

#ifndef _FRAMEV1PUBLIC_
#define _FRAMEV1PUBLIC_

#include "phenixOnline.h"

/*
**  Use C linkage for below structures
*/
#ifdef __cplusplus
extern "C" {
#endif

/*
** alignment block definitions (still not fully specified)
*/
  struct dcmAlignBlk {
    SWORD   timeStamp;
    SWORD   granuleEvtcnt; /* (lowest 16 bits) */
    PHDWORD   partitionVec;
  }; 
  typedef struct dcmAlignBlk DCMALIGNBLK;
  
  struct dcbAlignBlk {
    SWORD   timeStamp;                                          
    SWORD   granuleEvtcnt; /* (lowest 16 bits) */
    PHDWORD   partitionVec;
  }; 
  typedef struct dcbAlignBlk DCBALIGNBLK; 

  struct sebAlignBlk {
    PHDWORD   globalEventNum;
    PHDWORD   partitionVec;
  }; 
  typedef struct sebAlignBlk SEBALIGNBLK;
  
  struct atpAlignBlk {
    PHDWORD   globalEventNum;
    PHDWORD   partitionVec;
  }; 
  typedef struct atpAlignBlk ATPALIGNBLK; 
	
  /*
  **  For now handle the alignment block as a union.
  */
  typedef union alignBlkV1 {
    DCMALIGNBLK dcm;
    DCBALIGNBLK dcb;
    SEBALIGNBLK seb;
    ATPALIGNBLK atp;
    /*	ONCSALIGNBLK oncs; */
  } ALIGNBLKV1;
  
  /*
  **  We point to V1 frames using PHDWORD pointers
  */
  typedef PHDWORD* V1FRAME_ptr;


#ifdef __cplusplus
} 
/* end of extern "C"  */
#endif

#endif 
/* end of ifndef _FRAMESV1_ */
