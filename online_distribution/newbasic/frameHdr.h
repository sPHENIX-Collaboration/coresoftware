/* 
** frameHdr.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:13 $ 
** 
** $Log: frameHdr.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:13  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.4  1998/12/11 22:01:30  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**
*/

#ifndef _FRAMEHDR_
#define _FRAMEHDR_

#include "framePublic.h"
#include "frameHdrV1.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FRAME_LENGTH_OFFSET_OF_DWORD 0

#define FRAME_MARK_OFFSET_OF_DWORD 1

#define FRAME_HDR_VERSION_OFFSET_OF_DWORD 2
#define FRAME_HDR_VERSION_OFFSET_IN_DWORD 24
#define FRAME_HDR_VERSION_NUM_BITS 8
#define FRAME_HDR_VERSION_MASK 0xff000000

#define FRAME_HDR_LENGTH_OFFSET_OF_DWORD 2
#define FRAME_HDR_LENGTH_OFFSET_IN_DWORD 16
#define FRAME_HDR_LENGTH_NUM_BITS 8
#define FRAME_HDR_LENGTH_MASK 0x00ff0000

#define FRAME_LENGTH_OFFSET 0
#define FRAME_MARK_OFFSET 1

#define NUM_FRAME_VERSIONS 2        /* we have a dummy version 0 */
#define CURRENT_FRAME_HDR_VERSION 1
    
  CONSTANT UINT numFrameVersions = NUM_FRAME_VERSIONS;	
  CONSTANT UINT frameMarkV[NUM_FRAME_VERSIONS] = {0, V1_FRAMEMARK};
  CONSTANT UINT frameHdrLengthV[NUM_FRAME_VERSIONS] = {0, V1_FRAMEHDR_LENGTH};
  
  CONSTANT UINT frameLengthOffset = FRAME_LENGTH_OFFSET;
  CONSTANT UINT frameMarkOffset = FRAME_MARK_OFFSET;
  
  /*
  ** Define "current" header type for use in creating new frames
  */
  CONSTANT UINT currentFrameHdrVersion = CURRENT_FRAME_HDR_VERSION;
  CONSTANT UINT currentFrameHdrLength = V1_FRAMEHDR_LENGTH;
  CONSTANT UINT currentFrameMark = V1_FRAMEMARK;
  CONSTANT UINT currentAlignLength = V1_ALIGN_LENGTH;
  CONSTANT UINT currentFrameQuantum = V1_FRAME_QUANTUM;    

  CONSTANT UINT maxFrameHdrLength = V1_FRAMEHDR_LENGTH;
  
  /* *** temporarily in framePublic.h ***    typedef ALIGNBLKV1 ALIGNBLK; */
  
#ifdef __cplusplus
}
#endif

#endif




