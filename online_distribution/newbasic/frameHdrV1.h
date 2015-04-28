/* 
** frameHdrV1.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:13 $ 
** 
** $Log: frameHdrV1.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:13  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:33  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  Version 1 frame header definitions
**
*/

#ifndef _FRAMEHDRV1_
#define _FRAMEHDRV1_

#include "phenixOnline.h"
#include "frameV1Public.h"

/*
**  Use C linkage for below structures
*/
#ifdef __cplusplus
extern "C" {
#endif


#define STATUS_OFFSET_OF_DWORD 2
#define STATUS_OFFSET_IN_DWORD 0
#define STATUS_NUM_BITS 16
#define STATUS_MASK 0x0000ffff

#define SEQUENCE_NUMBER_OFFSET_OF_DWORD 3
#define SEQUENCE_NUMBER_OFFSET_IN_DWORD 24
#define SEQUENCE_NUMBER_NUM_BITS 8
#define SEQUENCE_NUMBER_MASK 0xff000000

#define SEQUENCE_CODE_OFFSET_OF_DWORD 3
#define SEQUENCE_CODE_OFFSET_IN_DWORD 16
#define SEQUENCE_CODE_NUM_BITS 8
#define SEQUENCE_CODE_MASK 0x00ff0000

#define SOURCE_ID_OFFSET_OF_DWORD 3
#define SOURCE_ID_OFFSET_IN_DWORD 0 
#define SOURCE_ID_NUM_BITS 16
#define SOURCE_ID_MASK 0x0000ffff

#define DATA_TYPE_OFFSET_OF_DWORD 4
#define DATA_TYPE_OFFSET_IN_DWORD 24
#define DATA_TYPE_NUM_BITS 8
#define DATA_TYPE_MASK 0xff000000

#define FRAME_TYPE_OFFSET_OF_DWORD 4
#define FRAME_TYPE_OFFSET_IN_DWORD 16 
#define FRAME_TYPE_NUM_BITS 8
#define FRAME_TYPE_MASK 0x00ff0000

#define ERROR_LENGTH_OFFSET_OF_DWORD 4
#define ERROR_LENGTH_OFFSET_IN_DWORD 0
#define ERROR_LENGTH_NUM_BITS 16
#define ERROR_LENGTH_MASK 0x0000ffff

#define HISTORY_LENGTH_OFFSET_OF_DWORD 5
#define HISTORY_LENGTH_OFFSET_IN_DWORD 16
#define HISTORY_LENGTH_NUM_BITS 16
#define HISTORY_LENGTH_MASK 0xffff0000

#define ALIGN_LENGTH_OFFSET_OF_DWORD 5
#define ALIGN_LENGTH_OFFSET_IN_DWORD 8
#define ALIGN_LENGTH_NUM_BITS 8
#define ALIGN_LENGTH_MASK 0x0000ff00

#define PADDING_OFFSET_OF_DWORD 5
#define PADDING_OFFSET_IN_DWORD 0
#define PADDING_NUM_BITS 8
#define PADDING_MASK 0x000000ff

#define ALIGN_OFFSET_OF_DWORD 6

/*
**  Parameters specific to V1 frames
*/
#define V1_FRAMEHDR_VERSION 1
#define V1_FRAMEHDR_LENGTH  8
#define V1_FRAMEMARK 0xFFFFFF00
#define V1_ALIGN_LENGTH 2
#define V1_FRAME_QUANTUM 2

/*
** Some constants for V1 frame definition
*/

CONSTANT UINT v1FrameHdrVersion = V1_FRAMEHDR_VERSION;
CONSTANT PHDWORD v1FrameHdrLength = V1_FRAMEHDR_LENGTH;
CONSTANT PHDWORD v1FrameMark = V1_FRAMEMARK;
CONSTANT UINT v1AlignLength = V1_ALIGN_LENGTH;
CONSTANT UINT v1FrameQuantum = V1_FRAME_QUANTUM;

#ifdef __cplusplus
} 
/* end of extern "C"  */
#endif

#endif 
/* end of ifndef _FRAMEHDRV1_ */








