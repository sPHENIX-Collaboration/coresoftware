/* 
** Cframe.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:10 $ 
** 
** $Log: Cframe.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:10  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.11  1999/05/17 20:47:24  markacs
** (stephen) last bit of addition of findNextError
**
** Revision 1.10  1999/04/13 18:27:36  markacs
** (stephen) put history block access in
**
** Revision 1.9  1998/12/17 15:23:01  phoncs
** (stephen markacs) OOPS, had a C++ type comment.
**
** Revision 1.8  1998/12/17 15:18:16  phoncs
** (stephen markacs) more bounds checking in Cpacket.C (one instance) and removal of debugging comments from checkFrame.C
**
** Revision 1.7  1998/12/16 15:40:51  markacs
** (stephen markacs) changes to rawFormat code to stop using bool,true,false
**
** Revision 1.6  1998/12/11 22:01:09  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
** Cframe.h
**
** This file contains functions to access header fields.
**
** The functions are named as follows:
**
**    get: return the value in the field
**   find: return a pointer to the start of a block
**    set: set the value in the field to a given value
** adjust: change some parameter (such as length) by some 
**           given value
**  check: return true if something is valid
**
** These are only general guides to what the functions do.
**   There are small variations, such as getAlignBlock,
**   which is passed a pointer to an area of memory to
**   which it copies the alignment block.
**
**
** The fields are accessed through direct reading and
**   writing of bits using offsets which are #defined.
**   Fields which are a full PHDWORD (4 bytes) long are 
**   defined by a single value, the offset of the PHDWORD
**   that comprises them.  The smaller fields have four
**   #defined values:
**
** 1.The offset of the PHDWORD the field is in
**    (which word it is in)
** 2.The offset of the field in that PHDWORD
**    (where it starts in that word)
** 3.The number of bits in the field
**    (how big it is)
** 4.The mask used to pull out the field from its PHDWORD
**    (put in by hand for optimization)
**
**
** The layout of this file is as follows:
**
** (1) definition of field parameters (the offsets, etc.)
** (2) includes of phenixOnline.h, etc.
** (3) prototypes of functions defined later in this file
**      and #defines of version-dependent function names to
**      V1 function names (if in DCM mode)
** (4) include of CframeV1.h
** (5) generic functions
** (6) definition of indirection pointers and specific functions 
**      (if not in DCM mode)
**
**
** There are two types of routines in this file, those that
**   work with fields whose configurations (offsets and size)
**   are version-independent and those whose configurations
**   are allowed to vary with version number.  Fields of the
**   first type are operated on directly by routines in this
**   file.  Version-dependent fields are accessed through
**   an array of pointers to functions, those functions lying
**   in their respective version files, such as "CframeV1.h".
**   In DCM code, DCM is #defined, which causes the
**   routines to access the version-dependent fields to be
**   #defined directly to the corresponding function for the
**   current version.  
**
**
** There are a few compilation modes that are implemented in this
**   file through #defining the appropriate name:
**
**         DCM: for compilation to the DCM; causes direct 
**              redefinition of functions names that access 
**              version-dependent fields to the names of the
**              current version's corresponding functions
**   DCM_CHECK: causes certain checks in low-level access routines
**              to be skipped for optimization
** STRONGCHECK: causes some extra-careful checks to be added
**
**
** Nomenclature note:  The "end" of a frame, a data block, etc.,
**   is taken to be the address of the LAST dword in the object,
**   NOT one past the last dword in the object.
**
**
*/

#include <stdio.h>  /* DEBUG */

#ifndef _FRAMES_
#define _FRAMES_

#ifdef __cplusplus
extern "C" {
#endif

#include "phenixOnline.h"
#include "framePublic.h"
#include "formatError.h"
#include "frameHdr.h"

  /* Version-INDEPENDENT functions' prototypes */

  INLINE_P VALUE_ret getFrameHdrLength (FRAME_ptr);   
  INLINE_P VALUE_ret getFrameHdrVersion (FRAME_ptr);
  INLINE_P VALUE_ret getFrameLength (FRAME_ptr); 
  INLINE_P VALUE_ret getFrameMark (FRAME_ptr);

  INLINE_P LOGIC_ret setFrameHdrLength(FRAME_ptr, UINT);    
  INLINE_P LOGIC_ret setFrameHdrVersion(FRAME_ptr, UINT);
  INLINE_P LOGIC_ret setFrameLength(FRAME_ptr, UINT);
  INLINE_P LOGIC_ret setFrameMark(FRAME_ptr, UINT);

  VALUE_ret checkFrameHdrVersion (FRAME_ptr);
    
  INLINE_P LOGIC_ret validFrameHdr (FRAME_ptr);
  INLINE_P LOGIC_ret currentFrameHdr (FRAME_ptr);
  INLINE_P LOGIC_ret validFrameMark (FRAME_ptr);
  INLINE_P LOGIC_ret emptyFrame(FRAME_ptr);

  VALUE_ret checkFrameEndianism(FRAME_ptr frame_ptr);
  INLINE_P void byteSwapFrame (FRAME_ptr frame_ptr);
   
  PTR_ret findFrameEnd (FRAME_ptr);
  PTR_ret findFrameDataStart (FRAME_ptr);
  PTR_ret findFrameDataEnd (FRAME_ptr);    

  VALUE_ret adjustFrameLength (FRAME_ptr, UINT, UINT, LOGIC_ret);
  VALUE_ret removeFramePadding (FRAME_ptr);

  /* Version-DEPENDENT functions' macro redirections or prototypes */

#ifdef DCM

#define makeFrameHdr makeFrameHdrV1
#define getFrameDataLength getFrameDataLengthV1
#define getFrameHistoryLength getFrameHistoryLengthV1
#define getFrameErrorLength getFrameErrorLengthV1
#define getFrameAlignLength getFrameAlignLengthV1
#define getFrameSourceId getFrameSourceIdV1
#define getFrameDataType getFrameDataTypeV1
#define getFrameType getFrameTypeV1
#define getFrameStatus getFrameStatusV1
#define getFramePadding getFramePaddingV1
#define setFramePadding setFramePaddingV1
#define adjustFrameDataLength adjustFrameDataLengthV1
#define adjustFrameHistoryLength adjustFrameHistoryLengthV1
#define adjustframeErrorLength adjustFrameErrorlengthV1
#define orFrameStatus orFrameStatusV1
#define setDataType setDataTypeV1
#define setFrameType setFrameTypeV1
#define setSourceId setSourceIdV1
#define setFrameHistoryLength setFrameHistoryLengthV1
#define setFrameErrorLength setFrameErrorLengthV1
#define setFrameAlignLength setFrameAlignLengthV1
#define setFrameStatus setFrameStatusV1
#define findFrameAlignBlock findFrameAlignBlockV1
#define findFrameHistoryStart findFrameHistoryStartV1
#define findFrameErrorStart findFrameErrorStartV1
#define getAlignBlock getAlignBlockV1
#define setAlignBlock setAlignBlockV1
#define getHistoryEntry getHistoryEntryV1
#define getHistoryStage getHistoryStageV1
#define getHistorySourceIndex getHistorySourceIndexV1
#define getHistoryStatus getHistoryStatusV1
#define findNextError findNextErrorV1

#else
  VALUE_ret makeFrameHdr (FRAME_ptr, UINT, UINT, UINT, UINT);
  VALUE_ret getFrameDataLength (FRAME_ptr);
  VALUE_ret getFrameHistoryLength (FRAME_ptr);
  VALUE_ret getFrameErrorLength (FRAME_ptr);
  VALUE_ret getFrameAlignLength (FRAME_ptr);
    
  VALUE_ret getFrameSourceId (FRAME_ptr);
  VALUE_ret getFrameDataType (FRAME_ptr);
  VALUE_ret getFrameType (FRAME_ptr);
  VALUE_ret getFrameStatus (FRAME_ptr);
  VALUE_ret getFramePadding (FRAME_ptr);

  VALUE_ret setFramePadding (FRAME_ptr, UINT);
  VALUE_ret adjustFrameDataLength (FRAME_ptr, UINT);
  VALUE_ret adjustFrameHistoryLength (FRAME_ptr, UINT);
  VALUE_ret adjustFrameErrorLength (FRAME_ptr, UINT);
  VALUE_ret orFrameStatus (FRAME_ptr, UINT);

  LOGIC_ret setDataType(FRAME_ptr, UINT);
  LOGIC_ret setFrameType(FRAME_ptr, UINT);
  LOGIC_ret setSourceId(FRAME_ptr, UINT);
  LOGIC_ret setFrameHistoryLength(FRAME_ptr, UINT);
  LOGIC_ret setFrameErrorLength(FRAME_ptr, UINT);
  LOGIC_ret setFrameAlignLength(FRAME_ptr, UINT);
  LOGIC_ret setFrameStatus(FRAME_ptr, UINT);

  PTR_ret findFrameAlignBlock (FRAME_ptr);
  PTR_ret findFrameErrorStart (FRAME_ptr);
  PTR_ret findFrameHistoryStart (FRAME_ptr);

  VALUE_ret getAlignBlock (FRAME_ptr, PHDWORD*, UINT);
  LOGIC_ret setAlignBlock (FRAME_ptr, PHDWORD*, UINT);

  VALUE_ret getHistoryEntry(FRAME_ptr, UINT);
  VALUE_ret getHistoryStage(FRAME_ptr, UINT);
  VALUE_ret getHistorySourceIndex(FRAME_ptr, UINT);
  VALUE_ret getHistoryStatus(FRAME_ptr, UINT);
  PTR_ret findNextError(FRAME_ptr, FRAME_ptr);

#endif
#ifdef __cplusplus
}                  
/* end of extern "C" */
#endif


/* Here begins the inlined code for the version-INDEPENDENT functions */
/* ================================================================== */

/*  Return the length of the frame header
*/
INLINE_D VALUE_ret getFrameHdrLength (FRAME_ptr frame_ptr) 
{
  UINT hdrLength = getBitsMACRO(frame_ptr,
				FRAME_HDR_LENGTH_OFFSET_OF_DWORD,
				FRAME_HDR_LENGTH_OFFSET_IN_DWORD,
				FRAME_HDR_LENGTH_MASK        );
#ifdef STRONGCHECK
  if (hdrLength == 0 || hdrLength > getFrameLength(frame_ptr)) {
    setFrameError(FORMAT_ERR_INVALID_HDRLENGTH, frame_ptr, length);
    return valueFailure;
  }
#endif
  return hdrLength;
}

/*  Return the version number of the frame. If enabled, 
**   check for invalid version.
*/
INLINE_D VALUE_ret getFrameHdrVersion (FRAME_ptr frame_ptr) 
{
  UINT version = getBitsMACRO(frame_ptr,
			      FRAME_HDR_VERSION_OFFSET_OF_DWORD,
			      FRAME_HDR_VERSION_OFFSET_IN_DWORD,
			      FRAME_HDR_VERSION_MASK        );
	/*#ifdef STRONGCHECK */
  if ((version >= numFrameVersions)||(version==0)) 
  {
    setFrameError(FORMAT_ERR_INVALID_HDRVERSION, frame_ptr, version);
    return valueFailure; 
  }  
	/*#endif*/
  return version;
}


/*  Return the length of the frame. If enabled check for overflow.
*/
INLINE_D VALUE_ret getFrameLength (FRAME_ptr frame_ptr) 
{
  PHDWORD length = getWordMACRO(frame_ptr, FRAME_LENGTH_OFFSET_OF_DWORD);
  
#ifdef STRONGCHECK
  if ((length&(1<<32))!=0) {
    setFrameError(FORMAT_ERR_LENGTH_OVERFLOW, frame_ptr, length);
    return valueFailure;
  }
#endif
  return length;
}

/*  Return the frame marker.
*/
INLINE_D VALUE_ret getFrameMark (FRAME_ptr frame_ptr) 
{
  PHDWORD frameMark = getWordMACRO(frame_ptr, FRAME_MARK_OFFSET_OF_DWORD);
  return frameMark;
} 

/*  Store the length of the header in the header itself.
*/ 
INLINE_D LOGIC_ret setFrameHdrLength (FRAME_ptr frame_ptr, UINT hdrLength) 
{
  setBitsMACRO(frame_ptr, FRAME_HDR_LENGTH_OFFSET_OF_DWORD,
	       FRAME_HDR_LENGTH_OFFSET_IN_DWORD, FRAME_HDR_LENGTH_MASK,
	       hdrLength);
  return TRUE;
}

/*  Store the version of the header in the header itself
*/
INLINE_D LOGIC_ret setFrameHdrVersion (FRAME_ptr frame_ptr, UINT hdrVersion) 
{
  setBitsMACRO(frame_ptr, FRAME_HDR_VERSION_OFFSET_OF_DWORD,
	       FRAME_HDR_VERSION_OFFSET_IN_DWORD, FRAME_HDR_VERSION_MASK,
	       hdrVersion);
  return TRUE; 
}

/*  Store the frame length in the header
*/
INLINE_D LOGIC_ret setFrameLength (FRAME_ptr frame_ptr, UINT length) {
  setWordMACRO(frame_ptr, 
	       FRAME_LENGTH_OFFSET_OF_DWORD,
	       length);
  return TRUE;
}

/*  Store the frame marker in the header
*/
INLINE_D LOGIC_ret setFrameMark (FRAME_ptr frame_ptr, UINT frameMark) 
{
  setWordMACRO(frame_ptr, 
	       FRAME_MARK_OFFSET_OF_DWORD,
	       frameMark);
  return TRUE;
}


/*  Check for a valid frame header.
*/
INLINE_D LOGIC_ret validFrameHdr (FRAME_ptr frame_ptr) 
{
  return (checkFrameHdrVersion(frame_ptr) != valueFailure);
}

/*  Check to see if the frame header is of the current version.
*/
INLINE_D LOGIC_ret currentFrameHdr (FRAME_ptr frame_ptr) 
{
  UINT version = getFrameHdrVersion (frame_ptr);
  return (version == currentFrameHdrVersion);
}

/*  Check to see if the frame has a valid marker.
*/
INLINE_D LOGIC_ret validFrameMark (FRAME_ptr frame_ptr) {
  UINT version = getFrameHdrVersion(frame_ptr);
  if (version != valueFailure)
    return (getFrameMark(frame_ptr) == frameMarkV[version]);
  else return FALSE;
}
  

/* Check to see if the frame is empty.
**   Note that an "empty frame" consists of a header only. (No tailer)
*/               
INLINE_D LOGIC_ret emptyFrame (FRAME_ptr frame_ptr) {
  return (getFrameHdrLength(frame_ptr) == getFrameLength(frame_ptr));
}

/*
**  Byte-swap an entire frame using the (swapped) length stored in the first PHDWORD.
**
**  This routine assumes that the endianism has been determined to be wrong AND 
**    that the byte-swapped header satisfies a set of minimal validity tests.
*/
INLINE_D void byteSwapFrame (FRAME_ptr frame_ptr)
{
  PHDWORD correctLength = singleDwordByteSwap(*frame_ptr);
  dwordByteSwap(frame_ptr, frame_ptr, correctLength);
}

/* ======================================================== */
/* Here ends the code for the version-INDEPENDENT functions */


#endif /* end of ifdef _FRAMES_ */

