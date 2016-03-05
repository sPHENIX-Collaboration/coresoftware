/* 
** dataBlock.h
** 
** Author: $Author: phoncs $  
**   Date: $Date: 2010/09/21 19:37:49 $ 
** 
** $Log: dataBlock.h,v $
** Revision 1.2  2010/09/21 19:37:49  phoncs
** DLW: change name of DWORD to PHDWORD
**
** Revision 1.1.1.1  2000/07/21 01:51:11  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.6  2000/03/22 22:58:37  phoncs
** mlp -- added DWORD
**
** Revision 1.5  1998/12/11 22:01:15  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**      dataBlock.h
**
**      This routine provides acces to information in the data descriptor.
**      The data descriptor has a format which changes depending on the 
**      data structure employed.  So there are as many routines for 
**      creating data descriptors as there are data structures.  
**      Currently, the types of structures are:
**      
**      - Unstructured: sets three fields, the hitFormat, dataPadding
**        and wordSize, in a one dword descriptor.  
**
*/

#ifndef _DATABLOCK_ 
#define _DATABLOCK_

#include "phenixOnline.h"
#include "packetPublic.h"
#include "Cpacket.h"
#include "formatError.h"

/*
**      Use C linkage.
*/

#ifdef __cplusplus
extern "C" {
#endif
 

  /*
  **  These are the offsets for unstructured data.  
  **  For unstructured data the descriptor has a length
  **  of one dword.
  */

#define UNSTRUCT_DATA_PADDING_OFFSET_OF_DWORD 0
#define UNSTRUCT_DATA_PADDING_OFFSET_IN_DWORD 24 
#define UNSTRUCT_DATA_PADDING_NUM_BITS 8
#define UNSTRUCT_DATA_PADDING_MASK 0xff000000

#define UNSTRUCT_WORD_SIZE_OFFSET_OF_DWORD 0
#define UNSTRUCT_WORD_SIZE_OFFSET_IN_DWORD 16 
#define UNSTRUCT_WORD_SIZE_NUM_BITS 8
#define UNSTRUCT_WORD_SIZE_MASK 0x00ff0000

#define UNSTRUCT_HIT_FORMAT_OFFSET_OF_DWORD 0
#define UNSTRUCT_HIT_FORMAT_OFFSET_IN_DWORD 0
#define UNSTRUCT_HIT_FORMAT_NUM_BITS 16
#define UNSTRUCT_HIT_FORMAT_MASK 0x0000ffff
 
CONSTANT PHDWORD UnstructDataDescrLength = 1;



  /* 
  ** These routines set/get unstructured descriptor fields given a pointer
  **  to a packet.
  */

INLINE_P LOGIC_ret makeUnstructPacketDataDescr (PACKET_ptr, UINT, UINT);

INLINE_P VALUE_ret getUnstructPacketWordSize (PACKET_ptr);
INLINE_P VALUE_ret getUnstructPacketDataPadding (PACKET_ptr);
INLINE_P VALUE_ret getUnstructPacketHitFormat (PACKET_ptr);

INLINE_P LOGIC_ret setUnstructPacketWordSize (PACKET_ptr, UINT);
INLINE_P LOGIC_ret setUnstructPacketDataPadding (PACKET_ptr, UINT);
INLINE_P LOGIC_ret setUnstructPacketHitFormat (PACKET_ptr, UINT);

VALUE_ret getUnstructPacketDataLengthBytes(PACKET_ptr);
VALUE_ret getUnstructPacketDataLengthWords(PACKET_ptr);
VALUE_ret extendUnstructPacketDataBlock(PACKET_ptr, UINT);


  /*
  **  These routines set/get fields from the descriptor given a pointer
  **    to the descriptor.
  */

INLINE_P LOGIC_ret makeUnstructDataDescr (PHDWORD*, UINT, UINT);

INLINE_P LOGIC_ret setUnstructDescrWordSize (PHDWORD*, UINT);
INLINE_P LOGIC_ret setUnstructDescrDataPadding (PHDWORD*, UINT);
INLINE_P LOGIC_ret setUnstructDescrHitFormat (PHDWORD*, UINT);

INLINE_P VALUE_ret getUnstructDescrWordSize (PHDWORD*);
INLINE_P VALUE_ret getUnstructDescrDataPadding (PHDWORD*);
INLINE_P VALUE_ret getUnstructDescrHitFormat (PHDWORD*);


INLINE_D LOGIC_ret makeUnstructPacketDataDescr (PACKET_ptr packet_ptr, 
					      UINT inWordSize, 
					      UINT inHitFormat)
{
  PHDWORD *descr_ptr;

  /*
  **  Set the descriptor length in the packet header
  */
  setPacketDataDescrLength(packet_ptr, UnstructDataDescrLength);

  descr_ptr = findPacketDataDescr(packet_ptr);
  if (descr_ptr == ptrFailure) return logicFailure;
  else return makeUnstructDataDescr(descr_ptr, inWordSize, inHitFormat);
}

INLINE_D VALUE_ret getUnstructPacketWordSize (PACKET_ptr packet_ptr)
{
  PHDWORD* descr_ptr = findPacketDataDescr(packet_ptr);
  if (descr_ptr == ptrFailure) return valueFailure;
  else return getUnstructDescrWordSize(descr_ptr);
}

INLINE_D VALUE_ret getUnstructPacketHitFormat (PACKET_ptr packet_ptr)
{
  PHDWORD* descr_ptr = findPacketDataDescr(packet_ptr);
  if (descr_ptr == ptrFailure) return valueFailure;
  else return getUnstructDescrHitFormat(descr_ptr);
}

INLINE_D VALUE_ret getUnstructPacketDataPadding (PACKET_ptr packet_ptr)
{
  PHDWORD* descr_ptr = findPacketDataDescr(packet_ptr);
  if (descr_ptr == ptrFailure) return valueFailure;
  else return getUnstructDescrDataPadding(descr_ptr);
}

INLINE_D LOGIC_ret setUnstructPacketHitFormat (PACKET_ptr packet_ptr, UINT unHitFormat)
{
  PHDWORD* descr_ptr = findPacketDataDescr(packet_ptr);
  if (descr_ptr == ptrFailure) return valueFailure;
  else {
    return setUnstructDescrHitFormat(descr_ptr, unHitFormat);
  }
}

INLINE_D LOGIC_ret setUnstructPacketWordSize (PACKET_ptr packet_ptr, UINT unWordSize)
{
  PHDWORD* descr_ptr = findPacketDataDescr(packet_ptr);
  if (descr_ptr == ptrFailure) return valueFailure;
  else {
    return setUnstructDescrWordSize(descr_ptr, unWordSize);
  }
}

INLINE_D LOGIC_ret setUnstructPacketDataPadding (PACKET_ptr packet_ptr, UINT paddingBytes)
{
  PHDWORD* descr_ptr = findPacketDataDescr(packet_ptr);
  if (descr_ptr == ptrFailure) return valueFailure;
  else {
    return setUnstructDescrDataPadding(descr_ptr, paddingBytes);
  }
}

  /*
  **  This routine creates a one dword data descriptor for unstructured data.  The 
  **  wordSize and hitFormat are user specified values, the data padding entry is
  **  initialized to zero.
  */

INLINE_D LOGIC_ret makeUnstructDataDescr (PHDWORD* dataDescr_ptr, UINT inWordSize, 
					UINT inHitFormat)
{
  *dataDescr_ptr = 0;
  setUnstructDescrWordSize (dataDescr_ptr, inWordSize);
  setUnstructDescrHitFormat (dataDescr_ptr, inHitFormat);

  return logicSuccess;
}

/*  Return the word size from the unstructured descriptor
 */
INLINE_D VALUE_ret getUnstructDescrWordSize(PHDWORD* descr_ptr)
{
  return getBitsMACRO(descr_ptr, UNSTRUCT_WORD_SIZE_OFFSET_OF_DWORD,
		      UNSTRUCT_WORD_SIZE_OFFSET_IN_DWORD, 
		      UNSTRUCT_WORD_SIZE_MASK);
}

/*  Return the hit format from the unstructured descriptor
 */
INLINE_D VALUE_ret getUnstructDescrHitFormat(PHDWORD* descr_ptr)
{
  return getBitsMACRO(descr_ptr, UNSTRUCT_HIT_FORMAT_OFFSET_OF_DWORD,
		      UNSTRUCT_HIT_FORMAT_OFFSET_IN_DWORD, 
		      UNSTRUCT_HIT_FORMAT_MASK);
}

/*  Return the data padding field from the unstructured descriptor
 */
INLINE_D VALUE_ret getUnstructDescrDataPadding(PHDWORD* descr_ptr)
{
  return getBitsMACRO(descr_ptr, UNSTRUCT_DATA_PADDING_OFFSET_OF_DWORD, 
		      UNSTRUCT_DATA_PADDING_OFFSET_IN_DWORD, 
		      UNSTRUCT_DATA_PADDING_MASK);
}

/*  Set the word size in the unstructured descriptor
 */
INLINE_D LOGIC_ret setUnstructDescrWordSize(PHDWORD* descr_ptr, UINT wordSize)
{
  setBitsMACRO(descr_ptr, UNSTRUCT_WORD_SIZE_OFFSET_OF_DWORD,
	       UNSTRUCT_WORD_SIZE_OFFSET_IN_DWORD, 
	       UNSTRUCT_WORD_SIZE_MASK, wordSize);
  return TRUE;
}

/*  Set the hit format in the unstructured descriptor
 */
INLINE_D LOGIC_ret setUnstructDescrHitFormat(PHDWORD* descr_ptr, UINT hitFormat)
{
  setBitsMACRO(descr_ptr, UNSTRUCT_HIT_FORMAT_OFFSET_OF_DWORD,
	       UNSTRUCT_HIT_FORMAT_OFFSET_IN_DWORD, 
	       UNSTRUCT_HIT_FORMAT_MASK, hitFormat);
  return TRUE;
}

/*  Set the padding field in the unstructured descriptor
 */
INLINE_D LOGIC_ret setUnstructDescrDataPadding(PHDWORD* descr_ptr, 
					     UINT paddingBytes)
{
  setBitsMACRO(descr_ptr, UNSTRUCT_DATA_PADDING_OFFSET_OF_DWORD, 
		      UNSTRUCT_DATA_PADDING_OFFSET_IN_DWORD, 
		      UNSTRUCT_DATA_PADDING_MASK, paddingBytes);
  return TRUE;
}

#ifdef __cplusplus
}
/* end of external C */
#endif

#endif
/* end of indef _DATABLOCK_ */




