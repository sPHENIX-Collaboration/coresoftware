/* 
** CpacketV1.h
** 
** Author: $Author: phnxbld $  
**   Date: $Date: 2009/09/19 14:34:09 $ 
** 
** $Log: CpacketV1.h,v $
** Revision 1.2  2009/09/19 14:34:09  phnxbld
** add parenthesis to take care of operator precedence
**
** Revision 1.1.1.1  2000/07/21 01:51:10  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.5  1998/12/11 22:01:13  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  CpacketV1.h
**
**  These are the version one routines for accessing and modifying
**  the version dependent data in the packet header.  
*/

#ifndef _CPACKETV1_
#define _CPACKETV1_
/*
**  Use C linkage
*/
#ifdef __cplusplus
extern "C" {
#endif

#include "phenixOnline.h"
#include "packetPublic.h"
#include "packetV1Public.h"
#include "packetHdrV1.h"
#include "formatError.h"
#include "dataBlock.h"

VALUE_ret makePacketV1Hdr (PACKET_ptr, UINT);

PTR_ret findPacketV1DataStart (PACKET_ptr);
PTR_ret findPacketV1DataEnd (PACKET_ptr);
PTR_ret findPacketV1DebugStart (PACKET_ptr);
PTR_ret findPacketV1ErrorStart (PACKET_ptr);

INLINE_P PTR_ret findPacketV1DataDescr (PACKET_ptr);

INLINE_P VALUE_ret getPacketV1DataLength (PACKET_ptr);
INLINE_P VALUE_ret getPacketV1DebugLength (PACKET_ptr);
INLINE_P VALUE_ret getPacketV1ErrorLength (PACKET_ptr);
INLINE_P VALUE_ret getPacketV1Padding (PACKET_ptr);
INLINE_P VALUE_ret getPacketV1NumErrors (PACKET_ptr);
INLINE_P VALUE_ret getPacketV1Id (PACKET_ptr);
INLINE_P VALUE_ret getPacketV1Structure (PACKET_ptr);
INLINE_P VALUE_ret getPacketV1Endianism (PACKET_ptr);
INLINE_P VALUE_ret getPacketV1Status (PACKET_ptr);
INLINE_P VALUE_ret getPacketV1DataDescrLength (PACKET_ptr);

INLINE_P VALUE_ret orPacketV1Status (PACKET_ptr, UINT);
INLINE_P VALUE_ret adjustPacketV1DataLength (PACKET_ptr, int);
INLINE_P VALUE_ret adjustPacketV1DebugLength (PACKET_ptr, int);
INLINE_P VALUE_ret adjustPacketV1ErrorLength (PACKET_ptr, int);

LOGIC_ret validPacketV1Hdr (PACKET_ptr);

INLINE_P LOGIC_ret emptyPacketV1 (PACKET_ptr);
INLINE_P LOGIC_ret setPacketV1Id (PACKET_ptr, UINT);
INLINE_P LOGIC_ret setPacketV1Endianism (PACKET_ptr, UINT);
INLINE_P LOGIC_ret setPacketV1Structure (PACKET_ptr, UINT);
INLINE_P LOGIC_ret setPacketV1Padding (PACKET_ptr, UINT);
INLINE_P LOGIC_ret setPacketV1Status (PACKET_ptr, UINT);
INLINE_P LOGIC_ret setPacketV1ErrorLength (PACKET_ptr, UINT);
INLINE_P LOGIC_ret setPacketV1DebugLength (PACKET_ptr, UINT);
INLINE_P LOGIC_ret setPacketV1NumErrors (PACKET_ptr, UINT);
INLINE_P LOGIC_ret setPacketV1DataDescrLength (PACKET_ptr, UINT);

/*
**  Check to see if this is an empty packet.
*/
INLINE_D LOGIC_ret emptyPacketV1 (PACKET_ptr packet_ptr)
{
  return getPacketHdrLength(packet_ptr) == getPacketLength (packet_ptr);
}

INLINE_D VALUE_ret getPacketV1DataLength (PACKET_ptr packet_ptr)
{

  PHDWORD dataLength;
  if (getPacketStructure(packet_ptr) == Unstructured)
    {

      int factor = 4 / getUnstructPacketWordSize(packet_ptr);

      dataLength = factor * ( getPacketLength(packet_ptr)
			      - getPacketHdrLength(packet_ptr)
			      - getPacketErrorLength(packet_ptr)
			      - getPacketDebugLength(packet_ptr) )
	- getPacketPadding(packet_ptr);

      return dataLength;

    }
  
  dataLength = getPacketLength(packet_ptr) - getPacketHdrLength(packet_ptr) -
                     getPacketErrorLength(packet_ptr) - getPacketDebugLength(packet_ptr) -
                     getPacketPadding(packet_ptr);
  return dataLength;
}

INLINE_D VALUE_ret getPacketV1DebugLength (PACKET_ptr packet_ptr)
{
  return getBitsMACRO(packet_ptr, DEBUG_LENGTH_OFFSET_OF_DWORD,
		      DEBUG_LENGTH_OFFSET_IN_DWORD, DEBUG_LENGTH_MASK);
}

INLINE_D VALUE_ret getPacketV1ErrorLength (PACKET_ptr packet_ptr)
{
  return  getBitsMACRO(packet_ptr, PACKET_ERROR_LENGTH_OFFSET_OF_DWORD,
		      PACKET_ERROR_LENGTH_OFFSET_IN_DWORD, 
		      PACKET_ERROR_LENGTH_MASK);
}

INLINE_D VALUE_ret getPacketV1Padding (PACKET_ptr packet_ptr)
{
  return  getBitsMACRO(packet_ptr, PACKET_PADDING_OFFSET_OF_DWORD,
		       PACKET_PADDING_OFFSET_IN_DWORD, PACKET_PADDING_MASK);
}

INLINE_D VALUE_ret getPacketV1NumErrors (PACKET_ptr packet_ptr)
{
  UINT errorBlockLength = getPacketV1ErrorLength(packet_ptr);
  return calcNumErrorsV1 (errorBlockLength);
}

INLINE_D VALUE_ret getPacketV1Id (PACKET_ptr packet_ptr)
{
  return  getWordMACRO(packet_ptr, ID_OFFSET_OF_DWORD);
}

INLINE_D VALUE_ret getPacketV1Structure (PACKET_ptr packet_ptr)
{
  return  getBitsMACRO(packet_ptr, STRUCTURE_OFFSET_OF_DWORD,
		      STRUCTURE_OFFSET_IN_DWORD, STRUCTURE_MASK);
}

INLINE_D VALUE_ret getPacketV1Status (PACKET_ptr packet_ptr)
{
  return  getBitsMACRO(packet_ptr, PACKET_STATUS_OFFSET_OF_DWORD, 
		       PACKET_STATUS_OFFSET_IN_DWORD, PACKET_STATUS_MASK);
} 

INLINE_D VALUE_ret getPacketV1DataDescrLength (PACKET_ptr packet_ptr)
{
  return  getBitsMACRO(packet_ptr, DESCR_LENGTH_OFFSET_OF_DWORD, 
		       DESCR_LENGTH_OFFSET_IN_DWORD, DESCR_LENGTH_MASK);
} 

INLINE_D VALUE_ret getPacketV1Endianism (PACKET_ptr packet_ptr)
{
  return  getBitsMACRO(packet_ptr, ENDIANISM_OFFSET_OF_DWORD,
		      ENDIANISM_OFFSET_IN_DWORD, ENDIANISM_MASK);
}

INLINE_D VALUE_ret orPacketV1Status (PACKET_ptr packet_ptr, UINT inBits)
{
  UINT status;
  if ((inBits & ((1<<PACKET_STATUS_NUM_BITS)-1)) == inBits) {
    status = getBitsMACRO(packet_ptr, PACKET_STATUS_OFFSET_OF_DWORD,
			  PACKET_STATUS_OFFSET_IN_DWORD, PACKET_STATUS_MASK);
    status|=inBits;
    setBitsMACRO(packet_ptr, PACKET_STATUS_OFFSET_OF_DWORD, 
		 PACKET_STATUS_OFFSET_IN_DWORD,
		 PACKET_STATUS_MASK, status);
    return status;
  }
  else {
    setUserError (FORMAT_ERR_INVALID_DATA, inBits);
    return valueFailure;
  }
}  

/*
**  Update the length of the data block in the header. In V1 packets, this is a NOOP
*/
INLINE_D VALUE_ret adjustPacketV1DataLength( PACKET_ptr packet_ptr, int addDwords)
{
  return getPacketV1DataLength(packet_ptr);
}

INLINE_D VALUE_ret adjustPacketV1DebugLength(PACKET_ptr packet_ptr, int addLength)
{
  UINT debugLength = getPacketV1DebugLength (packet_ptr);
  if (debugLength == valueFailure)
    return valueFailure;
  else {
    UINT newLength = debugLength + addLength;
    setPacketV1DebugLength(packet_ptr, debugLength);
    return newLength;
  }
}

INLINE_D VALUE_ret adjustPacketV1ErrorLength(PACKET_ptr packet_ptr, int addLength)
{
  UINT errorLength = getPacketV1ErrorLength (packet_ptr);
  if (errorLength == valueFailure)
    return valueFailure;
  else {
    UINT newLength = errorLength + addLength;
    setPacketV1ErrorLength(packet_ptr, newLength);
    return newLength;
  }
}

INLINE_D LOGIC_ret setPacketV1Id (PACKET_ptr packet_ptr, UINT inId)
{
  setWordMACRO(packet_ptr, ID_OFFSET_OF_DWORD, inId);
  return TRUE;
}

INLINE_D LOGIC_ret setPacketV1Endianism (PACKET_ptr packet_ptr, UINT inEndianism)
{
  setBitsMACRO(packet_ptr, ENDIANISM_OFFSET_OF_DWORD, 
	       ENDIANISM_OFFSET_IN_DWORD, ENDIANISM_MASK, 
	       inEndianism);
  return TRUE;
}

INLINE_D LOGIC_ret setPacketV1Structure (PACKET_ptr packet_ptr, UINT inStruct)
{
  setBitsMACRO(packet_ptr, STRUCTURE_OFFSET_OF_DWORD, 
	       STRUCTURE_OFFSET_IN_DWORD, STRUCTURE_MASK, inStruct);
  return TRUE;
}

INLINE_D LOGIC_ret setPacketV1Padding (PACKET_ptr packet_ptr, UINT inPadding)
{
  setBitsMACRO(packet_ptr, PACKET_PADDING_OFFSET_OF_DWORD, 
	       PACKET_PADDING_OFFSET_IN_DWORD, PACKET_PADDING_MASK, inPadding);
  return TRUE;
}

INLINE_D LOGIC_ret setPacketV1Status (PACKET_ptr packet_ptr, UINT inStatus)
{
  setBitsMACRO(packet_ptr, PACKET_STATUS_OFFSET_OF_DWORD, 
	       PACKET_STATUS_OFFSET_IN_DWORD, PACKET_STATUS_MASK, inStatus);
  return TRUE;
}

INLINE_D LOGIC_ret setPacketV1ErrorLength (PACKET_ptr packet_ptr, UINT inLength)
{
  setBitsMACRO(packet_ptr, PACKET_ERROR_LENGTH_OFFSET_OF_DWORD, 
	       PACKET_ERROR_LENGTH_OFFSET_IN_DWORD, PACKET_ERROR_LENGTH_MASK, inLength);
  return TRUE;
}

INLINE_D LOGIC_ret setPacketV1DebugLength (PACKET_ptr packet_ptr, UINT inLength)
{
  setBitsMACRO(packet_ptr, DEBUG_LENGTH_OFFSET_OF_DWORD, 
	       DEBUG_LENGTH_OFFSET_IN_DWORD, DEBUG_LENGTH_MASK, inLength);
  return TRUE;
}

INLINE_D LOGIC_ret setPacketV1DataDescrLength (PACKET_ptr packet_ptr, UINT inLength)
{
  setBitsMACRO(packet_ptr, DESCR_LENGTH_OFFSET_OF_DWORD, 
	       DESCR_LENGTH_OFFSET_IN_DWORD, DESCR_LENGTH_MASK, inLength);
  return TRUE;
}


/*
** Return a pointer to the data descriptor
*/
INLINE_D PTR_ret findPacketV1DataDescr (PACKET_ptr packet_ptr)
{
  return packet_ptr + DATADESCR_OFFSET_OF_DWORD;
}


#ifdef __cplusplus
} 
/* end of extern "C" */
#endif

#endif 
/* end of ifndef _CPACKET_ */
