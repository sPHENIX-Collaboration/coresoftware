/* 
** Cpacket.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:10 $ 
** 
** $Log: Cpacket.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:10  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.6  1998/12/17 21:57:06  phoncs
** (stephen markacs) better bounds checking in version gets
**
** Revision 1.5  1998/12/11 22:01:13  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  Cpacket.h
**
**    Defines prototypes for all public functions in Cpacket.C.
**
**    These routines provide access to information in the header, allow 
**    modification of data in the header and provide various logical
**    tests of packet integrity, etc.  The version independent
**    routines contained here are inlined.  Version independent 
**    routines which are not inlined are defined in Cpacket.C.
**
**    In addition to the packet header, there is a small data
**    descriptor of one or two dwords which immediately follows the
**    packet header at the head of the data block. The data descriptor
**    contains information specific to the data block and has a format which  
**    varies depending on the type of data structure used in the data
**    block. The descriptor is defined in dataBlock.h and is placed
**    after the packet in packetRoutines.C. 
**
*/
#ifndef _CPACKET_
#define _CPACKET_

#ifdef __cplusplus
extern "C" {
#endif

#include "phenixOnline.h"
#include "packetPublic.h"
#include "packetHdr.h"
#include "formatError.h"

  /*
  **  The following routines are GENERIC and do not depend on packet 
  **  header version. They are either inlined in this file or coded in 
  **  Cpacket.C
  */  

INLINE_P LOGIC_ret currentPacket (PACKET_ptr);

INLINE_P VALUE_ret getPacketLength (PACKET_ptr);
INLINE_P VALUE_ret getPacketHdrLength (PACKET_ptr);
INLINE_P VALUE_ret getPacketHdrVersion (PACKET_ptr);

INLINE_P LOGIC_ret setPacketHdrVersion(PACKET_ptr, UINT);
INLINE_P LOGIC_ret setPacketHdrLength(PACKET_ptr, UINT);
INLINE_P LOGIC_ret setPacketLength(PACKET_ptr, PHDWORD);

PTR_ret findPacketEnd (PACKET_ptr);
PTR_ret findPacketDataStart (PACKET_ptr);

VALUE_ret extendPacketDataBlock(const PACKET_ptr, UINT, int);
VALUE_ret extendPacketDebugBlock(const PACKET_ptr, UINT, int);
VALUE_ret extendPacketErrorBlock(const PACKET_ptr, UINT, int);
VALUE_ret extendPacketLength(const PACKET_ptr, UINT, int);

LOGIC_ret removePacketPadding(PACKET_ptr);

INLINE_P void endianSwapError (PACKETERROR* outError, PACKETERROR* inError);

  /*
  **  The following routines execute functions that will depend on the
  **  version of the packet header. Outside the DCM's these routines
  **  use pointer tables defined below to execute a version-specific
  **  function with identical definition. For use in the DCM's the
  **  generic function names are #defined to be the version-specific 
  **  for the "current" packet header version being created in the DAQ.
  */

#ifdef DCM 

#define makePacketHdr makePacketHdrV1

#define getPacketDebugLength getPacketV1DebugLength
#define getPacketErrorLength getPacketV1ErrorLength
#define getPacketNumErrors getPacketV1NumErrors
#define getPacketEndianism getPacketV1Endianism
#define getPacketId getPacketV1Id
#define getPacketPadding getPacketV1Padding
#define getPacketStructure getPacketV1Structure
#define getPacketStatus getPacketV1Status
#define getPacketDataDescrLength getPacketDataDescrLengthV1

#define setPacketDebugLength setPacketV1DebugLength
#define setPacketNumErrors setPacketV1NumErrors
#define setPacketEndianism setPacketV1Endianism
#define setPacketId setPacketV1Id
#define setPacketPadding setPacketV1Padding
#define setPacketStructure  setPacketV1Structure
#define setPacketStatus setPacketV1Status
#define setPacketErrorLength setPacketV1ErrorLength

#define orPacketStatus orPacketV1Status
#define adjustPacketDataLength adjustPacketV1DataLength
#define adjustPacketDebugLength adjustPacketV1DebugLength
#define adjustPacketErrorLength adjustPacketV1ErrorLength

#define findPacketErrorStart findPacketV1ErrorStart
#define findPacketDebugStart findPacketV1DebugStart
#define findPacketDataStart findPacketV1DataStart
#define findPacketDataEnd findPacketV1DataEnd

#define validPacketHdr validPacketV1Hdr
#define emptyPacket emptyPacketV1	

#else

VALUE_ret makePacketHdr(const PACKET_ptr, UINT);
VALUE_ret getPacketDataLength (const PACKET_ptr);
VALUE_ret getPacketDebugLength (const PACKET_ptr);
VALUE_ret getPacketErrorLength (const PACKET_ptr);
VALUE_ret getPacketNumErrors (const PACKET_ptr);
VALUE_ret getPacketPadding (const PACKET_ptr);
VALUE_ret getPacketId (const PACKET_ptr);
VALUE_ret getPacketStatus (PACKET_ptr);
VALUE_ret getPacketStructure (const PACKET_ptr);
VALUE_ret getPacketEndianism (const PACKET_ptr);
VALUE_ret getPacketDataDescrLength (const PACKET_ptr);

LOGIC_ret setPacketDebugLength (PACKET_ptr, UINT); 
LOGIC_ret setPacketNumErrors (PACKET_ptr, UINT);
LOGIC_ret setPacketEndianism (PACKET_ptr, UINT);
LOGIC_ret setPacketId (PACKET_ptr, UINT);
LOGIC_ret setPacketPadding (PACKET_ptr, UINT);
LOGIC_ret setPacketStructure (PACKET_ptr, UINT); 
LOGIC_ret setPacketStatus (PACKET_ptr, UINT);
LOGIC_ret setPacketErrorLength  (PACKET_ptr, UINT);
LOGIC_ret setPacketDataDescrLength (PACKET_ptr, UINT);

VALUE_ret orPacketStatus (PACKET_ptr, UINT);   
LOGIC_ret validPacketHdr (PACKET_ptr);
LOGIC_ret emptyPacket (PACKET_ptr);

PTR_ret findPacketDataDescr(PACKET_ptr);
PTR_ret findPacketDataStart(PACKET_ptr);
PTR_ret findPacketDataEnd (PACKET_ptr);
PTR_ret findPacketDebugStart (PACKET_ptr);
PTR_ret findPacketErrorStart (PACKET_ptr);

VALUE_ret adjustPacketDataLength (PACKET_ptr, int);
VALUE_ret adjustPacketErrorLength (PACKET_ptr, int);
VALUE_ret adjustPacketDebugLength (PACKET_ptr, int);

#endif


/* 
**  ==================================================================
**  Inlined GENERIC routines that do not need the indirection arrays
**  ==================================================================
*/

INLINE_D VALUE_ret getPacketHdrLength (PACKET_ptr packet_ptr)
{
  return getBitsMACRO(packet_ptr, 
		      PACKET_HDR_LENGTH_OFFSET_OF_DWORD,
		      PACKET_HDR_LENGTH_OFFSET_IN_DWORD,
		      PACKET_HDR_LENGTH_MASK);
}


INLINE_D VALUE_ret getPacketLength (PACKET_ptr packet_ptr)
{
  PHDWORD length = getWordMACRO(packet_ptr, PACKET_LENGTH_OFFSET_OF_DWORD);

#ifdef STRONGCHECK
  if ((length&(1<<32))!=0) return valueFailure;
#endif

  return length;
}

INLINE_D VALUE_ret getPacketHdrVersion (PACKET_ptr packet_ptr)
{
 BYTE hdrVersion = (BYTE)getBitsMACRO(packet_ptr, 
				      PACKET_HDR_VERSION_OFFSET_OF_DWORD,
				      PACKET_HDR_VERSION_OFFSET_IN_DWORD,
				      PACKET_HDR_VERSION_MASK);

  if ((hdrVersion >= numPacketVersions)||(hdrVersion==0)) {
    setPacketError(FORMAT_ERR_INVALID_PACKET_HDRVERSION, 
		   packet_ptr, hdrVersion);
    return valueFailure;
  }  

  return hdrVersion;
}

INLINE_D LOGIC_ret setPacketHdrVersion (PACKET_ptr packet_ptr, UINT hdrVersion)
{
#ifdef CHECK
  if ((hdrVersion >= numPacketVersions)||(hdrVersion==0)) {
    setPacketError(FORMAT_ERR_INVALID_PACKET_HDRVERSION, 
		   packet_ptr, hdrVersion);
    return valueFailure;
  }  
#endif 
  setBitsMACRO(packet_ptr, 
	       PACKET_HDR_VERSION_OFFSET_OF_DWORD,
	       PACKET_HDR_VERSION_OFFSET_IN_DWORD, 
	       PACKET_HDR_VERSION_MASK, hdrVersion);
  return TRUE;
}

INLINE_D LOGIC_ret setPacketLength (PACKET_ptr packet_ptr, PHDWORD length)
{
  /*
  **  Currently we do no checking -- return "true"
  */
#ifdef CHECK
  if (length >= MAX_DWORD_VALUE) {
    setPacketError(FORMAT_ERR_LENGTH_OVERFLOW, packet_ptr, length);
    return FALSE;
  }
#endif
  setWordMACRO(packet_ptr, PACKET_LENGTH_OFFSET_OF_DWORD, length);
  return TRUE;
}


INLINE_D LOGIC_ret setPacketHdrLength (PACKET_ptr packet_ptr, UINT hdrLength)
{
#ifdef CHECK
  if (hdrLength >= ) {
    setPacketError(FORMAT_ERR_INVALID_HDRLENGTH, packet_ptr, hdrLength);
    return FALSE;
  }
#endif

  setBitsMACRO(packet_ptr, PACKET_HDR_LENGTH_OFFSET_OF_DWORD,
	       PACKET_HDR_LENGTH_OFFSET_IN_DWORD, 
	       PACKET_HDR_LENGTH_MASK, hdrLength);
  return TRUE;
}


INLINE_D LOGIC_ret currentPacket (PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion(packet_ptr);
  return (version == currentPacketHdrVersion);
}

INLINE_D void endianSwapError (PACKETERROR* outError, PACKETERROR* inError)
{
  /*
  **  Execute the swap function provided in the V1 error routines
  */
  endianSwapErrorV1 (outError, inError);
}

#ifdef __cplusplus
}
#endif 
/* end of extern "C" */

#endif
/* end of ifndef _CPACKET_ */










