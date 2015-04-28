/* 
** packetRoutines.C
** 
** Author: $Author: pinkenbu $  
**   Date: $Date: 2004/10/19 19:21:05 $ 
** 
** $Log: packetRoutines.C,v $
** Revision 1.4  2004/10/19 19:21:05  pinkenbu
** fix insure warnings
**
** Revision 1.3  2004/07/15 13:06:07  phoncs
** DLW: add extendPacketDataBlock call to extendUnstructDataBlock
**
** Revision 1.8  2001/03/12 17:13:00  kelly
** committed packetRoutines bugfix
**
** Revision 1.7  1998/12/16 15:40:52  markacs
** (stephen markacs) changes to rawFormat code to stop using bool,true,false
**
** Revision 1.6  1998/12/11 22:02:06  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  Contains routines which manipulate packets in a fixed-length contiguous
**    array in memory. Uses generic routines in Cpacket to extract fields
**    from the packet headersand uses routines in dataBlock to maintain the
**    relevant fields in the data block descriptor.
*/

#ifndef VXWORKS
#include "malloc.h"
#include <stdlib.h>
#endif

#include "phenixOnline.h"
#include "packetPublic.h"
#include "Cpacket.h"
#include "dataBlock.h"
#include "packetRoutines.h"

/*
** makeEmptyPacket
**
**   Routine to make a new packet header in a buffer pointed to by "newPacketPtr". 
**   The header is created with "empty" data, debug and error blocks.
**
**  A pointer to the PHDWORD immediately following the "packet" (the header) is returned
*/

PTR_ret makeEmptyPacket (PACKET_ptr packet_ptr, UINT maxPacketLen, UINT packetId)
{	
  UINT packetLength = makePacketHdr (packet_ptr, maxPacketLen);
  if (packetLength == valueFailure) return ptrFailure;
 
  /*
  **	Now set user-specified fields
  */
  setPacketId(packet_ptr, packetId);

  /*
  ** Success
  */
  return findPacketEnd(packet_ptr)+1;
}

/*
** setPacketUnstructured
*/
LOGIC_ret setPacketUnstructured (PACKET_ptr packet_ptr, UINT inWordSize,  
				 UINT inHitFormat)
{
  /*
  **  Check for valid input
  */
  if (inWordSize > 4 ) return logicFailure;

  /*
  **  Set the packet structure field
  */
  setPacketStructure(packet_ptr, Unstructured);

  /*
  **  Make the data descriptor.
  */

  makeUnstructPacketDataDescr (packet_ptr, inWordSize,  inHitFormat);

  return logicSuccess;
}

/*
**  Make an unstructured empty packet. a pointer to the PHDWORD following the empty
**  packet is returned
*/
PTR_ret makeUnstructPacket (PACKET_ptr packet_ptr, UINT maxPacketLength,  
			    UINT packetId, UINT inWordSize,  UINT inHitFormat)
{
  VALUE_ret packetLength;
  LOGIC_ret logicalResult;

  packetLength = makePacketHdr (packet_ptr, maxPacketLength);
  if (packetLength == valueFailure) return ptrFailure;
		
  /*
  **  Now set user-specified fields
  */
  setPacketId(packet_ptr, packetId);
  setPacketStructure(packet_ptr, Unstructured);

  logicalResult = makeUnstructPacketDataDescr(packet_ptr, inWordSize, inHitFormat);
  if (!logicalResult) return ptrFailure;
  else return findPacketEnd(packet_ptr)+1;
}

/*
**  Store a new error entry in the packet. Return the length of the 
**    resulting error block
*/
LOGIC_ret appendPacketError (PACKET_ptr packet_ptr,  UINT maxPacketLength, 
			     ERRORENTRYV1_ptr errorEntry_ptr)
{
  PHDWORD newLength;
  PHDWORD* newError_ptr = findPacketErrorStart(packet_ptr) + 
    getPacketErrorLength(packet_ptr);

  /*
  ** Extend the length of the error block (and the packet)
  */
  newLength = extendPacketErrorBlock(packet_ptr, maxPacketLength, 
				     errorEntryV1Length);
  if (newLength == valueFailure) return logicFailure;

  /*
  **  Now copy in the error entry
  */
  dwordCopy(newError_ptr, (PHDWORD*) errorEntry_ptr, errorEntryV1Length);
  return logicSuccess;
}

/*
**  Reserve space for numDwords in the packet debug block. Return a pointer to the 
**  start of the debug block. The assumption is that the user will completely
**  fill the requested space.
*/
PTR_ret reservePacketDebugData (PACKET_ptr packet_ptr, UINT maxPacketLength, 
				UINT numDwords)
{
  PHDWORD newLength;
  
  /*
  **  We can only do this when there's no error block
  */
  if (getPacketErrorLength(packet_ptr) > 0) {
    setPacketError(FORMAT_ERR_INVALID_APPEND,packet_ptr,0);
  }

  /*
  **  Extend the length of the debug block.
  */
  newLength = extendPacketDebugBlock(packet_ptr, maxPacketLength, numDwords);
  if (newLength == valueFailure) return ptrFailure;
  else return findPacketDebugStart(packet_ptr);
}

/*
**  Initiate a write operation into a pristine unstructured packet.
**
**    The packet must be initially empty. It is extended to allow space
**    fir maxNumWords of unstructured data. A pointer to the start of
**    the packet is returned.
*/
PTR_ret startUnstructDataWrite (PACKET_ptr packet_ptr,  UINT maxPacketLength, 
				PHDWORD maxNumWords)
{
  PHDWORD* write_ptr;
  PHDWORD dataBlockLength;
 
  /*
  **  
  **    in the packet. Also make sure there's no debug or error blocks yet because
  **    we assume that these will be written after the data
  */
  if (!emptyPacket(packet_ptr)) {
    setUserError(FORMAT_ERR_INVALID_APPEND, getPacketLength(packet_ptr));
    return ptrFailure;
  }

  /*
  **  Check for unstructured packet.
  */
  if (getPacketStructure(packet_ptr) != Unstructured) {
    setUserError (FORMAT_ERR_WRONG_STRUCTURE, getPacketStructure(packet_ptr));
    return ptrFailure;
  }

  /*
  **  Try to extend the unstructured data block
  */
  dataBlockLength = extendUnstructPacketDataBlock(packet_ptr, maxNumWords);
  if (dataBlockLength == valueFailure) return ptrFailure;

  /*
  ** Now update the packet length
  */
  extendPacketDataBlock(packet_ptr, maxPacketLength, dataBlockLength);

  write_ptr = findPacketDataStart (packet_ptr);

  /*
  **  Now return to the user the pointer to where data is to be written
  */
  return write_ptr;
}

/*
**  Finish writing unstructured block
*/
PTR_ret finishUnstructDataWrite (PACKET_ptr packet_ptr, UINT maxPacketLength, 
				   PHDWORD actualWords)
{
  PHDWORD deltaWords;
  PHDWORD newLength;
  PHDWORD maxNumWords = getUnstructPacketDataLengthWords(packet_ptr);


  /*
  **  Check to make sure that user hasn't written beyond the pre-set maximum
  */
  if (actualWords > maxNumWords) {
    /*
    **  This is a serious error, data will be corrupted
    */
    return ptrFailure;
  }

  /*
  ** Now update the length
  */
  deltaWords = maxNumWords - actualWords;
  newLength = extendUnstructPacketDataBlock(packet_ptr, -((int) deltaWords) );
  if (newLength == valueFailure) return ptrFailure;

  return findPacketEnd(packet_ptr)+1;
}

/*
**  storePacketHits
**
**    General routine to store data in a packet. Works for unstructured and (hopefully)
**    all structured packets. Currently only implementationis for unstructured packets.
**
**    For unstructured data, numHits = number of words of size "wordSize"
**                           data_arr should point to the input data.
**                           address_arr is ignored.
**
*/
VALUE_ret storePacketHits (PACKET_ptr packet_ptr,  UINT maxPacketLen, 
			   UINT* address_arr, BYTE* data_arr,  
			   UINT numHits,  UINT maxNumHits) 
{
  /*
  **  Check to make sure packet has valid header and is empty
  */
  if (!validPacketHdr(packet_ptr))
  {  
    setPacketError(FORMAT_ERR_INVALID_HEADER,packet_ptr,0);
    return valueFailure;
  }
  if (!emptyPacket(packet_ptr))
  {
    setPacketError(FORMAT_ERR_NONEMPTY_PACKET,packet_ptr,0);
    return valueFailure;
  }

  /*
  **  Handle different structure types
  */
  switch (getPacketStructure(packet_ptr)) {
  case Unstructured: {
    PHDWORD newLength, result;

    /*
    **  Get the size of the data block.
    */
    newLength = extendUnstructPacketDataBlock(packet_ptr, numHits);
    
    result = extendPacketDataBlock(packet_ptr, maxPacketLen, newLength);
    if (result == valueFailure) return valueFailure;
    else {
      BYTE *data_ptr;
      PHDWORD numBytes = numHits*getUnstructPacketWordSize(packet_ptr);

      /*
      **  Now copy the data in.
      */ 
      data_ptr = (BYTE*) findPacketDataStart(packet_ptr);
      byteCopy (data_ptr, data_arr, numBytes);
      data_ptr += numBytes;
      byteClear (data_ptr, getUnstructPacketDataPadding(packet_ptr));

      return getPacketLength(packet_ptr);
    }
    break;
  }
  default:
    setPacketError(FORMAT_ERR_WRONG_STRUCTURE,packet_ptr,0);
    return valueFailure;
  }
  
  /*
  ** Success
  */
  return getPacketLength(packet_ptr);
}

/*
**  fetchPacketHits
*/
VALUE_ret fetchPacketHits (PACKET_ptr packet_ptr, UINT** address_arr, BYTE** hits_arr, UINT* hitLength) 
{
  /*
  **  Make sure we got a pointer to a real packet
  */
  if (!validPacketHdr (packet_ptr))
    return valueFailure;
  
  /*
  **  Handle different structure types
  */
  switch (getPacketStructure (packet_ptr)) {
  case Unstructured:
    {
      PHDWORD numWords = getUnstructPacketDataLengthWords(packet_ptr);
      UINT wordSize = getUnstructPacketWordSize(packet_ptr);
      PHDWORD* unformBlock_ptr = findPacketDataStart (packet_ptr);
      PHDWORD numBytes = numWords*wordSize;
      
      /*
      **	Allocate space for data array
      */
      BYTE* hits_ptr	= (BYTE *) malloc (numBytes);
      byteCopy (hits_ptr, (BYTE *) unformBlock_ptr, numBytes);
      
      /*
      **  Alllocate space for one address (0)
      */
      *address_arr = (UINT*) malloc (sizeof(UINT));
      **address_arr = 0;
      
      *hits_arr = hits_ptr;
      *hitLength = wordSize;
      return numWords;
    }
}
  
  /*
  ** Success
  */
  
  return 0;
}

















