/* 
** dataBlock.C
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:11 $ 
** 
** $Log: dataBlock.C,v $
** Revision 1.1.1.1  2000/07/21 01:51:11  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.4  1998/12/11 22:02:00  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  These are the routines for acces to and modification 
**  of the data block descriptor which are not inlined within
**  dataBlock.h.  
**/
#include "dataBlock.h"

/*
**  Return the number of "words" stored in the unstructured block
**    where a word is defined by the length stored in the descriptor
*/
VALUE_ret getUnstructPacketDataLengthBytes (PACKET_ptr packet_ptr)
{
  PHDWORD dataBlockLengthDwords, dataBlockLengthBytes, paddingBytes;

  dataBlockLengthDwords = getPacketDataLength(packet_ptr);
  if (dataBlockLengthDwords == valueFailure) return valueFailure;

  dataBlockLengthBytes = dataBlockLengthDwords*PHDWORD_SIZE;
  paddingBytes = getUnstructPacketDataPadding(packet_ptr);
  return dataBlockLengthBytes - paddingBytes;
}

/*
**  Return the number of "words" stored in the unstructured block
**    where a word is defined by the length stored in the descriptor
*/
VALUE_ret getUnstructPacketDataLengthWords (PACKET_ptr packet_ptr)
{
  PHDWORD dataLengthBytes = getUnstructPacketDataLengthBytes(packet_ptr);
  UINT wordSize = getUnstructPacketWordSize(packet_ptr);
  PHDWORD dataLengthWords = dataLengthBytes/wordSize;

  if (dataLengthBytes % wordSize != 0) {
    setPacketError(FORMAT_ERR_DATA_INCONSISTENCY,packet_ptr,dataLengthBytes);
    return valueFailure;
  }
  else return dataLengthWords;
}

/*
** Returns the length in PHDWORDS of the data block after addWords
**   "words" of size wordSize are added to the data block in the
**   packet. Updates the padding field in the packet.
*/
VALUE_ret extendUnstructPacketDataBlock (PACKET_ptr packet_ptr, 
					 UINT addWords)
{
  PHDWORD currentLength;
  PHDWORD newLength, newLengthBytes, newLengthDwords;
  UINT newPadding;
  UINT modulo;

  /*
  ** Get the current data block length
  */
  currentLength = getUnstructPacketDataLengthWords(packet_ptr);
  if (currentLength == valueFailure) return valueFailure;

  /*
  ** Calculate the new length 
  */
  newLength = currentLength + addWords;
  newLengthBytes = newLength*getUnstructPacketWordSize(packet_ptr);

  /*
  **  Calculate the new length in dwords
  */
  modulo = newLengthBytes % PHDWORD_SIZE;
  if (modulo > 0) newPadding = PHDWORD_SIZE - modulo;
  else newPadding = 0;

  newLengthDwords = (newLengthBytes + newPadding)/PHDWORD_SIZE;

  /*
  **  Now set the padding field
  */
  setUnstructPacketDataPadding(packet_ptr, newPadding);

  return newLengthDwords;
}
