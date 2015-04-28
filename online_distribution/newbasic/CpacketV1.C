/* 
** CpacketV1.C
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:10 $ 
** 
** $Log: CpacketV1.C,v $
** Revision 1.1.1.1  2000/07/21 01:51:10  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.4  1999/07/12 14:23:37  phoncs
** (stephen markacs) added endianism setting to makePacketV1Hdr() in CpacketV1.C
**
** Revision 1.3  1998/12/11 22:01:59  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  These are the Version 1 routines to provide
**  information from and acces to the fields in the 
**  packet header.  The inlined routines are 
**  contained in CpacketV1.h.  
**
*/
#include "phenixOnline.h"
#include "Cpacket.h"
#include "CpacketV1.h"

VALUE_ret makePacketV1Hdr (PHDWORD* packet_ptr, UINT maxPacketLength)
{
  if (maxPacketLength < packetV1HdrLength) 
    return valueFailure;

  dwordClear(packet_ptr, packetV1HdrLength);
  setPacketLength(packet_ptr, packetV1HdrLength);
  setPacketHdrVersion (packet_ptr, packetV1HdrVersion);
  setPacketHdrLength (packet_ptr, packetV1HdrLength);

  unsigned long n = 1;
  char * cp = (char*)&n;
  if (*cp==1) setPacketEndianism(packet_ptr, 1);
         else setPacketEndianism(packet_ptr, 2);

  return packetV1HdrLength;
}

/*
** Return a pointer to the debug block.
*/
PTR_ret findPacketV1DebugStart (PACKET_ptr packet_ptr)
{
  return (packet_ptr + getPacketLength(packet_ptr) - 
          getPacketV1DebugLength(packet_ptr) - getPacketV1ErrorLength(packet_ptr) );
}

/*
** Return a pointer to the error block.
*/
PTR_ret findPacketV1ErrorStart (PACKET_ptr packet_ptr)
{
  return (packet_ptr + getPacketLength(packet_ptr) - 
	  getPacketV1ErrorLength(packet_ptr));
}

PTR_ret findPacketV1DataStart (PACKET_ptr packet_ptr)
{
  return packet_ptr + getPacketHdrLength(packet_ptr);
} 

PTR_ret findPacketV1DataEnd (PACKET_ptr packet_ptr)
{
  return findPacketV1DebugStart(packet_ptr) - 1;
} 


/*
**  Check for a valid V1 header. Currently this test only checks the
**    length to make sure its consistent with the expected length.
**
*/
LOGIC_ret validPacketV1Hdr (PACKET_ptr packet_ptr)
{
  if (getPacketHdrLength(packet_ptr) == packetV1HdrLength) {
    setPacketSuccess();
    return TRUE;
  }
  else {
    setPacketError (FORMAT_ERR_INVALID_HEADER, packet_ptr, 0);
    return FALSE;
  }
}















