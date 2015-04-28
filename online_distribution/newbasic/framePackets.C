/* 
** framePackets.C
** 
** Author: $Author: phnxbld $  
**   Date: $Date: 2009/08/19 12:30:51 $ 
** 
** $Log: framePackets.C,v $
** Revision 1.3  2009/08/19 12:30:51  phnxbld
** fix check on failure
**
** Revision 1.2  2002/05/31 22:15:42  purschke
** mlp -- went through the insure report and tried to fix
** every little problem there is, unused variables, dead statements,
** all. It'll probably take another round to complete, but it should get
** rid of most warnings.
**
** Revision 1.1.1.1  2000/07/21 01:51:13  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.6  1998/12/23 20:34:57  markacs
** fixed findFramePacketId so it returns ptrFailure if called for an index of 0 in the case that there are no packets in the frame
**
** Revision 1.5  1998/12/11 22:02:04  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  Routines that perform functions on packets with knowledge of the 
**    frame that the packets reside in.
**
**/

#include "framePackets.h"
#include "Cframe.h"
#include "Cpacket.h"
#include "packetPublic.h"
#include "framePublic.h"
#include "frameRoutines.h"
#include "packetRoutines.h"
#include "formatError.h"


PACKET_ptr appendEmptyFramePacket (FRAME_ptr frame_ptr, PHDWORD maxFrameLength, UINT packetId)
{
  UINT extendDwords;
  UINT maxPacketLength;
  PACKET_ptr newPacket_ptr;

  /*
  **  Check for valid frame header
  */
  if (!validFrameHdr (frame_ptr))
    return ptrFailure;

  /*
  **  Now extend the frame to hold the new packet header
  */
  extendDwords = extendFrameData(frame_ptr, maxFrameLength, packetV1HdrLength);
  if (extendDwords == valueFailure) {
    setUserError(FORMAT_ERR_BUFFER_OVERFLOW, maxFrameLength);
    return ptrFailure;
  }

  /*
  **  Get pointer to data section
  */
  newPacket_ptr = findFrameDataStart(frame_ptr);

  /*
  **  Now make the empty packet
  */
  maxPacketLength = maxFrameLength - (UINT) (newPacket_ptr - frame_ptr);
  makeEmptyPacket(newPacket_ptr, maxPacketLength, packetId);
	
  setPacketSuccess();
  return newPacket_ptr;
}

LOGIC_ret isLastFramePacket (FRAME_ptr frame_ptr, PACKET_ptr packet_ptr) 
{
  PACKET_ptr lastPacket_ptr;

  /*
  **  Check for valid frame, packet
  */
  if (!validFrameHdr (frame_ptr)) return FALSE;
  if (!validPacketHdr (frame_ptr)) return FALSE;

  /*
  ** Now find the last packet
  */
  lastPacket_ptr = findLastFramePacket (frame_ptr);

  return (lastPacket_ptr == packet_ptr);

}

/*
**  Return the index'th packet in the frame where the indexing starts from 0.
**    (i.e. findFramePacketIndex(frame_ptr, 0) returns the first packet)
*/
PACKET_ptr findFramePacketIndex(FRAME_ptr frame_ptr, UINT index)
{
  PACKET_ptr testPacket_ptr;

  /*
  **  Check first for valid frame header
  */
  if (!validFrameHdr (frame_ptr)) {
    return ptrFailure;
  }

  /*
  **  Point to first packet
  */
  testPacket_ptr = (PACKET_ptr) findFrameDataStart (frame_ptr);

  if (index > 0) {
    UINT iPacket;

    /*
    **  Loop the necessary number of times but check to make sure we have valid pointer
    */
    for (iPacket=0; iPacket<index ; iPacket++) {
      testPacket_ptr = findNextFramePacket (frame_ptr, testPacket_ptr);
      if (testPacket_ptr == ptrFailure) {
	return ptrFailure;
      }
    }
  }
  else  /* if index=0 make sure there's a packet there */
    {
      if (!validPacketHdr(testPacket_ptr)) return ptrFailure;
    }

  /*
  **  If we get here we have the requested packet pointer
  */
  return testPacket_ptr;
}

/*
**  Return a pointer to the first packet in the frame
*/
PACKET_ptr findFirstFramePacket(FRAME_ptr frame_ptr)
{
  return findFramePacketIndex(frame_ptr, 0);
}

/*
**  Return a pointer to the first packet in the frame
*/
PACKET_ptr findLastFramePacket(FRAME_ptr frame_ptr)
{
  PACKET_ptr thisPacket_ptr, nextPacket_ptr;

  /*
  **  Find the first packet
  */
  if (!(nextPacket_ptr = findFirstFramePacket (frame_ptr))) return ptrFailure;

  /*
  **  Now iterate to the end
  */
  while (nextPacket_ptr != ptrFailure) {
    thisPacket_ptr = nextPacket_ptr;
    nextPacket_ptr = findNextFramePacket (frame_ptr, nextPacket_ptr);
  }

  return thisPacket_ptr;
}

/*
**  findFramePacketID 
**
**		Find a packet with a specified packet ID in a frame starting at the 
**		  packet in the frame. If no packet with the ID is found ptrFailure is returned.
**
*/
PACKET_ptr findFramePacketId (FRAME_ptr frame_ptr, UINT packetId)
{
  PACKET_ptr testPacket_ptr;

  /*
  **  Check to make sure we have a good frame
  */
  if (!validFrameHdr (frame_ptr)) {
    return ptrFailure;
  }

  /*
  **  Point to first packet
  */
  testPacket_ptr = (PACKET_ptr) findFrameDataStart (frame_ptr);

  /*
  **  Now iterate until we fail to get enxt packet or find the one we want
  */
  while ( (testPacket_ptr != ptrFailure) && (getPacketId(testPacket_ptr) != packetId) )
    testPacket_ptr = findNextFramePacket (frame_ptr, testPacket_ptr);

  return testPacket_ptr;
}


PACKET_ptr findNextFramePacketId (FRAME_ptr frame_ptr, PACKET_ptr previousPacket_ptr, UINT packetId)
{
  /*
  **  BEWARE !!! No check on valid frame once we're past the first packet.
  */
	
	/*
**  Skip the "previous" packet
*/	
  PACKET_ptr testPacket_ptr = findNextFramePacket(frame_ptr, previousPacket_ptr);

  /*
	**  Now iterate until we fail to get enxt packet or find the one we want
	*/
  while ( (testPacket_ptr != ptrFailure) && (getPacketId(testPacket_ptr) != packetId) )
    testPacket_ptr = findNextFramePacket (frame_ptr, testPacket_ptr);

  return testPacket_ptr;
}

/*
**  Private routine to search for a packet in a unbroken data block of packets
*/
PACKET_ptr findNextFramePacket (FRAME_ptr frame_ptr, PACKET_ptr currentPacket_ptr)
{
  PACKET_ptr nextPacket_ptr;

  /*
  ** Figure out where the frame ends
  */
  PHDWORD* frameDataStart_ptr = findFrameDataStart(frame_ptr);
  PHDWORD* frameDataEnd_ptr = findFrameDataEnd(frame_ptr);
	
  /*
  **  Check for valid current packet pointer
  */
  if ((currentPacket_ptr < frameDataStart_ptr) || (currentPacket_ptr > frameDataEnd_ptr)) {
    return ptrFailure;
  }

  /*
  **  Now skip past "current" packet
  */
  nextPacket_ptr = findPacketEnd(currentPacket_ptr) + 1;

  if (nextPacket_ptr > frameDataEnd_ptr) return ptrFailure;
  else return nextPacket_ptr;
}


/*
** routine to find the number of packets in a frame
*/
VALUE_ret frameNumPackets (FRAME_ptr frame_ptr)
{
  UINT numPackets = 1;
  PACKET_ptr packet_ptr = findFirstFramePacket(frame_ptr);
  if (packet_ptr == ptrFailure) return 0;
  while ( (packet_ptr = findNextFramePacket(frame_ptr,packet_ptr)) != ptrFailure )
    numPackets++;
  return numPackets;
}
