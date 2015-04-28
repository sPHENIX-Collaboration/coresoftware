/* 
** framePackets.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:13 $ 
** 
** $Log: framePackets.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:13  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.4  1998/12/11 22:01:36  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  Routines that work with both frames and packets
**
**
*/
#ifndef _FRAME_PACKETS_
#define _FRAME_PACKETS_

#include "phenixOnline.h"
#include "framePublic.h"
#include "packetPublic.h"

/*
**  Use C linkage for below structures
*/
#ifdef __cplusplus
extern "C" {
#endif
  
  PTR_ret appendEmptyFramePacket (FRAME_ptr, PHDWORD, UINT);
  
  PTR_ret findFirstFramePacket (FRAME_ptr);
  PTR_ret findLastFramePacket (FRAME_ptr);
  LOGIC_ret isLastFramePacket (FRAME_ptr, PACKET_ptr);
  PTR_ret findFramePacketIndex (FRAME_ptr, UINT);
  PTR_ret findFramePacketId (FRAME_ptr, UINT);
  PTR_ret findNextFramePacketId (FRAME_ptr, PACKET_ptr, UINT);
  PTR_ret findNextFramePacket (FRAME_ptr, PACKET_ptr);

  VALUE_ret frameNumPackets (FRAME_ptr);

#ifdef __cplusplus
					 }
#endif

#endif




