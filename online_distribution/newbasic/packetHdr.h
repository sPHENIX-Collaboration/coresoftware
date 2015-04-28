/* 
** packetHdr.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:16 $ 
** 
** $Log: packetHdr.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:16  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:45  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**   packetHdr.h
**
**
**   Include file which defines the header for packets and 
**   some constants used in accessing, interpreting or 
**   storing data in the packet header.  All constants 
**   defined here are for fields with version independent 
**   values.  Also defined here are the default values for 
**   all of the version independent fields.  Version dependent
**   values for the Version 1 packet header are found in  
**   packetHdrV1.h.  Value for the fields defined in the
**   data descriptor are found in dataBlockHdr.h.  All fields
**   have default values listed, even those which are 
**   typically user specified.
**/

#include "packetHdrV1.h"
#ifndef _PACKETHDR__
#define _PACKETHDR__

/*
**  Use C linkage
*/

#ifdef __cplusplus
extern "C" {
#endif
	
#define PACKET_LENGTH_OFFSET_OF_DWORD 0

#define PACKET_HDR_VERSION_OFFSET_OF_DWORD 1
#define PACKET_HDR_VERSION_OFFSET_IN_DWORD 24
#define PACKET_HDR_VERSION_NUM_BITS 8
#define PACKET_HDR_VERSION_MASK 0xff000000
  
#define PACKET_HDR_LENGTH_OFFSET_OF_DWORD 1
#define PACKET_HDR_LENGTH_OFFSET_IN_DWORD 16
#define PACKET_HDR_LENGTH_NUM_BITS 8
#define PACKET_HDR_LENGTH_MASK 0x00ff0000

  /*
  **  Define standard packet to be "current"
  */
#define CURRENT_PACKET_VERSION 1
#define CURRENT_PACKETHDR_LENGTH PACKETV1_HDR_LENGTH

  CONSTANT UINT currentPacketHdrLength = CURRENT_PACKETHDR_LENGTH;

#define PACKETV1_QUANTUM
  const UINT packetQuantum = 2;

#ifdef __cplusplus
} 
/* end of extern "C" */
#endif

#endif
/* end the ifndef _PACKETHDR_ block */












