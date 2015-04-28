/* 
** packetHdrV1.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:16 $ 
** 
** $Log: packetHdrV1.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:16  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:46  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**   packetHdrV1.h
**
**   Include file which defines the header for Version 1 packets 
**   and some constants used in accessing, interpreting or storing 
**   data in the header.  It also provides default settings for the 
**   fields manipulated in the Version 1 packet headers. Version 
**   independent constants and defaults are defined in packetHdr.h and 
**   constants and defaults used in the data descriptor fields are
**   set in dataBlockHdr.h.  
*/


#ifndef _PACKETHDRV1__
#define _PACKETHDRV1__

/*
**  Use C linkage
*/
#ifdef __cplusplus
extern "C" {
#endif	

#define PACKET_LENGTH_OFFSET_OF_DWORD 0

#define PACKET_STATUS_OFFSET_OF_DWORD 1 
#define PACKET_STATUS_OFFSET_IN_DWORD 0
#define PACKET_STATUS_NUM_BITS 16
#define PACKET_STATUS_MASK 0x0000ffff

#define ID_OFFSET_OF_DWORD 2

#define DEBUG_LENGTH_OFFSET_OF_DWORD 3
#define DEBUG_LENGTH_OFFSET_IN_DWORD 16
#define DEBUG_LENGTH_NUM_BITS 16
#define DEBUG_LENGTH_MASK 0xffff0000

#define PACKET_ERROR_LENGTH_OFFSET_OF_DWORD 3
#define PACKET_ERROR_LENGTH_OFFSET_IN_DWORD 0
#define PACKET_ERROR_LENGTH_NUM_BITS 16
#define PACKET_ERROR_LENGTH_MASK 0x0000ffff

#define STRUCTURE_OFFSET_OF_DWORD 4
#define STRUCTURE_OFFSET_IN_DWORD 24 
#define STRUCTURE_NUM_BITS 8
#define STRUCTURE_MASK 0xff000000

#define DESCR_LENGTH_OFFSET_OF_DWORD 4
#define DESCR_LENGTH_OFFSET_IN_DWORD 16
#define DESCR_LENGTH_NUM_BITS 8
#define DESCR_LENGTH_MASK 0x00ff0000

#define ENDIANISM_OFFSET_OF_DWORD 4
#define ENDIANISM_OFFSET_IN_DWORD 8
#define ENDIANISM_NUM_BITS 8
#define ENDIANISM_MASK 0x0000ff00

#define PACKET_PADDING_OFFSET_OF_DWORD 4
#define PACKET_PADDING_OFFSET_IN_DWORD 0
#define PACKET_PADDING_NUM_BITS 8
#define PACKET_PADDING_MASK 0x000000ff
 
#define DATADESCR_OFFSET_OF_DWORD 5

#define PACKETV1_HDR_LENGTH 6
  const UINT packetV1HdrLength = PACKETV1_HDR_LENGTH;

#define PACKETV1_HDR_VERSION 1
  const UINT packetV1HdrVersion = PACKETV1_HDR_VERSION;

#ifdef __cplusplus
} 
/* end of extern "C" */
#endif

#endif
/* end the ifndef _PACKETHDR_ block */








