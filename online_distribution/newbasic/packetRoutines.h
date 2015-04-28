/* 
** packetRoutines.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:17 $ 
** 
** $Log: packetRoutines.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:17  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.4  1998/12/11 22:01:48  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**
**	packetRoutines.h
**
**
**		Function prototypes for routines that perform more sophisticated
**		manipulation of packets than provided for by routines in packets.cpp.
*/

#ifndef _PACKET_ROUTINES_
#define _PACKET_ROUTINES_

#include "phenixOnline.h"
#include "packetPublic.h"

/*
**  Use C linkage for below structures
*/
#ifdef __cplusplus
extern "C" {
#endif

PTR_ret makeEmptyPacket (PACKET_ptr, UINT, UINT);

PTR_ret makeUnstructPacket (PACKET_ptr, UINT, UINT, UINT, UINT);

LOGIC_ret appendPacketError (PACKET_ptr, UINT, ERRORENTRYV1_ptr);

PTR_ret reservePacketDebugData (PACKET_ptr, UINT, UINT);

PTR_ret startUnstructDataWrite (PACKET_ptr, UINT, PHDWORD);

PTR_ret finishUnstructDataWrite (PACKET_ptr, UINT, PHDWORD);

VALUE_ret storePacketHits (PACKET_ptr, UINT, UINT*, BYTE*, UINT, UINT);

VALUE_ret fetchPacketHits (PACKET_ptr, UINT**, BYTE**, UINT*);

#ifdef __cplusplus
} 
#endif
/* end of ifdef __cplusplus */

#endif
/* end of ifndef _PACKET_ROUTINES_ */














