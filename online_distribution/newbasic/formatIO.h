/* 
** formatIO.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:13 $ 
** 
** $Log: formatIO.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:13  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:19  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**	formatIO.h
**
**
**		Defines function prototypes and useful macros for 
**		dumping packets and frames.
**
*/

#ifndef _FORMATIO_
#define _FORMATIO_

#include <stdio.h>
#include "Cpacket.h"
#include "framePackets.h"

/*
**  Use C linkage for below structures
*/
#ifdef __cplusplus
extern "C" {
#endif

#include "phenixOnline.h"
#include "framePublic.h"
#include "packetPublic.h"

/*
**  Function prototypes
*/

VALUE_ret dumpFrameHdr (FRAME_ptr);

VALUE_ret dumpFrame (FRAME_ptr);

VALUE_ret dumpFramePackets (FRAME_ptr);

VALUE_ret dumpPacket (PACKET_ptr);

#ifdef __cplusplus
}
/* End of "extern C" */
#endif

#endif






