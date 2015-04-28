/* 
** packetPublic.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:17 $ 
** 
** $Log: packetPublic.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:17  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:47  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
** packetPublic.h
**
**   This file contains enumerations, typedefs etc. that public users of 
**   packet C and C++ routines might need. 
**
*/

#include "phenixOnline.h"
#include "packetV1Public.h"
#include "errorBlock.h"

#ifndef _PACKETPUBLIC_
#define _PACKETPUBLIC_

/*
**  Use C linkage in C++ code
*/
#ifdef __cplusplus
extern "C" {
#endif

  /*
  **  We point to packets using a PHDWORD pointer.
  */
  typedef PHDWORD* PACKET_ptr ;

  /*
  **  Typedef for the error entry structure which users may interact with
  **  directly.
  **
  **  *** currently the packet and frame errors share the same structure 
  **      defined in errorBlock.h
  */
  typedef struct errorEntryV1 PACKETERROR;
  
  /*
  **  How many different packet header versions there are
  */
#define NUM_PACKET_VERSIONS 2
#define	CURRENT_PACKETHDR_VERSION 1

  CONSTANT UINT numPacketVersions = NUM_PACKET_VERSIONS;
  CONSTANT UINT currentPacketHdrVersion = CURRENT_PACKETHDR_VERSION;
  
  /*
  **  Enumerations for packet structure types.
  */
  enum PacketStructEnum {
    defaultStructure = 0,
    Unstructured = 0,
    Packets = 1,
    hitArray = 2,
    hitList = 3
  } ;
  typedef enum PacketStructEnum PACKETSTRUCTENUM;
  
  /*
  **  Status code enumerations
  */
  enum PacketStatusBits {
    errorBlockInvalid = 0
  } ;
  
#ifdef __cplusplus
}  /* end extern C block */
#endif

#endif 
/* end ifdef _PACKETPUBLIC_ */
