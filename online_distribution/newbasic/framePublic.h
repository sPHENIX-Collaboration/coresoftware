/* 
** framePublic.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:13 $ 
** 
** $Log: framePublic.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:13  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:42  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
** framePublic.h
**
**   Defines "public" enumerations, typedefs, etc. used in both
**   C and C++ frame routines.
**
*/

#ifndef _FRAMEPUBLIC_ 
#define _FRAMEPUBLIC_ 

#include "phenixOnline.h"
#include "frameV1Public.h"

/*
**  Use C linkage for below structures
*/
#ifdef __cplusplus
extern "C" {
#endif

  /*
  **  Here's a current failure in the frame versioning scheme as implemented
  **  Right now the alignment block is treated as being defined with a particular
  **  version. That means the contents are version specific. Until we figure out
  **  how to handle this problem, explicitly define the alignment block to be from V1
  */
  typedef ALIGNBLKV1 ALIGNBLK;

  /*
  ** Define the type of pointer used for a frame. For now a frame
  **   will be pointed to as if it is a PHDWORD array. When headers
  **   are accessed the pointer to the frame will be cast to that
  **   of the proper header type.
  */
  
  typedef PHDWORD* FRAME_ptr;
  
  /*
  **  Enumeration for the "dataType" field
  */
  enum DataTypes {
    normalData = 1,
    calibrationData = 2,
    monitorData = 3,
    rawData = 4
    /*
    ** more to be added
    */
  } ;
  
  /*
  **  Enumeration for the "frameType" field
  */
  enum FrameTypes {
    dcmFrame  = 1,
    dcbFrame = 2,
    sebFrame = 3,
    atpFrame = 4,
    oncsFrame = 5
    /*
    ** more to be added
    */
  } ;
  
#ifdef __cplusplus
}
#endif

#endif








