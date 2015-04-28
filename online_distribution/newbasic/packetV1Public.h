/* 
** packetV1Public.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:17 $ 
** 
** $Log: packetV1Public.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:17  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:49  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
** packetV1Public.h
**
**   This file contains enumerations, typedefs etc. that public 
**   users of packet C and C++ routines might need. 
**
*/

#include "phenixOnline.h"

#ifndef _PACKETV1PUBLIC_
#define _PACKETV1PUBLIC_

/*
**  Use C linkage in C++ code
*/
#ifdef __cplusplus
extern "C" {
#endif

  /*
  **  We point to packets using a PHDWORD pointer.
  */
  typedef PHDWORD* PACKETV1_ptr;
  
#ifdef __cplusplus
}  /* end extern C block */
#endif

#endif 
/* end ifdef _PACKETV1PUBLIC_ */








