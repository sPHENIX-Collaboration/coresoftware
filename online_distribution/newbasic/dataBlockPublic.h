/* 
** dataBlockPublic.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:11 $ 
** 
** $Log: dataBlockPublic.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:11  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:16  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
** dataBlockPublic.h
**
**   This file contains enumerations, typedefs etc. that public users of 
**   data descriptor C and C++ routines might need. 
**
*/

#include "phenixOnline.h"
#include "dataBlockHdr.h"
#include "packetPublic.h"

#ifndef _DATABLOCKPUBLIC_
#define _DATABLOCKPUBLIC_

/*
**  Use C linkage in C++ code
*/
#ifdef __cplusplus
extern "C" {
#endif

  /*
  **  We point to the data descriptor using a PHDWORD pointer.
  */
  typedef PHDWORD* DATABLOCK_ptr;
  
  /*
  **  The error entry structure for the data descriptor should be the 
  **  same as for the packet routines.  The typedef is given in 
  **  packetPublic.h and the error waring is "PACKETERROR", rather 
  **  than a seperate command for the descriptor.
  */

  /*
  **  Presumably, when the various sorts of data sturctures are 
  **  fully enumerated, there will be a set of enumerations corresponding
  **  to the various types of data descriptors used to describe
  **  those structures.  These should follow the packet structure types 
  **  enumerated in packetPublic.h.  As of now there is only the one dword 
  **  unstructured type of descriptor and a two dword type which has all 
  **  possible fields set. The latter of these will probably never be 
  **  actually used.
  */

#ifdef __cplusplus
}  /* end extern C block */
#endif

#endif
/* end ifdef _DATABLOCKPUBLIC_ */





