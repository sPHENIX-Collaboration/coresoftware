/* 
** fakeFrame.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:12 $ 
** 
** $Log: fakeFrame.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:12  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.4  1999/09/29 22:16:07  steinber
** mods to bring afs to nevis1 version
**
** Revision 1.3  1998/12/11 22:01:17  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
#ifndef _FAKEFRAME_
#define _FAKEFRAME_

/*
**  Use C linkage for below structures
*/
#ifdef __cplusplus
extern "C" {
#endif

#include "phenixOnline.h"
#include "framePublic.h"


int fakeFrame (FRAME_ptr,int,UINT,ALIGNBLK,int,int*,int*,int,int);
int splitFakeFrame (FRAME_ptr,int,int,PHDWORD* []);

#ifdef __cplusplus
}
#endif

#endif //_FAKEFRAME_



