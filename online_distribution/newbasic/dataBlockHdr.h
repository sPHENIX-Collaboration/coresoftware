/* 
** dataBlockHdr.h
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:11 $ 
** 
** $Log: dataBlockHdr.h,v $
** Revision 1.1.1.1  2000/07/21 01:51:11  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.3  1998/12/11 22:01:15  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**
**  dataBlockHdr.h
**
**  This file defines the descriptor for the data block and some constants
**  used for interpreting, accessing or storing data in the descriptor.
**  It also defines the default setting for all of the data descriptor
**  fields, even those set by the user.
**/


#ifndef _DATABLOCKHDR_
#define _DATABLOCKHDR_

/*
**  Use C linkage
*/
#ifdef __cplusplus
extern "C" {
#endif

  typedef
    struct dataBlockHdr {

      BYTE		wordSize		; /* "Word" size used to store packet data   */
      BYTE		addrLength		; /* number of bytes used for channel address*/
      BYTE		hitLength		; /* Length of a single "hit" in bytes       */
      BYTE              dataPadding             ; /* Padding of the data block               */
	
      SWORD		hitFormat		; /* Format of a single hit                  */
      SWORD		numEntries		; /* Number of "objects" stored in packet    *
	      				   * (definition depends on packetStruct)    */
    } DATABLOCKHDR, * DATABLOCKHDR_ptr;
  /*
  **   These are the default settings for all of the fields 
  **   in the data descriptor.  Some of the fields are user specified,
  **   but still have default values listed here.
  */

#define CURRENT_NUM_ENTRIES 0 
#define CURRENT_WORD_SIZE 0
#define CURRENT_HIT_LENGTH 0
#define CURRENT_HIT_FORMAT 0
#define CURRENT_ADDR_LENGTH 0 
#define CURRENT_DATA_PADDING 0
#define CURRENT_DESC_LENGTH 2
#define CURRENT_UNSTRUCT_DESC_LENGTH 1

CONSTANT UINT currentNumEntries = CURRENT_NUM_ENTRIES; 
CONSTANT UINT currentWordSize = CURRENT_WORD_SIZE; 
CONSTANT UINT currentHitLength = CURRENT_HIT_LENGTH;
CONSTANT UINT currentHitFormat = CURRENT_HIT_FORMAT;
CONSTANT UINT currentAddrLength = CURRENT_ADDR_LENGTH;
CONSTANT UINT currentDataPadding = CURRENT_DATA_PADDING;
CONSTANT UINT currentDescLength = CURRENT_DESC_LENGTH;
CONSTANT UINT currentUnstructDescLength = CURRENT_UNSTRUCT_DESC_LENGTH;
#ifdef __cplusplus
}
/* end of extern "C" */
#endif

#endif
/* end the ifndef in _DATABLOCKHDR_ block */
