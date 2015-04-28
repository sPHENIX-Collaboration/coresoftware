/* 
** Cpacket.C
** 
** Author: $Author: phoncs $  
**   Date: $Date: 2004/02/14 00:17:03 $ 
** 
** $Log: Cpacket.C,v $
** Revision 1.2  2004/02/14 00:17:03  phoncs
** J.J change by brian to fix the padding problem
**
** Revision 1.1.1.1  2000/07/21 01:51:10  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.6  1998/12/17 15:18:16  phoncs
** (stephen markacs) more bounds checking in Cpacket.C (one instance) and removal of debugging comments from checkFrame.C
**
** Revision 1.5  1998/12/11 22:01:59  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**  These are the version independent functions which 
**  provide access to and mocification of information
**  in the header, but which are not inlined.
**
*/

#include "phenixOnline.h"
#include "packetPublic.h"
#include "packetHdr.h"
#include "Cpacket.h"
#include "CpacketV1.h"

PTR_ret findPacketEnd (PACKET_ptr packet_ptr)
{
#ifdef STRONGCHECK
  if (!validPacketHdr(packet_ptr)) {
    setPacketError(FORMAT_ERR_INVALID_HDR, packet_ptr, 0);
    return logicFailure;
  }
#endif 

  PHDWORD packetLength = getPacketLength (packet_ptr);
  if (packetLength != valueFailure) {
    setPacketSuccess ();
    return ((PHDWORD *) packet_ptr + packetLength - 1);
  }
  else 
    return ptrFailure;
}

PTR_ret findPacketDataStart (PACKET_ptr packet_ptr)
{
  if (validPacketHdr(packet_ptr))
    return ((PHDWORD*)packet_ptr + getPacketHdrLength(packet_ptr) + 
	    getPacketDataDescrLength(packet_ptr) - 1);
  else
    return ptrFailure;
}

/*
** Adjust the length of the data block in the packet. 
**
**  Note: a negative extension is allowed
*/
VALUE_ret extendPacketDataBlock(const PACKET_ptr packet_ptr, UINT maxPacketLength, 
			      int addDwords)
{
  /*
  ** Extend the packet itself
  */
  PHDWORD newLength = extendPacketLength (packet_ptr, maxPacketLength, addDwords);

  /*
  ** Update the length of the data block itself (NOOP for V1)
  */
  if (newLength != valueFailure) {
    return adjustPacketDataLength(packet_ptr, addDwords);
  }
  else
    return valueFailure;
}

/*
** Adjust the length of the debug block in the packet. 
**
**  Note: a negative extension is allowed
*/
VALUE_ret extendPacketDebugBlock(const PACKET_ptr packet_ptr,  UINT maxPacketLength, 
				 int addDwords)
{
  /*
  ** Extend the packet itself
  */
  PHDWORD newLength = extendPacketLength (packet_ptr, maxPacketLength, addDwords);

  /*
  ** Update the length of the debug blcok itself (NOOP) for V1
  */
  if (newLength != valueFailure) 
    return adjustPacketDebugLength(packet_ptr, addDwords);
  else
    return valueFailure;
}


/*
** Adjust the length of the error block in the packet. 
**
**  Note: a negative extension is allowed
*/
VALUE_ret extendPacketErrorBlock(const PACKET_ptr packet_ptr,  UINT maxPacketLength, 
				 int addDwords)
{
  /*
  ** Extend the packet itself
  */
  PHDWORD newLength = extendPacketLength (packet_ptr, maxPacketLength, addDwords);

  /*
  ** Update the length of the debug blcok itself (NOOP) for V1
  */
  if (newLength != valueFailure) 
    return adjustPacketErrorLength(packet_ptr, addDwords);
  else
    return valueFailure;
}


/*
** Adjust the length of the packet and update the padding count. 
**
**  Note: a negative extension is allowed
*/
VALUE_ret extendPacketLength(const PACKET_ptr packet_ptr,  UINT maxPacketLength, 
			     int addDwords)
{
  UINT nPadDwords;
  PHDWORD newLength;

  /*
  ** Drop any existing padding.
  */
  if (!removePacketPadding(packet_ptr)) {
    return valueFailure;
  }

  /*
  **  Check for buffer overflow
  */
  newLength = getPacketLength(packet_ptr) + addDwords;

  if (newLength > maxPacketLength) {
    setPacketError(FORMAT_ERR_BUFFER_OVERFLOW,packet_ptr,newLength);
    return valueFailure;
  }

  {
    /*
    **  Now calculate the padding
    */
    UINT modulo = newLength % packetQuantum;
    if (modulo > 0)  nPadDwords = packetQuantum - modulo;
    else nPadDwords = 0;
  }

  newLength += nPadDwords;

  /*
  ** Store the number of padding words in the header
  */
  setPacketPadding(packet_ptr, nPadDwords);

  /*
  ** Update the packet length
  */
  setPacketLength(packet_ptr, newLength);

  return newLength;
}

LOGIC_ret removePacketPadding(PACKET_ptr packet_ptr)
{
  /*
  ** Get the current number of padding dwords
  */
  UINT nPadDwords = getPacketPadding(packet_ptr);

  /*
  ** decrease the packet length by this amount
  */
  if (nPadDwords > 0) {
    PHDWORD newPacketLength = getPacketLength(packet_ptr) - nPadDwords;

    /*
    ** Check for reasonable result
    */
    if (newPacketLength < packetV1HdrLength) {
      /*
      **  Oops, there's a problem. Set error result and return failure.
      */
      setPacketError(FORMAT_ERR_HDR_INCONSISTENCY,packet_ptr,0);
      return logicFailure;
    }

    /*
    **  Update the header length
    */
    setPacketLength(packet_ptr, newPacketLength);
    setPacketPadding(packet_ptr, 0);
  }

  return logicSuccess;
}


#ifndef DCM
/*
** Define indirection pointer arrays to various frame access functions
** To avoid conflict with frame indirection pointer arrays, a "P" has
** been placed in each definition before the word "FUNCTION".
*/
typedef LOGIC_ret CHECKPFUNC (PACKET_ptr);
typedef VALUE_ret ACCESSPFUNC (PACKET_ptr);
typedef VALUE_ret ADJUSTPFUNC (PACKET_ptr, int);
typedef VALUE_ret ORPFUNC (PACKET_ptr, UINT);
typedef PTR_ret   PTRACCESSPFUNC (PACKET_ptr);
typedef LOGIC_ret MODIFYPFUNC (PACKET_ptr, UINT);


typedef ACCESSPFUNC* ACCESSPFUNCPTR_arr[NUM_PACKET_VERSIONS]; 
typedef PTRACCESSPFUNC* PTRACCESSPFUNCPTR_arr[NUM_PACKET_VERSIONS];
typedef CHECKPFUNC*  CHECKPFUNCPTR_arr[NUM_PACKET_VERSIONS];
typedef MODIFYPFUNC* MODIFYPFUNCPTR_arr[NUM_PACKET_VERSIONS]; 
typedef ADJUSTPFUNC* ADJUSTPFUNCPTR_arr[NUM_PACKET_VERSIONS];
typedef ORPFUNC* ORPFUNCPTR_arr[NUM_PACKET_VERSIONS];

#ifdef __cplusplus
ACCESSPFUNC* const NULLACCESSPFUNC = 0;
PTRACCESSPFUNC* const NULLPTRACCESSPFUNC = 0;
CHECKPFUNC*  const NULLCHECKPFUNC = 0;
MODIFYPFUNC* const NULLMODIFYPFUNC = 0;
ADJUSTPFUNC* const NULLADJUSTPFUNC = 0;
ORPFUNC* const NULLORPFUNC = 0;
#else
/*
** Unfortunately for some C compilers we have to lose the benefits of 
**   type-checking because const's are not considered valid initializers
*/
#define NULLACCESSPFUNC 0
#define NULLPTRACCESSPFUNC 0
#define NULLCHECKPFUNC 0
#define NULLMODIFYPFUNC 0
#define NULLADJUSTPFUNC 0
#define NULLADJUSTDATAPFUNC 0
#define NULLORPFUNC 0
#endif


ACCESSPFUNCPTR_arr CONSTANT getPacketDataLengthV = {NULLACCESSPFUNC, &getPacketV1DataLength};
ACCESSPFUNCPTR_arr CONSTANT getPacketDebugLengthV = {NULLACCESSPFUNC, &getPacketV1DebugLength};
ACCESSPFUNCPTR_arr CONSTANT getPacketErrorLengthV = {NULLACCESSPFUNC, &getPacketV1ErrorLength};
ACCESSPFUNCPTR_arr CONSTANT getPacketNumErrorsV = {NULLACCESSPFUNC, &getPacketV1NumErrors};
ACCESSPFUNCPTR_arr CONSTANT getPacketEndianismV = {NULLACCESSPFUNC, &getPacketV1Endianism};
ACCESSPFUNCPTR_arr CONSTANT getPacketIdV = {NULLACCESSPFUNC, &getPacketV1Id};
ACCESSPFUNCPTR_arr CONSTANT getPacketPaddingV = {NULLACCESSPFUNC, &getPacketV1Padding};
ACCESSPFUNCPTR_arr CONSTANT getPacketStructureV = {NULLACCESSPFUNC, &getPacketV1Structure};
ACCESSPFUNCPTR_arr CONSTANT getPacketStatusV = {NULLACCESSPFUNC, &getPacketV1Status};
ACCESSPFUNCPTR_arr CONSTANT getPacketDataDescrLengthV = {NULLACCESSPFUNC, &getPacketV1DataDescrLength};

ORPFUNCPTR_arr CONSTANT orPacketStatusV = {NULLORPFUNC, &orPacketV1Status};
ADJUSTPFUNCPTR_arr CONSTANT adjustPacketDataLengthV = {NULLADJUSTPFUNC, 
						       &adjustPacketV1DataLength};
ADJUSTPFUNCPTR_arr CONSTANT adjustPacketDebugLengthV = {NULLADJUSTPFUNC, 
							&adjustPacketV1DebugLength};
ADJUSTPFUNCPTR_arr CONSTANT adjustPacketErrorLengthV = {NULLADJUSTPFUNC, 
							&adjustPacketV1ErrorLength};

ORPFUNCPTR_arr CONSTANT makePacketHdrV = {NULLORPFUNC, &makePacketV1Hdr};

MODIFYPFUNCPTR_arr CONSTANT setPacketDebugLengthV = {NULLMODIFYPFUNC, &setPacketV1DebugLength};
MODIFYPFUNCPTR_arr CONSTANT setPacketDataDescrLengthV = {NULLMODIFYPFUNC, 
							 &setPacketV1DataDescrLength};

MODIFYPFUNCPTR_arr CONSTANT setPacketEndianismV = {NULLMODIFYPFUNC, &setPacketV1Endianism};
MODIFYPFUNCPTR_arr CONSTANT setPacketIdV = {NULLMODIFYPFUNC, &setPacketV1Id};
MODIFYPFUNCPTR_arr CONSTANT setPacketPaddingV = {NULLMODIFYPFUNC, &setPacketV1Padding};
MODIFYPFUNCPTR_arr CONSTANT setPacketErrorLengthV = {NULLMODIFYPFUNC, &setPacketV1ErrorLength};
MODIFYPFUNCPTR_arr CONSTANT setPacketStructureV = {NULLMODIFYPFUNC, &setPacketV1Structure};
MODIFYPFUNCPTR_arr CONSTANT setPacketStatusV = {NULLMODIFYPFUNC, &setPacketV1Status};

PTRACCESSPFUNCPTR_arr CONSTANT findPacketErrorStartV = {NULLPTRACCESSPFUNC, 
							&findPacketV1ErrorStart};
PTRACCESSPFUNCPTR_arr CONSTANT findPacketDebugStartV = {NULLPTRACCESSPFUNC, 
							&findPacketV1DebugStart};
PTRACCESSPFUNCPTR_arr CONSTANT findPacketDataStartV = {NULLPTRACCESSPFUNC, 
						       &findPacketV1DataStart};
PTRACCESSPFUNCPTR_arr CONSTANT findPacketDataEndV = {NULLPTRACCESSPFUNC, &findPacketV1DataEnd};
PTRACCESSPFUNCPTR_arr CONSTANT findPacketDataDescrV = {NULLPTRACCESSPFUNC, 
						       &findPacketV1DataDescr};

CHECKPFUNCPTR_arr  CONSTANT validPacketHdrV = {NULLCHECKPFUNC, &validPacketV1Hdr};
CHECKPFUNCPTR_arr  CONSTANT emptyPacketV = {NULLCHECKPFUNC, &emptyPacketV1};

/*
**   "indirection" routines.
**
*/

VALUE_ret makePacketHdr (const PACKET_ptr packet_ptr, UINT maxPacketLength)
{
  UINT version = currentPacketHdrVersion;
  if (version != valueFailure)
    return (*makePacketHdrV[version])(packet_ptr, maxPacketLength);
  else
    return valueFailure;
}


VALUE_ret getPacketStatus (PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion(packet_ptr);
  if (version != valueFailure)
    return (*getPacketStatusV[version])(packet_ptr);
  else
    return valueFailure;
}

VALUE_ret orPacketStatus (PACKET_ptr packet_ptr, UINT inBits)
{
  UINT version = getPacketHdrVersion(packet_ptr);
  if (version != valueFailure)
    return (*orPacketStatusV[version])(packet_ptr, inBits);
  else
    return valueFailure;
}

VALUE_ret getPacketDataLength (const PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return getPacketDataLengthV[version](packet_ptr);
  else

    return valueFailure;
}

VALUE_ret getPacketErrorLength (const PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return getPacketErrorLengthV[version](packet_ptr);
  else

    return valueFailure;
}

VALUE_ret getPacketDebugLength (const PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return getPacketDebugLengthV[version](packet_ptr);
  else
    return valueFailure;
} 

VALUE_ret getPacketNumErrors (const PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return getPacketNumErrorsV[version](packet_ptr);
  else
    return valueFailure;
}

VALUE_ret getPacketId (const PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return getPacketIdV[version](packet_ptr);
  else
    return valueFailure;
} 

VALUE_ret getPacketStructure (const PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return getPacketStructureV[version](packet_ptr);
  else return valueFailure;
} 

VALUE_ret getPacketEndianism (const PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return getPacketEndianismV[version](packet_ptr);
  else
    return valueFailure;
}

VALUE_ret getPacketPadding (const PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return getPacketPaddingV[version](packet_ptr);
  else
    return valueFailure;
}

VALUE_ret getPacketDataDescrLength (const PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return getPacketDataDescrLengthV[version](packet_ptr);
  else
    return valueFailure;
}


LOGIC_ret setPacketDebugLenth(PACKET_ptr packet_ptr, UINT inLength)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return setPacketDebugLengthV[version](packet_ptr, inLength);
  else
    return logicFailure;
}

/*
**  This was commented out up top.
**
**
**  LOGIC_ret setPacketNumErrors(PACKET_ptr packet_ptr, UINT inLength)
**  {
**    UINT version = getPacketHdrVersion (packet_ptr);
**    if (version != valueFailure)
**      return setPacketNumErrorsV[version](packet_ptr, numErrors);
**    else
**      return logicFailure;
**  }
*/

LOGIC_ret setPacketId(PACKET_ptr packet_ptr, UINT inId)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return setPacketIdV[version](packet_ptr, inId);
  else
    return valueFailure;
}

LOGIC_ret setPacketStructure(PACKET_ptr packet_ptr, UINT structure)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return setPacketStructureV[version](packet_ptr, structure);
  else
    return logicFailure;
}

LOGIC_ret setPacketEndianism(PACKET_ptr packet_ptr, UINT endianism)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return setPacketEndianismV[version](packet_ptr, endianism);
  else
    return logicFailure;
}

LOGIC_ret setPacketPadding(PACKET_ptr packet_ptr, UINT padding)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return setPacketPaddingV[version](packet_ptr, padding);
  else
    return logicFailure;
}

LOGIC_ret setPacketStatus(PACKET_ptr packet_ptr, UINT status)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return setPacketStatusV[version](packet_ptr, status);
  else
    return logicFailure;
}

LOGIC_ret setPacketErrorLength(PACKET_ptr packet_ptr, UINT errorLength)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return setPacketErrorLengthV[version](packet_ptr, errorLength);
  else
    return logicFailure;
}

LOGIC_ret setPacketDataDescrLength(PACKET_ptr packet_ptr, UINT descrLength)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return setPacketDataDescrLengthV[version](packet_ptr, descrLength);
  else
    return logicFailure;
}

PTR_ret findPacketDataEnd (PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return findPacketDataEndV[version](packet_ptr);
  else
    return ptrFailure;
}

PTR_ret findPacketDebugStart (PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return findPacketDebugStartV[version](packet_ptr);
  else
    return ptrFailure;
}

PTR_ret findPacketErrorStart (PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return (*findPacketErrorStartV[version])(packet_ptr);
  else
    return ptrFailure;
}

PTR_ret findPacketDataDescr (PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return (*findPacketDataDescrV[version])(packet_ptr);
  else
    return ptrFailure;
}


LOGIC_ret validPacketHdr (PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if ( (version > 0) && (version < NUM_PACKET_VERSIONS) )
    return validPacketHdrV[version](packet_ptr);
  else
    return logicFailure;
}

LOGIC_ret emptyPacket (PACKET_ptr packet_ptr)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return emptyPacketV[version](packet_ptr);
  else
    return logicFailure;
}

VALUE_ret adjustPacketDataLength(PACKET_ptr packet_ptr, int addDwords)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return adjustPacketDataLengthV[version](packet_ptr, addDwords);
  else
    return valueFailure;
}

VALUE_ret adjustPacketDebugLength(PACKET_ptr packet_ptr, int addLength)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return adjustPacketDebugLengthV[version](packet_ptr, addLength);
  else
    return valueFailure;
}

VALUE_ret adjustPacketErrorLength(PACKET_ptr packet_ptr, int addLength)
{
  UINT version = getPacketHdrVersion (packet_ptr);
  if (version != valueFailure)
    return adjustPacketErrorLengthV[version](packet_ptr, addLength);
  else
    return valueFailure;
}

#endif










