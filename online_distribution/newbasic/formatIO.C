/* 
** formatIO.C
** 
** Author: $Author: phnxbld $  
**   Date: $Date: 2009/08/19 12:31:53 $ 
** 
** $Log: formatIO.C,v $
** Revision 1.4  2009/08/19 12:31:53  phnxbld
** fix compiler warning
**
** Revision 1.3  2005/12/19 16:29:24  pinkenbu
** fix insure bad format warnings
**
** Revision 1.2  2002/05/31 22:15:42  purschke
** mlp -- went through the insure report and tried to fix
** every little problem there is, unused variables, dead statements,
** all. It'll probably take another round to complete, but it should get
** rid of most warnings.
**
** Revision 1.1.1.1  2000/07/21 01:51:13  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.4  1998/12/11 22:02:03  markacs
** (stephen markacs) adding log into cvs tags
** 
*/
/*
**   formatIO.C
**
**
**   Routines to dump frames and packets
**
*/

#include <stdio.h>
#include "phenixOnline.h"
#include "Cframe.h"
#include "framePublic.h"
#include "packetPublic.h"
#include "Cpacket.h"
#include "dataBlock.h"
#include "formatIO.h"

/*
**  dumpFrame 
*/

VALUE_ret dumpFrame (FRAME_ptr frame_ptr)
{
  PHDWORD *dump_ptr;
  UINT  nDwords;
  VALUE_ret result;

  /*
  ** Dump the header
  */
  result = dumpFrameHdr (frame_ptr);
  if (result == valueFailure)
    {
      printf ("Unknown or invalid Frame header \n");
      return valueFailure;
    }

  /*
  **  Get number of Dwords in frame contents
  */
  nDwords = getFrameDataLength (frame_ptr);
  dump_ptr = findFrameDataStart (frame_ptr);	

  /*
  **  For now just dump data + history in Hex format
  */
  printf ("Dump of frame data: \n");

  {
    /*
    ** Loop through entire frame dumping each word
    */
    UINT iDword;

    for (iDword = 0; iDword < nDwords; iDword++)
      printf ("Frame Word %u = %#8x \n", iDword, *dump_ptr++);
  }

  return 0;
}

VALUE_ret dumpFramePackets (FRAME_ptr frame_ptr)
{
  PHDWORD *dump_ptr;
  VALUE_ret result;

  /*
  ** Dump the frame header
  */
  result = dumpFrameHdr (frame_ptr);
  if (result == valueFailure) 
    {
      printf ("Unknown or invalid Frame header \n");
      return valueFailure;
    }


  /*
  **  Now loop through the frame dumping the packets
  */
  dump_ptr = findFirstFramePacket(frame_ptr);
  while (dump_ptr != ptrFailure)
    {
      VALUE_ret result;

      /*
      **  Dump the packet.
      */
      result = dumpPacket((PACKET_ptr) dump_ptr);
      if (result == valueFailure)
	{
	  printf ("Error in frame data -- unknown or invalid packet \n");
	  return valueFailure;
	}

      dump_ptr = findNextFramePacket(frame_ptr, dump_ptr);
    }

  return 0;
}


VALUE_ret dumpFrameHdr (FRAME_ptr frame_ptr)
{
  /*
  **  Check for valid frame header and correct version
  */
  if (!validFrameHdr (frame_ptr))
    {
      /*
      **  We have received an "invalid" frame, report error
      */
      return valueFailure;
    }

  /*
  ** Now that we are satisfied that we have a valid header
  **	 define a proper pointer to the current header type
  */

  printf ("Frame Mark = %#.8x \n", getFrameMark(frame_ptr));

  printf ("Frame version = %u, Frame Hdr length = %u \n", 
	  getFrameHdrVersion(frame_ptr), getFrameHdrLength(frame_ptr));


  printf ("data type = %u, Frame type = %u, Source Id = %u \n", 
	  getFrameDataType(frame_ptr), getFrameType(frame_ptr), 
	  getFrameSourceId(frame_ptr));

  printf ("Frame length = %u, Error Length = %u, History Length = %u \n", 
	  getFrameLength(frame_ptr), getFrameErrorLength(frame_ptr), 
	  getFrameHistoryLength(frame_ptr));

  printf ("Frame status = %#.4x, Frame padding = %u, Align length = %u \n", 
	  getFrameStatus(frame_ptr), getFramePadding(frame_ptr), 
	  getFrameAlignLength(frame_ptr));
  {
    UINT iAlign;
    /*
    ** Dump the alignment block
    */

    PHDWORD* align_ptr = findFrameAlignBlock(frame_ptr);

    for (iAlign = 0; iAlign < getFrameAlignLength(frame_ptr); iAlign++)
      printf ("Alignment block word %u = %#.8x \n", iAlign, *align_ptr++);
  }

  return 0;
}

/*
**  Routine to dump the header and contents of a packet
*/
VALUE_ret dumpPacket (PACKET_ptr packet_ptr)
{
  /*
  **  Check for valid packet header and correct version
  */  
  if (!validPacketHdr (packet_ptr)) {
    /*
    **  We have received an "invalid" packet, report error
    */
    return valueFailure;
  }
  
  /*
  ** Dump the header
  */
  printf ("Packet Length: %u", getPacketLength(packet_ptr));

  printf ("Packet HDR version: %u, Packet HDR length:%u",
	  getPacketHdrVersion(packet_ptr), getPacketHdrLength(packet_ptr));

  printf ("Packet Status bits: %#.4x \n", getPacketStatus(packet_ptr));

  printf ("Packet Ident: %u", getPacketId(packet_ptr));

  printf ("Endianism: %u, Padding: %u \n", getPacketEndianism(packet_ptr), 
	  getPacketPadding(packet_ptr));
  
  printf ("Error Length: %u, Debug Length: %u \n", 
	  getPacketErrorLength(packet_ptr), getPacketDebugLength(packet_ptr));
  
  /*
  ** Print-out structure-specific information
  */
  switch (getPacketStructure(packet_ptr)) {
  case Unstructured: 
    {
      UINT wordSize = getUnstructPacketWordSize(packet_ptr);
      PHDWORD numWords = getUnstructPacketDataLengthWords(packet_ptr);

      printf ("Unstructured packet, Format: %u, Word Size: %u, Length (words): %u",
	      getUnstructPacketHitFormat(packet_ptr), wordSize, numWords);
	      
      /*
      **  Now dump the contents of the packet
      */
      {
        UINT iWord;
	PHDWORD dataWord = 0; // init not needed but it suppresses compiler warning
        BYTE *dump_ptr;

	printf ("Dump of packet data: \n");
	/*
	**  Do a word-by-word dump of the packet
	*/
	
	dump_ptr = (BYTE*) findPacketDataStart(packet_ptr);
	
	for (iWord = 0; iWord < numWords; iWord++) {
	  switch(wordSize) 
	    {
	    case 1: { dataWord = *dump_ptr; break; }
	    case 2: {dataWord = *((SWORD*) dump_ptr); break;}
	    case 4: {dataWord = *((PHDWORD*) dump_ptr); break;}
	    default: break; 
	    }
	  printf ("Packet Word %u = %#8x \n", iWord, dataWord);
	  dump_ptr+=wordSize;
	}
      }
    }
  }

  /*
  ** Return success
  */
  return 0;
}


