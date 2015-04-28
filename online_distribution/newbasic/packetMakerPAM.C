/* 
** packetMakerPAM.C
** 
** Author: $Author: purschke $  
**   Date: $Date: 2000/07/21 01:51:17 $ 
** 
** $Log: packetMakerPAM.C,v $
** Revision 1.1.1.1  2000/07/21 01:51:17  purschke
** mlp -- adding the new automakified "basic" module to CVS.
**
**
** Revision 1.5  1998/12/11 22:02:06  markacs
** (stephen markacs) adding log into cvs tags
** 
*/

#include "packetMakerPAM.h"

//
//  Source code for packetMakerPAM
//

//
// n is an upper limit on the number of tables to be used in
// a single event
//
PacketMakerTableList::PacketMakerTableList(int ntables)
  : maxNumberOfTables(ntables), numberOfTables(-1), totalLength(0)
{
  //  printf("Ntables = %d\n",ntables);
  packetList = new PacketDesc[ntables];
  for(int i=0;i<ntables;i++)
    {
      //      printf("PacketMakerTableList: Initializing packet %d\n",i);
      packetList[i].Location = NULL;
      packetList[i].Format = 0;
      packetList[i].Length = 0;
      packetList[i].Id = 0;
    }
};

PacketMakerTableList::~PacketMakerTableList()
{
  delete packetList;
}

PacketMakerTableList::addTableToList(PHDWORD* location, PHDWORD length, 
				     PHDWORD format, PHDWORD id)
{
  if (++numberOfTables < maxNumberOfTables)
    {
      packetList[numberOfTables].Location=location;
      packetList[numberOfTables].Length=length;
      packetList[numberOfTables].Format=format;
      packetList[numberOfTables].Id=id;
      totalLength += length;
      return(0);
    }
  else
    {
      printf("packetMakerTableList: Cannot add more tables!\n");
      return(-1);
    }
}

void PacketMakerTableList::dump()
{
  printf("PacketMakerTableList(%x)::dump():\n",this); 
  printf("Table#  Location  Length  Format    Id\n");
  printf("--------------------------------------\n");
  for (int i = 0;i<=numberOfTables;i++)
    {
      printf("%6d  %x  %6d  %6d  %4d\n",
	     i,
	     packetList[i].Location,
	     packetList[i].Length,
	     packetList[i].Format,
	     packetList[i].Id);
    }
}

PacketMaker::PacketMaker(int frameId, 
			 int eventNumber,
			 PHDWORD partitionVector,
			 PacketMakerTableList* pmtl)
{

  packetMakerTableList = pmtl;
  
  PHDWORD addr_arr[100000];
  
  // Estimate total length
  
  PHDWORD estimatedLength = 
    currentFrameHdrLength +
    packetMakerTableList->getTotalLength() +
    packetMakerTableList->getNumberOfTables() * currentPacketHdrLength;

  PHDWORD extraLength = (PHDWORD) ((float) estimatedLength * 1.1);

  estimatedLength += extraLength;  // give 10% more space for now...

  PHDWORD totalFrameSize = 0;
  
  frameBuffer = new PHDWORD[estimatedLength];
  dwordClear(frameBuffer,estimatedLength);

  PHDWORD numberOfPackets = packetMakerTableList->getMaxNumberOfTables();

  FRAME_ptr frame_ptr = frameBuffer;

  PHDWORD maxFrameLen = estimatedLength;
  ALIGNBLK alignBlk;
  
  alignBlk.atp.globalEventNum = eventNumber;
  alignBlk.atp.partitionVec = partitionVector;

  PHDWORD err = makeFrameHdr(frame_ptr,maxFrameLen,normalData,atpFrame,0x1);

  setAlignBlock(frame_ptr,(PHDWORD*) &alignBlk, 2);
  
  if (err != 0) {
    printf("PacketMaker: makeFrameHdr failed!\n");
    //		throw frameHdrFailed();
  }

  PHDWORD bufferSize = estimatedLength * sizeof(PHDWORD);
  PHDWORD* write_ptr;
  PHDWORD packetLength;

  for (PHDWORD i = 0;i<numberOfPackets;i++)
    {
      write_ptr = findFrameEnd(frame_ptr)+1;
      PHDWORD remainingSpace = bufferSize - (int) (write_ptr - (PHDWORD*) frame_ptr);

      makeUnstructPacket(
			 write_ptr,
			 remainingSpace,
			 packetMakerTableList->getPacketId(i), //packetId
			 4, //wordSize
			 packetMakerTableList->getPacketFormat(i) //hit format
			 );

      packetLength = storePacketHits(
				     write_ptr,
				     remainingSpace,
				     (UINT*) &addr_arr,
				     (BYTE*) packetMakerTableList->getPacketLocation(i),
				     packetMakerTableList->getPacketLength(i),
				     0
				     );
      printf("PacketMaker: packetLength =%d\n",packetLength);
      extendFrameDataNopad(frame_ptr,bufferSize,packetLength);
      totalFrameSize = getFrameLength(frame_ptr);

    }

  if (totalFrameSize > estimatedLength)
    {
      printf("PacketMaker: Frame too large. %d > %d\n",
	     totalFrameSize,estimatedLength);
      frameLength = 0;
    }
  frameLength = totalFrameSize;
}

PacketMaker::~PacketMaker()
{
  delete frameBuffer;
}

