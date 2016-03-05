//
// testEventIterator   mlp 9/11/96
//
// this iterator delivers test-type Event objects which can
// be used to debug the analysis software which uses the events
// delivered. The events contain 1 frame and one packet


#include <time.h>

#include "testEventiterator.h"
#include "oBuffer.h"
#include <stdlib.h>

#include "A_Event.h"
#include "Cframe.h"
#include "frameRoutines.h"
#include "packetRoutines.h"
#include "packetConstants.h"
#define EVTLENGTH 500

class oBuffer;

testEventiterator::~testEventiterator()
{
  delete R;
}

testEventiterator::testEventiterator()
{
  R = new simpleRandom(876565);
  current_event = 0;
}


void  
testEventiterator::identify (OSTREAM &os) const
{ 
  os << getIdTag()  << std::endl;

};

const char * 
testEventiterator::getIdTag () const
{ 
  return " -- testEventiterator (standard) ";
};



Event *
testEventiterator::getNextEvent()
{

  int i;
  int packetlength;

  evtdata_ptr Eptr;
  PHDWORD tdata[EVTLENGTH];
  ALIGNBLK alignBlk;

  for (i = 0; i<EVTLENGTH; ) tdata[i++]=0;
  Eptr = ( evtdata_ptr ) tdata;

  // construct the Event header
  Eptr->evt_length = EVTHEADERLENGTH;
  Eptr->evt_type=1;
  Eptr->evt_sequence=++current_event;
  Eptr->run_number=1331;
  Eptr->date =time(0);
  Eptr->time =-1;

  alignBlk.dcb.timeStamp     =   0;    
  alignBlk.dcb.granuleEvtcnt =   0;    
  alignBlk.dcb.partitionVec  =   0;

  int idata[20];
  short sdata[20];
  int rdata[4];

  for (i = 0; i<20; i++) idata[i] = i;
  for (i = 0; i<20; i++) sdata[i] = i;

  //preset the empty event header.
  tdata[0] = EVTHEADERLENGTH;

  PHDWORD *frame = &tdata[ tdata[0] ];
  makeFrameHdr(frame,EVTLENGTH-EVTHEADERLENGTH
		     ,rawData,oncsFrame,101);

  PHDWORD *packetstart;

  // --------- packet 1001 -------------------

  packetstart = findFrameDataEnd (frame) +1;

  makeUnstructPacket (packetstart, EVTLENGTH, 1001, 4, ID4EVT);

  packetlength = storePacketHits (packetstart, EVTLENGTH, 
                   0, (BYTE*) idata, 20, 0);

  adjustFrameLength (frame, EVTLENGTH , packetlength, 1);

  // --------- packet 1002 -------------------
  packetstart = findFrameDataEnd (frame) + 1;

  makeUnstructPacket (packetstart, EVTLENGTH, 1002, 2, ID2EVT);

  packetlength = storePacketHits (packetstart, EVTLENGTH, 
                   0, (BYTE*) sdata, 20, 0);

  adjustFrameLength (frame, EVTLENGTH , packetlength, 1);

  // --------- packet 1003 -------------------

  packetstart = findFrameDataEnd (frame) +1;

  makeUnstructPacket (packetstart, EVTLENGTH, 1003, 4, ID4EVT);

#ifndef WIN32
  rdata[0] = int(R->gauss(0.,10.));
  rdata[1] = int(R->gauss(0.,100.));
  rdata[2] = int(R->gauss(0.,1000.));
  rdata[3] = int(R->gauss(0.,10000.));
#else
  rdata[0] = 10 * rand();
  rdata[1] = 100 * rand();
  rdata[2] = 1000 * rand();
  rdata[3] = 10000 * rand();
#endif

  packetlength = storePacketHits (packetstart, EVTLENGTH, 
                   0, (BYTE*) rdata, 4, 0);

  adjustFrameLength (frame, EVTLENGTH , packetlength, 1);

  tdata[0] += getFrameLength (frame);

  Event *e = new  A_Event(tdata);
  e->convert();

  return e;
}
