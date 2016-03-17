#include "oEvent.h"
#include "EventTypes.h"

#include "A_Event.h"
#include "Cframe.h"
#include "frameRoutines.h"
#include "packetRoutines.h"
#include "packetConstants.h"

#include <time.h>

// the constructor first ----------------
oEvent::oEvent (PHDWORD * where, const int length
		, const int irun, const int etype, const int evtseq)
{
  event_base = where;
  evthdr = ( evtdata_ptr ) where;
  evthdr->evt_type = etype;
  max_length = length;
  maxSizeOfThisFrame = 0;
  prepare_next (evtseq, irun);
}

void oEvent::set_event_type(const int etype)
{
   evthdr->evt_type = etype;
}

int oEvent::prepare_next()
{
  // re-initialize the event header length
  evthdr->evt_length =  EVTHEADERLENGTH;

  // if < 0, just increment the current seq. number
  evthdr->evt_sequence++;

  // reset the current data index, and the leftover counter
  current = 0;
  in_frame = 0; // indicate that we are just now assembling a frame
  left=max_length - EVTHEADERLENGTH ;
  evthdr->date = time(0);
  evthdr->time = -1;

  return 0;

}

int oEvent::prepare_next( const int evtseq 
			  , const int irun )
{
  // re-initialize the event header length
  evthdr->evt_length =  EVTHEADERLENGTH;

  // if < 0, just increment the current seq. number

  evthdr->evt_sequence = evtseq;

  // if > 0, adjust the run number, else just keep it.
  evthdr->run_number=irun;
  
  // reset the current data index, and the leftover counter
  current = 0;
  in_frame = 0; // indicate that we are just now assembling a frame
  left=max_length - EVTHEADERLENGTH ;
  evthdr->date = time(0);
  evthdr->time = -1;

  return 0;
}
  

int oEvent::addFrame(PHDWORD *frame)
{
  int len,i;
  PHDWORD *to = &(evthdr->data[current]);
  PHDWORD *from = frame;

  len = getFrameLength(frame);
  for (i=0; i< len; i++) *to++ = *from++;

  evthdr->evt_length += len;
  current += len;
  left -= len;
  in_frame = 0;
  return len;
}



int  oEvent::addPacket( const Packet *p)
{

  int additionalFrameLength = 0;
  if (!in_frame) 
    {		
      currentFramePtr = &evthdr->data[current];
      makeFrameHdr(currentFramePtr  , left,
                     rawData,oncsFrame,101);
      in_frame = 1;
      left -= currentFrameHdrLength;
      maxSizeOfThisFrame = left;
      additionalFrameLength = currentFrameHdrLength;
      current += currentFrameHdrLength;
      evthdr->evt_length += currentFrameHdrLength;
    }

  PHDWORD* packetstart;
  int packetlength  = p->getLength(); 

  
  packetstart = findFrameDataEnd(currentFramePtr ) +1;

  p->copyMe( (int *)packetstart, packetlength );

  if (packetlength >0)  
    { 
      evthdr->evt_length += packetlength;
      current +=  packetlength;
      adjustFrameLength ( currentFramePtr, maxSizeOfThisFrame , packetlength, 1);
      left -= packetlength;
      return packetlength + additionalFrameLength;
    }
  else return  -1;
}



int  oEvent::addUnstructPacketData(PHDWORD * data, 
		    const int length,
		    const int id,
		    const int wordsize,
		    const int hitformat)
{

  int additionalFrameLength = 0;
  if (!in_frame) 
    {		
      currentFramePtr = &evthdr->data[current];
      makeFrameHdr(currentFramePtr  , left,
                     rawData,oncsFrame,101);
      in_frame = 1;
      left -= currentFrameHdrLength;
      maxSizeOfThisFrame = left;
      additionalFrameLength = currentFrameHdrLength;
      current += currentFrameHdrLength;
      evthdr->evt_length += currentFrameHdrLength;
    }

  PHDWORD* packetstart;
  int packetlength;
  
  packetstart = findFrameDataEnd(currentFramePtr ) +1;

  makeUnstructPacket (packetstart, left, id, wordsize, hitformat);

  packetlength = storePacketHits (packetstart, left, 
                   0, (BYTE*) data, length, 0);

  if (packetlength >0)  
    { 
      evthdr->evt_length += packetlength;
      current +=  packetlength;
      adjustFrameLength ( currentFramePtr, maxSizeOfThisFrame , packetlength, 1);
      left -= packetlength;
      return packetlength + additionalFrameLength;
    }
  else return  -1;
}

