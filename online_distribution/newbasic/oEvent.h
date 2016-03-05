#ifndef __oEVENT_BASE_H
#define __oEVENT_BASE_H

#include "EvtConstants.h"
#include "phenixTypes.h"
#include "EvtStructures.h"
#include "packet.h"

#define WINDOWSEXPORT

// virtual base class for an "event"

#ifndef __CINT__
class WINDOWSEXPORT oEvent {
#else
class  oEvent {
#endif

public:

  //** Constructors
  oEvent(PHDWORD *, const int maxlength
	 , const int irun, const int itype, const int eseq);

  virtual ~oEvent() {}
  virtual int prepare_next();
  virtual int prepare_next( const int, const int);
  virtual void set_event_type(const int);

  // frame and packet adding
  virtual int addFrame(PHDWORD *);

  virtual int addPacket ( const Packet * p);

  virtual int addUnstructPacketData(PHDWORD * data, 
		    const int length,
		    const int id,
		    const int wordsize,
		    const int hitformat);
protected:
  
  PHDWORD *event_base;
  int current;
  int in_frame;
  PHDWORD * currentFramePtr;
  evtdata_ptr evthdr;
  int max_length;
  int left;
  int maxSizeOfThisFrame;
};

#endif
