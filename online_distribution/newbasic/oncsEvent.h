#ifndef __ONCS_EVENT_H
#define __ONCS_EVENT_H

#include "Event.h"
#include "phenixTypes.h"
#include "oncsEvtConstants.h"
#include "oncsEvtStructures.h"
#include <map>

#ifndef __CINT__
class WINDOWSEXPORT oncsEvent : public Event {
#else
class  oncsEvent : public Event {
#endif

public:
  // constructors and destructors
  oncsEvent(int *);
  ~oncsEvent();

  virtual unsigned int getEvtLength();
  virtual int getEvtType();
  virtual int getEvtSequence();
  virtual int getRunNumber();
  //virtual PHTimeStamp * getTimeStamp() const;
  
  virtual void identify(std::ostream& os = std::cout) const;


  virtual Packet* getPacket(const int);
  virtual Packet* getPacket(const int, const int hitFormat);

  virtual int getPacketList(Packet*[], const int);

  virtual int existPacket (const int packetid);

  virtual int Copy ( int *destination, const unsigned int length, int *nw, const char *what="");

  virtual int is_pointer_type() const;
  virtual int convert ();

  virtual int getDate() { return 0;};
  virtual time_t getTime() const;
  virtual Packet * makePacket(PHDWORD *pp, const int hitFormat = 0);

protected:
  int is_data_type;  // 0 is pointer based --  1 is data based

  oncsevtdata_ptr EventData;

  int hasMap;
  int errorcode;
  virtual int createMap();
  std::map <int, PHDWORD *> pmap;
};

#endif
