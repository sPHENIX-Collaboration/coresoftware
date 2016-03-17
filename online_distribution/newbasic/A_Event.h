#ifndef __A_EVENT_H
#define __A_EVENT_H

#include "Event.h"
#include "EvtConstants.h"
#include "phenixOnline.h"
#include "EvtStructures.h"

#if !defined(SunOS) && !defined(OSF1)
#include <map>
#endif

#ifndef __CINT__
class WINDOWSEXPORT A_Event : public Event {
#else
class  A_Event : public Event {
#endif

public:
  // constructors and destructors
  A_Event( PHDWORD *);
  A_Event( int *);
  virtual ~A_Event();

  // info & debug utils

  virtual unsigned int getEvtLength();
  virtual int getEvtType();
  virtual int getEvtSequence();
  virtual int getRunNumber();
  // virtual PHTimeStamp * getTimeStamp() const;

  virtual void identify(std::ostream& os = std::cout) const;

  virtual void listFrame( const int id=0, OSTREAM& os=COUT ) const;

  virtual void listHistory( const int id=0, OSTREAM& os=COUT ) const;

  virtual void listError( const int id=0, OSTREAM& os=COUT ) const;

  unsigned int getFrameEntry(  const char *what, const int id=0, const int index =0) const;

  // packet handling
  virtual Packet* getPacket(const int);

  virtual Packet* getPacket(const int, const int hitFormat);

  virtual int getPacketList(Packet*[], const int);

  virtual int existPacket (const int);

  //event copying
  virtual int Copy ( int *, const unsigned int, int *, const char *what ="");



  int getErrorCode();
  unsigned int getTagWord( const int i=0) const 
    { 
      if ( i) return EventData->reserved[1]; 
      return EventData->reserved[0]; 
    };

  virtual int is_pointer_type() const;
  virtual int convert ();

  static void dumpFrame(PHDWORD *fp, OSTREAM &os = COUT);
  static void dumpErrorBlock(PHDWORD *fp, OSTREAM &os = COUT);

  static void dumpBlock(PHDWORD *p, UINT len, OSTREAM &os = COUT, const int how=EVT_HEXADECIMAL);

  virtual int getDate() { return 0;};
  virtual time_t getTime() const;

  static Packet *makePacket(PHDWORD *pp, const int hitformat=0);

protected:
  virtual int updateFramelist();

  virtual unsigned int getFrameValue(const char *what, PHDWORD *fp, const unsigned int index =0) const;


#if !defined(SunOS) && !defined(OSF1)
  virtual int createMap();
#endif
  int is_data_type;  // 0 is pointer based --  1 is data based
  evtdata_ptr EventData;
  PHDWORD **framelist;
  int NumberFrames;
  int hasMap;
  int errorcode;

#if !defined(SunOS) && !defined(OSF1)
  std::map <int, PHDWORD *> pmap;

#endif

};

#endif



