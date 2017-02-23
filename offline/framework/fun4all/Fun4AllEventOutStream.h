#ifndef __FUN4ALLEVENTOUTSTREAM_H__
#define __FUN4ALLEVENTOUTSTREAM_H__

// base class for output streams writing Events in
// one or the other form

#include "Fun4AllBase.h"
#include <Event/phenixTypes.h>

#ifndef __CINT__
#include <boost/version.hpp> // to get BOOST_VERSION
#if (__GNUC__ == 4 && __GNUC_MINOR__ == 8 && (BOOST_VERSION == 105700 || BOOST_VERSION == 106000 || BOOST_VERSION == 106300) )
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <boost/numeric/interval.hpp>
#pragma GCC diagnostic warning "-Wunused-local-typedefs"
#else
#include <boost/numeric/interval.hpp>
#endif
#endif

#include <map>
#include <set>
#include <string>


class Event;
class Packet;
class Fun4AllEventOutputManager;

class Fun4AllEventOutStream: public Fun4AllBase
{
 public:
  virtual ~Fun4AllEventOutStream();
  virtual int StreamStatus() {return 0;}
  virtual int WriteEvent(Event *evt);
  virtual int WriteEventOut(Event* /*evt*/) {return 0;}
  virtual int CloseOutStream() {return 0;}

  int  AddPacket(const int ipkt); 
  int DropPacket(const int ipkt);
  int AddPacketRange(const int minpacket, const int maxpacket);
  int DropPacketRange(const int minpacket, const int maxpacket);
  void SetManager(Fun4AllEventOutputManager *myman) {mymanager = myman;}
 
 protected:

  Fun4AllEventOutStream(const std::string &name= "OUTSTREAM");
  int resize_evtbuf(const unsigned int newsize);
  PHDWORD *evtbuf;
  unsigned int evtbuf_size;  
// flag to stear behavior, if 1 only add packets (drop all others), if 0 no filtering, 
// if -1 accept all, drop selected and afterwards add back selected ones
  int add_or_remove; 
  Packet **plist;
  int max_npackets;
  int npackets;
  int default_addall;
#ifndef __CINT__
  std::map<int,boost::numeric::interval<int> > addpktrange; 
  std::map<int,boost::numeric::interval<int> > droppktrange; 
#endif
   Fun4AllEventOutputManager *mymanager; // pointer to my master
};

#endif /* __FUN4ALLEVENTOUTSTREAM_H__ */
