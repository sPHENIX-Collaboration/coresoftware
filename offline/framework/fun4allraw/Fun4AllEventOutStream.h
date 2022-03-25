// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLEVENTOUTSTREAM_H
#define FUN4ALLRAW_FUN4ALLEVENTOUTSTREAM_H

// base class for output streams writing Events in
// one or the other form

#include <fun4all/Fun4AllBase.h>

#include <Event/phenixTypes.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/numeric/interval.hpp>
#pragma GCC diagnostic pop

#include <map>
#include <string>

class Event;
class Packet;
class Fun4AllEventOutputManager;

class Fun4AllEventOutStream : public Fun4AllBase
{
 public:
  virtual ~Fun4AllEventOutStream();
  virtual int StreamStatus() { return 0; }
  virtual int WriteEvent(Event *evt);
  virtual int WriteEventOut(Event * /*evt*/) { return 0; }
  virtual int CloseOutStream() { return 0; }

  int AddPacket(const int ipkt);
  int DropPacket(const int ipkt);
  int AddPacketRange(const int minpacket, const int maxpacket);
  int DropPacketRange(const int minpacket, const int maxpacket);
  void SetManager(Fun4AllEventOutputManager *myman) { m_MyManager = myman; }

 protected:
  Fun4AllEventOutStream(const std::string &name = "OUTSTREAM");
  int resize_evtbuf(const unsigned int newsize);
  Fun4AllEventOutputManager *MyManager() { return m_MyManager; }

 private:
  PHDWORD *evtbuf = nullptr;
  Fun4AllEventOutputManager *m_MyManager = nullptr;  // pointer to my master
  unsigned int evtbuf_size = 0;
  // flag to stear behavior, if 1 only add packets (drop all others), if 0 no filtering,
  // if -1 accept all, drop selected and afterwards add back selected ones
  int add_or_remove = 0;
  Packet **plist = nullptr;
  int max_npackets = 1000;
  int npackets = 0;
  int default_addall = 0;
  std::map<int, boost::numeric::interval<int> > addpktrange;
  std::map<int, boost::numeric::interval<int> > droppktrange;
};

#endif
