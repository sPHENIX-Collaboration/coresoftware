#ifndef __PACKET_GL1_EVCLOCKS__
#define __PACKET_GL1_EVCLOCKS__

#include <packet_w124.h>

/**
   This is the packet which deals with data in GL1 format.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_gl1_evclocks : public Packet_w4 {
#else
class  Packet_gl1_evclocks : public Packet_w4 {
#endif
public:
  Packet_gl1_evclocks(PACKET_ptr);
  ~Packet_gl1_evclocks();
 virtual int  iValue(const int channel, const char *what);
 virtual int  iValue(const int channel, const int what);
 virtual void dump ( OSTREAM& );

protected:
 enum GL1_Evclock_Types {EVCLOCK,PARVECT};
 virtual int *decode (int *);
};

#endif /* __PACKET_GL1_EVCLOCKS__ */

















