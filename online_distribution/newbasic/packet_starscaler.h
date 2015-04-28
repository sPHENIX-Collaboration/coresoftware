#ifndef __PACKET_STARSCALER_H__
#define __PACKET_STARSCALER_H__


#include "packet_w124.h"
#if !defined(SunOS) && !defined(OSF1)
#include <map>
#endif

/**
   This is the packet which deals with data in STARSCALER format.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_starscaler : public Packet_w4 {
#else
class  Packet_starscaler : public Packet_w4 {
#endif


public:
  Packet_starscaler(PACKET_ptr);
  ~Packet_starscaler();
  long long  lValue(const int channel);
  int    iValue(const int channel);
  int    iValue(const int channel,const char *what);

  void  dump ( OSTREAM& ) ;

protected:

  virtual int  *decode (int *);
  long long *s_vector;
  unsigned int s_vectorlength;

#if !defined(SunOS) && !defined(OSF1)
  std::map <int, int> smap;
#endif


};

#endif /* __PACKET_STARSCALER_H__ */
