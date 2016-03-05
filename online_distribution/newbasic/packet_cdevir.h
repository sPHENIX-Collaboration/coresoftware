#ifndef __PACKET_CDEVIR_H__
#define __PACKET_CDEVIR_H__


#include <packet_w124.h>

/**
   This is the packet which deals with data in CDEVIR format.
   It inherits from Packet\_w4 because the data are 32bit entities.


*/



#ifndef __CINT__
class WINDOWSEXPORT Packet_cdevir : public Packet_w4 {
#else
class  Packet_cdevir : public Packet_w4 {
#endif

public:
  Packet_cdevir(PACKET_ptr);
  virtual void  dump ( OSTREAM& ) ;
  virtual double  dValue(const int channel,const char *what);

protected:
  virtual int *decode (int *);
  struct cdevIrData *ps;
};

#endif /* __PACKET_CDEVIR_H__ */
