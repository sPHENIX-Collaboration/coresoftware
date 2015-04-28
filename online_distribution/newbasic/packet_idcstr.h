#ifndef __PACKET_IDCSTR_H__
#define __PACKET_IDCSTR_H__


#include "packet_w124.h"

/**
   This is the packet which deals with data in IDCSTR format.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_idcstr : public Packet_w1 {
#else
class  Packet_idcstr : public Packet_w1 {
#endif

public:
  Packet_idcstr(PACKET_ptr);
  ~Packet_idcstr();
  int iValue(const int channel);
  void  dump ( OSTREAM& os=COUT) ;

protected:
  virtual int *decode (int *);
  unsigned char *sarray;
  int allocated_length;


};

#endif /* __PACKET_IDCSTR_H__ */
