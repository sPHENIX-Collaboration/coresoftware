#ifndef __PACKET_ID4EVT_H__
#define __PACKET_ID4EVT_H__


#include "packet_w124.h"

/**
   This is the packet which deals with data in ID4EVT format.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet_id4evt : public Packet_w4 {
#else
class  Packet_id4evt : public Packet_w4 {
#endif

public:
  Packet_id4evt(PACKET_ptr);

protected:
  virtual int *decode (int *);
};

#endif /* __PACKET_ID4EVT_H__ */
