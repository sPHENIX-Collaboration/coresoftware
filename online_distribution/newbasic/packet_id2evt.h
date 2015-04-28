#ifndef __PACKET_ID2EVT_H__
#define __PACKET_ID2EVT_H__


#include "packet_w124.h"
/**
   This is the packet which deals with data in ID2EVT format.
   It inherits from Packet\_w2 because the data are 16bit entities.
*/

#ifndef __CINT__
class WINDOWSEXPORT Packet_id2evt : public Packet_w2 {
#else
class  Packet_id2evt : public Packet_w2 {
#endif

public:
  Packet_id2evt(PACKET_ptr);

protected:
  virtual int *decode (int *);
};

#endif /* __PACKET_ID2EVT_H__ */
