#ifndef __PACKET_hbd_fpgashort_H__
#define __PACKET_hbd_fpgashort_H__


#include "packet_w124.h"

/**
   This is the packet which deals with data in hbd_fpgashort format.
   It inherits from Packet\_w4 because the data are 32bit entities.
*/

#ifndef __CINT__
class WINDOWSEXPORT Packet_hbd_fpgashort : public Packet_w4{ 
#else
class  Packet_hbd_fpgashort : public Packet_w4 {
#endif

public:
  Packet_hbd_fpgashort(PACKET_ptr);

  /** with the "what" parameter you can decide which aspect of
 the data is made available. This class is one of those which have
 several different "kinds" of data; we use this to bring up the AMU
 cell information and all the misc. items in the FEM headers and
 trailers.


  */



  virtual int    iValue(const int channel,const char *what);
  virtual int    iValue(const int channel,const int y);

  void setNumSamples(const int ns) { HBD_NSAMPLES = ns; }
  
  virtual void   dump ( OSTREAM& );

 protected:
  virtual int *decode (int *);
  int nr_modules;
  int HBD_NSAMPLES;
  int hbd_parity;
};

#endif /* __PACKET_hbd_fpgashort_H__ */


