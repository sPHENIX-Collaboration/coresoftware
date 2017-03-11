#ifndef __PACKET_IDDIGITIZERV2_H__
#define __PACKET_IDDIGITIZERV2_H__


#include "packet_w124.h"

#ifndef __CINT__
class WINDOWSEXPORT Packet_iddigitizerv2 : public  Packet_w4 {
#else
class  Packet_iddigitizerv2 : public  Packet_w4 {
#endif

public:
  Packet_iddigitizerv2( PACKET_ptr);
  ~Packet_iddigitizerv2();

  int    iValue(const int sample, const int ch);
  int    iValue(const int ,const char * what);
  void  dump ( OSTREAM& os = COUT) ;


protected:
  int decode ();


  int _evtnr;
  int _detid;
  int _module_address;
  int _clock;
  int _fem_slot;
  int _fem_evtnr;
  int _fem_clock;

  int _nsamples;
  int _nr_modules;


  int _nchannels;
  int _is_decoded;

  int array[32][128];

};


#endif /* __PACKET_IDDIGITIZERV2_H__ */
