#ifndef __ONCSSUB_IDDIGITIZERV1_H__
#define __ONCSSUB_IDDIGITIZERV1_H__

#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_iddigitizerv1 : public  oncsSubevent_w4 {
#else
class  oncsSub_iddigitizerv1 : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_iddigitizerv1( subevtdata_ptr);
  ~oncsSub_iddigitizerv1();

  int    iValue(const int sample, const int ch);
  int    iValue(const int ,const char * what);
  void  dump ( OSTREAM& os = COUT) ;


protected:
  int decode ();


  int _nsamples;
  int _l1_delay;
  int _slot_nr;
  int _nr_modules;
  int _clock;
  int _evtnr;
  int _nchannels;
  int _slot_nr_from_data;

  int _is_decoded;

  int array[32][128];

};


#endif /* __ONCSSUB_IDDIGITIZERV1_H__ */
