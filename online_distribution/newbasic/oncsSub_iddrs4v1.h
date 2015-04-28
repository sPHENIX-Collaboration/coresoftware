#ifndef __ONCSSUB_IDDRS4V1_H__
#define __ONCSSUB_IDDRS4V1_H__

#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_iddrs4v1 : public  oncsSubevent_w4 {
#else
class  oncsSub_iddrs4v1 : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_iddrs4v1( subevtdata_ptr);
  ~oncsSub_iddrs4v1();

  int      iValue(const int,const char *);
  float    rValue(const int,const char *);
  float    rValue(const int,const int );
  void  dump ( OSTREAM& os = COUT) ;

protected:
  int *decode (int *);
  float *wave;
  int samples;
  int enabled_channelmask;
  int channel_offset[4];

};


#endif /* __ONCSSUB_IDDRS4V1_H__ */
