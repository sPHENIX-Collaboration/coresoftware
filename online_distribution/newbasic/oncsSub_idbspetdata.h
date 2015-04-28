#ifndef __ONCSSUB_IDBSPETDATA_H__
#define __ONCSSUB_IDBSPETDATA_H__


#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idbspetdata : public  oncsSubevent_w4 {
#else
class  oncsSub_idbspetdata : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_idbspetdata( subevtdata_ptr);
  ~oncsSub_idbspetdata();

  //  int    iValue(const int);
  int    iValue(const int,const char *);
  double dValue(const int channel);
  long long  lValue(const int channel);
  void  dump ( OSTREAM& os = COUT) ;

protected:
  int *decode (int *);
  int samples;
  int evnr;
  double *d_time;
  long long *lval;
  int d_timelength;
};



#endif /* __ONCSSUB_IDBSPETDATA_H__ */
