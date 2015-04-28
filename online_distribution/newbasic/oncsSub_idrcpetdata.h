#ifndef __ONCSSUB_IDRCPETDATA_H__
#define __ONCSSUB_IDRCPETDATA_H__


#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idrcpetdata : public  oncsSubevent_w4 {
#else
class  oncsSub_idrcpetdata : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_idrcpetdata( subevtdata_ptr);
  ~oncsSub_idrcpetdata();

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



#endif /* __ONCSSUB_IDRCPETDATA_H__ */
