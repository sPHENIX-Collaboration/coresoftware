#ifndef __ONCSSUB_IDFIFOBOARD_H__
#define __ONCSSUB_IDFIFOBOARD_H__


#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idfifoboard : public  oncsSubevent_w4 {
#else
class  oncsSub_idfifoboard : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_idfifoboard( subevtdata_ptr);
  ~oncsSub_idfifoboard();

  //  int    iValue(const int);
  int    iValue(const int,const char *);
  double dValue(const int channel);
  void  dump ( OSTREAM& os = COUT) ;

protected:
  int *decode (int *);
  int samples;
  int evnr;
  double *d_time;
  int d_timelength;
};



#endif /* __ONCSSUB_IDFIFOBOARD_H__ */
