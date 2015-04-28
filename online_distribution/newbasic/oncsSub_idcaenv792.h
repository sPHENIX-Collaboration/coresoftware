#ifndef __ONCSSUB_IDCAENV792_H__
#define __ONCSSUB_IDCAENV792_H__


#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idcaenv792 : public  oncsSubevent_w4 {
#else
class  oncsSub_idcaenv792 : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_idcaenv792( subevtdata_ptr);

  int    iValue(const int);
  int    iValue(const int,const char *);
  void  dump ( OSTREAM& os = COUT) ;

protected:
  int *decode (int *);
  int samples;
  int evnr;
};



#endif /* __ONCSSUB_IDCAENV792_H__ */
