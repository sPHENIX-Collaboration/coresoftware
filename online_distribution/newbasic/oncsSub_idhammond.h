#ifndef __ONCSSUB_IDHAMMOND_H__
#define __ONCSSUB_IDHAMMOND_H__

#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idhammond : public oncsSubevent_w4 {
#else
class  oncsSub_idhammond : public oncsSubevent_w4 {
#endif

public:
  oncsSub_idhammond( subevtdata_ptr);
  int    iValue(const int,const int);
  float  rValue(const int,const int);

protected:
  int *decode (int *);
};



#endif /* __ONCSSUB_IDHAMMOND_H__ */
