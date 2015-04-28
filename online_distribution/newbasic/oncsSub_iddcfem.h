#ifndef __ONCSSUB_IDDCFEM_H__
#define __ONCSSUB_IDDCFEM_H__

#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_iddcfem : public oncsSubevent_w4 {
#else
class  oncsSub_iddcfem : public oncsSubevent_w4 {
#endif

public:
  oncsSub_iddcfem( subevtdata_ptr);
  virtual int    iValue(const int,const int);
  virtual float  rValue(const int,const int);

protected:
  int *decode (int *);
};


#endif /* __ONCSSUB_IDDCFEM_H__ */
