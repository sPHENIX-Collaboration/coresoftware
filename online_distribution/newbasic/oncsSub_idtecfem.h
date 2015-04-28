#ifndef __ONCSSUB_IDTECFEM_H__
#define __ONCSSUB_IDTECFEM_H__

#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idtecfem : public oncsSubevent_w4 {
#else
class  oncsSub_idtecfem : public oncsSubevent_w4 {
#endif

public:
  oncsSub_idtecfem( subevtdata_ptr);
  virtual int    iValue(const int,const int);
  virtual float  rValue(const int,const int);

protected:
  int *decode (int *);
};


#endif /* __ONCSSUB_IDTECFEM_H__ */
