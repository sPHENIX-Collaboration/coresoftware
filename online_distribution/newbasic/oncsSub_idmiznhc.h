#ifndef __ONCSSUB_IDMIZNHC_H__
#define __ONCSSUB_IDMIZNHC_H__

#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idmiznhc : public oncsSubevent_w4{
#else
class  oncsSub_idmiznhc : public oncsSubevent_w4{
#endif

public:
  oncsSub_idmiznhc( subevtdata_ptr);

protected:
  int *decode (int *);
};


#endif /* __ONCSSUB_IDMIZNHC_H__ */
