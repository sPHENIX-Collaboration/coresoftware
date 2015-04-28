#ifndef __ONCSSUB_IDSAM_H__
#define __ONCSSUB_IDSAM_H__

#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idsam : public oncsSubevent_w4 {
#else
class  oncsSub_idsam : public oncsSubevent_w4 {
#endif

public:
  oncsSub_idsam( subevtdata_ptr);

protected:
  int *decode (int *);
};


#endif /* __ONCSSUB_IDSAM_H__ */
