#ifndef __ONCSSUB_ID2EVT_H__
#define __ONCSSUB_ID2EVT_H__

#include "oncsSubevent.h"


#ifndef __CINT__
class WINDOWSEXPORT oncsSub_id2evt : public oncsSubevent_w2 {
#else
class  oncsSub_id2evt : public oncsSubevent_w2 {
#endif

public:
  oncsSub_id2evt( subevtdata_ptr);

protected:
  int *decode (int *);
};


#endif /* __ONCSSUB_ID2EVT_H__ */
