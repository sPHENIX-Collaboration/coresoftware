#ifndef __ONCSSUB_ID4EVT_H__
#define __ONCSSUB_ID4EVT_H__

#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_id4evt : public  oncsSubevent_w4 {
#else
class  oncsSub_id4evt : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_id4evt( subevtdata_ptr);

protected:
  int *decode (int *);
};


#endif /* __ONCSSUB_ID4EVT_H__ */
