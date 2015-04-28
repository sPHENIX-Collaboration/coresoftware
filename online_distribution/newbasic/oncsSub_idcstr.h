#ifndef __ONCSSUB_IDCSTR_H__
#define __ONCSSUB_IDCSTR_H__

#include "oncsSubevent.h"


#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idcstr : public oncsSubevent_w1 {
#else
class  oncsSub_idcstr : public oncsSubevent_w1 {
#endif

public:
  oncsSub_idcstr( subevtdata_ptr);
  ~oncsSub_idcstr();
  int iValue (const int i);
  void  dump ( OSTREAM& os = COUT) ;


protected:
  int *decode (int *);
  unsigned char *sarray;
  int allocated_length;
};


#endif /* __ONCSSUB_IDCSTR_H__ */
