#ifndef __ONCSSUB_IDSIS3300_H__
#define __ONCSSUB_IDSIS3300_H__


#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idsis3300 : public  oncsSubevent_w4 {
#else
class  oncsSub_idsis3300 : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_idsis3300( subevtdata_ptr);

  int    iValue(const int,const int);
  int    iValue(const int,const char *);
  void  dump ( OSTREAM& os = COUT) ;

protected:
  int *decode (int *);
  int samples;
  int wraparound;
};



#endif /* __ONCSSUB_IDSIS3300_H__ */
