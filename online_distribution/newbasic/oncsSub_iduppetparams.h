#ifndef __ONCSSUB_IDUPPETPARAMS_H__
#define __ONCSSUB_IDUPPETPARAMS_H__


#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_iduppetparams : public  oncsSubevent_w4 {
#else
class  oncsSub_iduppetparams : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_iduppetparams( subevtdata_ptr);
  ~oncsSub_iduppetparams();

  //  int    iValue(const int);
  int    iValue(const int,const char *);
  int    iValue(const int);
  void  dump ( OSTREAM& os = COUT) ;

protected:
  int *decode (int *);

};



#endif /* __ONCSSUB_IDUPPETPARAMS_H__ */
