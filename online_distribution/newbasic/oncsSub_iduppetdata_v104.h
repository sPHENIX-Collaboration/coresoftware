#ifndef __ONCSSUB_IDUPPETDATA_V104_H__
#define __ONCSSUB_IDUPPETDATA_V104_H__


#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_iduppetdata_v104 : public  oncsSubevent_w4 {
#else
class  oncsSub_iduppetdata_v104 : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_iduppetdata_v104( subevtdata_ptr);
  ~oncsSub_iduppetdata_v104();

  //  int    iValue(const int);
  int    iValue(const int,const char *);
  double dValue(const int channel);
  long long  lValue(const int channel);
  void  dump ( OSTREAM& os = COUT) ;

protected:
  int *decode (int *);
  int samples;
  int evnr;
  double *d_time;
  long long *lval;
  long long *tval;
  unsigned int serialnumber;
  unsigned int udpheader1;
  unsigned int udpheader2;
  unsigned int sysword;
  unsigned int eventrate;
  unsigned int ratcap_upenn;
  unsigned int clock_sel;
  
  int d_timelength;
};



#endif /* __ONCSSUB_IDUPPETDATA_V104_H__ */
