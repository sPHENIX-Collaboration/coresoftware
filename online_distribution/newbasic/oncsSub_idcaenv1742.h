#ifndef __ONCSSUB_IDCAENV1742_H__
#define __ONCSSUB_IDCAENV1742_H__


#include "oncsSubevent.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idcaenv1742 : public  oncsSubevent_w4 {
#else
class  oncsSub_idcaenv1742 : public  oncsSubevent_w4 {
#endif

public:
  oncsSub_idcaenv1742( subevtdata_ptr);

  int    iValue(const int ch);
  int    iValue(const int sample, const int ch);
  int    iValue(const int,const char *);
  void  dump ( OSTREAM& os = COUT) ;

protected:
  int *decode (int *);
  int samples;
  int dlength;
  int evnr;
  int freq;
  int group_mask;
  int index_cell[4];
  int tr_present[4];

};



#endif /* __ONCSSUB_IDCAENV1742_H__ */
