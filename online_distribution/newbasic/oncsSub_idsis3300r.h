#ifndef __ONCSSUB_IDSIS3300R_H__
#define __ONCSSUB_IDSIS3300R_H__


#include "oncsSub_idsis3300.h"
#include <vector>

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idsis3300r : public oncsSub_idsis3300 {
#else
class  oncsSub_idsis3300r : public oncsSub_idsis3300 {
#endif

public:
  oncsSub_idsis3300r( subevtdata_ptr);

  int    iValue(const int,const int);
  //int    iValue(const int,const char *);
  //void  dump ( OSTREAM& os = COUT) ;

protected:
  int *decode (int *);
  //  int samples;
  int decoded;
  std::vector<int> v[8];

};


#endif /* __ONCSSUB_IDSIS3300R_H__ */
