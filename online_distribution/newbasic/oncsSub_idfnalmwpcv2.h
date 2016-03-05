#ifndef __ONCSSUB_IDFNALMWPCV2_H__
#define __ONCSSUB_IDFNALMWPCV2_H__

#include <vector>

#include "oncsSub_idfnalmwpc.h"

#ifndef __CINT__
class WINDOWSEXPORT oncsSub_idfnalmwpcv2 :  public oncsSub_idfnalmwpc {
#else
class oncsSub_idfnalmwpcv2 :  public oncsSub_idfnalmwpc {
#endif

public:
  oncsSub_idfnalmwpcv2( subevtdata_ptr);
  ~oncsSub_idfnalmwpcv2();


protected:
  int *decode (int *);

};

#endif 
