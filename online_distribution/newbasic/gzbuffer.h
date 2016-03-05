#ifndef __GZBUFFER_H
#define __GZBUFFER_H

#include "prdfBuffer.h"

#ifndef __CINT__
class WINDOWSEXPORT gzbuffer : public prdfBuffer{
#else
class  gzbuffer : public prdfBuffer{
#endif

public:

  //** Constructors

  gzbuffer( PHDWORD *array, const int length);
  ~gzbuffer();

  Event * getEvent();


protected:

  PHDWORD  *bufferarray;
  buffer *theBuffer;


};

#endif
