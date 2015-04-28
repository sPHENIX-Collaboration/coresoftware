#ifndef __GZBUFFER_H
#define __GZBUFFER_H

#include "buffer.h"

#ifndef __CINT__
class WINDOWSEXPORT gzbuffer : public buffer{
#else
class  gzbuffer : public buffer{
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
