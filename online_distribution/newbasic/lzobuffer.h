#ifndef __LZOBUFFER_H
#define __LZOBUFFER_H

#include "prdfBuffer.h"
#include <lzo/lzo1x.h>


#ifndef __CINT__
class WINDOWSEXPORT lzobuffer : public prdfBuffer{
#else
class  lzobuffer : public prdfBuffer{
#endif

public:

  //** Constructors

  lzobuffer( PHDWORD *array, const int length);
  ~lzobuffer();

  Event * getEvent();


protected:
  static int lzo_initialized;

  PHDWORD  *bufferarray;
  buffer *theBuffer;
  lzo_byte *wrkmem;

  int _broken;

};

#endif
