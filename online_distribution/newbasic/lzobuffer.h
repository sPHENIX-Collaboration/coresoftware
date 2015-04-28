#ifndef __LZOBUFFER_H
#define __LZOBUFFER_H

#include "buffer.h"
#include <lzo/lzo1x.h>


#ifndef __CINT__
class WINDOWSEXPORT lzobuffer : public buffer{
#else
class  lzobuffer : public buffer{
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
