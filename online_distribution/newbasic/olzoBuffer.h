#ifndef __OLZOBUFFER_H__
#define __OLZOBUFFER_H__

#include <lzo/lzo1x.h>

#include "oBuffer.h"



#ifndef __CINT__
class WINDOWSEXPORT olzoBuffer : public oBuffer{
#else
class  olzoBuffer : public oBuffer{
#endif

public:

  //** Constructors

#ifndef WIN32
  olzoBuffer (int fdin, PHDWORD * where, 
	     const int length,
	     const int irun=1, 
	     const int iseq=0 );
#else
  olzoBuffer (const char *fpp, PHDWORD * where, 
	     const int length,
             int &status,
	     const int irun=1, 
	     const int iseq=0 );
#endif
  virtual  ~olzoBuffer();


  virtual int writeout ();


protected:

  static int lzo_initialized;

  int _broken;

  lzo_byte *wrkmem;

  PHDWORD  *outputarray;
  lzo_uint  outputarraylength;
 

};



#endif /* __OLZOBUFFER_H__ */

