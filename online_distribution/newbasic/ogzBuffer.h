#ifndef __OGZBUFFER_H
#define __OGZBUFFER_H

#include <zlib.h>
#include "oBuffer.h"



#ifndef __CINT__
class WINDOWSEXPORT ogzBuffer : public oBuffer{
#else
class  ogzBuffer : public oBuffer{
#endif

public:

  //** Constructors
#ifndef WIN32
  ogzBuffer (int fd , PHDWORD * where, 
	     const int length,
	     const int level =3,
	     const int irun=1, 
	     const int iseq=0 );
#else
  ogzBuffer (const char *fpp, PHDWORD * where, 
	     const int length,
             int &status,
	     const int level =3,
	     const int irun=1, 
	     const int iseq=0 );
#endif
  virtual  ~ogzBuffer();


  virtual int writeout ();


protected:

  PHDWORD  *outputarray;
  uLongf  outputarraylength;
  int compressionlevel;

};

#endif

