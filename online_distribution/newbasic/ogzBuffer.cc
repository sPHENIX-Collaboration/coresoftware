
#include "ogzBuffer.h"
#include "BufferConstants.h"


// the constructor first ----------------
#ifndef WIN32
ogzBuffer::ogzBuffer (int fdin, PHDWORD * where, 
		      const int length, 
		      const int level,
		      const int irun, 
		      const int iseq): 
  oBuffer(fdin,where,length,irun,iseq)
#else
ogzBuffer::ogzBuffer (const char *fpp, PHDWORD * where, 
		      const int length, 
                      int & status,
		      const int level,
		      const int irun, 
		      const int iseq): 
  oBuffer(fpp,where,length,status,irun,iseq)
#endif
{
  // get a buffer for zlib

  compressionlevel = level;
  outputarraylength = (int)(length *1.1) + 2048;
  outputarray = new PHDWORD[outputarraylength];

}

// ----------------------------------------------------------
// returns the number of bytes written, including record wasted space.
//
int ogzBuffer::writeout()
{
  if (! dirty) return 0;

  if (! has_end) addEoB();

  uLongf outputlength_in_bytes = outputarraylength*4;
  uLong bytes_to_be_written = bptr->Length; 

  compress2 ( (Bytef*) &outputarray[4], &outputlength_in_bytes, (Bytef*) bptr, 
		bytes_to_be_written, compressionlevel);

  outputarray[0] = outputlength_in_bytes +4*BUFFERHEADERLENGTH;
  outputarray[1] = GZBUFFERMARKER; // -518;
  outputarray[2] = bptr->Bufseq;
  outputarray[3] = bptr->Length;

  unsigned int ip =0;
  char *cp = (char *) outputarray;

  while (ip<outputarray[0])
    {
      write ( fd, cp, BUFFERBLOCKSIZE);
      cp += BUFFERBLOCKSIZE;
      ip += BUFFERBLOCKSIZE;
    }
  dirty = 0;
  byteswritten += ip;
  return 0;
}


// ----------------------------------------------------------
ogzBuffer::~ogzBuffer()
{
  writeout();
  delete [] outputarray;

}

