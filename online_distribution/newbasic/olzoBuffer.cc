
#include "olzoBuffer.h"
#include "BufferConstants.h"

#include <lzo/lzoutil.h>
#include <cstring>
#include <stdlib.h>
#include <unistd.h>

int olzoBuffer::lzo_initialized = 0;


// the constructor first ----------------
#ifndef WIN32
olzoBuffer::olzoBuffer (int fdin, PHDWORD * where, 
		      const int length, 
		      const int irun, 
		      const int iseq): 
  oBuffer(fdin,where,length,irun,iseq)
#else
olzoBuffer::olzoBuffer (const char *fpp, PHDWORD * where, 
		      const int length, 
                      int &status,
		      const int irun, 
		      const int iseq): 
  oBuffer(fpp,where,length,status,irun,iseq)
#endif
{
  // get a buffer for zlib

  _broken = 0;

  if ( !  lzo_initialized )
    {
      if (lzo_init() != LZO_E_OK)
	{
	  COUT << "Could not initialize LZO" << std::endl;
	  _broken = 1;
	}
      
      lzo_initialized = 1;
    }


  wrkmem = (lzo_bytep) lzo_malloc(LZO1X_1_12_MEM_COMPRESS);
  if (wrkmem)
    {
      memset(wrkmem, 0, LZO1X_1_12_MEM_COMPRESS);
    }
  outputarraylength = (int)(length *1.1) + 2048;
  outputarray = new PHDWORD[outputarraylength];

}

// ----------------------------------------------------------
// returns the number of bytes written, including record wasted space.
//
int olzoBuffer::writeout()
{


  if (! dirty) return 0;

  if (! has_end) addEoB();

  lzo_uint outputlength_in_bytes = outputarraylength*4-16;
  lzo_uint in_len = bptr->Length; 

  lzo1x_1_12_compress( (lzo_byte *) bptr,
			in_len,  
		       (lzo_byte *)&outputarray[4],
			&outputlength_in_bytes,wrkmem);


  outputarray[0] = outputlength_in_bytes +4*BUFFERHEADERLENGTH;
  outputarray[1] =  LZO1XBUFFERMARKER;
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
olzoBuffer::~olzoBuffer()
{
  writeout();
  delete [] outputarray;
  lzo_free(wrkmem);

}

