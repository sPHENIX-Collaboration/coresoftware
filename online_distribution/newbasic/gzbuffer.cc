#include "gzbuffer.h"
#include <zlib.h>

// the constructor first ----------------
gzbuffer::gzbuffer (PHDWORD *array , const int length )

{
  is_good =1;
  bufferarray=0;
  theBuffer=0;

  uLong bytes; 
  uLongf outputlength_in_bytes;
  if (array[1] == GZBUFFERMARKER )
    {
      bytes = array[0]-4*BUFFERHEADERLENGTH;
      outputlength_in_bytes = array[3];
    }
  else if ( u4swap(array[1]) == GZBUFFERMARKER)
    {
      bytes = i4swap(array[0])-16;
      outputlength_in_bytes = i4swap(array[3]);
    }

  else
    {
 	COUT << " wrong buffer" << std::endl;
	is_good = 0;
	return;
    }


  int outputlength =  (outputlength_in_bytes+3)/4;
  bufferarray = new PHDWORD[outputlength];
  uLongf i =  outputlength_in_bytes;

  uncompress ( (Bytef*) bufferarray, &i,  (Bytef*) &array[4], bytes); 

  theBuffer = new prdfBuffer(bufferarray, outputlength);

}

// ---------------------------------------------------------
Event * gzbuffer::getEvent()
{
  return theBuffer->getEvent();
}

// ---------------------------------------------------------

gzbuffer::~gzbuffer()
{
  if ( bufferarray) delete [] bufferarray;
  if ( theBuffer) delete theBuffer;
}


