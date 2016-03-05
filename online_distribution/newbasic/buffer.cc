#include "buffer.h"
#include "gzbuffer.h"
#include "lzobuffer.h"
#include "prdfBuffer.h"
#include "oncsBuffer.h"


buffer::buffer ()
{

}


buffer::~buffer()
{}


int buffer::makeBuffer( PHDWORD *bp, const int allocatedsize, buffer **bptr)
{
  if ( bp[1]== GZBUFFERMARKER || buffer::u4swap(bp[1])== GZBUFFERMARKER )
    {
      *bptr = new gzbuffer(bp, allocatedsize );
      return 0;
    }

  else if ( bp[1]== LZO1XBUFFERMARKER || buffer::u4swap(bp[1])== LZO1XBUFFERMARKER )
    {
      *bptr = new lzobuffer ( bp, allocatedsize );
      return 0;
    }

  else if ( bp[1]== ONCSBUFFERMARKER || buffer::u4swap(bp[1])== ONCSBUFFERMARKER )
    {
      *bptr = new oncsBuffer ( bp, allocatedsize );
      return 0;
    }

  else if ( bp[1]== BUFFERMARKER || buffer::u4swap(bp[1])== BUFFERMARKER )
    {
      *bptr = new prdfBuffer ( bp, allocatedsize );
      return 0;
    }

  *bptr = 0;
  return -1;
}



// ---------------------------------------------------------
int buffer::i4swap(const int in)
{
  union 
  {
    int i4;
    char c[4];
  } i,o;

  i.i4 = in;
  o.c[0] = i.c[3];
  o.c[1] = i.c[2];
  o.c[2] = i.c[1];
  o.c[3] = i.c[0];
  return o.i4;
}

// ---------------------------------------------------------
unsigned int buffer::u4swap(const unsigned int in)
{
  union 
  {
    unsigned int i4;
    char c[4];
  } i,o;

  i.i4 = in;
  o.c[0] = i.c[3];
  o.c[1] = i.c[2];
  o.c[2] = i.c[1];
  o.c[3] = i.c[0];
  return o.i4;
}


// ---------------------------------------------------------
int buffer::i22swap(const int in)
{
  union 
  {
    int i4;
    char c[4];
  } i,o;

  i.i4 = in;
  o.c[0] = i.c[1];
  o.c[1] = i.c[0];
  o.c[2] = i.c[3];
  o.c[3] = i.c[2];
  return o.i4;
}

// ---------------------------------------------------------
short buffer::i2swap(const  short in)
{
  union 
  {
    short i2;
    char c[2];
  } i,o;

  i.i2 = in;
  o.c[0] = i.c[1];
  o.c[1] = i.c[0];

  return o.i2;
}


