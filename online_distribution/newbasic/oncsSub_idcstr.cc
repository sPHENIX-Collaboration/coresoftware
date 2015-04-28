#include "oncsSub_idcstr.h"
#include <string.h>

oncsSub_idcstr::oncsSub_idcstr(subevtdata_ptr data)
  :oncsSubevent_w1 (data)
{
  sarray = 0;
  allocated_length = 0;
}


oncsSub_idcstr::~oncsSub_idcstr()
{
  if ( sarray) delete [] sarray;
}

int oncsSub_idcstr::iValue ( const int i)
{
  int n;
  if ( ! sarray ) decode ( &n);

  if ( i < 0 || i >= allocated_length)
    {
      return 0;
    }

  int c = sarray[i];
  return c;

}

void oncsSub_idcstr::dump ( OSTREAM &os)
{
  //  identify(os);
  int n;
  if ( ! sarray ) decode ( &n);

  int i = 1;
  int c;

  os.write( (const char *) sarray, allocated_length);

  //  os << std::endl;
  
  // c = iValue(0);
  // while ( c) 
  //   {
  //     os << (char ) c;
  //     c = iValue(i++);
  //   }

  os << std::flush;

}


int *oncsSub_idcstr::decode ( int *nwout)

{
  int dlength = ( getLength()-4)*4 - getPadding();
  unsigned char *SubeventData = ( unsigned char * ) &SubeventHdr->data;
  sarray = new unsigned char[dlength+1];
  memcpy ( sarray, SubeventData,dlength); 
  sarray[dlength] = 0;
  allocated_length = dlength;

  return 0;
}

