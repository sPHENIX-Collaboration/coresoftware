#include "packet_idcstr.h"
#include <stdio.h>

Packet_idcstr::Packet_idcstr(PACKET_ptr data)
  : Packet_w1 (data)
{
  sarray = 0;
}
  
Packet_idcstr::~Packet_idcstr()
{
  if ( sarray) delete [] sarray;
}


int Packet_idcstr::iValue(const int i)
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



int *Packet_idcstr::decode ( int *nwout)
{

  int dlength = getDataLength();

  //   std::cout << __FILE__ << "  " << __LINE__ << " datalength: " 
  //	      << getDataLength() << " padding " << getPadding() << std::endl;

  unsigned char *SubeventData = ( unsigned char * ) findPacketDataStart(packet);
  sarray = new unsigned char[dlength+1];
  memcpy ( sarray, SubeventData,dlength); 
  sarray[dlength] = 0;
  allocated_length = dlength;

  return 0;

}

void Packet_idcstr::dump ( OSTREAM &os)
{

  int n;
  if ( ! sarray ) decode ( &n);

  os.write( (const char *) sarray, allocated_length);
  os << std::flush;

}

