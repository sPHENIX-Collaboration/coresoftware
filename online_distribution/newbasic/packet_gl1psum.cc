#include <packet_gl1psum.h>
#include <buffer.h>

Packet_gl1psum::Packet_gl1psum(PACKET_ptr data)
  : Packet_w4(data)
{
}

Packet_gl1psum::~Packet_gl1psum()
{
}

int *Packet_gl1psum::decode ( int *nwout)
{
  int *p,*k, *from;
  int i;
  p = new int[11];

  from = (int *)  findPacketDataStart(packet);
  
  k = p;
  if ( getHitFormat() == IDGL1PSUMOBS) 
    {
      for (i =0; i<11; i++) *k++ = buffer::i4swap(from[i]);
    }
  else
    {
      for (i =0; i<11; i++) *k++ = from[i];
    }

  p[0] = p[0] & 0xFF;
  p[1] = ( p[1] &0xF)+((p[1]>>4)&0x7)*15; 

  
  *nwout = 11;
  return p;
}

int Packet_gl1psum::iValue(const int ich)
{
  if (ich < 0 || ich >= 8) return 0;
  
  if ( decoded_data1 ==0 ) decoded_data1 = decode(&data1_length);
  
  return decoded_data1[ich+2];
}



int Packet_gl1psum::iValue(const int ich, const char *what)
{
  
  
  if ( decoded_data1 ==0 ) decoded_data1 = decode(&data1_length);
  
  if (!strcmp(what,"EVTNR"))
    {
	return decoded_data1[0] ;
    }

  else if (!strcmp(what,"BEAMCROSSID"))
    {
	return decoded_data1[1];
    }
  else if (!strcmp(what,"CROSSCTR"))
    {
	return decoded_data1[10];
    }
  else return 0;

}


void Packet_gl1psum::dump ( OSTREAM &os)
{

  int i;


  this->identify(os); 
  os << "Evt number:       " << (unsigned int) iValue(0,"EVTNR") << std::endl;
  os << "Crossing ID:      " << (unsigned int) iValue(0, "BEAMCROSSID") << std::endl;
  os << "Crossing Counter: " << (unsigned int) iValue(0, "CROSSCTR") << std::endl;
  
  for ( i = 0; i < 4; i++) 
    {
      os << std::setw(12) << iValue(i) << "  ";
    }
  os << std::endl;

  for ( i = 4; i < 8; i++) 
    {
      os << std::setw(12) << iValue(i) << "  ";
    }
  os << std::endl << std::endl;

}

