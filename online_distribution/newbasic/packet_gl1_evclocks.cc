#include <packet_gl1_evclocks.h>
#include <string.h>

Packet_gl1_evclocks::Packet_gl1_evclocks(PACKET_ptr data)
  : Packet_w4 (data){ }

Packet_gl1_evclocks::~Packet_gl1_evclocks()
{

}

int *Packet_gl1_evclocks::decode ( int *nwout)
{
  int *p,*k;
  int i;
  int dlength = getDataLength();

  unsigned int* buf = (unsigned int *) findPacketDataStart(packet);
  if (buf == 0) return 0;
 
  p = new int[dlength];
  k = p;
  for (i =0; i<dlength; i++) *k++ = buf[i];

  *nwout = dlength;
  return p;

}

int Packet_gl1_evclocks::iValue(const int ich, const char *what)
{

  
  if (strcmp(what,"EVCLOCK")==0)
    {

      if( (ich<0) || (ich>3) ) return -1;  
      unsigned int* buf = (unsigned int *) findPacketDataStart(packet);
      if (buf == 0) 
	return -1;
      else
	return buf[ich];
    }
  else if (strcmp(what,"PARVECT")==0)
    {

      if( (ich<0) || (ich>3) || (getDataLength()<8) ) return -1;  
      unsigned int* buf = (unsigned int *) findPacketDataStart(packet);
      if (buf == 0) 
	return -1;
      else
	return buf[4+ich];
    }

  else return 0;

}

int Packet_gl1_evclocks::iValue(const int ich, const  int what)
{

  switch (what) {

  case EVCLOCK:
    {
      
      if( (ich<0) || (ich>3) ) return -1;  
      unsigned int* buf = (unsigned int *) findPacketDataStart(packet);
      if (buf == 0) 
	return -1;
      else
	return buf[ich];

      break;
    }
  case PARVECT:
    {
      
      if( (ich<0) || (ich>3) || (getDataLength()<8) ) return -1;  
      unsigned int* buf = (unsigned int *) findPacketDataStart(packet);
      if (buf == 0) 
	return -1;
      else
	return buf[4+ich];
      
      break;
    }
  default:
    return -1; 

  }

}

void Packet_gl1_evclocks::dump ( OSTREAM &os)
{
  int j,l,m,dlength;
  unsigned int delta = 0; 

  dlength = getDataLength();
  this->identify(os); 

  unsigned int* buf = (unsigned int *) findPacketDataStart(packet);
  if (buf == 0) return; 

  os << std::endl;
  os << "GL1 Clock Counter Event Difference Data Packet:" << std::endl;
  os << " Clock counter for this event (N) = " << buf[0]; 
  if(dlength>=8) 
    os << ",  Partition = 0x" << std::hex << buf[4] << std::dec << std::endl;
  else 
    os << std::endl; 
 
  os << " Clock counter for N-1 event = " << buf[1]; 
  if(dlength>=8) os << ",  Partition = 0x" << std::hex << buf[5] << std::dec; 
  if(buf[0]>buf[1]) 
    delta = buf[0]-buf[1]; 
  else
    delta = 0xFFFFFFFF-buf[1] + buf[0] + 1; // Add 1 clock = 0x0
  os << std::endl << "     (delta from this event = " << delta << " clocks, " << delta*0.106 << " microseconds)" << std::endl; 

  os << " Clock counter for N-2 event = " << buf[2];
  if(dlength>=8) os << ",  Partition = 0x" << std::hex << buf[6] << std::dec; 
  if(buf[0]>buf[2]) 
    delta = buf[0]-buf[2]; 
  else
    delta = 0xFFFFFFFF-buf[2] + buf[0] + 1; // Add 1 clock = 0x0
  os << std::endl << "     (delta from this event = " << delta << " clocks, " << delta*0.106 << " microseconds)" << std::endl; 

  os << " Clock counter for N-3 event = " << buf[3];
  if(dlength>=8) os << ",  Partition = 0x" << std::hex << buf[7] << std::dec; 
  if(buf[0]>buf[3]) 
    delta = buf[0]-buf[3]; 
  else
    delta = 0xFFFFFFFF-buf[3] + buf[0] + 1; // Add 1 clock = 0x0
  os << std::endl << "     (delta from this event = " << delta << " clocks, " << delta*0.106 << " microseconds)" << std::endl; 
 
  os << std::endl; 
  
  // Finally, a generic HEX dump:

  char oldFill;
  oldFill=os.fill('0');

  j = 0;
  int* k=( int *) findPacketDataStart(packet);
  if ( k ==0 ) return;

  while (1)
  {
    os << std::endl << std::dec << SETW(5) << j << " |  ";
    for (l=0;l<4;l++)
    {
      os << std::hex << SETW(8) << k[j++] << " ";
      if (j>=dlength) break;
    }
    if (j>=dlength) break;
  }	
  os << std::endl;
  for(m=0;m<54;m++) os<<"=";
  os << std::dec << std::endl;
  os.fill(oldFill);


}

