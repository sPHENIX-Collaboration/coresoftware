#include "packet_starscaler.h"

Packet_starscaler::Packet_starscaler(PACKET_ptr data)
  : Packet_w4 (data)
{
  s_vectorlength= 0;
  s_vector= 0;
}
  
Packet_starscaler::~Packet_starscaler()
{
#if !defined(SunOS) && !defined(OSF1)
  smap.clear();
  if ( s_vector) delete []  s_vector;
#endif
}

int  *Packet_starscaler::decode ( int *nwout)
{

#if !defined(SunOS) && !defined(OSF1)
  int dlength = getDataLength();
  s_vectorlength= dlength/2;

  

  int *p = new int [s_vectorlength];
  s_vector = new long long [s_vectorlength];
  s_vector[0] = 0;
  int *m = (int *) findPacketDataStart(packet);

  unsigned int i;
  for ( i = 0; i < s_vectorlength; i++)
    {
      int l = i*2;
      p[i] = (( m[l+1]>>8) & 0xffffff);
      long long lsb = (unsigned int) m[l];
      long long msb = (unsigned int) m[l+1];
      s_vector[i] = lsb | (( msb &0xff) <<32);
      smap[p[i]] = i;

    }

  
  *nwout = s_vectorlength;
  return p;
#else
  return 0;
#endif

}

// ------------------------------------------------------

long long  Packet_starscaler::lValue (const int chan)
{
#if !defined(SunOS) && !defined(OSF1)
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }
  


  if (chan < 0 || chan  >= (1<<24)  ) return 0;
  
  std::map <int, int>::iterator it;
  if ( ( it = smap.find(chan) ) == smap.end() )
    {
      return 0;
    }

  int index = it->second;

  if ( index < 0 || index >=  s_vectorlength )
    {

      return 0;
    }
  return s_vector[index];

#else
  return 0;
#endif

}

int Packet_starscaler::iValue (const int index)
{
#if !defined(SunOS) && !defined(OSF1)
  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }
  


  if (index < 0 || index >= s_vectorlength ) return 0;
  
  return decoded_data1[index];

#else
  return 0;
#endif

}

int  Packet_starscaler::iValue(const int ich, const char *what)
{

  if (decoded_data1 == NULL )
    {
      if ( (decoded_data1 = decode(&data1_length))==NULL)
	return 0;
    }
  
 if ( strcmp(what,"CHANNELCOUNT")==0)
    {
      return s_vectorlength -1;
    }
 return 0;

}

void Packet_starscaler::dump ( OSTREAM &os) 
{

  
  this->identify(os);
  
  int i;

  for ( i= 0; i < iValue(0,"CHANNELCOUNT"); i++)
    {
      os << i << "  " << std::setw(10) << iValue(i) <<  " | " << std::setw(12) << lValue(iValue(i)) << std::endl;
    }
} 

