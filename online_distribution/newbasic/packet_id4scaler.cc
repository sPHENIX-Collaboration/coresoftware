#include "packet_id4scaler.h"

Packet_id4scaler::Packet_id4scaler(PACKET_ptr data)
  : Packet_w4 (data){}
  
int *Packet_id4scaler::decode ( int *nwout)
{
  int *p,*k;

  int dlength = getDataLength();

  //  int status = decode_id4scaler( temp
  //			      ,(int *)  findPacketDataStart(packet) 
  //			      ,dlngth
  //			      ,MAX_OUTLENGTH, &olength);


  
  int * from = (int *)  findPacketDataStart(packet);
  if (from == 0) 
    {
      *nwout = 0;
      return 0;
    }
  int i;
  p = new int[dlength];
  k = p;
  for (i =0; i<dlength; i++) *k++ = from[i];
  *nwout = dlength;
  return p;
}

int Packet_id4scaler::iValue(const int channel,const char *what)
{
  int *k = (int *) findPacketDataStart(packet);
  if (k == 0) 
    {
      return 0;
    }

  int triggermasklength = k[2];
  int scalerlength = k[3+triggermasklength];

  int rawIndex    = 3+triggermasklength +1;
  int lifeIndex   = rawIndex + scalerlength;
  int scaledIndex = lifeIndex + scalerlength;
  int stringIndex = scaledIndex + scalerlength;

  if (strcmp(what,"NUMBERSCALERS") == 0)  // 
    {			
      return scalerlength;
    }

  else if (strcmp(what,"NUMBERMASKS") == 0)  // 
    {			
      return triggermasklength;
    }

  else if (strcmp(what,"TRIGGERMASK") == 0)  // 
    {			
      if (channel <0 || channel >=triggermasklength) return 0;
      return k[3+channel];
    }

  else if (strcmp(what,"RAWSCALERS") == 0)  // 
    {			
      if (channel <0 || channel >=scalerlength) return 0;
      return k[rawIndex+channel];
    }

  else if (strcmp(what,"LIFESCALERS") == 0 ||  
	   strcmp(what,"LIVESCALERS") == 0)  // 
    {			
      if (channel <0 || channel >=scalerlength) return 0;
      return k[lifeIndex+channel];
    }

  else if (strcmp(what,"SCALEDSCALERS") == 0)  // 
    {			
      if (channel <0 || channel >=scalerlength) return 0;
      return k[scaledIndex+channel];
    }

  else if ( strcmp(what,"TIMESTRING")==0)
    {
      if  (channel <0 || channel >=  k[stringIndex] ) return EOF;
      char *c = (char *)  &k[stringIndex+1];
      return c[channel];
    }

  std::cout << "packet_id4scaler::iValue error unknown datum: " << what << std::endl;
  return 0;
}



void Packet_id4scaler::dump ( OSTREAM& os) 
{
  int i;
  this->identify(os);
  int *k = (int *) findPacketDataStart(packet);
  if (k == 0) 
    {
      return;
    }

#ifdef WIN32 
  ULONGLONG *beamclock = (ULONGLONG  * ) k;
  os << "Beamclock " << std::hex << (ULONG) *beamclock << std::dec << std::endl;
#else
  unsigned long long *beamclock = ( unsigned long long * ) k;
  os << "Beamclock " << std::hex <<  *beamclock << std::dec << std::endl;
#endif

  int triggermasklength = k[2];
  int scalerlength = k[3+triggermasklength];
  int rawIndex    = 3+triggermasklength +1;
  int lifeIndex   = rawIndex + scalerlength;
  int scaledIndex = lifeIndex + scalerlength;
  int stringIndex = scaledIndex + scalerlength;

  for (i=0; i< triggermasklength ; i++)
    {
      os << "Tiggermask " << i << "   " << std::hex << iValue(0,"TRIGGERMASK") << std::dec << std::endl;
    }

  os << "Scalers: raw   |   live   |   scaled" <<std::endl;

  for (i = 0; i< iValue(0,"NUMBERSCALERS");  i++)
    {
      os << SETW(3) << i << "   "  
	 << SETW(9) << iValue(i,"RAWSCALERS") << "  " 
	 << SETW(9) << iValue(i,"LIFESCALERS") << "  " 
	 << SETW(9) << iValue(i,"SCALEDSCALERS") << "  " 
	 << std::endl;
    }

  os << "time string: ";
  char *c = (char *)  &k[stringIndex+1];
  for ( i = 0; i< k[stringIndex]; i++)
    {
      os << *c++;
    }

  os << std::endl;
}

