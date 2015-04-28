#include "oncsSub_idsis3300.h"
#include <cstring>

oncsSub_idsis3300::oncsSub_idsis3300(subevtdata_ptr data)
  :oncsSubevent_w4 (data)
{
  samples = 0;
  wraparound = 0;

}
  
int *oncsSub_idsis3300::decode ( int *nwout)
{
  int *p;


  int i,j;
  int *SubeventData = &SubeventHdr->data;

  samples = (*SubeventData) & 0xffff;
  wraparound = (*SubeventData) >> 16;

  //  cout << "Samples: " << samples << std::endl;

  p = new int [samples*8];
  j = 0;
  for ( i = 0; i< samples*4; i++)
    {
      p[j++] = (SubeventData[i+1] >> 16) & 0x3fff;
      p[j++] = SubeventData[i+1] & 0x3fff;
    }

  *nwout = j-1;
  return p;
}

int oncsSub_idsis3300::iValue(const int ch ,const int s)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( ch < 0 || ch >7 ) return 0;
  if ( s < 0 || s >=samples ) return 0;

  return decoded_data1[ch + 8*s];

}

int oncsSub_idsis3300::iValue(const int,const char *what)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( strcmp(what,"SAMPLES") == 0 )
  {
    return samples;
  }
  // see if this is a "wraparound" event, where the samples cross
  // over the end of the memory and continue a low memory. This is
  // a debugging function since we think there might be a problem with those events. 
  if ( strcmp(what,"WRAPAROUND") == 0 )
  {
    return wraparound;
  }

  return 0;

}


void  oncsSub_idsis3300::dump ( OSTREAM& os )
{
  int i,j;

  os << "Samples: " << iValue(0,"SAMPLES");
  if ( iValue(0,"WRAPAROUND") )
    {
      os << " - wrap-around reaodut";
    }
  os << std::endl;

  for ( i = 0; i < iValue(0,"SAMPLES"); i++)
    {

      os << std::setw(6) << i << " |  "; 
      for ( j = 0; j < 8; j++)
	{
	  os << std::setw(8) << iValue(j,i) << "  ";
	}
      os << std::endl;
    }
  os << std::endl;
}

