#include "oncsSub_iddrs4v1.h"
#include <string.h>
#include <iostream>
#include <iomanip>

using namespace std;

oncsSub_iddrs4v1::oncsSub_iddrs4v1(subevtdata_ptr data)
  :oncsSubevent_w4 (data)
{
  wave = 0;
  //  int dummy;
  // decode(&dummy);
}

oncsSub_iddrs4v1::~oncsSub_iddrs4v1()
{
  if ( wave) delete [] wave;
}

int *oncsSub_iddrs4v1::decode ( int *nwout)
{

  int dlength = ( getLength()-4) - getPadding();


  samples = SubeventHdr->data & 0xffff;
  enabled_channelmask = (SubeventHdr->data >> 16) & 0xf;

  float *d = (float *) &SubeventHdr->data;
  d++;
  int len = samples;
  for ( int i =0; i<4; i++)
    {
      if ( enabled_channelmask & ( 1<<i) ) 
	{
	  channel_offset[i] = len;
	  len += samples;
	}
    }
  wave = new float[len];
  memcpy ( wave, d, len * sizeof(float));
    
  *nwout = 0;
  return 0;
}

float oncsSub_iddrs4v1::rValue ( const int isample,const char * what)
{
  if ( isample < 0 || isample > samples) return 0;
  
  if ( ! wave) 
    {
      int dummy;
      decode (&dummy);
    }


  if ( strcmp(what,"TIME") == 0 )
  {
    return wave[isample]; // the first 1024 sample are the time base
  }
  return 0;

}

int oncsSub_iddrs4v1::iValue ( const int i,const char * what)
{
  if ( ! wave) 
    {
      int dummy;
      decode (&dummy);
    }

  if ( strcmp(what,"SAMPLES") == 0 )
  {
    return samples; 
  }

  else if ( strcmp(what,"ENABLED") == 0 )
  {
    return ( (enabled_channelmask & ( 1 <<  i)) !=0);
  }

  return 0;
}

float oncsSub_iddrs4v1::rValue ( const int isample,const int ich)
{
  if ( isample < 0 || isample > 1023) return 0;
  if ( ich < 0 || ich > 4) return 0;

  if (ich == 4) return rValue(isample,"TIME");
  
  if ( ! wave) 
    {
      int dummy;
      decode (&dummy);
    }

  return wave[channel_offset[ich] +isample]; 

}


void  oncsSub_iddrs4v1::dump ( OSTREAM& os )
{
  identify(os);
  int i;
  
  os << " Samples " << iValue(0, "SAMPLES")  << "  enabled channnels: ";
  for ( i = 0; i< 4; i++)
    {
      if ( iValue (i, "ENABLED" ) ) 
	{
	  os << "  " << i; 
	}
    }
  os << endl;

  os << "   ch |       time";
  if ( enabled_channelmask & 1)  os << "          ch0";
  if ( enabled_channelmask & 2)  os << "          ch1";
  if ( enabled_channelmask & 4)  os << "          ch2";
  if ( enabled_channelmask & 8)  os << "          ch3";
  os << endl;

  for ( i = 0; i < samples; i++)
    {
      os << setw(5) << i << " | ";
      os << setw(10) << rValue(i, "TIME") << "   ";
      if ( enabled_channelmask & 1)  os << setw(10) << rValue(i, 0) << "   ";
      if ( enabled_channelmask & 2)  os << setw(10) << rValue(i, 1) << "   ";
      if ( enabled_channelmask & 4)  os << setw(10) << rValue(i, 2) << "   ";
      if ( enabled_channelmask & 8)  os << setw(10) << rValue(i, 3) << "   ";
      os << endl;
    }
}


