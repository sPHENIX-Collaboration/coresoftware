#include "oncsSub_iddigitizerv1.h"

#include <string.h>

using namespace std;

oncsSub_iddigitizerv1::oncsSub_iddigitizerv1(subevtdata_ptr data)
  :oncsSubevent_w4 (data)
{

 _nsamples = 0;
 _l1_delay = 0;
 _slot_nr = 0;
 _nr_modules = 0;
 _clock = 0;
 _evtnr = 0;
 _nchannels = 0;
 _slot_nr_from_data = 0;
 _is_decoded = 0;
 // array = 0;

}


oncsSub_iddigitizerv1::~oncsSub_iddigitizerv1()
{
  //  if (array) delete [][] array;
}



const int offset=3;
const int samples=12;
const int channels=64;

  
int oncsSub_iddigitizerv1::decode ()
{

  if (_is_decoded ) return 0;
  _is_decoded = 1;


  int *k;


  // check later  int dlength = ( getLength()-4) - getPadding();

  int *SubeventData = &SubeventHdr->data;

  _nsamples    = (SubeventData[0] >> 24 ) & 0xff;
  _l1_delay    = (SubeventData[0] >> 16 ) & 0xff;
  _slot_nr     = (SubeventData[0] >> 8  ) & 0xff;
  _nr_modules =  SubeventData[0] & 0xff;

  _nchannels = _nr_modules * 64;

  _slot_nr_from_data     = (SubeventData[1] >> 16  ) & 0xff;

  _clock                 = (SubeventData[2] >> 16  ) & 0xffff;
  _evtnr                 =  SubeventData[2]  & 0xffff;

  //sanity check
  if ( _nr_modules > 5 || _nsamples >31 || _slot_nr > 22) 
    {
      return 0;
    }


  k = &SubeventData[offset];

  int index=0;
  for ( int ch  = 0; ch < _nchannels/2 ; ch++)
    {
      for ( int sample = 0; sample < _nsamples ; sample++)
	{
	  array[sample][2*ch]    =   k[index] & 0xffff;
	  array[sample][2*ch+1]  =  (k[index]>>16) & 0xffff;
	  index++;
	}
    }

  return 0;
}


int oncsSub_iddigitizerv1::iValue(const int sample, const int ch)
{
     decode();

  if ( sample >= _nsamples || sample < 0 
       || ch >= _nchannels || ch < 0 ) return 0;

  return array[sample][ch];

}

int oncsSub_iddigitizerv1::iValue(const int n, const char *what)
{

  decode();

  if ( strcmp(what,"SAMPLES") == 0 )
  {
    return _nsamples;
  }

  if ( strcmp(what,"CHANNELS") == 0 )
  {
    return _nchannels;
  }

  if ( strcmp(what,"L1DELAY") == 0 )
  {
    return _l1_delay;
  }

  if ( strcmp(what,"SLOTNR") == 0 )
  {
    return _slot_nr;
  }

  if ( strcmp(what,"NRMODULES") == 0 )
  {
    return _nr_modules;
  }

  if ( strcmp(what,"CLOCK") == 0 )
  {
    return _clock;
  }

  if ( strcmp(what,"EVTNR") == 0 )
  {
    return _evtnr;
  }

  return 0;

}

void  oncsSub_iddigitizerv1::dump ( OSTREAM& os )
{
  identify(os);

  os << "Evt Nr:     " << iValue(0,"EVTNR") << std::endl;
  os << "Clock:      " << iValue(0,"CLOCK") << std::endl;
  os << "Channels:   " << iValue(0,"CHANNELS") << std::endl;
  os << "Samples:    " << iValue(0,"SAMPLES") << std::endl;
  os << "L1 Delay:   " << iValue(0,"L1DELAY") << std::endl;
  os << "Slot Nr:    " << iValue(0,"SLOTNR") << std::endl;
  os << "Nr Modules: " << iValue(0,"NRMODULES") << std::endl;
  os << endl;

  //   0 |    8140   8136   8133   8136   8134   8133   8139   8135   8134   8134   8135   8130
  //  ch       s01    s1     s2     s3     s4     s5     s6     s7     s8      
  for ( int c = 0; c < _nchannels; c++)
    {
      os << setw(4) << c << " | ";

      for ( int s = 0; s < _nsamples; s++)
	{
	  os << setw(8) << iValue(s,c);
	}
      os << endl;
    }

}
