#include "packet_iddigitizerv2.h"

#include <string.h>

using namespace std;

Packet_iddigitizerv2::Packet_iddigitizerv2(PACKET_ptr data)
  :Packet_w4 (data)
{

 _nsamples = 0;

 _evtnr = 0;
 _detid = 0;
 _module_address  = 0;
 _clock = 0;
 _fem_slot = 0;
 _fem_evtnr = 0;
 _fem_clock = 0;


 _nchannels = 0;
 _is_decoded = 0;


}


Packet_iddigitizerv2::~Packet_iddigitizerv2()
{
  //  if (array) delete [][] array;
}



const int offset=9;
const int samples=12;
const int channels=64;

  
int Packet_iddigitizerv2::decode ()
{

  if (_is_decoded ) return 0;
  _is_decoded = 1;


  int *k;


  // check later  int dlength = ( getLength()-4) - getPadding();

  int *SubeventData = (int *) findPacketDataStart(packet); 

  //  _nsamples    = (SubeventData[0] >> 24 ) & 0xff;

  _evtnr           =  SubeventData[0]  & 0xffff;
  _detid           =  SubeventData[2]  & 0xffff;

  _module_address  = SubeventData[3] & 0xffff;
  _clock           = SubeventData[4] & 0xffff;
  _nsamples      = SubeventData[5] & 0xff;

  // if we are looking at older data without the nr_samples encoded,
  // they are 31 samples wide. 
  if ( _nsamples == 0xff) _nsamples = 31;

  _fem_slot        = SubeventData[6] & 0xffff;
  _fem_evtnr       = SubeventData[7] & 0xffff;
  _fem_clock       = SubeventData[8] & 0xffff;

  //  _l1_delay    = (SubeventData[0] >> 16 ) & 0xff;
  _nr_modules =  1;

  _nchannels = _nr_modules * 64;


  k = &SubeventData[offset];

  int index=0;
  for ( int ch  = 0; ch < _nchannels ; ch++)
    {
      for ( int sample = 0; sample < _nsamples ; sample++)
	{
	  array[sample][2*ch]    =   k[index++] & 0xffff;
	  array[sample][2*ch+1]  =   k[index++] & 0xffff;

	}
    }

  return 0;
}


int Packet_iddigitizerv2::iValue(const int sample, const int ch)
{
  decode();

  if ( sample >= _nsamples || sample < 0 
       || ch >= _nchannels || ch < 0 ) return 0;

  return array[sample][ch];

}

int Packet_iddigitizerv2::iValue(const int n, const char *what)
{

  decode();

  if ( strcmp(what,"CLOCK") == 0 )
  {
    return _clock;
  }

  if ( strcmp(what,"EVTNR") == 0 )
  {
    return _evtnr;
  }

  if ( strcmp(what,"SAMPLES") == 0 )
  {
    return _nsamples;
  }

  if ( strcmp(what,"NRMODULES") == 0 )
  {
    return _nr_modules;
  }


  if ( strcmp(what,"CHANNELS") == 0 )
  {
    return _nchannels;
  }

  if ( strcmp(what,"DETID") == 0 )
  {
    return _detid;
  }

  if ( strcmp(what,"MODULEADDRESS") == 0 )
  {
    return _module_address;
  }


  if ( strcmp(what,"FEMSLOT") == 0 )
  {
    return _fem_slot;
  }

  if ( strcmp(what,"FEMEVTNR") == 0 )
  {
    return _fem_evtnr;
  }

  if ( strcmp(what,"FEMCLOCK") == 0 )
  {
    return _fem_clock;
  }


  return 0;

}

void  Packet_iddigitizerv2::dump ( OSTREAM& os )
{
  identify(os);

  os << "Evt Nr:     " << iValue(0,"EVTNR") << std::endl;
  os << "Clock:      " << iValue(0,"CLOCK") << std::endl;
  os << "Channels:   " << iValue(0,"CHANNELS") << std::endl;
  os << "Samples:    " << iValue(0,"SAMPLES") << std::endl;
  os << "Det. ID:    " << iValue(0,"DETID") << std::endl;
  os << "Mod. Addr:  " << iValue(0,"MODULEADDRESS") << std::endl;
  os << "Nr Modules: " << iValue(0,"NRMODULES") << std::endl;
  os << "FEM Slot:   " << iValue(0,"FEMSLOT") << std::endl;
  os << "FEM Evt nr: " << iValue(0,"FEMEVTNR") << std::endl;
  os << "FEM Clock:  " << iValue(0,"FEMCLOCK") << std::endl;
  os << endl;

  //   0 |    8140   8136   8133   8136   8134   8133   8139   8135   8134   8134   8135   8130
  //  ch       s01    s1     s2     s3     s4     s5     s6     s7     s8      
  for ( int c = 0; c < _nchannels; c++)
    {
      os << setw(4) << c << " | ";

      for ( int s = 0; s < _nsamples; s++)
	{
	  os << setw(6) << iValue(s,c);
	}
      os << endl;
    }

}
