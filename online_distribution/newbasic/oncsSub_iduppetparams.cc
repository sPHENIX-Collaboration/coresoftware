#include "oncsSub_iduppetparams.h"
#include <cstring>

#include <arpa/inet.h>

using namespace std;

oncsSub_iduppetparams::oncsSub_iduppetparams(subevtdata_ptr data)
  :oncsSubevent_w4 (data)
{

}
  
oncsSub_iduppetparams::~oncsSub_iduppetparams()
{

}


int *oncsSub_iduppetparams::decode ( int *nwout)
{

  unsigned int *SubeventData = (unsigned int *) &SubeventHdr->data;

  int *p = new int [42];

  for ( int i = 0; i < 41; i++ ) 
    {
      p[i] = SubeventData[i];

    }
  
  *nwout = 41;
  return p;
}
 

int oncsSub_iduppetparams::iValue(const int ch)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( ch < 0 || ch >40 ) return 0;

  return decoded_data1[ch];

}


 
int oncsSub_iduppetparams::iValue(const int ich,const char *what)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  // register 1 - asic enable
  if ( strcmp(what,"ASICENABLED") == 0 )
  {
    if ( ich <0 || ich >=24) return 0;
    if ( ( decoded_data1[1] >> ich) &1 ) return 0;
    return 1;
  }

  if ( strcmp(what,"UDPHEADER1") == 0 )
  {
    return decoded_data1[9];
  }

  if ( strcmp(what,"UDPHEADER2") == 0 )
  {

    return decoded_data1[10];
  }

  if ( strcmp(what,"VTHHIGH") == 0 )
  {
    if ( ich < 0 || ich >= 6*24) return 0;
    return ( decoded_data1[11+ich] & 0xfff);
  }

  if ( strcmp(what,"VTHLOW") == 0 )
  {
    if ( ich < 0 || ich >= 6*24) return 0;
    return ( (decoded_data1[11+ich]>>16) & 0xfff);
  }

  if ( strcmp(what,"CLOCKPHASE") == 0 )
  {
    return decoded_data1[35];
  }

  if ( strcmp(what,"UDPTIMEOUT") == 0 )
  {
    return decoded_data1[36];
  }

  if ( strcmp(what,"CLOCKCORRECTION") == 0 )
  {
    return (decoded_data1[38] & 0xffffff);
  }

  if ( strcmp(what,"FIRMWAREVERSION") == 0 )
  {
    return decoded_data1[40];
  }


  return 0;

}


void  oncsSub_iduppetparams::dump ( OSTREAM& os )
{

 

  os << " Firmware version: 0x" << hex << iValue(0,"FIRMWAREVERSION") << dec << std::endl;
  os << " Enabled Chips:    0x" << hex << iValue(1) << dec << std::endl;
  os << " UDP Header 1:     0x" << hex << iValue(0,"UDPHEADER1") << dec << "  " << setw(6) << iValue(0,"UDPHEADER1") << std::endl;
  os << " UDP Header 2:     0x" << hex << iValue(0,"UDPHEADER2") << dec  << "  " << setw(6) << iValue(0,"UDPHEADER2") << std::endl;
  os << " UDP Timeout:      " <<  iValue(0,"UDPTIMEOUT") << std::endl;
  os << " Clock Phase:      0x" << hex << iValue(0,"CLOCKPHASE") << dec << std::endl;
  os << " Clock Correction: 0x" << hex << iValue(0,"CLOCKCORRECTION") << dec << std::endl;
  
  int i, k;

  os << std::endl;
  
  os << "                 VTH Low (* = ASIC enabled) " << endl;
  os << "   row |     T5      T4      T3      T2      T1      T0" << endl;
  os << "   ----+------------------------------------------" << endl;
 
  int row = 3;

  for ( i = 23; i >= 20 ; i--)
    {

      os << std::setw(5) << row-- << "  |  "; 
      for ( k = i;  k >= 0 ; k-=4)
	{

	  os << std::setw(5) << iValue(k,"VTHLOW");
	  if ( iValue(k,"ASICENABLED") )
	    {
	      cout << "*  ";
	    }
	  else
	    {
	      cout << "   ";
	    }
	      
	}	
      os << std::endl;
    }

  os << std::endl;

  os << "                 VTH High (* = ASIC enabled) " << endl;

  row = 3;
  for ( i = 23; i >= 20 ; i--)
    {

      os << std::setw(5) << row-- << "  |  "; 
      for ( k = i;  k >= 0 ; k-=4)
	{

	  os << std::setw(5) << iValue(k,"VTHHIGH");
	  if ( iValue(k,"ASICENABLED") )
	    {
	      cout << "*  ";
	    }
	  else
	    {
	      cout << "   ";
	    }
	      
	}	
      os << std::endl;
    }

  os << "          ( as seen on the motherboard with towers up )" << std::endl << std::endl;

}

