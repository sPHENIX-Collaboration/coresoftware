#include "oncsSub_idcaenv792.h"
#include <cstring>

oncsSub_idcaenv792::oncsSub_idcaenv792(subevtdata_ptr data)
  :oncsSubevent_w4 (data)
{
  samples = 0;
}
  
int *oncsSub_idcaenv792::decode ( int *nwout)
{
  int *p;


  int i,channel;
  int *SubeventData = &SubeventHdr->data;

  samples = ( *SubeventData >> 8) & 0x3f; // get bits 8-13 
  
  //  std::cout << "Samples: " << samples << std::endl;

  p = new int [32]; // we always get 32 channels, even if some channels are 0.
  memset ( p, 0, 32*4);  // clear the array


  for ( i = 0; i< samples; i++)
    {
      if ( ( ( SubeventData[i+1] >> 24) & 7) == 0)
	{
	  if ( getHitFormat() == IDCAENV785N)  // the 785N has the channel info in bits 17...20
	    {
	      channel = (SubeventData[i+1] >> 17) & 0x1f;
	    }
	  else
	    {
	      channel = (SubeventData[i+1] >> 16) & 0x3f;
	    }

	  p[channel] = SubeventData[i+1] & 0xfff;

	  //	  std::cout << std::hex << SubeventData[i+1] << "  " << channel << "  " << p[channel] << std::dec << std::endl;
	}
    }
  evnr = SubeventData[samples +1] & 0xffffff;
  
  *nwout = 32;
  return p;
}

int oncsSub_idcaenv792::iValue(const int ch)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( ch < 0 || ch >31 ) return 0;

  return decoded_data1[ch];

}

int oncsSub_idcaenv792::iValue(const int,const char *what)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( strcmp(what,"SAMPLES") == 0 )
  {
    return samples;
  }

  if ( strcmp(what,"EVNR") == 0 )
  {
    return evnr;
  }

  return 0;

}


void  oncsSub_idcaenv792::dump ( OSTREAM& os )
{
  int i,j;
  //  int *SubeventData = &SubeventHdr->data;
  
  os << "Samples: " << iValue(0,"SAMPLES") << std::endl;
  os << "Evt Nr:  " << iValue(0,"EVNR") << std::endl;
  
  
  //  for ( i = 0; i <= samples ; i++)
  //   {
  //    os << std::hex << std::setw(3) << i << "  " << (( SubeventData[i] >> 24) & 7) << "   "
  //	 << (  SubeventData[i] & 0xffffff ) << std::dec << std::endl ;
  //  }
  
  
  for ( i = 0; i < 32 ; i+=4)
    {
      
      os << std::setw(6) << i << " |  "; 
      for ( j = 0; j < 4; j++)
	{
	  os << std::setw(8) << iValue(j+i) << "  ";
	}
      os << std::endl;
    }
  os << std::endl;
}

