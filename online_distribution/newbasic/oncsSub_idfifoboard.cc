#include "oncsSub_idfifoboard.h"
#include <cstring>

oncsSub_idfifoboard::oncsSub_idfifoboard(subevtdata_ptr data)
  :oncsSubevent_w4 (data)
{
  samples = 0;
  d_time = 0;
  d_timelength = 0;
}
  
oncsSub_idfifoboard::~oncsSub_idfifoboard()
{
  if ( d_timelength )
    {
      delete [] d_time;
      d_time = 0;
    }
}


int *oncsSub_idfifoboard::decode ( int *nwout)
{
  int *p;


  int i;
  unsigned int *SubeventData = (unsigned int *) &SubeventHdr->data;
  SubeventData++;
  int dlength = getLength()-4 - getPadding() -1;
  samples = dlength /2;

  p = new int [samples];  //block id

  decoded_data2 = new int [samples]; //  crystal id
  data2_length = samples;

  decoded_data3 = new int [samples];  // event type
  data3_length = samples;

  d_time = new double [samples];     // time
  d_timelength = samples;

  for ( i = 0; i < samples; i++) 
    {
      p[i] = ((SubeventData[2*i+1] >> 24) & 0xf);
      decoded_data2[i] = ((SubeventData[2*i+1] >> 19) & 0x1f);
      decoded_data3[i] = ((SubeventData[2*i+1] >> 11) & 0xff);
      unsigned long long xx =  ( SubeventData[2*i+1] & 0x7ff);
      xx = ( xx << 32);
      xx += SubeventData[2*i];
      //      std::cout << std::hex <<SubeventData[2*i+1] << "  " << SubeventData[2*i] << "  " <<  xx << std::dec << std::endl;
      d_time[i] = xx;
    }

  *nwout = samples;
  return p;
}

//int oncsSub_idfifoboard::iValue(const int ch)
  //{

//  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

//  if ( ch < 0 || ch >31 ) return 0;

//  return decoded_data1[ch];

//}

int oncsSub_idfifoboard::iValue(const int ich,const char *what)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( strcmp(what,"SAMPLES") == 0 )
  {
    return samples;
  }

  if ( strcmp(what,"BLOCKID") == 0 )
  {
    if ( ich < 0 || ich >= samples) return 0;
    return decoded_data1[ich];
  }

  if ( strcmp(what,"CRYSTALID") == 0 )
  {
    if ( ich < 0 || ich >= samples) return 0;
    return decoded_data2[ich];
  }

  if ( strcmp(what,"EVENTTYPE") == 0 )
  {
    if ( ich < 0 || ich >= samples) return 0;
    return decoded_data3[ich];
  }

  return 0;

}

double oncsSub_idfifoboard::dValue(const int ich)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( ich < 0 || ich >= samples) return 0;
    return d_time[ich];
}


void  oncsSub_idfifoboard::dump ( OSTREAM& os )
{
  int i;
  //  int *SubeventData = &SubeventHdr->data;
  int is =  iValue(0,"SAMPLES");

  //  os << "Samples: " << iValue(0,"SAMPLES") << std::endl;

  os << " sample  B id    C id   time   ( " << is << " Samples )" << std::endl;
  
  for ( i = 0; i < is ; i++)
    {
      int prec = os.precision();

      os << std::setw(5) << i << "  |  "; 
 
      os << std::setw(4) << iValue(i,"BLOCKID") << "  ";
      os << std::setw(4) << iValue(i,"CRYSTALID") << "  ";
      //      os << std::setw(6) << std::hex << "0x" << iValue(i,"EVENTTYPE") 
      //	 << std::dec<< "    ";
      os << std::setw(24) << std::setprecision(24) << dValue(i)
	 << std::setprecision(prec) ;
	
      os << std::endl;
    }
  os << std::endl;
}

