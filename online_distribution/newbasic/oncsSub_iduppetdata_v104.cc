#include "oncsSub_iduppetdata_v104.h"
#include <cstring>

#include <arpa/inet.h>

using namespace std;

oncsSub_iduppetdata_v104::oncsSub_iduppetdata_v104(subevtdata_ptr data)
  :oncsSubevent_w4 (data)
{
  samples = 0;
  d_time = 0;
  d_timelength = 0;
  lval = 0;
  tval = 0;
  serialnumber = 0;
  udpheader1 = 0;
  udpheader2 = 0;
  sysword = 0;
  eventrate = 0;
  ratcap_upenn = 0;
  clock_sel = 0;
}
  
oncsSub_iduppetdata_v104::~oncsSub_iduppetdata_v104()
{
  if ( d_timelength )
    {
      delete [] d_time;
      d_time = 0;
    }
  if ( lval)
    {
      delete [] lval;
      lval = 0;
    }
  if ( tval)
    {
      delete [] tval;
      tval = 0;
    }
}


int *oncsSub_iduppetdata_v104::decode ( int *nwout)
{
  int *p;

  unsigned long long fine_timestamp;


  int i;

  unsigned int *SubeventData = (unsigned int *) &SubeventHdr->data;

  unsigned short *sdata = ( unsigned short *) &SubeventData[4];

  int dlength = getLength()-4 - getPadding();

  samples = (dlength - 3) /2;

  p = new int [samples];             // block id

  decoded_data2 = new int [samples]; //  crystal id
  data2_length = samples;

  decoded_data3 = new int [samples];  // rfgate
  data3_length = samples;

  decoded_data4 = new int [samples];  // cardiac
  data4_length = samples;

  lval  = new long long [samples];  // raw samples
  tval  = new long long [samples];  // binary time

  d_time = new double [samples];     // time
  d_timelength = samples;

  serialnumber = ntohl(SubeventData[0]);
  udpheader1   = ntohl(SubeventData[1]);
  udpheader2   = ntohl(SubeventData[2]);
  sysword      = ntohl(SubeventData[3]);

  eventrate = sysword & 0xfffffff;
  ratcap_upenn = ( sysword >> 28 ) & 1;
  clock_sel = ( sysword >> 29 ) & 1;



  int k = 0;

  for ( i = 0; i < 2*(dlength-4); i+=4 ) 
    {
      
      unsigned long long y;
      unsigned long long x64 = ntohs( sdata[i+3]);
      x64 <<= 48;

      y =  ntohs( sdata[i+2]);
      y <<=32;
      x64 |= y;

      y =  ntohs( sdata[i+1]);
      y <<=16;
      x64 |= y;

      y =  ntohs( sdata[i]);
      x64 |= y;

      fine_timestamp = (x64 >> 46) & 0x3f;
      unsigned long long t = x64 & 0x3fffffffffffL;
      unsigned int crystalid = ( x64 >> 52) & 0x1f;
      unsigned int block = ( x64 >> 57) & 0x1f;
      unsigned int rfgate =   ( x64 >> 62) & 0x1;
      unsigned int cardiac =   ( x64 >> 63) & 0x1;

      p[k] = block;
      decoded_data2[k] = crystalid;
      lval[k] = x64;
      tval[k] = t;
      d_time[k] = t;
      decoded_data3[k] = rfgate;
      decoded_data4[k] = cardiac;
      k++;

      //      cout << k<< endl;

      // cout << hex << x << dec << "   "  
      // 	   << " time:  " << t 
      // 	   << " block: " << block 
      // 	   << " chan:  " << channel 
      // 	   << " rfgate " << rfgate
      // 	   << " Cardiac " << cardiac
      // 	   <<  endl;
    }
  
  //  samples = k;
  *nwout = k;
  return p;
}
  
int oncsSub_iduppetdata_v104::iValue(const int ich,const char *what)
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

  if ( strcmp(what,"RFGATE") == 0 )
  {
    if ( ich < 0 || ich >= samples) return 0;
    return decoded_data3[ich];
  }

  if ( strcmp(what,"CARDIAC") == 0 )
  {
    if ( ich < 0 || ich >= samples) return 0;
    return decoded_data4[ich];
  }

  if ( strcmp(what,"SERIAL") == 0 )
  {
    return serialnumber;
  }

  if ( strcmp(what,"UDPHEADER1") == 0 )
  {
    return udpheader1;
  }

  if ( strcmp(what,"UDPHEADER2") == 0 )
  {
    return udpheader2;
  }

  if ( strcmp(what,"EVENTRATE") == 0 )
  {
    return eventrate;
  }

  if ( strcmp(what,"RC_UP") == 0 )
  {
    return ratcap_upenn;
  }

  if ( strcmp(what,"CLOCKSELECT") == 0 )
  {
    return clock_sel;
  }

  return 0;

}

double oncsSub_iduppetdata_v104::dValue(const int ich)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( ich < 0 || ich >= samples) return 0;
    return d_time[ich];
}

long long  oncsSub_iduppetdata_v104::lValue(const int ich)
{

  if ( decoded_data1 == 0 ) decoded_data1 = decode(&data1_length);

  if ( ich < 0 || ich >= samples) return 0;
    return lval[ich];
}

void  oncsSub_iduppetdata_v104::dump ( OSTREAM& os )
{

  identify(os);

  int i;
  int is =  iValue(0,"SAMPLES");


  os << " Format: ";
  if ( iValue(0, "RC_UP") )
    {
      os << "Ratcap ";
    }
  else
    {
      os << "UPenn ";
    }
  os << endl;

  os << " Clock setting: ";
  if ( iValue(0, "CLOCKSELECT") )
    {
      os << "Fine ";
    }
  else
    {
      os << "Coarse ";
    }
  os << endl;

  os << " Event rate:   " << iValue(0,"EVENTRATE") << endl;
  os << " UDP Header 1, 2: " << iValue(0,"UDPHEADER1") 
     << "  " << iValue(0,"UDPHEADER2") << endl;
  os << " Serial Number :   " << iValue(0,"SERIAL") << endl;

  os << " Samples: " << iValue(0,"SAMPLES") << std::endl;
  os << " sample  B id    C id   time   ( " << is << " Samples )" << std::endl;
  
  for ( i = 0; i < is ; i++)
    {
      int prec = os.precision();

      os << std::setw(5) << i << "  |  "; 
 
      os << std::setw(4) << iValue(i,"BLOCKID") << "  ";
      os << std::setw(4) << iValue(i,"CRYSTALID") << "  ";
      os << std::setw(24) << std::setprecision(24) << hex 
	//<< tval[i]
	 << dValue(i)
	 << dec   	 << std::setprecision(prec) << "   " ;
      os << std::setw(24) << std::setprecision(24) << hex  << lValue(i) << dec
  	 << std::setprecision(prec) ;
	
      os << std::endl;
    }
  os << std::endl;
}

