#include <packet_cdevmadch.h>
#include <time.h>
#include <stdio.h>
Packet_cdevmadch::Packet_cdevmadch(PACKET_ptr data)
  : Packet_w4 (data)
{
  ps = 0;
}
  
int *Packet_cdevmadch::decode ( int *nwout)
{

  if (ps != 0) return 0; 
 int il = getDataLength();
  no_structures = 4* il / sizeof ( struct cdevMadchData );
  std::cout << "no_structures = " << no_structures << std::endl;

  int *k = (int *) findPacketDataStart(packet);
  if (k == 0) 
    {
      ps = 0;
      *nwout = 0;
      return 0;
    }
  ps = ( struct cdevMadchData *) k;


  //fix_endianess ( &ps->current);
  //fix_endianess ( ps->avgOrbPosition);
  //fix_endianess ( ps->avgOrbVariance);
  //fix_endianess ( ps->avgOrbStat);


  *nwout = 0;

  return 0;
}


// ------------------------------------------------------


void Packet_cdevmadch::dump ( OSTREAM &os) 
{

  int i;

  decode (&i);

  this->identify(os);
  os << "Number of readings: " << iValue(0,"NOREADINGS") << std::endl;

  os << "current " << std::endl;
  for ( i = 0; i < iValue(0,"NOREADINGS") ; i++)
    {
      /*
      os << std::setw(4 ) << i ;
      os << std::setw(16) << ps[i].current;
      os << std::endl;
      */
      
      printf("%d  hex current %x  current %lf\n",i,(unsigned)ps[i].current,ps[i].current);
      //printf("%d    current %lf\n",i,ps[i].current);
      
      /*
      os << std::setw(16) << dValue(i,"current");
      os << std::setw(16) << " int value " << iValue(i,"current");
      os << std::setw(16) << " time  " << iValue(i,"timestamp");
      os << std::endl;
      */

    }

  dumpErrorBlock(os);
  dumpDebugBlock(os);
}


int   Packet_cdevmadch::iValue(const int ich, const char *what)
{

  //  std::cout << "IN  Packet_cdevmadch::rValue " << std::endl;
  int i;
  decode (&i);
// Unix time
  if ( strcmp(what, "NOREADINGS") == 0 ) return   no_structures ;	  
  if ( strcmp(what, "timestamp") == 0 ) return ps[ich].cdevCaptureTimeStamp;

  std::cout << "packet_cdevmadc::iValue error unknown datum: " << what << std::endl;
  return 0;
}


double   Packet_cdevmadch::dValue(const int ich, const char *what)
{


  //  std::cout << "IN  Packet_cdevmadch::rValue " << std::endl;
  int i;
  decode (&i);

  if ( ich < 0 || ich >= no_structures ) return 0;

  if ( strcmp(what, "current") == 0 ) return   ps[ich].current ;


  std::cout << "packet_cdevmadc::dValue error unknown datum: " << what << std::endl;
  return 0;
  

  
}


