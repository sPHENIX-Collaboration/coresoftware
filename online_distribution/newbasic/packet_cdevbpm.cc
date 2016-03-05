#include <packet_cdevbpm.h>
#include <time.h>

Packet_cdevbpm::Packet_cdevbpm(PACKET_ptr data)
  : Packet_w4 (data)
{
  ps = 0;
}
  
int *Packet_cdevbpm::decode ( int *nwout)
{

  if (ps != 0) return 0; 

  int il = getDataLength();

  no_structures = 4* il / sizeof ( struct cdevBPMData );
  std::cout << "no_structures = " << no_structures << std::endl;
  int *k = (int *) findPacketDataStart(packet);
  if (k == 0) 
    {
      ps = 0;
      *nwout = 0;
      return 0;
    }

  ps = ( struct cdevBPMData *) k;



  // no byte swat for floats
  //fix_endianess ( ps->avgOrbTimeStamp);
  //fix_endianess ( ps->avgOrbPosition);
  //fix_endianess ( ps->avgOrbVariance);
  //fix_endianess ( ps->avgOrbStat);


  *nwout = 0;

  return 0;
}


// ------------------------------------------------------


void Packet_cdevbpm::dump ( OSTREAM &os) 
{

  int i;

  decode (&i);

  this->identify(os);

  os << "Number of readings: " << iValue(0,"NOREADINGS") << std::endl;
  
  os << "index      ";
  os << "avgOrbTimeStamp " ;
  os << "avgOrbPosition  " ;
  os << "argOrbVariance  " ;
  os << "argOrbStat      " << std::endl;

  for ( i = 0; i < iValue(0,"NOREADINGS") ; i++)
    {
      os << std::setw(4 ) << i ;
      os << std::setw(16) << iValue(i,"avgOrbTimeStamp") ;
      os << std::setw(16) << rValue(i,"avgOrbPosition") ;
      os << std::setw(16) << rValue(i,"avgOrbVariance") ;
      os << std::setw(16) << rValue(i,"avgOrbStat");
      os << std::endl;
    }



  dumpErrorBlock(os);
  dumpDebugBlock(os);
}
//-------------------------------------------------------------

int   Packet_cdevbpm::iValue(const int ich, const char *what)
{

  //  std::cout << "IN  Packet_cdevbpm::rValue " << std::endl;
  int i;
  decode (&i);

  if ( ich < 0 || ich >= no_structures ) return 0;

// Unix time
  if ( strcmp(what, "NOREADINGS") == 0 ) return   no_structures ;
  if ( strcmp(what, "avgOrbTimeStamp") == 0 ) return   ps[ich].avgOrbTimeStamp ;
  if ( strcmp(what, "datavalidMask") == 0 ) return   ps[ich].datavalidMask ;

  std::cout << "packet_cdevbpm::iValue error unknown datum: " << what << std::endl;
  return 0;
}


float   Packet_cdevbpm::rValue(const int ich, const char *what)
{


  //  std::cout << "IN  Packet_cdevbpm::rValue " << std::endl;
  int i;
  decode (&i);

  if ( ich < 0 || ich >= no_structures ) return 0;


  if ( strcmp(what, "avgOrbPosition")  == 0 ) return   ps[ich].avgOrbPosition ;    
  if ( strcmp(what, "avgOrbVariance")  == 0 ) return   ps[ich].avgOrbVariance ;	
  if ( strcmp(what, "avgOrbStat") == 0 ) return   ps[ich].avgOrbStat;

  std::cout << "packet_cdevbpm::rValue error unknown datum: " << what << std::endl;
  return 0;
  

  
}


