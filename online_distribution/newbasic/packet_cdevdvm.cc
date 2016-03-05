#include <packet_cdevdvm.h>

Packet_cdevdvm::Packet_cdevdvm(PACKET_ptr data)
  : Packet_w4 (data)
{
  ps = 0;
  decoded = 0;
}
  
int *Packet_cdevdvm::decode ( int *nwout)
{
  if (decoded) {
    *nwout = 0;
    return 0;
  }

  decoded=1;  // decode only once...

  int *k;
  k = (int *)  findPacketDataStart(packet);
  if (k == 0) 
    {
      ps = 0;
      *nwout = 0;
      return 0;
    }

  ps =  ( struct cdevDvmData *) k;
  
 
 
  *nwout = 0;
  return 0;
}

float   Packet_cdevdvm::fValue(const int ich, const char *what)
{
  int i;
  decode (&i);

  if ( strcmp(what, "beamLifeTime") == 0 ) return ps->beamLifeTime;
  
  std::cout << "packet_cdevdvm::fValue error unknown datum: " << what << std::endl;
  return 0;

}


float   Packet_cdevdvm::rValue(const int ich, const char *what)
{
  int i;
  decode (&i);

  if ( strcmp(what, "beamLifeTime") == 0 ) return ps->beamLifeTime;
  
  std::cout << "packet_cdevdvm::rValue error unknown datum: " << what << std::endl;
  return 0;

}


double Packet_cdevdvm::dValue(const int channel,const char *what)
{
  int i;
  decode (&i);

  if ( strcmp(what, "beamCurrent") == 0 ) return ps->beamCurrent;
  if ( strcmp(what, "beamLifeTime") == 0 ) return ps->beamLifeTime;
  

  std::cout << "packet_cdevdvm::dValue error unknown datum: " << what << std::endl;
  return 0;
}



void Packet_cdevdvm::dump ( OSTREAM& os) 
{
  int i;
  decode (&i);

  this->identify(os);

  os << std::setw(12)<< "beamCurrent   " <<  ps->beamCurrent  << std::endl;
  os << std::setw(12)<< "beamLifeTime  " <<  ps->beamLifeTime  << std::endl;

}





