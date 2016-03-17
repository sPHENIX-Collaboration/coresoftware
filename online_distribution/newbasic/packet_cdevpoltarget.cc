#include <packet_cdevpoltarget.h>
#include <stdio.h>

Packet_cdevpoltarget::Packet_cdevpoltarget(PACKET_ptr data)
  : Packet_w4 (data)
{
  ps = 0;
}
  
int *Packet_cdevpoltarget::decode ( int *nwout)
{

  if (ps != 0) return 0; 

  int *k = (int *) findPacketDataStart(packet);
  if (k == 0) 
    {
      ps = 0;
      *nwout = 0;
      return 0;
    }
  ps = ( struct cdevPolTargetData *) k;
	

  *nwout = 0;

  return 0;
}


// ------------------------------------------------------


void Packet_cdevpoltarget::dump ( OSTREAM &os) 
{

  int i;

  decode (&i);

  this->identify(os);

  os << "positionEncLinear   " <<  iValue(i,"positionEncLinear") << std::endl;
  os << "positionEncRot      " <<  iValue(i,"positionEncRot")  << std::endl;


  dumpErrorBlock(os);
  dumpDebugBlock(os);
}


int   Packet_cdevpoltarget::iValue(const int ich, const char *what)
{

  //  std::cout << "IN  Packet_cdevpoltarget::iValue " << std::endl;
  int i;
  decode (&i);

  
  if ( strcmp(what,"positionEncLinear") == 0)
    {
	return  ps->positionEncLinear;
    }
  
  if ( strcmp(what,"positionEncRot") == 0)
    {

	return  ps->positionEncRot;
    }


  std::cout << "packet_cdevpoltarget::iValue error unknown datum: " << what << std::endl;
  return 0;
}

