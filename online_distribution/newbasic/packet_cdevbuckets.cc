#include <packet_cdevbuckets.h>

Packet_cdevbuckets::Packet_cdevbuckets(PACKET_ptr data)
  : Packet_w4 (data)
{
  ps = 0;
  decoded = 0;  // only decode once...
}
  
int *Packet_cdevbuckets::decode ( int *nwout)
{
  if (decoded) {
    *nwout = 0;
    return 0;
  }

  decoded = 1;

  int *k;
  k = (int *)  findPacketDataStart(packet);
  if (k == 0) 
    {
      ps = 0;
      *nwout = 0;
      return 0;
    }

  if ( getHitFormat() == IDCDEVBUCKETS) 
    {
       ps =  ( struct cdevBucketsData *) k;
    }



  *nwout = 0;
  return 0;
}

int   Packet_cdevbuckets::iValue(const int ich, const char *what)
{
  int i;
  decode (&i);


  if ( strcmp(what, "intendedFillPattern") == 0 ) {
    if (ich >= 0 && ich < 360 ) {
      return ps->m_intendedFillPattern[ich];
    }
  }

  if ( strcmp(what, "measuredFillPattern") == 0 ) {
    if (ich >= 0 && ich < 360 ) {
      return ps->m_measuredFillPattern[ich];
    }
  }

  if ( strcmp(what, "polarizationFillPattern") == 0 ) {
    if (ich >= 0 && ich < 360 ) {
      return ps->m_polarizationFillPattern[ich];
    }
  }








  std::cout << "packet_cdevbuckets::iValue error unknown datum: " << what << std::endl;
  return 0;

}


double Packet_cdevbuckets::dValue(const int channel,const char *what)
{
  int i;
  decode (&i);

  if ( strcmp(what, "bunchLength") == 0 ) return ps->m_bunchLength;
  if ( strcmp(what, "fillPatternThreshold") == 0 ) return ps->m_fillPatternThreshold;
  if ( strcmp(what, "bunchOneRelativePhase") == 0 ) return ps->m_bunchOneRelativePhase;

  std::cout << "packet_cdevbuckets::dValue error unknown datum: " << what << std::endl;
  return 0;
}



void Packet_cdevbuckets::dump ( OSTREAM& os) 
{
  int i;
  decode (&i);

  this->identify(os);

  os << "m_bunchLength             " <<  ps->m_bunchLength  << std::endl;
  os << "m_fillPatternThreshol     " <<  ps->m_fillPatternThreshold  << std::endl;
  os << "m_bunchOneRelativePhase   " <<  ps->m_bunchOneRelativePhase  << std::endl;

      os << "index " <<" intendedFillPattern " << " measuredFillPattern "<< " polarizationFillPattern " << std::endl;
	for (i =  0 ; i < 360 ;i++) 
	  {
	    os  << std::setw(3) <<i << std::setw(3)<< "|" <<  std::setw(10) <<  iValue(i,"intendedFillPattern")<< std::setw( 20)<< iValue(i,"measuredFillPattern")<< std::setw( 23)<<iValue(i,"polarizationFillPattern");
      os << std::endl;
	}
    os << std::endl;

}

