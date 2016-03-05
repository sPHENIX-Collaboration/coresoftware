#include <packet_cdevring.h>

Packet_cdevring::Packet_cdevring(PACKET_ptr data)
  : Packet_w4 (data)
{
  ps = 0;
  haspoldata = 0;
  hasfilldata = 0;
  decoded = 0;  // only decode once...
}
  
int *Packet_cdevring::decode ( int *nwout)
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

  if ( getHitFormat() == IDCDEVRING) 
    {
      haspoldata = 0;
       ps =  ( struct cdevRingData *) k;
    }
  else if (getHitFormat() == IDCDEVRINGPOL )    
    {
      haspoldata = 1;
      ps =  ( struct cdevRingData *) k;
    }
  else if (getHitFormat() == IDCDEVRINGFILL )
    {
      hasfilldata = 1;
      ps = ( struct cdevRingData *) k;
    }



  *nwout = 0;
  return 0;
}

int   Packet_cdevring::iValue(const int ich, const char *what)
{
  int i;
  decode (&i);

  if ( strcmp(what, "ringState") == 0 ) {
    if (ich >=0 && ich < 256) {
      return ps->m_ringState[ich] ;
    }
  }

  if ( strcmp(what, "ionSpecies") == 0 ) {
    if (ich>=0 && ich<1024) {
      return ps->m_ionSpecies[ich] ;
    }
  }

  if ( strcmp(what, "stoneType") == 0 ) return ps->m_stoneType ;

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


  if ( strcmp(what, "timeOfFillStart") == 0 ) return ps->m_timeOfFillStart;
  if ( strcmp(what, "timeOfLuminosityStart") == 0 ) return ps->m_timeOfLuminosityStart;

  if (hasfilldata )
    {
      if ( strcmp(what, "fillNumber") == 0 ) return ps->m_fillNumber;
      if ( strcmp(what, "datavalidMask") == 0 ) return ps->m_datavalidMask;
      haspoldata = 1;  // fill data is a superset of poldata
    }

  if (haspoldata )
    {
      if ( strcmp(what, "measuredPolarizationUp") == 0 ) {
	if (ich >= 0 && ich < 360 ) {
	  return ps->m_measuredPolarizationUp[ich];
	}
      }
      if ( strcmp(what, "measuredPolarizationDown") == 0 ) {
	if (ich >= 0 && ich < 360 ) {
	  return ps->m_measuredPolarizationDown[ich];
	}
      }


    } // haspoldata test


  std::cout << "packet_cdevring::iValue error unknown datum: " << what << std::endl;
  return 0;

}


double Packet_cdevring::dValue(const int channel,const char *what)
{
  int i;
  decode (&i);

  if ( strcmp(what, "beamEnergy") == 0 ) return ps->m_beamEnergy;
  if ( strcmp(what, "gamma") == 0 ) return ps->m_gamma;
  if ( strcmp(what, "momentumSpread") == 0 ) return ps->m_momentumSpread;
  if ( strcmp(what, "bunchLength") == 0 ) return ps->m_bunchLength;
  if ( strcmp(what, "bunchOneRelativePhase") == 0 ) return ps->m_bunchOneRelativePhase;
  if ( strcmp(what, "synchrotronTune") == 0 ) return ps->m_synchrotronTune;
  if ( strcmp(what, "chromaticityVertical") == 0 ) return ps->m_chromaticityVertical;
  if ( strcmp(what, "chromaticityHorizontal") == 0 ) return ps->m_chromaticityHorizontal;
  if ( strcmp(what, "emittanceVertical") == 0 ) return ps->m_emittanceVertical;
  if ( strcmp(what, "emittanceHorizontal") == 0 ) return ps->m_emittanceHorizontal;
  if ( strcmp(what, "betaIPMHorizontal") == 0 ) return ps->m_betaIPMHorizontal;
  if ( strcmp(what, "betaIPMVertical") == 0 ) return ps->m_betaIPMVertical;

  std::cout << "packet_cdevring::dValue error unknown datum: " << what << std::endl;
  return 0;
}



void Packet_cdevring::dump ( OSTREAM& os) 
{
  int i;
  decode (&i);

  this->identify(os);

  os << "m_ringState               " <<  ps->m_ringState  << std::endl;
  os << "m_ionSpecies              " <<  ps->m_ionSpecies  << std::endl;
  os << "m_beamEnergy              " <<  ps->m_beamEnergy  << std::endl;
  os << "m_gamma                   " <<  ps->m_gamma  << std::endl;
  os << "m_stoneType               " <<  ps->m_stoneType  << std::endl;
  os << "m_momentumSpread          " <<  ps->m_momentumSpread  << std::endl;
  os << "m_bunchLength             " <<  ps->m_bunchLength  << std::endl;
  os << "m_chromaticityVertical    " <<  ps->m_chromaticityVertical  << std::endl;
  os << "m_chromaticityHorizontal  " <<  ps->m_chromaticityHorizontal  << std::endl;
  os << "m_emittanceVertical       " <<  ps->m_emittanceVertical  << std::endl;


  if (hasfilldata)
    {
      os << "fillNumber              " << iValue(i,"fillNumber") << std::endl;
      os << "datavalidMask           " << iValue(i,"datavalidMask") << std::endl;
      haspoldata = 1;
    }

  if (haspoldata )
    {
      os << "index " <<" intendedFillPattern " << " measuredFillPattern "<< " polarizationFillPattern "<< " measuredPolarizationUp" << " measuredPolarizationDown"<< std::endl;
	for (i =  0 ; i < 360 ;i++) 
	  {
	    os  << std::setw(3) <<i << std::setw(3)<< "|" <<  std::setw(10) <<  iValue(i,"intendedFillPattern")<< std::setw( 20)<< iValue(i,"measuredFillPattern")<< std::setw( 23)<<iValue(i,"polarizationFillPattern")<< std::setw(25) <<iValue(i,"measuredPolarizationUp")<< std::setw(23) << iValue(i,"measuredPolarizationDown");
      os << std::endl;
	}
    os << std::endl;
	

    } else // haspoldata test
      {

  os << "index " <<" intendedFillPattern " << " measuredFillPattern "<< " polarizationFillPattern "<< std::endl;
  for (i = 0; i< 360; i++) {
    os  << std::setw(3) <<i << std::setw(3)<< "|" <<  std::setw(10) <<  iValue(i,"intendedFillPattern")<< std::setw( 20)<< iValue(i,"measuredFillPattern")<< std::setw( 15)<<iValue(i,"polarizationFillPattern"); //ps->m_intendedFillPattern[i] ;
    os << std::endl;
  }
  os << std::endl;

      } // no polarity data

}

