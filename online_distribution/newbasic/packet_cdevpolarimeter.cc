#include <packet_cdevpolarimeter.h>
#include <stdio.h>

Packet_cdevpolarimeter::Packet_cdevpolarimeter(PACKET_ptr data)
  : Packet_w4 (data)
{
  ps = 0;
  haspoldata = 0;
}
  
int *Packet_cdevpolarimeter::decode ( int *nwout)
{

  if (ps != 0) return 0; 

  int *k = (int *) findPacketDataStart(packet);
  if (k == 0) 
    {
      ps = 0;
      *nwout = 0;
      return 0;
    }

  if ( getHitFormat() == IDCDEVPOLARIMETER ) 
    {
      haspoldata = 0;
      ps = ( struct cdevPolarimeterData *) &k[1];
    }
  else if (getHitFormat() == IDCDEVPOLARIMETERZ )
    {
      haspoldata = 1;
      ps =  ( struct cdevPolarimeterData *) &k[1];
    }

  //  fix_endianess (&ps->m_cdevCaptureTimeStamp);
  // ejd 4/26/03 strings are garbled do not byte swap
  //  fix_endianess (&ps->runIdS);
  //fix_endianess ( ps->daqVersionS, 80);
  //fix_endianess ( ps->cutIdS, 80);
  //fix_endianess ( ps->targetIdS, 80);
  //fix_endianess ( ps->statusStringS, 80);


  //  fix_endianess(&ps->startTimeS);	// Unix time
  //  fix_endianess(&ps->stopTimeS);	// Unix time
    
  //  for (i=0; i<2; i++) fix_endianess(&ps->encoderPositionS[i]);
  //  fix_endianess(&ps->statusS);	// bit pattern if <0 data is not usable
  //  fix_endianess(&ps->totalCountsS);
  //  fix_endianess(&ps->upCountsS);
  //  fix_endianess(&ps->downCountsS);
  //  fix_endianess(&ps->unpolCountsS);
		
//    for (i=0; i<360; i++) fix_endianess(&ps->countsUpLeftS[i]);
//    for (i=0; i<360; i++) fix_endianess(&ps->countsLeftS[i]);
//    for (i=0; i<360; i++)fix_endianess(&ps->countsDownLeftS[i]);
//    for (i=0; i<360; i++)fix_endianess(&ps->countsDownRightS[i]);
//    for (i=0; i<360; i++)fix_endianess(&ps->countsRightS[i]);
//    for (i=0; i<360; i++)fix_endianess(&ps->countsUpRightS[i]);



//    fix_endianess(&ps->numberEventsS);	// provided by MCR before measurement
//    fix_endianess(&ps->maxTimeS);	

  *nwout = 0;

  return 0;
}


// ------------------------------------------------------


void Packet_cdevpolarimeter::dump ( OSTREAM &os) 
{

  int i;

  decode (&i);

  this->identify(os);

  os << "daqVersionS   " <<  ps->daqVersionS  << std::endl;
  os << "cutIdS        " <<  ps->cutIdS  << std::endl;
  os << "targetIdS     " <<  ps->targetIdS << std::endl;
  os << "statusStringS " <<  ps->statusStringS << std::endl;

  os << "m_cdevCaptureTimeStamp   " <<  ps->m_cdevCaptureTimeStamp << std::endl;
  os << "runIdS                   " <<  ps->runIdS << std::endl;	// FILL.XXX --- where XXX is the run number 
  os << "startTimeS               " <<  ps->startTimeS << std::endl;	// Unix time
  os << "stopTimeS                " <<  ps->stopTimeS << std::endl;	// Unix time
  os << "encoderPositionS         " <<  ps->encoderPositionS[0] <<"     " <<  ps->encoderPositionS[1] << std::endl;

  //  os <<  ps->daqVersionS << std::endl;
  // os <<  ps->cutIdS << std::endl;
  //os <<  ps->targetIdS << std::endl;

  os << "statusS                  " <<  ps->statusS << std::endl;	// bit pattern if <0 data is not usable
  os << "totalCountsS             " <<  ps->totalCountsS << std::endl;
  os << "upCountsS                " <<  ps->upCountsS << std::endl;
  os << "downCountsS              " <<  ps->downCountsS << std::endl;
  os << "unpolCountsS             " <<  ps->unpolCountsS << std::endl;
  os << "avgAsymXS                " <<  ps->avgAsymXS << std::endl;
  os << "avgAsymX45S              " <<  ps->avgAsymX45S << std::endl;
  os << "avgAsymX90S              " <<  ps->avgAsymX90S << std::endl;    
  os << "avgAsymYS                " <<  ps->avgAsymYS << std::endl;
  os << "avgAsymErrorXS           " <<  ps->avgAsymErrorXS << std::endl;
  os << "avgAsymErrorX45S         " <<  ps->avgAsymErrorX45S << std::endl;
  os << "avgAsymErrorX90S         " <<  ps->avgAsymErrorX90S << std::endl;
  os << "avgAsymErrorYS           " <<  ps->avgAsymErrorYS << std::endl;

  os << "beamEnergyS              " <<  ps->beamEnergyS << std::endl;	// the same as ringSpec.color:beamEnergyM just for reference

  
  os << "analyzingPowerS          " <<  ps->analyzingPowerS << std::endl;
  os << "analyzingPowerErrorS     " <<  ps->analyzingPowerErrorS << std::endl;
  os << "numberEventsS            " <<  ps->numberEventsS << std::endl;	// provided by MCR before measurement
  os << "maxTimeS                 " <<  ps->maxTimeS << std::endl;	
  if (haspoldata )
    {
      os << "polarizationM        " << ps->polarizationM << std::endl;
    }




  os << " index        AsymXS      AsymYS AsymErrorXS  AsymErrorYS    UpLeftS       LeftS   DownLeftS  DownRightS      RightS    UpRightS" << std::endl;


  for (i=0; i< 360; i++)
  {
    
    os << std::setw(5) << i << " | "; 

    
    os << std::setw(12) << rValue(i,"bunchAsymXS"); //ps->bunchAsymXS[i];
    os << std::setw(12) << rValue(i,"bunchAsymYS"); //ps->bunchAsymYS[i];
    os << std::setw(12) << ps->bunchAsymErrorXS[i];
    os << std::setw(12) << ps->bunchAsymErrorYS[i];

    os << std::setw(12) << ps->countsUpLeftS[i];
    os << std::setw(12) << ps->countsLeftS[i];
    os << std::setw(12) << ps->countsDownLeftS[i];
    os << std::setw(12) << ps->countsDownRightS[i];
    os << std::setw(12) << ps->countsRightS[i];
    os << std::setw(12) << ps->countsUpRightS[i] << std::endl;
  }


  dumpErrorBlock(os);
  dumpDebugBlock(os);
}


int   Packet_cdevpolarimeter::iValue(const int ich, const char *what)
{

  int i;
  decode (&i);

  if ( strcmp(what, "cdevCaptureTimeStamp") == 0 ) return   ps->m_cdevCaptureTimeStamp ;  // is this not an int??
  if ( strcmp(what, "startTimeS") == 0 ) return   ps->startTimeS ;	// Unix time
  if ( strcmp(what, "stopTimeS") == 0 ) return   ps->stopTimeS ;	// Unix time

  if ( strcmp(what,"daqVersionS") == 0)
    {
      if ( ich >= 0 || ich < 80)
	return  ps->daqVersionS[ich];
    }

  if ( strcmp(what,"cutIdS") == 0)
    {
      if ( ich >= 0 || ich < 80)
	return  ps->cutIdS[ich];
    }

  if ( strcmp(what,"targetIdS") == 0)
    {
      if ( ich >= 0 || ich < 80)
	return  ps->targetIdS[ich];
    }

  if ( strcmp(what, "encoderPositionS") == 0 ) 
    {
      if ( ich >= 0 || ich < 2 ) 
	return ps->encoderPositionS[ich];
    }

  if ( strcmp(what, "statusS") == 0 ) return   ps->statusS ;

  if ( strcmp(what, "statusStringS") == 0 ) 
    {
      if ( ich >= 0 || ich < 2 ) 
	return ps->statusStringS[ich];
    }

  if ( strcmp(what, "totalCountsS") == 0 ) return   ps->totalCountsS ;
  if ( strcmp(what, "upCountsS") == 0 ) return   ps->upCountsS ;
  if ( strcmp(what, "downCountsS") == 0 ) return   ps->downCountsS ;
  if ( strcmp(what, "unpolCountsS") == 0 ) return   ps->unpolCountsS ;
  if ( strcmp(what,"countsUpLeftS") == 0)
    {
      if ( ich >= 0 || ich < 360)
	return  ps->countsUpLeftS[ich];
    }


  if ( strcmp(what,"countsLeftS") == 0)
    {
      if ( ich >= 0 || ich < 360)
	return ps->countsLeftS[ich];
    }
  
  if ( strcmp(what,"countsDownLeftS") == 0)
    {
      if ( ich >= 0 || ich < 360)
	return  ps->countsDownLeftS[ich];
    }
  
  if ( strcmp(what,"countsDownRightS") == 0)
    {
      if ( ich >= 0 || ich < 360)
	return  ps->countsDownRightS[ich];
    }
  
  if ( strcmp(what,"countsRightS") == 0)
    {
      if ( ich >= 0 || ich < 360)
	return  ps->countsRightS[ich];
    }
  
  if ( strcmp(what,"countsUpRightS") == 0)
    {
      if ( ich >= 0 || ich < 360)
	return  ps->countsUpRightS[ich];
    }

  if ( strcmp(what, "numberEventsS") == 0 ) return   ps->numberEventsS ;
  if ( strcmp(what, "maxTimeS")      == 0 ) return   ps->maxTimeS ;

  std::cout << "packet_cdevpolarimeter::iValue error unknown datum: " << what << std::endl;
  return 0;
}


float   Packet_cdevpolarimeter::rValue(const int ich, const char *what)
{


  //  std::cout << "IN  Packet_cdevpolarimeter::rValue " << std::endl;
  int i;
  decode (&i);

  if ( strcmp(what, "avgAsymXS") == 0 ) return   ps->avgAsymXS ;
  if ( strcmp(what, "avgAsymX45S") == 0 ) return   ps->avgAsymX45S ;
  if ( strcmp(what, "avgAsymX90S") == 0 ) return   ps->avgAsymX90S ;    
  if ( strcmp(what, "avgAsymYS") == 0 ) return   ps->avgAsymYS ;
  if ( strcmp(what, "avgAsymErrorXS") == 0 ) return   ps->avgAsymErrorXS ;
  if ( strcmp(what, "avgAsymErrorX45S") == 0 ) return   ps->avgAsymErrorX45S ;
  if ( strcmp(what, "avgAsymErrorX90S") == 0 ) return   ps->avgAsymErrorX90S ;
  if ( strcmp(what, "avgAsymErrorYS") == 0 ) return   ps->avgAsymErrorYS ;


  if ( strcmp(what,"bunchAsymXS") == 0)
    {
      if ( ich >= 0 || ich < 360)
	return ps->bunchAsymXS[ich];
    }
  
  if ( strcmp(what,"bunchAsymYS") == 0)
    {
      if ( ich >= 0 || ich < 360)
	return  ps->bunchAsymYS[ich];
    }
  
  if ( strcmp(what,"bunchAsymErrorXS") == 0)
    {
      if ( ich >= 0 || ich < 360)
	return  ps->bunchAsymErrorXS[ich];
    }
  
  if ( strcmp(what,"bunchAsymErrorYS") == 0)
    {
      if ( ich >= 0 || ich < 360)
	return  ps->bunchAsymErrorYS[ich];
    }
  
  
  if ( strcmp(what, "beamEnergyS") == 0 ) return   ps->beamEnergyS ;
  if ( strcmp(what, "analyzingPowerS") == 0 ) return   ps->analyzingPowerS ;
  if ( strcmp(what, "analyzingPowerErrorS") == 0 ) return   ps->analyzingPowerErrorS ;
  if (haspoldata )
    {
       if ( strcmp(what, "polarizationM") == 0 ) return ps->polarizationM ;
    }

  std::cout << "packet_cdevpolarimeter::iValue error unknown datum: " << what << std::endl;
  return 0.;  
}


double   Packet_cdevpolarimeter::dValue(const int ich, const char *what)
{


  //  std::cout << "IN  Packet_cdevpolarimeter::rValue " << std::endl;
  int i;
  decode (&i);

  if ( strcmp(what, "runIdS") == 0 ) return   ps->runIdS ;	// FILL.XXX --- where XXX is the run number 

  std::cout << "packet_cdevpolarimeter::iValue error unknown datum: " << what << std::endl;
  return 0;
}



