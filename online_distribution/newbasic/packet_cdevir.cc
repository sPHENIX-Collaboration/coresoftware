#include <packet_cdevir.h>

Packet_cdevir::Packet_cdevir(PACKET_ptr data)
  : Packet_w4 (data)
{
  ps = 0;
}
  
int *Packet_cdevir::decode ( int *nwout)
{
  int *k;
  k = (int *)  findPacketDataStart(packet);
  if (k == 0) 
    {
      ps = 0;
      *nwout = 0;
      return 0;
    }
  int i;


  ps =  ( struct cdevIrData *) k;
  
  for (i=0; i<12; i++)fix_endianess(&ps->m_tripletTrimCurrents[i]);
  for (i=0; i<360; i++)fix_endianess(&ps->m_polarPerBunchYellowX[i]);
  fix_endianess ( ps->m_irState, 256);

  *nwout = 0;
  return 0;
}

double Packet_cdevir::dValue(const int channel,const char *what)
{

  return 0;
}



void Packet_cdevir::dump ( OSTREAM& os) 
{
  int i;
  this->identify(os);

  decode (&i);


  os << "m_irState           " <<  ps->m_irState  << std::endl;

  for (i = 0; i< 12; i++)
    os << "tripletTrimCurrents " << std::setw(3) << i <<  std::setw(12) <<  ps->m_tripletTrimCurrents[i] << std::endl;

  os << "m_irVacuum              " << ps->m_irVacuum  << std::endl;
  os << "m_estBeamSizeYellowVert " << ps->m_estBeamSizeYellowVert  << std::endl;
  os << "m_estBeamSizeYellowHorz " << ps->m_estBeamSizeYellowHorz  << std::endl;
  os << "m_estBeamSizeBlueVert   " << ps->m_estBeamSizeBlueVert  << std::endl;
  os << "m_estBeamSizeBlueHorz   " << ps->m_estBeamSizeBlueHorz  << std::endl;
  os << "m_estimatedLuminosity   " << ps->m_estimatedLuminosity  << std::endl;
  os << "m_betaStarYellowHorz    " << ps->m_betaStarYellowHorz  << std::endl;
  os << "m_betaStarBlueHorz      " << ps->m_betaStarBlueHorz  << std::endl;
  os << "m_betaStarYellowVert    " << ps->m_betaStarYellowVert  << std::endl;
  os << "m_betaStarBlueVert      " << ps->m_betaStarBlueVert  << std::endl;


  os << "m_avgOrbitDXBpmYellowHorzOdd  " << ps->m_avgOrbitDXBpmYellowHorzOdd   << std::endl;
  os << "m_avgOrbitDXBpmYellowHorzEven " << ps->m_avgOrbitDXBpmYellowHorzEven   << std::endl;
  os << "m_avgOrbitDXBpmYellowVertOdd  " << ps->m_avgOrbitDXBpmYellowVertOdd   << std::endl;
  os << "m_avgOrbitDXBpmYellowVertEven " << ps->m_avgOrbitDXBpmYellowVertEven    << std::endl;
  os << "m_avgOrbitDXBpmBlueHorzOdd    " << ps->m_avgOrbitDXBpmBlueHorzOdd   << std::endl;
  os << "m_avgOrbitDXBpmBlueHorzEven   " << ps->m_avgOrbitDXBpmBlueHorzEven   << std::endl;
  os << "m_avgOrbitDXBpmBlueVertOdd    " << ps->m_avgOrbitDXBpmBlueVertOdd   << std::endl;
  os << "m_avgOrbitDXBpmBlueVertEven   " << ps->m_avgOrbitDXBpmBlueVertEven   << std::endl;
  os << "m_vertexStartTime             " << ps->m_vertexStartTime   << std::endl;
  os << "m_vertexEndTime               " << ps->m_vertexEndTime   << std::endl;
  os << "m_datavalidMask               " << ps->m_datavalidMask   << std::endl; 

  for (i = 0; i< 100; i++)
    {
      
      os << std::setw(3)  << i << " | "
	 << std::setw(12) << ps->m_experimentVertexX[i]
	 << std::setw(12) << ps->m_experimentVertexY[i]
	 << std::setw(12) << ps->m_experimentVertexZ[i]
	 << std::endl;

    }

  for (i = 0; i< 360; i++)
    {
      
      os << std::setw(3) << i <<  " | "  
	 <<  std::setw(12) << ps->m_polarPerBunchYellowX[i]
	 <<  std::setw(12) << ps->m_polarPerBunchYellowY[i]
	 <<  std::setw(12) << ps->m_polarPerBunchYellowZ[i]
	 <<  std::setw(12) << ps->m_polarPerBunchBlueX[i]
	 <<  std::setw(12) << ps->m_polarPerBunchBlueY[i]
	 <<  std::setw(12) << ps->m_polarPerBunchBlueZ[i]
	 << std::endl;
    }


}
