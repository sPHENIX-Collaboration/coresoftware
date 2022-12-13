#include "RawHitTpc.h"
#include <climits>

RawHitTpc::RawHitTpc()
{
}

void RawHitTpc::setAdc(const unsigned int adc)
 {
   if(adc > USHRT_MAX)
     m_adc = USHRT_MAX;
   else
     m_adc = (unsigned short) adc;
 }

unsigned int RawHitTpc::getAdc() { 
    return (unsigned int) m_adc;
  }

void RawHitTpc::setPhiBin(const unsigned int phibin) {
  std::cout << "TPC RAW HIT HAS NO PHIBIN " << phibin << std::endl;
}
unsigned int RawHitTpc::getPhiBin() { std::cout << "TPC RAW HIT HAS NO PHIBIN" << std::endl; return 0;}

void RawHitTpc::setTBin(const unsigned int tbin) {
  m_tbin = tbin;
}
unsigned int RawHitTpc::getTBin() { return m_tbin;}
