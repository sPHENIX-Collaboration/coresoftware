#include "RawHitv1.h"
#include <climits>

RawHitv1::RawHitv1()
{
}

void RawHitv1::setAdc(const unsigned int adc)
 {
   if(adc > USHRT_MAX)
     m_adc = USHRT_MAX;
   else
     m_adc = (unsigned short) adc;
 }

unsigned int RawHitv1::getAdc() { 
    return (unsigned int) m_adc;
  }

void RawHitv1::setPhiBin(const unsigned int phibin) {
  m_phibin = phibin;
}
unsigned int RawHitv1::getPhiBin() { return m_phibin;}

void RawHitv1::setTBin(const unsigned int tbin) {
  m_tbin = tbin;
}
unsigned int RawHitv1::getTBin() { return m_tbin;}
