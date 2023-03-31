#ifndef TPC_RAWHIT_H
#define TPC_RAWHIT_H

#include <trackbase/RawHit.h>

#include <trackbase/TrkrDefs.h>

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <iomanip>
#include <iostream>

/**
 * @brief TPC raw hit class
 *
 * mlp - first attempt to get something going
 */
class TPC_RawHit : public RawHit
{
 public:
  TPC_RawHit(const unsigned int phi, const unsigned int samplenr, const unsigned int adc)
  {
    m_phi = phi;
    m_samplenr = samplenr;
    m_adc = adc;
  }

  ~TPC_RawHit() override {}
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "phi: " << std::setw(4) << m_phi
       << " sample: " << std::setw(5) << m_samplenr
       << " ADC: " << std::setw(5) << m_adc
       << std::endl;
  }
  void Reset() override { m_adc = m_phi = m_samplenr = 0; }
  int isValid() const override { return 1; }

  // after digitization, these are the adc values
  virtual void setAdc(const unsigned int v) override { m_adc = v; }
  virtual unsigned int getAdc() override { return m_adc; }

  virtual void setPhiBin(const unsigned int v) override { m_phi = v; }
  virtual unsigned int getPhiBin() override { return m_phi; }

  virtual void setTBin(const unsigned int bv) override { m_samplenr = bv; }
  virtual unsigned int getTBin() override { return m_samplenr; }

 protected:
  unsigned int m_adc;
  unsigned int m_phi;
  unsigned int m_samplenr;

  ClassDefOverride(TPC_RawHit, 1);
};

#endif  // TRACKBASE_RAWHIT_H
