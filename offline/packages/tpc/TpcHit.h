/**
 * @file tpc/TpcHit.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Tpc hit object
 */
#ifndef TPC_TPCHIT_H
#define TPC_TPCHIT_H

#include <trackbase/TrkrHit.h>

#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif

/**
 * @brief Tpc hit object
 *
 * Container for Tpc hit object representing
 * the adc for a single phi/z bin
 */
class TpcHit : public TrkrHit
{
public:
  //! ctor
  TpcHit();
  //! dtor
  virtual ~TpcHit() {};

  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();
  virtual int isValid() const;
  
  /**
   * @brief Set the ADC information
   * @param[in] adc ADC value
   */
  void setAdc(const short adc) { m_adc = adc; }

  /**
   * @brief Get ADC value
   * @param[out] ADC value
   */
  short getAdc() const { return m_adc; }

private:
  short m_adc;
  ClassDef(TpcHit,1);
};

#endif //TPC_TPCHIT_H
