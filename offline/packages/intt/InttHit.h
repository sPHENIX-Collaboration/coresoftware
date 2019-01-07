/**
 * @file intt/InttHit.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Intt hit object
 */
#ifndef INTT_INTTHIT_H
#define INTT_INTTHIT_H

#include <trackbase/TrkrHit.h>

#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif

/**
 * @brief Intt hit object
 *
 * Container for Intt hit object representing
 * a single hit strip within a chip
 */
class InttHit : public TrkrHit
{
 public:
  //! ctor
  InttHit();
  //! dtor
  virtual ~InttHit(){};

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
  ClassDef(InttHit, 1);
};

#endif  //INTT_INTTHIT_H
