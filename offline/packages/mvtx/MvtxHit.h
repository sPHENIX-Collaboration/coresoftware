/**
 * @file mvtx/MvtxHit.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Mvtx hit object
 */
#ifndef MVTX_MVTXHIT_H
#define MVTX_MVTXHIT_H

#include <trackbase/TrkrHit.h>

#include <iostream>

/**
 * @brief Mvtx hit object
 *
 * Container for Mvtx hit object representing
 * a single hit pixel within a chip
 */
class MvtxHit : public TrkrHit
{
 public:
  //! ctor
  MvtxHit();
  //! dtor
  virtual ~MvtxHit(){};

  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();
  virtual int isValid() const;

  //void addEnergy(const double edep) {m_edep += edep;}
  //double getEnergy() const {return m_edep;}

  /**
   * @brief Set the ADC information
   * @param[in] adc ADC value
   */
  //void setAdc(const short adc) { m_adc = adc; }

  /**
   * @brief Get ADC value
   * @param[out] ADC value
   */
  //short getAdc() const { return m_adc; }
 private:
  //short m_adc;
  //double m_edep;
  ClassDef(MvtxHit, 1);
};

#endif  //MVTX_MVTXHIT_H
