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
 * a single hit strip within a chip
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
   * @brief Set hit information
   * @param[in] col Column index
   * @param[in] row Row index
   * 
   * This function sets the key of the TrkrHit base class
   * using the column and row information given. This is
   * meant as a convenience interface.
   */
  void setColumnRow(uint16_t col, uint16_t row);

  /**
   * @brief Set the ADC information
   * @param[in] adc ADC value
   */
  void setAdc(const int adc) { m_adc = adc; }

  /**
   * @brief Get column index
   * @param[out] Column index
   *
   * Get the column index from the hitkey value
   */
  uint16_t getColumn() const;

  /**
   * @brief Get row index
   * @param[out] Row index
   *
   * Get the row index from the hitkey value
   */
  uint16_t getRow() const;

  /**
   * @brief Get ADC value
   * @param[out] ADC value
   */
  int getAdc() const { return m_adc; }

private:
  int m_adc;
  ClassDef(TpcHit,1);
};

#endif //TPC_TPCHIT_H
