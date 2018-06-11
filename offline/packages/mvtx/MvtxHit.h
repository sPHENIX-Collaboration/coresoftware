/**
 * @file mvtx/MvtxHit.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Mvtx hit object
 */
#ifndef MVTX_MVTXHIT_H
#define MVTX_MVTXHIT_H

#include <trackbase/TrkrHit.h>

#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif

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
  virtual ~MvtxHit() {};

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

private:
  ClassDef(MvtxHit,1);
};

#endif //MVTX_MVTXHIT_H
