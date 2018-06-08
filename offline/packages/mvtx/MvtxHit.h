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

  /**
   * @brief Set hit information
   * @param[in] col Column index
   * @param[in] row Row index
   */
  void setColumnRow(uint16_t col, uint16_t row);

  /**
   * @brief Get column index
   * @param[out] Column index
   */
  uint16_t getColumn() { return m_col; }

  /**
   * @brief Get row index
   * @param[out] Row index
   */
  uint16_t getRow() { return m_row; }

private:
  uint16_t m_col; /// column index [,]
  uint16_t m_row; /// row index [,]
  ClassDef(MvtxHit,1);
};

#endif //MVTX_MVTXHIT_H
