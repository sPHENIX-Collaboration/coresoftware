/**
 * @file mvtx/MvtxHit.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of Mvtx hit object
 */
#include "MvtxHit.h"
#include "MvtxDefs.h"

MvtxHit::MvtxHit()
  : TrkrHit()
  , m_col(MvtxDefs::MAXCOL)
  , m_row(MvtxDefs::MAXROW)
{
}

void
MvtxHit::setColumnRow(uint16_t col, uint16_t row)
{
  if ( col > MvtxDefs::MAXCOL )
  {
  }
  if ( row > MvtxDefs::MAXROW )
  {
  }

  m_col = col;
  m_row = row;

}
