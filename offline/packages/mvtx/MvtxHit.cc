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
{
}

void 
MvtxHit::identify(std::ostream& os) const
{
  os << "MvtxHit with key:" << getKey() 
     << " and col:" << getColumn() 
     << " and row:" << getRow() 
     << std::endl;
}

void 
MvtxHit::Reset()
{
  TrkrHit::Reset();
}

int 
MvtxHit::isValid() const
{
  // valid if the key is not equal to the default value
  return getKey() != TrkrDefs::HITKEYMAX;
}

void
MvtxHit::setColumnRow(uint16_t col, uint16_t row)
{
  if ( col > MvtxDefs::MAXCOL || row > MvtxDefs::MAXROW )
  {
    std::cout << "MvtxHit::setColumnRow() - Error - Invalid input! " << std::endl
	      << "  col:" << col << " (MAXCOL:" << MvtxDefs::MAXCOL << ")" << std::endl
	      << "  row:" << row << " (MAXROW:" << MvtxDefs::MAXROW << ")" << std::endl;
    // reset the hit if we used invalid inputs
    exit(2);
  }

  setKey(MvtxDefs::genHitKey(col, row));
}

uint16_t
MvtxHit::getColumn() const
{
  return MvtxDefs::getCol(getKey());
}

uint16_t
MvtxHit::getRow() const
{
  return MvtxDefs::getRow(getKey());
}
