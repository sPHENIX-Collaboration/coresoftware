/**
 * @file tpc/TpcHit.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of Tpc hit object
 */
#include "TpcHit.h"
#include "TpcDefs.h"

TpcHit::TpcHit()
  : TrkrHit()
{
}

void 
TpcHit::identify(std::ostream& os) const
{
  os << "TpcHit with key:" << getKey() 
     << " and col:" << getColumn() 
     << " and row:" << getRow() 
     << std::endl;
}

void 
TpcHit::Reset()
{
  TrkrHit::Reset();
}

int 
TpcHit::isValid() const
{
  // valid if the key is not equal to the default value
  return getKey() != TrkrDefs::HITKEYMAX;
}

void
TpcHit::setColumnRow(uint16_t col, uint16_t row)
{
  if ( col > TpcDefs::MAXCOL || row > TpcDefs::MAXROW )
  {
    std::cout << "TpcHit::setColumnRow() - Error - Invalid input! " << std::endl
	      << "  col:" << col << " (MAXCOL:" << TpcDefs::MAXCOL << ")" << std::endl
	      << "  row:" << row << " (MAXROW:" << TpcDefs::MAXROW << ")" << std::endl;
    // reset the hit if we used invalid inputs
    exit(2);
  }

  setKey(TpcDefs::genHitKey(col, row));
}

uint16_t
TpcHit::getColumn() const
{
  return TpcDefs::getCol(getKey());
}

uint16_t
TpcHit::getRow() const
{
  return TpcDefs::getRow(getKey());
}
