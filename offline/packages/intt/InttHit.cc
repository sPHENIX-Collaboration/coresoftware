/**
 * @file intt/InttHit.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of Intt hit object
 */
#include "InttHit.h"
#include "InttDefs.h"

InttHit::InttHit()
  : TrkrHit()
{
}

void 
InttHit::identify(std::ostream& os) const
{
  os << "InttHit with key:" << getKey() 
     << " and col:" << getColumn() 
     << " and row:" << getRow() 
     << std::endl;
}

void 
InttHit::Reset()
{
  TrkrHit::Reset();
}

int 
InttHit::isValid() const
{
  // valid if the key is not equal to the default value
  return getKey() != TrkrDefs::HITKEYMAX;
}

void
InttHit::setColumnRow(uint16_t col, uint16_t row)
{
  if ( col > InttDefs::MAXCOL || row > InttDefs::MAXROW )
  {
    std::cout << "InttHit::setColumnRow() - Error - Invalid input! " << std::endl
	      << "  col:" << col << " (MAXCOL:" << InttDefs::MAXCOL << ")" << std::endl
	      << "  row:" << row << " (MAXROW:" << InttDefs::MAXROW << ")" << std::endl;
    // reset the hit if we used invalid inputs
    exit(2);
  }

  setKey(InttDefs::genHitKey(col, row));
}

uint16_t
InttHit::getColumn() const
{
  return InttDefs::getCol(getKey());
}

uint16_t
InttHit::getRow() const
{
  return InttDefs::getRow(getKey());
}
