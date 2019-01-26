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
  m_adc = 0;
  m_edep = 0.0;
}

void MvtxHit::identify(std::ostream& os) const
{
  os << "I am an MvtxHit" << std::endl;
}

void MvtxHit::Reset()
{
  TrkrHit::Reset();
}

int MvtxHit::isValid() const
{
  return 1;
}

