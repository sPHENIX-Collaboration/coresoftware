/**
 * @file mvtx/MvtxHit.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of Mvtx hit object
 */
#include "MvtxHit.h"

#include <trackbase/TrkrHit.h>

#include <ostream>              // for operator<<, endl, ostream, basic_ostream

MvtxHit::MvtxHit()
  : TrkrHit()
{
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

