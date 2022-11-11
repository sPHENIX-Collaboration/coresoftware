#include "TrkrClusterIterationMap.h"

#include <TSystem.h>
#include <iostream>

namespace
{
  TrkrClusterIterationMap::Map dummy;
}

void TrkrClusterIterationMap::Reset()
{
  std::cout << "TrkrClusterIterationMap: Reset() not implemented by daughter class" << std::endl;
  gSystem->Exit(1);
}

TrkrClusterIterationMap::ConstIter TrkrClusterIterationMap::begin() const
{
  return dummy.end();
}

TrkrClusterIterationMap::ConstIter TrkrClusterIterationMap::end() const
{
  return dummy.end();
}
