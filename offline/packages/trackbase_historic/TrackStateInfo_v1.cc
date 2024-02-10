#include "TrackStateInfo_v1.h"
#include <trackbase/TrkrDefs.h>

#include <cmath>
#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <utility>  // for pair

namespace
{
  // get unique index in cov. matrix array from i and j
  inline unsigned int covar_index(unsigned int i, unsigned int j)
  {
    if (i > j)
    {
      std::swap(i, j);
    }
    return i + 1 + (j + 1) * (j) / 2 - 1;
  }
}  // namespace
float TrackStateInfo_v1::get_covariance(int i, int j) const
{
  return m_Covariance[covar_index(i, j)];
}
void TrackStateInfo_v1::set_covariance(int i, int j, float value)
{
  m_Covariance[covar_index(i, j)] = value;
}