
#include "SvtxTrackInfo_v1.h"
#include "TrackStateInfo_v1.h"

#include <trackbase/TrkrDefs.h>

#include <cmath>
#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <utility>  // for pair

void SvtxTrackInfo_v1::CopyFrom(const SvtxTrackInfo& source)
{
  set_chisq(source.get_chisq());
  set_ndf(source.get_ndf());
  set_hitbitmap(source.get_hitbitmap());
  set_crossing(source.get_crossing());

  set_x(source.get_x());
  set_y(source.get_y());
  set_z(source.get_z());
  set_px(source.get_px());
  set_py(source.get_py());
  set_pz(source.get_pz());

  for (int i = 0; i < 6; i++)
  {
    for (int j = i; j < 6; j++)
    {
      set_covariance(i, j, source.get_covariance(i, j));
    }
  }
}

SvtxTrackInfo_v1& SvtxTrackInfo_v1::operator=(const SvtxTrackInfo_v1& source)
{
  if (this != &source)
  {
    CopyFrom(source);
  }
  return *this;
}