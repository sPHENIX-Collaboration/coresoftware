#include "SvtxTrackInfo_v2.h"
#include "TrackStateInfo_v1.h"

/**
 * SvtxTrackInfo_v2 is a class developed to project the track information from
 * the outermost TPC surface to the calorimeters
 **/

void SvtxTrackInfo_v2::CopyFrom(const SvtxTrackInfo& source)
{
  set_track_id(source.get_track_id());
  set_subsurfkey(source.get_subsurfkey());

  set_chisq(source.get_chisq());
  set_ndf(source.get_ndf());
  set_hitbitmap(source.get_hitbitmap());
  set_crossing(source.get_crossing());

  set_x(source.get_x());
  set_y(source.get_y());
  set_z(source.get_z());
  set_phi(source.get_phi());
  set_theta(source.get_theta());
  set_qOp(source.get_qOp());

  set_x_outer_tpc(source.get_x_outer_tpc());
  set_y_outer_tpc(source.get_y_outer_tpc());
  set_z_outer_tpc(source.get_z_outer_tpc());
  set_phi_outer_tpc(source.get_phi_outer_tpc());
  set_theta_outer_tpc(source.get_theta_outer_tpc());
  set_qOp_outer_tpc(source.get_qOp_outer_tpc());

  for (int i = 0; i < 5; i++)
  {
    for (int j = i; j < 5; j++)
    {
      set_covariance(i, j, source.get_covariance(i, j));
      set_covariance_outer_tpc(i, j, source.get_covariance_outer_tpc(i, j));
    }
  }
}

SvtxTrackInfo_v2& SvtxTrackInfo_v2::operator=(const SvtxTrackInfo_v2& source)
{
  if (this != &source)
  {
    CopyFrom(source);
  }
  return *this;
}
