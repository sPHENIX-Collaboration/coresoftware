#include "SvtxTrackInfo_v3.h"
#include "TrackStateInfo_v1.h"

/**
 * SvtxTrackInfo_v3 is a class developed to hold the information from the track
 * projection from the outermost TPC surface to the calorimeters
**/

void SvtxTrackInfo_v3::CopyFrom(const SvtxTrackInfo& source)
{
  set_track_id(source.get_track_id());
  set_subsurfkey(source.get_subsurfkey());
  set_chisq(source.get_chisq());
  set_ndf(source.get_ndf());
  set_hitbitmap(source.get_hitbitmap());
  set_crossing(source.get_crossing());

  set_x(STATE::VERTEX, source.get_x());
  set_y(STATE::VERTEX, source.get_y());
  set_z(STATE::VERTEX, source.get_z());
  set_phi(STATE::VERTEX, source.get_phi());
  set_theta(STATE::VERTEX, source.get_theta());
  set_qOp(STATE::VERTEX, source.get_qOp());

  for(int istate = STATE::OUTER_TPC; istate <= STATE::HCALOUT_BACKFACE; istate++)
  {
    set_x(istate, source.get_x(istate));
    set_y(istate, source.get_y(istate));
    set_z(istate, source.get_z(istate));
    set_phi(istate, source.get_phi(istate));
    set_theta(istate, source.get_theta(istate));
    set_qOp(istate, source.get_qOp(istate));

    for (int i = 0; i < 5; i++)
    {
      for (int j = i; j < 5; j++)
      {
        set_covariance(istate, i, j, source.get_covariance(istate, i, j));
      }
    }
  }

  if(std::isnan(source.get_phi_outer_tpc()) || std::isnan(source.get_theta_outer_tpc()) || std::isnan(source.get_qOp_outer_tpc()))
  {
    return;
  }

  set_x(STATE::OUTER_TPC, source.get_x_outer_tpc());
  set_y(STATE::OUTER_TPC, source.get_y_outer_tpc());
  set_z(STATE::OUTER_TPC, source.get_z_outer_tpc());
  set_phi(STATE::OUTER_TPC, source.get_phi_outer_tpc());
  set_theta(STATE::OUTER_TPC, source.get_theta_outer_tpc());
  set_qOp(STATE::OUTER_TPC, source.get_qOp_outer_tpc());

  for (int i = 0; i < 5; i++)
  {
    for (int j = i; j < 5; j++)
    {
      set_covariance(STATE::VERTEX, i, j, source.get_covariance(i, j));
      set_covariance(STATE::OUTER_TPC, i, j, source.get_covariance_outer_tpc(i, j));
    }
  }

}

SvtxTrackInfo_v3& SvtxTrackInfo_v3::operator=(const SvtxTrackInfo_v3& source)
{
  if(this != &source)
  {
    CopyFrom(source);
  }
  return *this;
}
