#include "TrackAnalysisUtils.h"

#include <cmath>
#include "SvtxTrack.h"

namespace TrackAnalysisUtils
{

  TrackAnalysisUtils::DCAPair get_dca(SvtxTrack* track,
                                      Acts::Vector3& vertex)
  {
    TrackAnalysisUtils::DCAPair pair;
    if (!track)
    {
      std::cout << "You passed a nullptr, can't calculate DCA" << std::endl;
      return pair;
    }

    Acts::Vector3 pos(track->get_x(),
                      track->get_y(),
                      track->get_z());
    Acts::Vector3 mom(track->get_px(),
                      track->get_py(),
                      track->get_pz());

    pos -= vertex;

    Acts::ActsSquareMatrix<3> posCov;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        posCov(i, j) = track->get_error(i, j);
      }
    }

    Acts::Vector3 r = mom.cross(Acts::Vector3(0., 0., 1.));
    float phi = atan2(r(1), r(0));
    phi *= -1;
    Acts::RotationMatrix3 rot;
    Acts::RotationMatrix3 rot_T;
    rot(0, 0) = std::cos(phi);
    rot(0, 1) = -std::sin(phi);
    rot(0, 2) = 0;
    rot(1, 0) = std::sin(phi);
    rot(1, 1) = std::cos(phi);
    rot(1, 2) = 0;
    rot(2, 0) = 0;
    rot(2, 1) = 0;
    rot(2, 2) = 1;

    rot_T = rot.transpose();

    Acts::Vector3 pos_R = rot * pos;
    Acts::ActsSquareMatrix<3> rotCov = rot * posCov * rot_T;
    //! First pair is DCA_xy +/- DCA_xy_err
    pair.first.first = pos_R(0);
    pair.first.second = sqrt(rotCov(0, 0));

    //! Second pair is DCA_z +/- DCA_z_err
    pair.second.first = pos_R(2);
    pair.second.second = sqrt(rotCov(2, 2));

    return pair;
  }

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track)
  {
    std::vector<TrkrDefs::cluskey> out;
    for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (seed)
      {
        std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
      }
    }
    return out;
  }

}  // namespace TrackAnalysisUtils
