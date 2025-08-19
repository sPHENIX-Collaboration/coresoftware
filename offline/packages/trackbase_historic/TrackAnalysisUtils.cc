#include "TrackAnalysisUtils.h"

#include "SvtxTrack.h"
#include "TrackSeed.h"
#include <trackbase/TrkrCluster.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterContainer.h>

#include <cmath>

namespace TrackAnalysisUtils
{
  float calc_dedx(TrackSeed* tpcseed, 
                  TrkrClusterContainer *cluster_map,
                  ActsGeometry* tgeometry,
                  float thickness_per_region[4])
  {

    std::vector<TrkrDefs::cluskey> clusterKeys;
    clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
                       tpcseed->end_cluster_keys());

    std::vector<float> dedxlist;
    for (unsigned long cluster_key : clusterKeys)
    {
      unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
      if (TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::TrkrId::tpcId)
      {
        continue;
      }
      TrkrCluster* cluster = cluster_map->findCluster(cluster_key);
      Acts::Vector3 cglob;
      cglob = tgeometry->getGlobalPosition(cluster_key, cluster);
      //      float x = cglob(0);
      // float y = cglob(1);
      float z = cglob(2);
      int csegment = 0;
      float thickness = 0;
      if (layer_local >= 7 + 32)
      {
        csegment = 2;
        thickness = thickness_per_region[3];
      }
      else if (layer_local >= 7 + 16)
      {
        csegment = 1;
        thickness = thickness_per_region[2];
      }
      if(csegment == 0)
      {
        if (layer_local % 2 == 0)
        {
          thickness = thickness_per_region[1];
        }
        else
        {
          thickness = thickness_per_region[0];
        }
      }
      float adc = cluster->getAdc();

      float r = std::sqrt(cglob(0) * cglob(0) + cglob(1) * cglob(1));
      float alpha = (r * r) / (2 * r * std::abs(1.0 / tpcseed->get_qOverR()));
      float beta = std::atan(tpcseed->get_slope());
      float alphacorr = std::cos(alpha);
      if (alphacorr < 0 || alphacorr > 4)
      {
        alphacorr = 4;
      }
      float betacorr = std::cos(beta);
      if (betacorr < 0 || betacorr > 4)
      {
        betacorr = 4;
      }
      adc /= (1 - ((105 - abs(z)) * 0.50 / 105));
      adc /= thickness;
      adc *= alphacorr;
      adc *= betacorr;
      dedxlist.push_back(adc);
      sort(dedxlist.begin(), dedxlist.end());
    }
    int trunc_min = 0;
    int trunc_max = (int) dedxlist.size() * 0.7;
    float sumdedx = 0;
    int ndedx = 0;
    for (int j = trunc_min; j <= trunc_max; j++)
    {
      sumdedx += dedxlist.at(j);
      ndedx++;
    }
    sumdedx /= ndedx;
    return sumdedx;
  }
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
