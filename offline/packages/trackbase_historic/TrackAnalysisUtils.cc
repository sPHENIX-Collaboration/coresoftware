#include "TrackAnalysisUtils.h"

#include <globalvertex/GlobalVertex.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <g4detectors/PHG4TpcGeomContainer.h>
#include "SvtxTrack.h"
#include "TrackSeed.h"

#include <cmath>

namespace TrackAnalysisUtils
{

  float calc_dedx(TrackSeed* tpcseed,
                  TrkrClusterContainer* clustermap,
                  ActsGeometry* tgeometry,
                  float thickness_per_region[4])
  {
    std::vector<TrkrDefs::cluskey> clusterKeys;
    clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
                       tpcseed->end_cluster_keys());

    std::vector<float> dedxlist;
    for (unsigned long cluster_key : clusterKeys)
    {
      auto detid = TrkrDefs::getTrkrId(cluster_key);
      if (detid != TrkrDefs::TrkrId::tpcId)
      {
        continue;  // the micromegas clusters are added to the TPC seeds
      }
      unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
      TrkrCluster* cluster = clustermap->findCluster(cluster_key);
      float adc = cluster->getAdc();
      float thick = 0;
      if (layer_local < 23)
      {
        if (layer_local % 2 == 0)
        {
          thick = thickness_per_region[1];
        }
        else
        {
          thick = thickness_per_region[0];
        }
      }
      else if (layer_local < 39)
      {
        thick = thickness_per_region[2];
      }
      else
      {
        thick = thickness_per_region[3];
      }
      auto cglob = tgeometry->getGlobalPosition(cluster_key, cluster);
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
      adc /= thick;
      adc *= alphacorr;
      adc *= betacorr;
      dedxlist.push_back(adc);
      sort(dedxlist.begin(), dedxlist.end());
    }
    int trunc_min = 0;
    if (dedxlist.empty())
    {
      return std::numeric_limits<float>::quiet_NaN();
    }
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

  float calc_dedx_calib(SvtxTrack* track,
                        TrkrClusterContainer* cluster_map,
                        ActsGeometry* tgeometry,
                        float thickness_per_region[4])
  {
    auto clusterKeys = get_cluster_keys(track->get_tpc_seed());

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
      if (csegment == 0)
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
      auto tpcseed = track->get_tpc_seed();
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
      if (track->get_crossing() < SHRT_MAX)
      {
        double z_crossing_corrected =
            TpcClusterZCrossingCorrection::correctZ(cglob.z(),
                                                    TpcDefs::getSide(cluster_key), track->get_crossing());

        double maxz = tgeometry->get_max_driftlength() + tgeometry->get_CM_halfwidth();
        adc /= (1 - ((maxz - abs(z_crossing_corrected)) * 0.50 / maxz));
      }
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
                                      GlobalVertex* vertex)
  {
    Acts::Vector3 vpos(vertex->get_x(),
                       vertex->get_y(),
                       vertex->get_z());
    Acts::Vector3 mom(track->get_px(),
                      track->get_py(),
                      track->get_pz());
    auto dca = get_dca(track, vpos);
    auto rot = rotationMatrixToLocal(mom);
    Acts::RotationMatrix3 rot_T = rot.transpose();

    Acts::ActsSquareMatrix<3> posCov;
    Acts::ActsSquareMatrix<3> vertexCov;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        posCov(i, j) = track->get_error(i, j);
        vertexCov(i, j) = vertex->get_error(i, j);
      }
    }

    Acts::ActsSquareMatrix<3> rotCov = rot * (posCov + vertexCov) * rot_T;
    dca.first.second = sqrt(rotCov(0, 0));
    dca.second.second = sqrt(rotCov(2, 2));

    return dca;
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

    auto rot = rotationMatrixToLocal(mom);
    Acts::RotationMatrix3 rot_T = rot.transpose();

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
  Acts::RotationMatrix3 rotationMatrixToLocal(const Acts::Vector3& mom)
  {
    Acts::Vector3 r = mom.cross(Acts::Vector3(0., 0., 1.));
    float phi = atan2(r(1), r(0));
    phi *= -1;
    Acts::RotationMatrix3 rot;

    rot(0, 0) = std::cos(phi);
    rot(0, 1) = -std::sin(phi);
    rot(0, 2) = 0;
    rot(1, 0) = std::sin(phi);
    rot(1, 1) = std::cos(phi);
    rot(1, 2) = 0;
    rot(2, 0) = 0;
    rot(2, 1) = 0;
    rot(2, 2) = 1;
    return rot;
  }
  std::vector<TrkrDefs::cluskey> get_cluster_keys(TrackSeed* seed)
  {
    std::vector<TrkrDefs::cluskey> out;

    if (seed)
    {
      std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
    }

    return out;
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

  std::pair<Acts::Vector2, Acts::Vector3>
  get_residual(TrkrDefs::cluskey& ckey, SvtxTrack* track, TrkrClusterContainer* clustermap,
                PHCompositeNode* topNode)
  {
    TpcGlobalPositionWrapper globalWrapper;
    globalWrapper.loadNodes(topNode);
    globalWrapper.set_suppressCrossing(true);
    TpcClusterMover mover;
    auto* tpccellgeo = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
    mover.initialize_geometry(tpccellgeo);
    mover.set_verbosity(0);
    auto* geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    auto* cluster = clustermap->findCluster(ckey);
    std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;
    for (const auto& key : get_cluster_keys(track))
    {
      auto* clus = clustermap->findCluster(key);

      // Fully correct the cluster positions for the crossing and all distortions
      Acts::Vector3 global = globalWrapper.getGlobalPositionDistortionCorrected(key, clus, track->get_crossing());
      // add the global positions to a vector to give to the cluster mover
      global_raw.emplace_back(key, global);
    }

    auto global_moved = mover.processTrack(global_raw);

    // loop over global vectors and get this cluster
    Acts::Vector3 clusglob(0, 0, 0);
    for (const auto& pair : global_raw)
    {
      auto thiskey = pair.first;
      clusglob = pair.second;
      if (thiskey == ckey)
      {
        break;
      }
    }

    Acts::Vector3 clusglob_moved(0, 0, 0);
    for (const auto& pair : global_moved)
    {
      auto thiskey = pair.first;
      clusglob_moved = pair.second;
      if (thiskey == ckey)
      {
        break;
      }
    }
    SvtxTrackState* state = nullptr;
    for (auto state_iter = track->begin_states();
         state_iter != track->end_states();
         ++state_iter)
    {
      SvtxTrackState* tstate = state_iter->second;
      auto stateckey = tstate->get_cluskey();
      if (stateckey == ckey)
      {
        state = tstate;
        break;
      }
    }
    Surface surf = geometry->maps().getSurface(ckey, cluster);
    Surface surf_ideal = geometry->maps().getSurface(ckey, cluster);  // Unchanged by distortion corrections
    // if this is a TPC cluster, the crossing correction may have moved it across the central membrane, check the surface
    auto trkrid = TrkrDefs::getTrkrId(ckey);
    if (trkrid == TrkrDefs::tpcId)
    {
      TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(ckey);
      TrkrDefs::subsurfkey new_subsurfkey = 0;
      surf = geometry->get_tpc_surface_from_coords(hitsetkey, clusglob_moved, new_subsurfkey);
    }

    auto loc = geometry->getLocalCoords(ckey, cluster, track->get_crossing());
    // in this case we get local coords from transform of corrected global coords
    clusglob_moved *= Acts::UnitConstants::cm;  // we want mm for transformations
    Acts::Vector3 normal = surf->normal(geometry->geometry().getGeoContext(),
                                        Acts::Vector3(1, 1, 1), Acts::Vector3(1, 1, 1));
    auto local = surf->globalToLocal(geometry->geometry().getGeoContext(),
                                     clusglob_moved, normal);
    if (local.ok())
    {
      loc = local.value() / Acts::UnitConstants::cm;
    }
    else
    {
      // otherwise take the manual calculation for the TPC
      // doing it this way just avoids the bounds check that occurs in the surface class method
      Acts::Vector3 loct = surf->transform(geometry->geometry().getGeoContext()).inverse() * clusglob_moved;  // global is in mm
      loct /= Acts::UnitConstants::cm;

      loc(0) = loct(0);
      loc(1) = loct(1);
    }
   clusglob_moved /= Acts::UnitConstants::cm;  // we want cm for the tree
    Acts::Vector2 stateloc(state->get_localX(), state->get_localY());
    Acts::Vector3 stateglob(state->get_x(), state->get_y(), state->get_z());

    return std::make_pair(stateloc - loc, stateglob - clusglob_moved);
  }

}  // namespace TrackAnalysisUtils
