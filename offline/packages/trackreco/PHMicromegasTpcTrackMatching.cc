#include "PHMicromegasTpcTrackMatching.h"

//#include "PHTrackPropagating.h"     // for PHTrackPropagating

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>

/// Tracking includes
#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrClusterv3.h>  // for TrkrCluster
#include <trackbase/TrkrDefs.h>       // for cluskey, getLayer, TrkrId
#include <trackbase_historic/TrackSeedHelper.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <tpc/TpcDistortionCorrectionContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TF1.h>
#include <TVector3.h>

#include <array>
#include <cmath>     // for sqrt, std::abs, atan2, cos
#include <iostream>  // for operator<<, basic_ostream
#include <map>       // for map
#include <set>       // for _Rb_tree_const_iterator
#include <utility>   // for pair, make_pair
#include <optional>  // for line-circle function

namespace
{

  //! convenience square method
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  //! get radius from x and y
  template <class T>
  inline constexpr T get_r(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }

  //! assemble list of cluster keys associated to seeds
  std::vector<TrkrDefs::cluskey> get_cluster_keys( const std::vector<TrackSeed*>& seeds )
  {
    std::vector<TrkrDefs::cluskey> out;
    for( const auto& seed: seeds )
    {
      if( seed )
      { std::copy( seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter( out ) ); }
    }
    return out;
  }

  /// calculate intersection from circle to line, in 2d. return true on success
  /**
  * circle is defined as (x-xc)**2 + (y-yc)**2 = r**2
  * line is defined as nx(x-x0) + ny(y-y0) = 0
  * to solve we substitute y by y0 - nx/ny*(x-x0) in the circle equation and solve the 2nd order polynom
  * there is the extra complication that ny can be 0 (vertical line) to prevent this, we multiply all terms of the polynom by ny**2
  * and account for this special case when calculating x from y
  */
  bool circle_line_intersection(
      double r, double xc, double yc,
      double x0, double y0, double nx, double ny,
      double& xplus, double& yplus, double& xminus, double& yminus)
  {
    if (ny == 0)
    {
      // vertical lines are defined by ny=0 and x = x0
      xplus = xminus = x0;

      // calculate y accordingly
      const double delta = square(r) - square(x0 - xc);
      if (delta < 0)
      {
        return false;
      }

      const double sqdelta = std::sqrt(delta);
      yplus = yc + sqdelta;
      yminus = yc - sqdelta;
    }
    else
    {
      const double a = square(nx) + square(ny);
      const double b = -2. * (square(ny) * xc + square(nx) * x0 + nx * ny * (y0 - yc));
      const double c = square(ny) * (square(xc) - square(r)) + square(ny * (y0 - yc) + nx * x0);
      const double delta = square(b) - 4. * a * c;
      if (delta < 0)
      {
        return false;
      }

      const double sqdelta = std::sqrt(delta);
      xplus = (-b + sqdelta) / (2. * a);
      xminus = (-b - sqdelta) / (2. * a);

      yplus = y0 - (nx / ny) * (xplus - x0);
      yminus = y0 - (nx / ny) * (xminus - x0);
    }

    return true;
  }

  bool line_line_intersection(
      double m, double b,
      double x0, double y0, double nx, double ny,
      double& xplus, double& yplus, double& xminus, double& yminus)
  {
    if (ny == 0)
    {
      // vertical lines are defined by ny=0 and x = x0
      xplus = xminus = x0;

      // calculate y accordingly
      yplus = yminus = m * x0 + b;
    }
    else
    {

      double denom = nx + ny*m;
      if(denom == 0) {
        return false; // lines are parallel and there is no intersection
      }

      double x = (nx*x0 + ny*y0 - ny*b)/denom;
      double y = m*x + b;
      // a straight line has a unique intersection point
      xplus = xminus = x;
      yplus = yminus = y;

    }

    return true;
  }

  // streamer of TVector3
  [[maybe_unused]] inline std::ostream& operator<<(std::ostream& out, const TVector3& vector)
  {
    out << "( " << vector.x() << "," << vector.y() << "," << vector.z() << ")";
    return out;
  }

}  // namespace

//____________________________________________________________________________..
PHMicromegasTpcTrackMatching::PHMicromegasTpcTrackMatching(const std::string& name)
  : SubsysReco(name)
{
}
//____________________________________________________________________________..
int PHMicromegasTpcTrackMatching::Init(PHCompositeNode* /* topNode */)
{
  std::cout
            << "PHMicromegasTpcTrackMatching::Init - "
            << " rphi_search_win inner layer " << _rphi_search_win[0]
            << " z_search_win inner layer " << _z_search_win[0]
            << " rphi_search_win outer layer " << _rphi_search_win[1]
            << " z_search_win outer layer " << _z_search_win[1]
            << std::endl;

  std::cout << "PHMicromegasTpcTrackMatching::Init - _use_silicon: " << _use_silicon << std::endl;
  std::cout << "PHMicromegasTpcTrackMatching::Init - _zero_field: " << _zero_field << std::endl;
  std::cout << "PHMicromegasTpcTrackMatching::Init - _min_tpc_layer: " << _min_tpc_layer << std::endl;
  std::cout << "PHMicromegasTpcTrackMatching::Init - _max_tpc_layer: " << _max_tpc_layer << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHMicromegasTpcTrackMatching::InitRun(PHCompositeNode* topNode)
{
  // load micromegas geometry
  _geomContainerMicromegas = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  if (!_geomContainerMicromegas)
  {
    std::cout << PHWHERE << "Could not find CYLINDERGEOM_MICROMEGAS_FULL." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // ensures there are at least two micromegas layers
  if (_geomContainerMicromegas->get_NLayers() != _n_mm_layers)
  {
    std::cout << PHWHERE << "Inconsistent number of micromegas layers." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get first micromegas layer
  _min_mm_layer = static_cast<CylinderGeomMicromegas*>(_geomContainerMicromegas->GetFirstLayerGeom())->get_layer();

  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  return ret;
}

//____________________________________________________________________________..
int PHMicromegasTpcTrackMatching::process_event(PHCompositeNode* topNode)
{
  if (_n_iteration > 0)
  {
    _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    if (!_iteration_map)
    {
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  // We will add the micromegas cluster to the TPC tracks already on the node tree

  _event++;

  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << " Event " << _event << " Seed track map size " << _svtx_seed_map->size() << std::endl;
  }

  // loop over the seed tracks - these are the seeds formed from matched tpc and silicon track seeds
  for (unsigned int seedID = 0;
       seedID != _svtx_seed_map->size(); ++seedID)
  {
    auto seed = _svtx_seed_map->get(seedID);
    auto siID = seed->get_silicon_seed_index();
    auto tracklet_si = _si_track_map->get(siID);

    short int crossing = 0;
    if (_pp_mode)
    {
      if (!tracklet_si)
      {
        continue;  // cannot use tracks not matched to silicon because crossing is unknown
      }

      crossing = tracklet_si->get_crossing();
      if (crossing == SHRT_MAX)
      {
        if (Verbosity() > 0)
        {
          std::cout << " svtx seed " << seedID << " with si seed " << siID
                    << " crossing not defined: crossing = " << crossing << " skip this track" << std::endl;
        }
        continue;
      }
    }

    auto tpcID = seed->get_tpc_seed_index();

    auto tracklet_tpc = _tpc_track_map->get(tpcID);
    if (!tracklet_tpc)
    {
      continue;
    }

    if (Verbosity() >= 1)
    {
      std::cout << std::endl
                << __LINE__
                << ": Processing TPC seed track: " << tpcID
                << ": crossing: " << crossing
                << ": nhits: " << tracklet_tpc->size_cluster_keys()
                << ": Total TPC tracks: " << _tpc_track_map->size()
                << ": phi: " << tracklet_tpc->get_phi()
                << std::endl;
    }

    // Get the outermost TPC clusters for this tracklet
    std::vector<Acts::Vector3> clusGlobPos;
    std::vector<Acts::Vector3> clusGlobPos_silicon;
    std::vector<Acts::Vector3> clusGlobPos_mvtx;

    // try extrapolate track to Micromegas and find corresponding tile
    /* this is all copied from PHMicromegasTpcTrackMatching */
    const auto cluster_keys = get_cluster_keys({tracklet_tpc,tracklet_si});
    for( const auto& cluster_key:cluster_keys )
    {

      // detector id and layer
      const auto detid = TrkrDefs::getTrkrId(cluster_key);
      switch( detid )
      {

        case TrkrDefs::tpcId:
        {
          // layer
          const unsigned int layer = TrkrDefs::getLayer(cluster_key);

          // check layer range
          if( layer < _min_tpc_layer || layer >= _max_tpc_layer )
          { continue; }

          // get matching
          const auto cluster = _cluster_map->findCluster(cluster_key);
          clusGlobPos.push_back( m_globalPositionWrapper.getGlobalPositionDistortionCorrected(cluster_key, cluster, crossing) );
          break;
        }

        case TrkrDefs::mvtxId:
        {
          // get matching
          const auto cluster = _cluster_map->findCluster(cluster_key);
          const auto global_position = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(cluster_key, cluster, crossing);
          clusGlobPos_silicon.push_back( global_position );
          clusGlobPos_mvtx.push_back( global_position );
          break;
        }

        case TrkrDefs::inttId:
        {
          // get matching
          const auto cluster = _cluster_map->findCluster(cluster_key);
          const auto global_position = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(cluster_key, cluster, crossing);
          clusGlobPos_silicon.push_back( global_position );
          break;
        }

        default:
        break;
      }

    }

    // check number of clusters
    if( _use_silicon )
    {

      if( clusGlobPos_mvtx.size()<3 ) { continue; }

    } else {

      if( clusGlobPos.size()<3 ) { continue; }

    }

    // r,z linear fit
    /* applies to both field ON and field OFF configurations */
    const auto [slope_rz, intersect_rz] = _use_silicon ?
      TrackFitUtils::line_fit(clusGlobPos_mvtx):
      TrackFitUtils::line_fit(clusGlobPos);

    if (Verbosity() > 10)
    {
      std::cout << " r,z fitted line has slope_rz " << slope_rz << " intersect_rz " << intersect_rz << std::endl;
    }

    // x,y straight fit paramers
    double slope_xy = 0, intersect_xy = 0;

    // circle fit parameters
    double R = 0, X0 = 0, Y0 = 0;

    if(_zero_field) {

      if (Verbosity() > 10)
      {
        std::cout << "zero field is ON, starting TPC clusters linear fit" << std::endl;
      }

      // x,y straight fit
      std::tie( slope_xy, intersect_xy ) = _use_silicon ?
        TrackFitUtils::line_fit_xy(clusGlobPos_silicon):
        TrackFitUtils::line_fit_xy(clusGlobPos);

      if (Verbosity() > 10)
      {
        std::cout << " zero field x,y fit has slope_xy " << slope_xy << " intersect_xy " << intersect_xy << std::endl;
      }

    } else {

      if(Verbosity() > 10)
      {
        std::cout << "zero field is OFF, starting TPC clusters circle fit" << std::endl;
      }

      // x,y circle
      std::tie( R, X0, Y0 ) = _use_silicon ?
        TrackFitUtils::circle_fit_by_taubin(clusGlobPos_silicon):
        TrackFitUtils::circle_fit_by_taubin(clusGlobPos);

      // toss tracks for which the fitted circle could not have come from the vertex
      if (R < 40.0)
      {
        continue;
      }

    }

    // loop over micromegas layer
    for (unsigned int imm = 0; imm < _n_mm_layers; ++imm)
    {
      // get micromegas geometry object
      const unsigned int layer = _min_mm_layer + imm;
      const auto layergeom = static_cast<CylinderGeomMicromegas*>(_geomContainerMicromegas->GetLayerGeom(layer));
      const auto layer_radius = layergeom->get_radius();

      // get intersection to track
      auto [xplus, yplus, xminus, yminus] =
        _zero_field ?
        TrackFitUtils::line_circle_intersection(layer_radius, slope_xy, intersect_xy):
        TrackFitUtils::circle_circle_intersection(layer_radius, R, X0, Y0);

      if (Verbosity() > 10)
      {
        std::cout << "xplus: " << xplus << " yplus " << yplus << " xminus " << xminus << " yminus " << std::endl;
      }

      if (!std::isfinite(xplus))
       {
         if (Verbosity() > 10)
         {
           std::cout << PHWHERE << " circle/circle intersection calculation failed, skip this case" << std::endl;
           std::cout << PHWHERE << " mm_radius " << layer_radius << " fitted R " << R << " fitted X0 " << X0 << " fitted Y0 " << Y0 << std::endl;
         }

         continue;
      }
      // we can figure out which solution is correct based on the last cluster position in the TPC
      const double last_clus_phi = std::atan2(clusGlobPos.back()(1), clusGlobPos.back()(0));
      double phi_plus = std::atan2(yplus, xplus);
      double phi_minus = std::atan2(yminus, xminus);

      // calculate z
      double r = layer_radius;
      double z = intersect_rz + slope_rz * r;

      // select the angle that is the closest to last cluster
      // store phi, apply coarse space charge corrections in calibration mode
      double phi = std::abs(last_clus_phi - phi_plus) < std::abs(last_clus_phi - phi_minus) ? phi_plus : phi_minus;

      // create cylinder intersection point in world coordinates
      const TVector3 world_intersection_cylindrical(r * std::cos(phi), r * std::sin(phi), z);

      // find matching tile
      int tileid = layergeom->find_tile_cylindrical(world_intersection_cylindrical);
      if (tileid < 0)
      {
        continue;
      }

      // get tile center and norm vector
      const auto tile_center = layergeom->get_world_from_local_coords(tileid, _tGeometry, {0, 0});
      const double x0 = tile_center.x();
      const double y0 = tile_center.y();

      const auto tile_norm = layergeom->get_world_from_local_vect(tileid, _tGeometry, {0, 0, 1});
      const double nx = tile_norm.x();
      const double ny = tile_norm.y();

      if(_zero_field) {

        // calculate intersection to tile
        if (!line_line_intersection(slope_xy, intersect_xy, x0, y0, nx, ny, xplus, yplus, xminus, yminus))
        {
          if (Verbosity() > 10)
          {
            std::cout << PHWHERE << "line_line_intersection - failed" << std::endl;
          }
          continue;
        }

      } else {

        // calculate intersection to tile
        if (!circle_line_intersection(R, X0, Y0, x0, y0, nx, ny, xplus, yplus, xminus, yminus))
        {
          if (Verbosity() > 10)
          {
            std::cout << PHWHERE << "circle_line_intersection - failed" << std::endl;
          }
          continue;
        }

      }

      // select again angle closest to last cluster
      phi_plus = std::atan2(yplus, xplus);
      phi_minus = std::atan2(yminus, xminus);
      const bool is_plus = (std::abs(last_clus_phi - phi_plus) < std::abs(last_clus_phi - phi_minus));

      // calculate x, y and z
      const double x = (is_plus ? xplus : xminus);
      const double y = (is_plus ? yplus : yminus);
      r = get_r(x, y);
      z = intersect_rz + slope_rz * r;

      /*
       * create planar intersection point in world coordinates
       * this is the position to be compared to the clusters
       */
      const TVector3 world_intersection_planar(x, y, z);

      // convert to tile local reference frame, apply SC correction
      const auto local_intersection_planar = layergeom->get_local_from_world_coords(tileid, _tGeometry, world_intersection_planar);

      // store segmentation type
      const auto segmentation_type = layergeom->get_segmentation_type();

      // generate tilesetid and get corresponding clusters
      const auto tilesetid = MicromegasDefs::genHitSetKey(layer, segmentation_type, tileid);
      const auto mm_clusrange = _cluster_map->getClusters(tilesetid);

      // do nothing if cluster range is empty
      if( mm_clusrange.first == mm_clusrange.second )
      { continue; }

      // keep track of cluster with smallest distance to local intersection
      double drphi_min = 0;
      double dz_min = 0;
      TrkrDefs::cluskey ckey_min = 0;
      bool first = true;
      for (auto clusiter = mm_clusrange.first; clusiter != mm_clusrange.second; ++clusiter)
      {
        const auto& [ckey, cluster] = *clusiter;
        if (_iteration_map)
        {
          if (_iteration_map->getIteration(ckey) > 0)
          {
            continue;
          }
        }

        // compute residuals and store
        /* in local tile coordinate, x is along rphi, and z is along y) */
        const double drphi = local_intersection_planar.x() - cluster->getLocalX();
        const double dz = local_intersection_planar.y() - cluster->getLocalY();
        switch( segmentation_type )
        {
          case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
          {
            // reject if outside of strip boundary
            if( std::abs(dz)>_z_search_win[imm] )
            { continue; }

            // keep as best if closer to projection
            if( first || std::abs(drphi) < std::abs(drphi_min) )
            {
              first = false;
              drphi_min = drphi;
              dz_min = dz;
              ckey_min = ckey;
            }
            break;
          }

          case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
          {
            // reject if outside of strip boundary
            if( std::abs(drphi)>_rphi_search_win[imm] )
            { continue; }

            // keep as best if closer to projection
            if( first || std::abs(dz) < std::abs(dz_min) )
            {
              first = false;
              drphi_min = drphi;
              dz_min = dz;
              ckey_min = ckey;
            }
            break;
          }
        }

        // prints out a line that can be grep-ed from the output file to feed to a display macro
        // compare to cuts and add to track if matching
        if( _test_windows && std::abs(drphi) < _rphi_search_win[imm] && std::abs(dz) < _z_search_win[imm])
        {
          // cluster rphi and z
          const auto glob = _tGeometry->getGlobalPosition(ckey, cluster);
          const double mm_clus_rphi = get_r(glob.x(), glob.y()) * std::atan2(glob.y(), glob.x());
          const double mm_clus_z = glob.z();

          // projection phi and z, without correction
          const double rphi_proj = get_r(world_intersection_planar.x(), world_intersection_planar.y()) * std::atan2(world_intersection_planar.y(), world_intersection_planar.x());
          const double z_proj = world_intersection_planar.z();

          /*
           * Note: drphi and dz might not match the difference of the rphi and z quoted values. This is because
           * 1/ drphi and dz are actually calculated in Tile's local reference frame, not in world coordinates
           * 2/ drphi also includes SC distortion correction, which the world coordinates don't
          */
          std::cout
            << "  Try_mms: " << (int) layer
            << " drphi " << drphi
            << " dz " << dz
            << " mm_clus_rphi " << mm_clus_rphi << " mm_clus_z " << mm_clus_z
            << " rphi_proj " << rphi_proj << " z_proj " << z_proj
            << " pt " << tracklet_tpc->get_pt()
            << " charge " << tracklet_tpc->get_charge()
            << std::endl;
        }
      }  // end loop over clusters

      // compare to cuts and add to track if matching
      if( (!first) && ckey_min > 0 && std::abs(drphi_min) < _rphi_search_win[imm] && std::abs(dz_min) < _z_search_win[imm])
      {
        tracklet_tpc->insert_cluster_key(ckey_min);
        if (Verbosity() > 0)
        {
          std::cout << " Match to MM's found for seedID " << seedID << " tpcID " << tpcID << " siID " << siID << std::endl;
        }
      }

    }  // end loop over Micromegas layers

    if (Verbosity() > 3)
    {
      tracklet_tpc->identify();
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << " Final seed map size " << _svtx_seed_map->size() << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_________________________________________________________________________________________________
int PHMicromegasTpcTrackMatching::End(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//_________________________________________________________________________________________________
int PHMicromegasTpcTrackMatching::GetNodes(PHCompositeNode* topNode)
{
  // all clusters
  if (_use_truth_clusters)
  {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_TRUTH");
  }
  else
  {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }

  if (!_cluster_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!_tGeometry)
  {
    std::cerr << PHWHERE << "No acts tracking geometry, can't continue."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _svtx_seed_map = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  if (!_svtx_seed_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find "
              << "SvtxTrackSeedContainer" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tpc_track_map = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_tpc_track_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find "
              << "TpcTrackSeedContainer" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _si_track_map = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!_si_track_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find "
              << "SiliconTrackSeedContainer" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // micromegas geometry
  _geomContainerMicromegas = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  if (!_geomContainerMicromegas)
  {
    std::cout << PHWHERE << "Could not find CYLINDERGEOM_MICROMEGAS_FULL." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // global position wrapper
  m_globalPositionWrapper.loadNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<TrkrDefs::cluskey> PHMicromegasTpcTrackMatching::getTrackletClusterList(TrackSeed* tracklet)
{
  std::vector<TrkrDefs::cluskey> cluskey_vec;
  for (auto clusIter = tracklet->begin_cluster_keys();
       clusIter != tracklet->end_cluster_keys();
       ++clusIter)
  {
    auto key = *clusIter;
    auto cluster = _cluster_map->findCluster(key);
    if (!cluster)
    {
      if(Verbosity() > 1)
      {
        std::cout << PHWHERE << "Failed to get cluster with key " << key << std::endl;
      }
      continue;
    }

    /// Make a safety check for clusters that couldn't be attached to a surface
    auto surf = _tGeometry->maps().getSurface(key, cluster);
    if (!surf)
    {
      continue;
    }

    // drop some bad layers in the TPC completely
    unsigned int layer = TrkrDefs::getLayer(key);
    if (layer == 7 || layer == 22 || layer == 23 || layer == 38 || layer == 39)
    {
      continue;
    }

    /* if (layer > 2 && layer < 7) */
    /* { */
      /* continue; */
    /* } */



    cluskey_vec.push_back(key);
  }  // end loop over clusters for this track
  return cluskey_vec;
}
