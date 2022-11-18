#include "MakeMilleFiles.h"

#include "Mille.h"

/// Tracking includes

#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <trackbase/TpcDefs.h>  // for side

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>

#include <fun4all/Fun4AllReturnCodes.h>

#include <Acts/Definitions/Algebra.hpp>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TF1.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <utility>

namespace
{
  /// square
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

//____________________________________________________________________________..
MakeMilleFiles::MakeMilleFiles(const std::string& name)
  : SubsysReco(name)
  , _mille(nullptr)
  , _trajectories(nullptr)
{
}

//____________________________________________________________________________..
int MakeMilleFiles::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  // Instantiate Mille and open output data file
  //  _mille = new Mille(data_outfilename.c_str(), false);   // write text in data files, rather than binary, for debugging only
  _mille = new Mille(data_outfilename.c_str(), _binary);

  // Write the steering file here, and add the data file path to it
  std::ofstream steering_file(steering_outfilename.c_str());
  steering_file << data_outfilename.c_str() << std::endl;
  steering_file.close();

  // print grouping setup to log file:
  std::cout << "MakeMilleFiles::InitRun: Surface groupings are silicon " << si_group << " tpc " << tpc_group << " mms " << mms_group << std::endl;

  return ret;
}

//____________________________________________________________________________..
int MakeMilleFiles::process_event(PHCompositeNode* /*topNode*/)
{
  // Outline:
  //
  // loop over tracks
  //   Make any track cuts here to skip undesirable tracks (maybe low pT?)
  //   loop over track states+measurements for each track
  //      for each measurement
  //         Get measurement value and error (global, what to use for error?)
  //         Calculate derivatives and residuals from Acts jacobians
  //         Rotate residual-derivative matrices to global coordinates
  //         These are stored in a map and unpacked for mille
  //   Call _mille->mille() with arguments obtained from previous iteration:
  //     local pars
  //     array of local derivatives
  //     global pars
  //     array of global derivatives
  //     array of integer global par labels
  //     residual value (float) z = measurement - track state
  //     sigma of measurement
  //   After processing all measurements for this track, call _mille->end() to add buffer to file and reset buffer
  // After all tracks are processed, file is closed when Mille destructor is called
  // Note: all units are in the Acts units of mm and GeV to avoid converting matrices

  if (Verbosity() > 0)
    std::cout << PHWHERE << " track map size " << _track_map->size() << std::endl;

  for(auto [key, track] : *_track_map)
  {
    auto crossing = track->get_silicon_seed()->get_crossing();

    if (Verbosity() > 0)
    {
      std::cout << std::endl
                << __LINE__ << ": Processing track itrack: " << key << ": nhits: " << track->size_cluster_keys()
                << ": Total tracks: " << _track_map->size() << ": phi: " << track->get_phi() << std::endl;
    }

    // Make any desired track cuts here
    // Maybe set a lower pT limit - low pT tracks are not very sensitive to alignment

    /// Get the corresponding acts trajectory to look up the state info
    const auto traj = _trajectories->find(key)->second;

    AlignmentStateMap alignStates = getAlignmentStates(traj, track, crossing);

    addTrackToMilleFile(alignStates, traj);

    /// Finish this track
    _mille->end();
 
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MakeMilleFiles::End(PHCompositeNode* /*topNode*/)
{
  delete _mille;

  return Fun4AllReturnCodes::EVENT_OK;
}

int MakeMilleFiles::GetNodes(PHCompositeNode* topNode)
{
  _trajectories = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsTrajectories");
  if (!_trajectories)
  {
    std::cout << PHWHERE << "ERROR: Can't find Acts trajectories" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!_tGeometry)
  {
    std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc distortion corrections
  _dcc_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
  if (_dcc_static)
  {
    std::cout << PHWHERE << "  found static TPC distortion correction container" << std::endl;
  }
  _dcc_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerAverage");
  if (_dcc_average)
  {
    std::cout << PHWHERE << "  found average TPC distortion correction container" << std::endl;
  }
  _dcc_fluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerFluctuation");
  if (_dcc_fluctuation)
  {
    std::cout << PHWHERE << "  found fluctuation TPC distortion correction container" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 MakeMilleFiles::getPCALinePoint(Acts::Vector3 global, SvtxTrackState* state)
{
  // Approximate track with a straight line consisting of the state position and the vector (px,py,pz)

  Acts::Vector3 track_dir(state->get_px(), state->get_py(), state->get_pz());
  track_dir = track_dir / track_dir.norm();
  Acts::Vector3 track_base(state->get_x(), state->get_y(), state->get_z());

  // The position of the closest point on the line is:
  // track_base + projection of difference between the point and track_base on the line vector
  Acts::Vector3 pca = track_base + ((global - track_base).dot(track_dir)) * track_dir;

  return pca;
}

std::vector<Acts::Vector3> MakeMilleFiles::getDerivativesAlignmentAngles(Acts::Vector3& global, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster, Surface surface, int crossing)
{
  // The value of global is from the geocontext transformation
  // we add to that transformation a small rotation around the relevant axis in the surface frame

  std::vector<Acts::Vector3> derivs_vector;

  // get the transformation from the geocontext
  Acts::Transform3 transform = surface->transform(_tGeometry->geometry().getGeoContext());

  // Make an additional transform that applies a small rotation angle in the surface frame
  for (unsigned int iangle = 0; iangle < 3; ++iangle)
  {
    // creates transform that adds a perturbation rotation around one axis, uses it to estimate derivative wrt perturbation rotation

    unsigned int trkrId = TrkrDefs::getTrkrId(cluster_key);
    unsigned int layer = TrkrDefs::getLayer(cluster_key);

    Acts::Vector3 derivs(0, 0, 0);
    Eigen::Vector3d theseAngles(0, 0, 0);
    theseAngles[iangle] = sensorAngles[iangle];  // set the one we want to be non-zero

    Acts::Vector3 keeper(0, 0, 0);
    for (int ip = 0; ip < 2; ++ip)
    {
      if (ip == 1)
      {
        theseAngles[iangle] *= -1;
      }  // test both sides of zero

      if (Verbosity() > 1)
      {
        std::cout << "     trkrId " << trkrId << " layer " << layer << " cluster_key " << cluster_key
                  << " sensorAngles " << theseAngles[0] << "  " << theseAngles[1] << "  " << theseAngles[2] << std::endl;
      }

      Acts::Transform3 perturbationTransformation = makePerturbationTransformation(theseAngles);
      Acts::Transform3 overallTransformation = transform * perturbationTransformation;

      // transform the cluster local position to global coords with this additional rotation added
      auto x = cluster->getLocalX() * 10;  // mm
      auto y = cluster->getLocalY() * 10;
      if (trkrId == TrkrDefs::tpcId)
      {
        y = convertTimeToZ(cluster_key, cluster);
      }

      Eigen::Vector3d clusterLocalPosition(x, y, 0);                                      // follows the convention for the acts transform of local = (x,z,y)
      Eigen::Vector3d finalCoords = overallTransformation * clusterLocalPosition / 10.0;  // convert mm back to cm

      // have to add corrections for TPC clusters after transformation to global
      if (trkrId == TrkrDefs::tpcId)
      {
        makeTpcGlobalCorrections(cluster_key, crossing, global);
      }

      if (ip == 0)
      {
        keeper(0) = (finalCoords(0) - global(0));
        keeper(1) = (finalCoords(1) - global(1));
        keeper(2) = (finalCoords(2) - global(2));
      }
      else
      {
        keeper(0) -= (finalCoords(0) - global(0));
        keeper(1) -= (finalCoords(1) - global(1));
        keeper(2) -= (finalCoords(2) - global(2));
      }

      if (Verbosity() > 1)
      {
        std::cout << "        finalCoords(0) " << finalCoords(0) << " global(0) " << global(0) << " finalCoords(1) " << finalCoords(1)
                  << " global(1) " << global(1) << " finalCoords(2) " << finalCoords(2) << " global(2) " << global(2) << std::endl;
        std::cout << "        keeper now:  keeper(0) " << keeper(0) << " keeper(1) " << keeper(1) << " keeper(2) " << keeper(2) << std::endl;
      }
    }

    // derivs vector contains:
    //   (dx/dalpha,     dy/dalpha,     dz/dalpha)     (for iangle = 0)
    //   (dx/dbeta,        dy/dbeta,        dz/dbeta)        (for iangle = 1)
    //   (dx/dgamma, dy/dgamma, dz/dgamma) (for iangle = 2)

    // Average the changes to get the estimate of the derivative
    // Check what the sign of this should be !!!!
    derivs(0) = keeper(0) / (2.0 * fabs(theseAngles[iangle]));
    derivs(1) = keeper(1) / (2.0 * fabs(theseAngles[iangle]));
    derivs(2) = keeper(2) / (2.0 * fabs(theseAngles[iangle]));
    derivs_vector.push_back(derivs);

    if (Verbosity() > 1)
    {
      std::cout << "        derivs(0) " << derivs(0) << " derivs(1) " << derivs(1) << " derivs(2) " << derivs(2) << std::endl;
    }
  }

  return derivs_vector;
}

Acts::Transform3 MakeMilleFiles::makePerturbationTransformation(Acts::Vector3 angles)
{
  // Note: Here beta is apllied to the z axis and gamma is applied to the y axis because the geocontext transform
  // will flip those axes when transforming to global coordinates
  Eigen::AngleAxisd alpha(angles(0), Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd beta(angles(2), Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd gamma(angles(1), Eigen::Vector3d::UnitZ());
  Eigen::Quaternion<double> q = gamma * beta * alpha;
  Eigen::Matrix3d perturbationRotation = q.matrix();

  // combine rotation and translation into an affine matrix
  Eigen::Vector3d nullTranslation(0, 0, 0);
  Acts::Transform3 perturbationTransformation;
  perturbationTransformation.linear() = perturbationRotation;
  perturbationTransformation.translation() = nullTranslation;

  return perturbationTransformation;
}

float MakeMilleFiles::convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster* cluster)
{
  // must convert local Y from cluster average time of arival to local cluster z position
  double drift_velocity = _tGeometry->get_drift_velocity();
  double zdriftlength = cluster->getLocalY() * drift_velocity;
  double surfCenterZ = 52.89;                // 52.89 is where G4 thinks the surface center is
  double zloc = surfCenterZ - zdriftlength;  // converts z drift length to local z position in the TPC in north
  unsigned int side = TpcDefs::getSide(cluster_key);
  if (side == 0) zloc = -zloc;
  float z = zloc * 10.0;

  return z;
}

void MakeMilleFiles::makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global)
{
  // make all corrections to global position of TPC cluster
  unsigned int side = TpcDefs::getSide(cluster_key);
  float z = m_clusterCrossingCorrection.correctZ(global[2], side, crossing);
  global[2] = z;

  // apply distortion corrections
  if (_dcc_static)
  {
    global = _distortionCorrection.get_corrected_position(global, _dcc_static);
  }
  if (_dcc_average)
  {
    global = _distortionCorrection.get_corrected_position(global, _dcc_average);
  }
  if (_dcc_fluctuation)
  {
    global = _distortionCorrection.get_corrected_position(global, _dcc_fluctuation);
  }
}

SvtxTrack::StateIter MakeMilleFiles::getStateIter(Acts::Vector3& global, SvtxTrack* track)
{
  float clus_radius = sqrt(global[0] * global[0] + global[1] * global[1]);
  auto state_iter = track->begin_states();
  float dr_min = -1;
  //for( auto iter = state_iter; iter != track->end_states(); ++iter )
  for (auto iter = track->begin_states(); iter != track->end_states(); ++iter)
  {
    const auto dr = std::abs(clus_radius - sqrt(iter->second->get_x() * iter->second->get_x() + iter->second->get_y() * iter->second->get_y()));
    if (dr_min < 0 || dr < dr_min)
    {
      state_iter = iter;
      dr_min = dr;
    }
    else
    {
      break;
    }
  }

  if (Verbosity() > 2)
  {
    std::cout << "track state for pathlength " << clus_radius << " dr_min " << dr_min << " position "
              << state_iter->second->get_x() << "  "
              << state_iter->second->get_y() << "  "
              << state_iter->second->get_z() << std::endl;
  }
  return state_iter;
}

int MakeMilleFiles::getTpcRegion(int layer)
{
  int region = 0;
  if(layer > 23 && layer < 39)
    region = 1;
  if(layer > 38 && layer < 55)
    region = 2;

  return region;  
}

int MakeMilleFiles::getLabelBase(Acts::GeometryIdentifier id)
{
  unsigned int volume = id.volume();
  unsigned int acts_layer = id.layer();
  unsigned int layer = base_layer_map.find(volume)->second + acts_layer / 2 - 1;
  unsigned int sensor = id.sensitive() - 1;  // Acts starts at 1

  int label_base = 1;  // Mille wants to start at 1

  // decide what level of grouping we want
  if (layer < 7)
    {
      if(si_group == siliconGroup::sensor)
	{
	  // every sensor has a different label
	  int stave = sensor / nsensors_stave[layer];
	  label_base += layer*1000000 + stave*10000 + sensor*10;
	  return label_base;
	}
      if(si_group == siliconGroup::stave)
	{
	  // layer and stave, assign all sensors to the stave number
	  int stave = sensor / nsensors_stave[layer];
	  label_base += layer*1000000 + stave*10000;
	  return label_base;
	}
      if(si_group == siliconGroup::barrel)
	{
	  // layer only, assign all sensors to sensor 0 
	  label_base += layer*1000000 + 0;
      
	  return label_base;
	}
    }
  else if (layer > 6 && layer < 55)
    {
      if(tpc_group == tpcGroup::hitset)
	{
	  // want every hitset (layer, sector, side) to have a separate label
	  // each group of 12 subsurfaces (sensors) is in a single hitset
	  int hitset = sensor/12; // hitsets 0-11 on side 0, 12-23 on side 1
	  label_base += layer*1000000 + hitset*10000;
	  return label_base;
	}
      if(tpc_group == tpcGroup::sector)
	{
	  // group all tpc layers in each region and sector, assign layer 7 and side and sector number to all layers and hitsets
	  int side = sensor / 144; // 0-143 on side 0, 144-287 on side 1
	  int sector = (sensor - side *144) / 12; 
	  // for a given layer there are only 12 sectors x 2 sides
	  // The following gives the sectors in the inner, mid, outer regions unique group labels
	  int region = getTpcRegion(layer);  // inner, mid, outer
	  label_base += 7*1000000 + (region * 24 + side*12 + sector) *10000; 
	  std::cout << " layer " << layer << " sensor " << sensor << " region " << region << " side " << side << " sector " << sector << " label_base " << label_base << std::endl;
	  return label_base;
	}
      if(tpc_group == tpcGroup::tpc)
	{
	  // all tpc layers and all sectors, assign layer 7 and sensor 0 to all layers and sensors
	  label_base += 7*1000000 + 0;
	  return label_base;
	}
    }
  else
    {
      if(mms_group == mmsGroup::tile)
	{
	  // every tile has different label
	  int tile = sensor;
	  label_base += layer*1000000 + tile*10000 + sensor*10;
	  return label_base;
	}
      if(mms_group == mmsGroup::mms)
	{
	  // assign layer 55 and tile 0 to all
	  label_base += 55*1000000 + 0;	  
	  return label_base;
	}
    }

  return -1;
}
  
AlignmentStateMap MakeMilleFiles::getAlignmentStates(const Trajectory& traj,
                                                     SvtxTrack* track,
                                                     short int crossing)
{
  const auto mj = traj.multiTrajectory();
  const auto& tips = traj.tips();
  const auto& trackTip = tips.front();
  std::map<TrkrDefs::cluskey, AlignmentState> alignStates;

  /// Collect the track states
  mj.visitBackwards(trackTip, [&](const auto& state) {
    /// Collect only track states which were used in smoothing of KF and are measurements
    if (not state.hasSmoothed() or
        not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag))
    {
      return true;
    }

    const auto& surface = state.referenceSurface();
    const auto& sl = static_cast<const ActsExamples::IndexSourceLink&>(state.uncalibrated());
    auto ckey = sl.cluskey();

    if (Verbosity() > 2)
    {
      std::cout << "sl index and ckey " << sl.index() << ", "
                << sl.cluskey() << std::endl;
    }

    auto clus = _cluster_map->findCluster(ckey);
    const auto trkrId = TrkrDefs::getTrkrId(ckey);

    if (!clus)
    {
      /// For some reason, the SL index and cluster key of the first state in the Acts trajectory
      /// lose scope outside of the track fit result when the fitter is first called. They are
      /// reset to either 0 or the largest possible value. So we have
      /// to look up the first cluster manually
      if (sl.index() == 0 or sl.index() > 57)
      {
        auto siseed = track->get_silicon_seed();
        ckey = *(siseed->begin_cluster_keys());
        clus = _cluster_map->findCluster(ckey);
      }

      /// if we still couldn't find it, skip it
      if (!clus)
      {
        if (Verbosity() > 2)
        {
          std::cout << "no cluster with key " << ckey << " and index " << sl.index() << std::endl;
        }

        return true;
      }
    }

    /// Gets the global parameters from the state
    const Acts::FreeVector freeParams =
        Acts::MultiTrajectoryHelpers::freeSmoothed(_tGeometry->geometry().getGeoContext(), state);

    const Acts::ActsDynamicMatrix measCovariance =
        state.effectiveCalibratedCovariance();

    /// Calculate the residual in global coordinates
    Acts::Vector3 clusGlobal = _tGeometry->getGlobalPosition(ckey, clus);
    if (trkrId == TrkrDefs::tpcId)
    {
      makeTpcGlobalCorrections(ckey, crossing, clusGlobal);
    }

    /// convert to acts units
    clusGlobal *= Acts::UnitConstants::cm;

    const Acts::FreeVector globalStateParams = Acts::detail::transformBoundToFreeParameters(surface, _tGeometry->geometry().getGeoContext(), state.smoothed());
    Acts::Vector3 stateGlobal = globalStateParams.segment<3>(Acts::eFreePos0);
    Acts::Vector3 residual = clusGlobal - stateGlobal;

    if (Verbosity() > 2)
    {
      Acts::Vector3 clus_sigma(0, 0, 0);
      if (_cluster_version == 3)
      {
        clus_sigma(2) = clus->getZError() * Acts::UnitConstants::cm;
        clus_sigma(0) = clus->getRPhiError() / sqrt(2) * Acts::UnitConstants::cm;
        clus_sigma(1) = clus->getRPhiError() / sqrt(2) * Acts::UnitConstants::cm;
      }
      else if (_cluster_version == 4)
      {
        double clusRadius = sqrt(clusGlobal(0) * clusGlobal(0) + clusGlobal(1) * clusGlobal(1));
        auto para_errors = _ClusErrPara.get_simple_cluster_error(clus, clusRadius, ckey);
        float exy2 = para_errors.first * Acts::UnitConstants::cm2;
        float ez2 = para_errors.second * Acts::UnitConstants::cm2;
        clus_sigma(2) = sqrt(ez2);
        clus_sigma(0) = sqrt(exy2 / 2.0);
        clus_sigma(1) = sqrt(exy2 / 2.0);
      }
      std::cout << "clus global is " << clusGlobal.transpose() << std::endl
                << "state global is " << stateGlobal.transpose() << std::endl
                << "Residual is " << residual.transpose() << std::endl;
      std::cout << "   clus errors are " << clus_sigma.transpose() << std::endl;
    }

    // Get the derivative of alignment (global) parameters w.r.t. measurement or residual
    const Acts::Vector3 direction = freeParams.segment<3>(Acts::eFreeDir0);
    // The derivative of free parameters w.r.t. path length. @note Here, we
    // assumes a linear track model, i.e. negecting the change of track
    // direction. Otherwise, we need to know the magnetic field at the free
    // parameters
    Acts::FreeVector pathDerivative = Acts::FreeVector::Zero();
    pathDerivative.head<3>() = direction;

    const Acts::ActsDynamicMatrix H = state.effectiveProjector();

    /// Acts residual, in local coordinates
    auto actslocres = state.effectiveCalibrated() - H * state.smoothed();

    // Get the derivative of bound parameters w.r.t. alignment parameters
    Acts::AlignmentToBoundMatrix d =
        surface.alignmentToBoundDerivative(_tGeometry->geometry().getGeoContext(), freeParams, pathDerivative);
    // Get the derivative of bound parameters wrt track parameters
    Acts::FreeToBoundMatrix j = surface.freeToBoundJacobian(_tGeometry->geometry().getGeoContext(), freeParams);

    // derivative of residual wrt track parameters
    auto dLocResTrack = -H * j;
    // derivative of residual wrt alignment parameters
    auto dLocResAlignment = -H * d;

    if (Verbosity() > 4)
    {
      std::cout << "local resids " << actslocres.transpose() << std::endl
                << " derivative of resiudal wrt track params " << std::endl
                << dLocResTrack << std::endl
                << " derivative of residual wrt alignment params " << std::endl
                << dLocResAlignment << std::endl;
    }

    /// The above matrices are in Acts local coordinates, so we need to convert them to
    /// global coordinates with a rotation d(l0,l1)/d(x,y,z)

    const double cosTheta = direction.z();
    const double sinTheta = std::sqrt(square(direction.x()) + square(direction.y()));
    const double invSinTheta = 1. / sinTheta;
    const double cosPhi = direction.x() * invSinTheta;
    const double sinPhi = direction.y() * invSinTheta;
    Acts::ActsMatrix<3, 2> rot = Acts::ActsMatrix<3, 2>::Zero();
    rot(0, 0) = -sinPhi;
    rot(0, 1) = -cosPhi * cosTheta;
    rot(1, 0) = cosPhi;
    rot(1, 1) = -sinPhi * cosTheta;
    rot(2, 1) = sinTheta;

    /// this is a 3x6 matrix now of d(x_res,y_res,z_res)/d(dx,dy,dz,alpha,beta,gamma)
    const auto dGlobResAlignment = rot * dLocResAlignment;

    /// Acts global coordinates are (x,y,z,t,hatx,haty,hatz,q/p)
    /// So now rotate to x,y,z, px,py,pz
    const double p = 1. / abs(globalStateParams[Acts::eFreeQOverP]);
    const double p2 = square(p);
    Acts::ActsMatrix<6, 8> sphenixRot = Acts::ActsMatrix<6, 8>::Zero();
    sphenixRot(0, 0) = 1;
    sphenixRot(1, 1) = 1;
    sphenixRot(2, 2) = 1;
    sphenixRot(3, 4) = p;
    sphenixRot(4, 5) = p;
    sphenixRot(5, 6) = p;
    sphenixRot(3, 7) = direction.x() * p2;
    sphenixRot(4, 7) = direction.y() * p2;
    sphenixRot(5, 7) = direction.z() * p2;

    const auto dGlobResTrack = rot * dLocResTrack * sphenixRot.transpose();

    if (Verbosity() > 3)
    {
      std::cout << "derivative of residual wrt alignment parameters glob " << std::endl
                << dGlobResAlignment << std::endl;
      std::cout << "derivative of residual wrt trakc parameters glob " << std::endl
                << dGlobResTrack << std::endl;
    }

   

    AlignmentState astate(state.index(), residual, dGlobResAlignment, dGlobResTrack, clusGlobal);

    /// To switch to creating things WRT local coordinates, change 
    /// AlignmentState::NLC to 8 and AlignmentState::NRES to 2 and 
    /// uncomment the following
    /// AlignmentState astate(state.index(actslocres, dLocResAlignment, dLocResTrack, clusGlobal);

    if(m_useAnalytic)
      {
	auto surf = _tGeometry->maps().getSurface(ckey, clus);
	auto anglederivs = getDerivativesAlignmentAngles(clusGlobal, ckey, 
							 clus, surf, 
							 crossing);
	AlignmentState::GlobalMatrix analytic = AlignmentState::GlobalMatrix::Zero();
        analytic(0,0) = 1;
        analytic(1,1) = 1;
        analytic(2,2) = 1;
	for(int i=0; i<AlignmentState::NRES; ++i) 
	  {
	    for(int j=3; j<AlignmentState::NGL; ++j)
	      {
		/// convert to mm
		analytic(i,j) = anglederivs.at(i)(j-3) * Acts::UnitConstants::cm;
	      }
	  }

	astate.set_dResAlignmentPar(analytic);
      }
  
    alignStates.insert(std::make_pair(ckey, astate));
    
    return true;
    });

  return alignStates;
}

void MakeMilleFiles::addTrackToMilleFile(AlignmentStateMap& alignStates, const Trajectory& traj)
{
  for (auto& [ckey, astate] : alignStates)
  {
    if (Verbosity() > 2)
      {
	std::cout << "adding state for ckey " << ckey << std::endl;
      }
    // The global alignment parameters are given initial values of zero by default, we do not specify them
    // We identify the global alignment parameters for this surface
    const auto cluster = _cluster_map->findCluster(ckey);

    const auto residual = astate.get_residual();
    const auto& global = astate.get_clusglob();

    // need standard deviation of measurements
    AlignmentState::ResidualVector clus_sigma = AlignmentState::ResidualVector::Zero();
    if (AlignmentState::NRES == 3)
    {
      if (_cluster_version == 3)
      {
        clus_sigma(2) = cluster->getZError() * Acts::UnitConstants::cm;
        clus_sigma(0) = cluster->getRPhiError() / sqrt(2) * Acts::UnitConstants::cm;
        clus_sigma(1) = cluster->getRPhiError() / sqrt(2) * Acts::UnitConstants::cm;
      }
      else if (_cluster_version == 4)
      {
        double clusRadius = sqrt(global[0] * global[0] + global[1] * global[1]);
        auto para_errors = _ClusErrPara.get_simple_cluster_error(cluster, clusRadius, ckey);
        float exy2 = para_errors.first * Acts::UnitConstants::cm2;
        float ez2 = para_errors.second * Acts::UnitConstants::cm2;
        clus_sigma(2) = sqrt(ez2);
        clus_sigma(0) = sqrt(exy2 / 2.0);
        clus_sigma(1) = sqrt(exy2 / 2.0);
      }

      if (std::isnan(clus_sigma(0)) ||
          std::isnan(clus_sigma(1)) ||
          std::isnan(clus_sigma(2)))
      {
        continue;
      }
    }
    else if (AlignmentState::NRES == 2)
    {
      if (_cluster_version == 3)
      {
        clus_sigma(1) = cluster->getZError() * Acts::UnitConstants::cm;
        clus_sigma(0) = cluster->getRPhiError() * Acts::UnitConstants::cm;
      }
      else if (_cluster_version == 4)
      {
        double clusRadius = sqrt(global[0] * global[0] + global[1] * global[1]);
        auto para_errors = _ClusErrPara.get_simple_cluster_error(cluster, clusRadius, ckey);
        float exy2 = para_errors.first * Acts::UnitConstants::cm2;
        float ez2 = para_errors.second * Acts::UnitConstants::cm2;
        clus_sigma(1) = sqrt(ez2);
        clus_sigma(0) = sqrt(exy2);
      }
    }

    const auto& multiTraj = traj.multiTrajectory();
    const auto& state = multiTraj.getTrackState(astate.get_tsIndex());
    Acts::GeometryIdentifier id = state.referenceSurface().geometryId();
    int label_base = getLabelBase(id);  // This value depends on how the surfaces are grouped

    int glbl_label[AlignmentState::NGL];
    for (int i = 0; i < AlignmentState::NGL; ++i)
    {
      glbl_label[i] = label_base + i;
      if (Verbosity() > 1)
      {
        std::cout << "  glbl " << i << " label " << glbl_label[i] << " ";
      }
    }

    if (Verbosity() > 1)
    {
      std::cout << std::endl;
    }

    /// For N residual coordinates x,y,z or local x,z
    for (int i = 0; i < AlignmentState::NRES; ++i)
    {
      // Add the measurement separately for each coordinate direction to Mille
      float glbl_derivative[AlignmentState::NGL];
      for (int j = 0; j < AlignmentState::NGL; ++j)
      {
        glbl_derivative[j] = astate.get_dResAlignmentPar()(i, j);
      }

      float lcl_derivative[AlignmentState::NLC];
      for (int j = 0; j < AlignmentState::NLC; ++j)
      {
        lcl_derivative[j] = astate.get_dResTrackPar()(i, j);
      }
      if (Verbosity() > 2)
      {
        std::cout << "coordinate " << i << " has residual " << residual(i) << " and clus_sigma " << clus_sigma(i) << std::endl
                  << "global deriv " << std::endl;
        
        for (int k = 0; k < AlignmentState::NGL; k++)
        {
          std::cout << glbl_derivative[k] << ", ";
        }
        std::cout << std::endl
                  << "local deriv " << std::endl;
        for (int k = 0; k < AlignmentState::NLC; k++)
        {
          std::cout << lcl_derivative[k] << ", ";
        }
        std::cout << std::endl;
      }

      if (clus_sigma(i) < 1.0)  // discards crazy clusters
      {
        _mille->mille(AlignmentState::NLC, lcl_derivative, AlignmentState::NGL, glbl_derivative, glbl_label, residual(i), clus_sigma(i));
      }
    }
  }
 
  return;
}
