#include "MakeMilleFiles.h"

#include "Mille.h"

/// Tracking includes

#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>  // for side
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxAlignmentState.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/TrackAnalysisUtils.h>

#include <trackreco/ActsPropagator.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

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
{
}

//____________________________________________________________________________..
int MakeMilleFiles::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  // Instantiate Mille and open output data file
  //  _mille = new Mille(data_outfilename.c_str(), false);   // write text in data files, rather than binary, for debugging only
  _mille = new Mille(data_outfilename.c_str(), _binary);

  // Write the steering file here, and add the data file path to it
  std::ofstream steering_file(steering_outfilename);
  steering_file << data_outfilename << std::endl;
  steering_file << m_constraintFileName << std::endl;
  steering_file.close();

  m_constraintFile.open(m_constraintFileName);
  if (m_useEventVertex)
  {
    for (int i : AlignmentDefs::glbl_vtx_label)
    {
      m_constraintFile << " Constraint   0.0" << std::endl;
      m_constraintFile << "       " << i << "   1" << std::endl;
    }
  }
  // print grouping setup to log file:
  std::cout << "MakeMilleFiles::InitRun: Surface groupings are mvtx " << mvtx_group << " intt " << intt_group << " tpc " << tpc_group << " mms " << mms_group << std::endl;

  return ret;
}

//____________________________________________________________________________..
int MakeMilleFiles::process_event(PHCompositeNode* /*topNode*/)
{
  // Outline:
  //
  // loop over track alignment states
  //   Make any track cuts here to skip undesirable tracks (maybe low pT?)
  //   loop over track states+measurements for each track
  //      for each measurement, performed in trackreco/ActsAlignmentStates.cc
  //         Get measurement value and error
  //         Calculate derivatives and residuals from Acts jacobians
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
  {
    std::cout << PHWHERE << " track map size " << _track_map->size() << std::endl;
    std::cout << "state map size " << _state_map->size() << std::endl;
  }

  Acts::Vector3 eventVertex = Acts::Vector3::Zero();
  if (m_useEventVertex)
  {
    eventVertex = getEventVertex();
  }

  ActsPropagator propagator(_tGeometry);

  for (auto [key, statevec] : *_state_map)
  {
    // Check if track was removed from cleaner
    auto iter = _track_map->find(key);
    if (iter == _track_map->end())
    {
      continue;
    }

    SvtxTrack* track = iter->second;

    if (Verbosity() > 0)
    {
      std::cout << std::endl
                << __LINE__ << ": Processing track itrack: " << key << ": nhits: " << track->size_cluster_keys()
                << ": Total tracks: " << _track_map->size() << ": phi: " << track->get_phi() << std::endl;
    }

    //! Make any desired track cuts here
    //! Maybe set a lower pT limit - low pT tracks are not very sensitive to alignment
    addTrackToMilleFile(statevec);

    //! Only take tracks that have 2 mm within event vertex
    if (m_useEventVertex &&
        fabs(track->get_z() - eventVertex.z()) < 0.2 &&
        fabs(track->get_x()) < 0.2 &&
        fabs(track->get_y()) < 0.2)
    {
      //! set x and y to 0 since we are constraining to the x-y origin
      //! and add constraints to pede later
      eventVertex(0) = 0;
      eventVertex(1) = 0;

      auto dcapair = TrackAnalysisUtils::get_dca(track, eventVertex);
      Acts::Vector2 vtx_residual(-dcapair.first.first, -dcapair.second.first);
      vtx_residual *= Acts::UnitConstants::cm;

      float lclvtx_derivative[SvtxAlignmentState::NRES][SvtxAlignmentState::NLOC];
      bool success = getLocalVtxDerivativesXY(track, propagator,
                                              eventVertex, lclvtx_derivative);

      // The global derivs dimensions are [alpha/beta/gamma](x/y/z)
      float glblvtx_derivative[SvtxAlignmentState::NRES][3];
      getGlobalVtxDerivativesXY(track, eventVertex, glblvtx_derivative);

      if (Verbosity() > 2)
      {
        std::cout << "vertex info for trakc " << track->get_id() << " with charge " << track->get_charge() << std::endl;
        std::cout << "vertex is " << eventVertex.transpose() << std::endl;
        std::cout << "vertex residuals " << vtx_residual.transpose()
                  << std::endl;
        std::cout << "global vtx derivatives " << std::endl;
        for (auto& i : glblvtx_derivative)
        {
          for (float j : i)
          {
            std::cout << j << ", ";
          }
          std::cout << std::endl;
        }
      }
      if (success)
      {
        for (int i = 0; i < 2; i++)
        {
          if (!isnan(vtx_residual(i)))
          {
            _mille->mille(SvtxAlignmentState::NLOC, lclvtx_derivative[i],
                          AlignmentDefs::NGLVTX, glblvtx_derivative[i],
                          AlignmentDefs::glbl_vtx_label, vtx_residual(i),
                          m_vtxSigma(i));
          }
        }
      }
    }

    //! Finish this track
    _mille->end();
  }

  if (Verbosity() > 0)
  {
    std::cout << "Finished processing mille file " << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MakeMilleFiles::End(PHCompositeNode* /*unused*/)
{
  delete _mille;
  m_constraintFile.close();

  return Fun4AllReturnCodes::EVENT_OK;
}

int MakeMilleFiles::GetNodes(PHCompositeNode* topNode)
{
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

  _state_map = findNode::getClass<SvtxAlignmentStateMap>(topNode, "SvtxAlignmentStateMap");
  if (!_state_map)
  {
    std::cout << PHWHERE << "Error, can't find alignment state map" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool MakeMilleFiles::getLocalVtxDerivativesXY(SvtxTrack* track,
                                              ActsPropagator& propagator,
                                              const Acts::Vector3& vertex,
                                              float lclvtx_derivative[SvtxAlignmentState::NRES][SvtxAlignmentState::NLOC])
{
  //! Get the first track state beyond the vertex, which will be the
  //! innermost track state and propagate it to the vertex surface to
  //! get the jacobian at the vertex
  SvtxTrackState* firststate = (*std::next(track->begin_states(), 1)).second;

  TrkrDefs::cluskey ckey = firststate->get_cluskey();
  auto cluster = _cluster_map->findCluster(ckey);
  auto surf = _tGeometry->maps().getSurface(ckey, cluster);

  auto param = propagator.makeTrackParams(firststate, track->get_charge(), surf).value();
  auto perigee = propagator.makeVertexSurface(vertex);
  auto actspropagator = propagator.makePropagator();

  Acts::PropagatorOptions<> options(_tGeometry->geometry().getGeoContext(),
                                    _tGeometry->geometry().magFieldContext);

  auto result = actspropagator.propagate(param, *perigee, options);

  if (result.ok())
  {
    auto jacobian = *result.value().transportJacobian;

    Eigen::Matrix<double, 2, 6> projector = Eigen::Matrix<double, 2, 6>::Zero();
    projector(0, 0) = 1;
    projector(1, 1) = 1;
    auto deriv = projector * jacobian;
    if (Verbosity() > 2)
    {
      std::cout << "local vtxderiv " << std::endl
                << deriv << std::endl;
    }

    for (int i = 0; i < deriv.rows(); i++)
    {
      for (int j = 0; j < deriv.cols(); j++)
      {
        lclvtx_derivative[i][j] = deriv(i, j);
      }
    }

    return true;
  }

  return false;
}

void MakeMilleFiles::getGlobalVtxDerivativesXY(SvtxTrack* track,
                                               const Acts::Vector3& vertex,
                                               float glblvtx_derivative[SvtxAlignmentState::NRES][3])
{
  Acts::SquareMatrix3 identity = Acts::SquareMatrix3::Identity();

  Acts::Vector3 projx = Acts::Vector3::Zero();
  Acts::Vector3 projy = Acts::Vector3::Zero();
  getProjectionVtxXY(track, vertex, projx, projy);

  glblvtx_derivative[0][0] = identity.col(0).dot(projx);
  glblvtx_derivative[0][1] = identity.col(1).dot(projx);
  glblvtx_derivative[0][2] = identity.col(2).dot(projx);
  glblvtx_derivative[1][0] = identity.col(0).dot(projy);
  glblvtx_derivative[1][1] = identity.col(1).dot(projy);
  glblvtx_derivative[1][2] = identity.col(2).dot(projy);

  /// swap sign to fit expectation of pede to have derivative of fit rather than
  /// residual
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      glblvtx_derivative[i][j] *= -1;
    }
  }
}
void MakeMilleFiles::getProjectionVtxXY(SvtxTrack* track,
                                        const Acts::Vector3& vertex,
                                        Acts::Vector3& projx,
                                        Acts::Vector3& projy)
{
  Acts::Vector3 tangent(track->get_px(), track->get_py(), track->get_pz());
  Acts::Vector3 normal(track->get_px(), track->get_py(), 0);

  tangent /= tangent.norm();
  normal /= normal.norm();

  Acts::Vector3 localx(1, 0, 0);
  Acts::Vector3 localz(0, 0, 1);

  Acts::Vector3 xglob = localToGlobalVertex(track, vertex, localx);
  Acts::Vector3 yglob = localz + vertex;

  Acts::Vector3 X = (xglob - vertex) / (xglob - vertex).norm();
  Acts::Vector3 Y = (yglob - vertex) / (yglob - vertex).norm();

  // see equation 31 of the ATLAS paper (and discussion) for this
  projx = X - (tangent.dot(X) / tangent.dot(normal)) * normal;
  projy = Y - (tangent.dot(Y) / tangent.dot(normal)) * normal;

  return;
}
Acts::Vector3 MakeMilleFiles::localToGlobalVertex(SvtxTrack* track,
                                                  const Acts::Vector3& vertex,
                                                  const Acts::Vector3& localx) const
{
  Acts::Vector3 mom(track->get_px(),
                    track->get_py(),
                    track->get_pz());

  Acts::Vector3 r = mom.cross(Acts::Vector3(0., 0., 1.));
  float phi = atan2(r(1), r(0));
  Acts::RotationMatrix3 rot;
  Acts::RotationMatrix3 rot_T;

  rot(0, 0) = cos(phi);
  rot(0, 1) = -sin(phi);
  rot(0, 2) = 0;
  rot(1, 0) = sin(phi);
  rot(1, 1) = cos(phi);
  rot(1, 2) = 0;
  rot(2, 0) = 0;
  rot(2, 1) = 0;
  rot(2, 2) = 1;

  rot_T = rot.transpose();

  Acts::Vector3 pos_R = rot * localx;

  pos_R += vertex;

  return pos_R;
}

Acts::Vector3 MakeMilleFiles::getEventVertex()
{
  /**
   * Returns event vertex in cm as averaged track positions
   */
  float xsum = 0;
  float ysum = 0;
  float zsum = 0;
  int nacceptedtracks = 0;

  for (auto [key, statevec] : *_state_map)
  {
    // Check if track was removed from cleaner
    auto iter = _track_map->find(key);
    if (iter == _track_map->end())
    {
      continue;
    }

    SvtxTrack* track = iter->second;

    /// The track vertex is given by the fit as the PCA to the beamline
    xsum += track->get_x();
    ysum += track->get_y();
    zsum += track->get_z();

    nacceptedtracks++;
  }

  return Acts::Vector3(xsum / nacceptedtracks,
                       ysum / nacceptedtracks,
                       zsum / nacceptedtracks);
}

void MakeMilleFiles::addTrackToMilleFile(SvtxAlignmentStateMap::StateVec& statevec)
{
  for (auto state : statevec)
  {
    TrkrDefs::cluskey ckey = state->get_cluster_key();

    if (Verbosity() > 2)
    {
      std::cout << "adding state for ckey " << ckey << " with hitsetkey "
                << (int) TrkrDefs::getHitSetKeyFromClusKey(ckey) << std::endl;
    }
    // The global alignment parameters are given initial values of zero by default, we do not specify them
    // We identify the global alignment parameters for this surface

    TrkrCluster* cluster = _cluster_map->findCluster(ckey);
    const unsigned int layer = TrkrDefs::getLayer(ckey);
    const unsigned int trkrid = TrkrDefs::getTrkrId(ckey);
    const SvtxAlignmentState::ResidualVector residual = state->get_residual();
    const Acts::Vector3 global = _tGeometry->getGlobalPosition(ckey, cluster);

    // need standard deviation of measurements
    SvtxAlignmentState::ResidualVector clus_sigma = SvtxAlignmentState::ResidualVector::Zero();

    double clusRadius = sqrt(global[0] * global[0] + global[1] * global[1]);
    auto para_errors = _ClusErrPara.get_clusterv5_modified_error(cluster, clusRadius, ckey);
    double phierror = sqrt(para_errors.first);
    double zerror = sqrt(para_errors.second);
    clus_sigma(1) = zerror * Acts::UnitConstants::cm;
    clus_sigma(0) = phierror * Acts::UnitConstants::cm;

    if (std::isnan(clus_sigma(0)) ||
        std::isnan(clus_sigma(1)))
    {
      continue;
    }

    auto surf = _tGeometry->maps().getSurface(ckey, cluster);

    int glbl_label[SvtxAlignmentState::NGL];
    if (layer < 3)
    {
      AlignmentDefs::getMvtxGlobalLabels(surf, glbl_label, mvtx_group);
    }
    else if (layer > 2 && layer < 7)
    {
      AlignmentDefs::getInttGlobalLabels(surf, glbl_label, intt_group);
    }
    else if (layer < 55)
    {
      AlignmentDefs::getTpcGlobalLabels(surf, ckey, glbl_label, tpc_group);
    }
    else if (layer < 57)
    {
      AlignmentDefs::getMMGlobalLabels(surf, glbl_label, mms_group);
    }

    if (Verbosity() > 1)
    {
      std::cout << std::endl;
    }

    /// For N residual local coordinates x, z
    for (int i = 0; i < SvtxAlignmentState::NRES; ++i)
    {
      // Add the measurement separately for each coordinate direction to Mille
      float glbl_derivative[SvtxAlignmentState::NGL];
      for (int j = 0; j < SvtxAlignmentState::NGL; ++j)
      {
        glbl_derivative[j] = state->get_global_derivative_matrix()(i, j);

        if (is_layer_fixed(layer) ||
            is_layer_param_fixed(layer, j, fixed_layer_gparams))
        {
          glbl_derivative[j] = 0.;
        }

        if (trkrid == TrkrDefs::mvtxId)
        {
          // need stave to get clamshell
          auto stave = MvtxDefs::getStaveId(ckey);
          auto clamshell = AlignmentDefs::getMvtxClamshell(layer, stave);
          if (is_layer_param_fixed(layer, j, fixed_layer_gparams) ||
              is_mvtx_layer_fixed(layer, clamshell))
          {
            glbl_derivative[j] = 0;
          }
        }
        else if (trkrid == TrkrDefs::inttId)
        {
          if (is_layer_param_fixed(layer, j, fixed_layer_gparams) ||
              is_layer_fixed(layer))
          {
            glbl_derivative[j] = 0;
          }
        }
        if (trkrid == TrkrDefs::tpcId)
        {
          auto sector = TpcDefs::getSectorId(ckey);
          auto side = TpcDefs::getSide(ckey);
          if (is_layer_param_fixed(layer, j, fixed_layer_gparams) ||
              is_tpc_sector_fixed(layer, sector, side) ||
              is_layer_fixed(layer))
          {
            glbl_derivative[j] = 0.0;
          }
        }
      }

      float lcl_derivative[SvtxAlignmentState::NLOC];
      for (int j = 0; j < SvtxAlignmentState::NLOC; ++j)
      {
        lcl_derivative[j] = state->get_local_derivative_matrix()(i, j);

        if (is_layer_param_fixed(layer, j, fixed_layer_lparams))
        {
          lcl_derivative[j] = 0.;
        }
      }
      if (Verbosity() > 2)
      {
        std::cout << "coordinate " << i << " has residual " << residual(i) << " and clus_sigma " << clus_sigma(i) << std::endl
                  << "global deriv " << std::endl;

        for (float k : glbl_derivative)
        {
          if (k > 0 || k < 0)
          {
            std::cout << "NONZERO GLOBAL DERIVATIVE" << std::endl;
          }
          std::cout << k << ", ";
        }
        std::cout << std::endl
                  << "local deriv " << std::endl;
        for (float k : lcl_derivative)
        {
          std::cout << k << ", ";
        }
        std::cout << std::endl;
      }

      if (clus_sigma(i) < 1.0)  // discards crazy clusters
      {
        if (Verbosity() > 3)
        {
          std::cout << "ckey " << ckey << " and layer " << layer << " buffers:" << std::endl;
          AlignmentDefs::printBuffers(i, residual, clus_sigma, lcl_derivative, glbl_derivative, glbl_label);
        }
        float errinf = 1.0;
        if (m_layerMisalignment.find(layer) != m_layerMisalignment.end())
        {
          errinf = m_layerMisalignment.find(layer)->second;
        }

        _mille->mille(SvtxAlignmentState::NLOC, lcl_derivative, SvtxAlignmentState::NGL, glbl_derivative, glbl_label, residual(i), errinf * clus_sigma(i));
      }
    }
  }

  return;
}

bool MakeMilleFiles::is_layer_fixed(unsigned int layer)
{
  bool ret = false;
  auto it = fixed_layers.find(layer);
  if (it != fixed_layers.end())
  {
    ret = true;
  }

  return ret;
}
void MakeMilleFiles::set_layers_fixed(unsigned int minlayer,
                                      unsigned int maxlayer)
{
  for (unsigned int i = minlayer; i < maxlayer; i++)
  {
    fixed_layers.insert(i);
  }
}
void MakeMilleFiles::set_layer_fixed(unsigned int layer)
{
  fixed_layers.insert(layer);
}

bool MakeMilleFiles::is_mvtx_layer_fixed(unsigned int layer, unsigned int clamshell)
{
  bool ret = false;

  std::pair pair = std::make_pair(layer, clamshell);
  auto it = fixed_mvtx_layers.find(pair);
  if (it != fixed_mvtx_layers.end())
  {
    ret = true;
  }

  return ret;
}
void MakeMilleFiles::set_mvtx_layer_fixed(unsigned int layer, unsigned int clamshell)
{
  fixed_mvtx_layers.insert(std::make_pair(layer, clamshell));
}
bool MakeMilleFiles::is_layer_param_fixed(unsigned int layer, unsigned int param, std::set<std::pair<unsigned int, unsigned int>>& param_fixed)
{
  bool ret = false;
  std::pair<unsigned int, unsigned int> pair = std::make_pair(layer, param);
  auto it = param_fixed.find(pair);
  if (it != param_fixed.end())
  {
    ret = true;
  }

  return ret;
}

void MakeMilleFiles::set_layer_gparam_fixed(unsigned int layer, unsigned int param)
{
  fixed_layer_gparams.insert(std::make_pair(layer, param));
}
void MakeMilleFiles::set_layer_lparam_fixed(unsigned int layer, unsigned int param)
{
  fixed_layer_lparams.insert(std::make_pair(layer, param));
}
bool MakeMilleFiles::is_tpc_sector_fixed(unsigned int layer, unsigned int sector, unsigned int side)
{
  bool ret = false;
  unsigned int region = AlignmentDefs::getTpcRegion(layer);
  unsigned int subsector = region * 24 + side * 12 + sector;
  auto it = fixed_sectors.find(subsector);
  if (it != fixed_sectors.end())
  {
    ret = true;
  }

  return ret;
}
