#include "MakeMilleFiles.h"

#include "Mille.h"

/// Tracking includes

#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxAlignmentState.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <trackbase/TpcDefs.h>  // for side

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
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  // Instantiate Mille and open output data file
  //  _mille = new Mille(data_outfilename.c_str(), false);   // write text in data files, rather than binary, for debugging only
  _mille = new Mille(data_outfilename.c_str(), _binary);

  // Write the steering file here, and add the data file path to it
  std::ofstream steering_file(steering_outfilename);
  steering_file << data_outfilename << std::endl;
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

    // Make any desired track cuts here
    // Maybe set a lower pT limit - low pT tracks are not very sensitive to alignment

    addTrackToMilleFile(statevec);

    /// Finish this track
    _mille->end();
  }

  if (Verbosity() > 0)
  {
    std::cout << "Finished processing mille file " << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MakeMilleFiles::End(PHCompositeNode*)
{
  delete _mille;

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

void MakeMilleFiles::addTrackToMilleFile(SvtxAlignmentStateMap::StateVec statevec)
{
  for (auto state : statevec)
  {
    TrkrDefs::cluskey ckey = state->get_cluster_key();

    if (Verbosity() > 2)
    {
      std::cout << "adding state for ckey " << ckey << std::endl;
    }
    // The global alignment parameters are given initial values of zero by default, we do not specify them
    // We identify the global alignment parameters for this surface

    TrkrCluster* cluster = _cluster_map->findCluster(ckey);
    const unsigned int layer = TrkrDefs::getLayer(ckey);

    const SvtxAlignmentState::ResidualVector residual = state->get_residual();
    const Acts::Vector3 global = _tGeometry->getGlobalPosition(ckey, cluster);

    // need standard deviation of measurements
    SvtxAlignmentState::ResidualVector clus_sigma = SvtxAlignmentState::ResidualVector::Zero();

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
    else if (_cluster_version == 5)
    {
      double clusRadius = sqrt(global[0] * global[0] + global[1] * global[1]);
      TrkrClusterv5* clusterv5 = dynamic_cast<TrkrClusterv5*>(cluster);
      auto para_errors = _ClusErrPara.get_clusterv5_modified_error(clusterv5, clusRadius, ckey);
      double phierror = sqrt(para_errors.first);
      double zerror = sqrt(para_errors.second);
      clus_sigma(1) = zerror * Acts::UnitConstants::cm;
      clus_sigma(0) = phierror * Acts::UnitConstants::cm;
    }

    if (std::isnan(clus_sigma(0)) ||
        std::isnan(clus_sigma(1)))
    {
      continue;
    }

    auto surf = _tGeometry->maps().getSurface(ckey, cluster);

    int glbl_label[SvtxAlignmentState::NGL];
    if (layer < 7)
    {
      AlignmentDefs::getSiliconGlobalLabels(surf, glbl_label, si_group);
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

        if (is_layer_fixed(layer) || is_layer_param_fixed(layer, j))
        {
          glbl_derivative[j] = 0.0;
        }
	if(TrkrDefs::getTrkrId(ckey) == TrkrDefs::tpcId)
	  {
	    auto sector = TpcDefs::getSectorId(ckey);
	    auto side = TpcDefs::getSide(ckey);
	    if(is_tpc_sector_fixed(layer, sector, side))
	      {
		glbl_derivative[j] = 0.0;
	      }
	  }
      }

      float lcl_derivative[SvtxAlignmentState::NLOC];
      for (int j = 0; j < SvtxAlignmentState::NLOC; ++j)
      {
        lcl_derivative[j] = state->get_local_derivative_matrix()(i, j);
      }
      if (Verbosity() > 2)
      {
        std::cout << "coordinate " << i << " has residual " << residual(i) << " and clus_sigma " << clus_sigma(i) << std::endl
                  << "global deriv " << std::endl;

        for (int k = 0; k < SvtxAlignmentState::NGL; k++)
        {
          if (glbl_derivative[k] > 0 || glbl_derivative[k] < 0)
            std::cout << "NONZERO GLOBAL DERIVATIVE" << std::endl;
          std::cout << glbl_derivative[k] << ", ";
        }
        std::cout << std::endl
                  << "local deriv " << std::endl;
        for (int k = 0; k < SvtxAlignmentState::NLOC; k++)
        {
          std::cout << lcl_derivative[k] << ", ";
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
    ret = true;

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

bool MakeMilleFiles::is_layer_param_fixed(unsigned int layer, unsigned int param)
{
  bool ret = false;
  std::pair<unsigned int, unsigned int> pair = std::make_pair(layer, param);
  auto it = fixed_layer_params.find(pair);
  if (it != fixed_layer_params.end())
    ret = true;

  return ret;
}

void MakeMilleFiles::set_layer_param_fixed(unsigned int layer, unsigned int param)
{
  std::pair<unsigned int, unsigned int> pair = std::make_pair(layer, param);
  fixed_layer_params.insert(pair);
}
bool MakeMilleFiles::is_tpc_sector_fixed(unsigned int layer, unsigned int sector, unsigned int side)
 {
   bool ret = false;
   unsigned int region = AlignmentDefs::getTpcRegion(layer);
   unsigned int subsector = region * 24 + side * 12 + sector;
   auto it = fixed_sectors.find(subsector);
   if(it != fixed_sectors.end()) 
     ret = true;

   return ret;
 }
