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
    {
      std::cout << PHWHERE << " track map size " << _track_map->size() << std::endl;
      std::cout << "state map size " << _state_map->size() << std::endl;
    }

  for (auto [key, statevec] : *_state_map)
  {
    
    // Check if track was removed from cleaner
    auto iter = _track_map->find(key);
    if(iter == _track_map->end())
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

  if(Verbosity() > 0 )
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

int MakeMilleFiles::getTpcRegion(int layer)
{
  int region = 0;
  if (layer > 23 && layer < 39)
    region = 1;
  if (layer > 38 && layer < 55)
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
    if (si_group == siliconGroup::sensor)
    {
      // every sensor has a different label
      int stave = sensor / nsensors_stave[layer];
      label_base += layer * 1000000 + stave * 10000 + sensor * 10;
      return label_base;
    }
    if (si_group == siliconGroup::stave)
    {
      // layer and stave, assign all sensors to the stave number
      int stave = sensor / nsensors_stave[layer];
      label_base += layer * 1000000 + stave * 10000;
      return label_base;
    }
    if (si_group == siliconGroup::barrel)
    {
      // layer only, assign all sensors to sensor 0
      label_base += layer * 1000000 + 0;

      return label_base;
    }
  }
  else if (layer > 6 && layer < 55)
  {
    if (tpc_group == tpcGroup::hitset)
    {
      // want every hitset (layer, sector, side) to have a separate label
      // each group of 12 subsurfaces (sensors) is in a single hitset
      int hitset = sensor / 12;  // hitsets 0-11 on side 0, 12-23 on side 1
      label_base += layer * 1000000 + hitset * 10000;
      return label_base;
    }
    if (tpc_group == tpcGroup::sector)
    {
      // group all tpc layers in each region and sector, assign layer 7 and side and sector number to all layers and hitsets
      int side = sensor / 144;  // 0-143 on side 0, 144-287 on side 1
      int sector = (sensor - side * 144) / 12;
      // for a given layer there are only 12 sectors x 2 sides
      // The following gives the sectors in the inner, mid, outer regions unique group labels
      int region = getTpcRegion(layer);  // inner, mid, outer
      label_base += 7 * 1000000 + (region * 24 + side * 12 + sector) * 10000;
      // std::cout << " layer " << layer << " sensor " << sensor << " region " << region << " side " << side << " sector " << sector << " label_base " << label_base << std::endl;
      return label_base;
    }
    if (tpc_group == tpcGroup::tpc)
    {
      // all tpc layers and all sectors, assign layer 7 and sensor 0 to all layers and sensors
      label_base += 7 * 1000000 + 0;
      return label_base;
    }
  }
  else
  {
    if (mms_group == mmsGroup::tile)
    {
      // every tile has different label
      int tile = sensor;
      label_base += layer * 1000000 + tile * 10000 + sensor * 10;
      return label_base;
    }
    if (mms_group == mmsGroup::mms)
    {
      // assign layer 55 and tile 0 to all
      label_base += 55 * 1000000 + 0;
      return label_base;
    }
  }

  return -1;
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

    const auto cluster = _cluster_map->findCluster(ckey);
    const auto layer = TrkrDefs::getLayer(ckey);

    const auto residual = state->get_residual();
    const auto& global = _tGeometry->getGlobalPosition(ckey, cluster);

    // need standard deviation of measurements
    SvtxAlignmentState::ResidualVector clus_sigma = SvtxAlignmentState::ResidualVector::Zero();
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

    Acts::GeometryIdentifier id = _tGeometry->maps().getSurface(ckey, cluster)->geometryId();
    int label_base = getLabelBase(id);  // This value depends on how the surfaces are grouped

    int glbl_label[SvtxAlignmentState::NGL];
    for (int i = 0; i < SvtxAlignmentState::NGL; ++i)
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

    /// For N residual coordinates x,y,z
    for (int i = 0; i < SvtxAlignmentState::NRES; ++i)
    {
      // Add the measurement separately for each coordinate direction to Mille
      float glbl_derivative[SvtxAlignmentState::NGL];
      for (int j = 0; j < SvtxAlignmentState::NGL; ++j)
      {
	/// swap the order to match what is expected from the workflow
        glbl_derivative[j] = state->get_global_derivative_matrix()(i, (j+3)%SvtxAlignmentState::NGL);

        if (is_layer_fixed(layer) || is_layer_param_fixed(layer, j))
        {
          glbl_derivative[j] = 0.0;
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
	 if(Verbosity() > 3)
	    { 
	      std::cout << "ckey " << ckey << " and layer " << layer << " buffers:" << std::endl; 
	      printBuffers(i, residual, clus_sigma, lcl_derivative, glbl_derivative, glbl_label); 
	    }

        _mille->mille(SvtxAlignmentState::NLOC, lcl_derivative, SvtxAlignmentState::NGL, glbl_derivative, glbl_label, residual(i), clus_sigma(i));
      }
    }
  }

  return;
}

void MakeMilleFiles::printBuffers(int index, Acts::Vector3 residual, Acts::Vector3 clus_sigma, float lcl_derivative[], float glbl_derivative[], int glbl_label[])
{
  std::cout << " float buffer: " << " residual " << "  " << residual(index);
  for (int il=0;il<SvtxAlignmentState::NLOC;++il) { if(lcl_derivative[il] != 0) std::cout << " lcl_deriv["<< il << "] " << lcl_derivative[il] << "  ";  }
  std::cout  << " sigma " << "  " << clus_sigma(index) << "  ";
  for (int ig=0;ig<SvtxAlignmentState::NGL;++ig) { if(glbl_derivative[ig] != 0)  std::cout << " glbl_deriv["<< ig << "] " << glbl_derivative[ig] << "  ";  }
  std::cout << " int buffer: " << " 0 " << "  ";
  for (int il=0;il<SvtxAlignmentState::NLOC;++il) { if(lcl_derivative[il] != 0) std::cout << " lcl_label["<< il << "] " << il << "  ";  }
  std::cout << " 0 " << "  ";
  for (int ig=0;ig<SvtxAlignmentState::NGL;++ig) { if(glbl_derivative[ig] != 0) std::cout << " glbl_label["<< ig << "] " << glbl_label[ig] << "  ";  }
  std::cout << " end of meas " << std::endl;		    
}
bool MakeMilleFiles::is_layer_fixed(unsigned int layer)
{
  bool ret = false;
  auto it = fixed_layers.find(layer);
  if (it != fixed_layers.end())
    ret = true;

  return ret;
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
