#include "KFParticle_caloInfo.h"
#include "KFParticle_util.h"

#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::...
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap, SvtxTrac...

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <phool/getClass.h>  // for getClass

#include <KFParticle.h>  // for KFParticle

#include <TSystem.h>  // for gSystem->Exit()
#include <TTree.h>    // for TTree

#include <algorithm>  // for max, find
#include <cmath>      // for pow, sqrt
#include <cstdlib>    // for NULL, abs
#include <iostream>   // for operator<<, endl, basi...
#include <iterator>   // for end, begin
#include <map>        // for _Rb_tree_iterator, map
#include <memory>     // for allocator_traits<>::va...
#include <utility>    // for pair

void KFParticle_caloInfo::initializeBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number)
{
  m_tree->Branch((daughter_number + "_EMCAL_DeltaPhi").c_str(), &detector_emcal_deltaphi[daughter_id], (daughter_number + "_EMCAL_DeltaPhi/F").c_str());
  m_tree->Branch((daughter_number + "_EMCAL_DeltaEta").c_str(), &detector_emcal_deltaeta[daughter_id], (daughter_number + "_EMCAL_DeltaEta/F").c_str());
  m_tree->Branch((daughter_number + "_EMCAL_DeltaZ").c_str(), &detector_emcal_deltaz[daughter_id], (daughter_number + "_EMCAL_DeltaZ/F").c_str());
  m_tree->Branch((daughter_number + "_EMCAL_energy_3x3").c_str(), &detector_emcal_energy_3x3[daughter_id], (daughter_number + "_EMCAL_energy_3x3/F").c_str());
  m_tree->Branch((daughter_number + "_EMCAL_energy_5x5").c_str(), &detector_emcal_energy_5x5[daughter_id], (daughter_number + "_EMCAL_energy_5x5/F").c_str());
  m_tree->Branch((daughter_number + "_EMCAL_energy_cluster").c_str(), &detector_emcal_cluster_energy[daughter_id], (daughter_number + "_EMCAL_energy_cluster/F").c_str());
  m_tree->Branch((daughter_number + "_EMCAL_cluster_prob").c_str(), &detector_emcal_cluster_prob[daughter_id], (daughter_number + "_EMCAL_cluster_prob/F").c_str());

  // Detailed calo info
  if (m_config.detailedCalorimetry)
  {
    m_tree->Branch((daughter_number + "_EMCAL_projection_z").c_str(),   &detector_emcal_projection_z[daughter_id], (daughter_number + "_EMCAL_projection_z/F").c_str());
    m_tree->Branch((daughter_number + "_EMCAL_projection_eta").c_str(), &detector_emcal_projection_eta[daughter_id], (daughter_number + "_EMCAL_projection_eta/F").c_str());
    m_tree->Branch((daughter_number + "_EMCAL_projection_phi").c_str(), &detector_emcal_projection_phi[daughter_id], (daughter_number + "_EMCAL_projection_phi/F").c_str());
    m_tree->Branch((daughter_number + "_EMCAL_cluster_z").c_str(),   &detector_emcal_cluster_z[daughter_id], (daughter_number + "_EMCAL_cluster_z/F").c_str());
    m_tree->Branch((daughter_number + "_EMCAL_cluster_eta").c_str(), &detector_emcal_cluster_eta[daughter_id], (daughter_number + "_EMCAL_cluster_eta/F").c_str());
    m_tree->Branch((daughter_number + "_EMCAL_cluster_phi").c_str(), &detector_emcal_cluster_phi[daughter_id], (daughter_number + "_EMCAL_cluster_phi/F").c_str());
    m_tree->Branch((daughter_number + "_EMCAL_Chi2").c_str(), &detector_emcal_chi2[daughter_id], (daughter_number + "_EMCAL_Chi2/F").c_str());
    m_tree->Branch((daughter_number + "_EMCAL_nTowers").c_str(), &detector_emcal_ntowers[daughter_id], (daughter_number + "_EMCAL_nTowers/i").c_str());
  }

  // 5x5 cell info
  if (m_config.calo5x5CellInfo)
  {
    m_tree->Branch((daughter_number + "_EMCAL_5x5Cell_binPhi").c_str(), &detector_emcal_5x5Cell_Phi[daughter_id]);
    m_tree->Branch((daughter_number + "_EMCAL_5x5Cell_binEta").c_str(), &detector_emcal_5x5Cell_Eta[daughter_id]);
    m_tree->Branch((daughter_number + "_EMCAL_5x5Cell_E").c_str(), &detector_emcal_5x5Cell_E[daughter_id]);
  }

  // Inner and outer HCal -- not yet functional; fills with NAN by default until software is developed
  m_tree->Branch((daughter_number + "_IHCAL_DeltaPhi").c_str(), &detector_ihcal_deltaphi[daughter_id], (daughter_number + "_IHCAL_DeltaPhi/F").c_str());
  m_tree->Branch((daughter_number + "_IHCAL_DeltaEta").c_str(), &detector_ihcal_deltaeta[daughter_id], (daughter_number + "_IHCAL_DeltaEta/F").c_str());
  m_tree->Branch((daughter_number + "_IHCAL_energy_3x3").c_str(), &detector_ihcal_energy_3x3[daughter_id], (daughter_number + "_IHCAL_energy_3x3/F").c_str());
  m_tree->Branch((daughter_number + "_IHCAL_energy_5x5").c_str(), &detector_ihcal_energy_5x5[daughter_id], (daughter_number + "_IHCAL_energy_5x5/F").c_str());
  m_tree->Branch((daughter_number + "_IHCAL_energy_cluster").c_str(), &detector_ihcal_cluster_energy[daughter_id], (daughter_number + "_IHCAL_energy_cluster/F").c_str());
  m_tree->Branch((daughter_number + "_OHCAL_DeltaPhi").c_str(), &detector_ohcal_deltaphi[daughter_id], (daughter_number + "_OHCAL_DeltaEta/F").c_str());
  m_tree->Branch((daughter_number + "_OHCAL_DeltaEta").c_str(), &detector_ohcal_deltaeta[daughter_id], (daughter_number + "_OHCAL_DeltaEta/F").c_str());
  m_tree->Branch((daughter_number + "_OHCAL_energy_3x3").c_str(), &detector_ohcal_energy_3x3[daughter_id], (daughter_number + "_OHCAL_energy_3x3/F").c_str());
  m_tree->Branch((daughter_number + "_OHCAL_energy_5x5").c_str(), &detector_ohcal_energy_5x5[daughter_id], (daughter_number + "_OHCAL_energy_5x5/F").c_str());
  m_tree->Branch((daughter_number + "_OHCAL_energy_cluster").c_str(), &detector_ohcal_cluster_energy[daughter_id], (daughter_number + "_OHCAL_energy_cluster/F").c_str());
}

/*
The following function matches tracks to calo clusters. As of 8/18/2025, this only extends to the EMCal. HCal matching is in development.
*/
void KFParticle_caloInfo::fill(PHCompositeNode *topNode, const KFParticle &daughter, int daughter_id, bool &isTrackEMCalmatch,
                               const KFParticle &vertex_in)
{

  // Get track map
  auto *dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_config.trackMapNodeName);
  if (!dst_trackmap)
  {
    std::cout << "KFParticle truth matching: " << m_config.trackMapNodeName << " does not exist" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  // Get track
  SvtxTrack *track = KFParticle_util::getTrack(daughter.Id(), dst_trackmap);

  // EMCal info
  if (!clustersEM)
  {
    clustersEM = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
    if (!clustersEM)
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "CLUSTERINFO_CEMC" << std::endl;
    }
  }
  if (!EMCalGeo)
  {
    EMCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if (!EMCalGeo)
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWERGEOM_CEMC" << std::endl;
    }
  }
  if (!_towersEM)
  {
    _towersEM = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
    if (!_towersEM)
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWER_CALIB_CEMC" << std::endl;
    }
  }
  if (!clustersIH)
  {
    clustersIH = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALIN");
  }
  if (!IHCalGeo)
  {
    IHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  }
  if (!clustersOH)
  {
    clustersOH = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALOUT");
  }
  if (!OHCalGeo)
  {
    OHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  }

  // Radii for track projections
  const double caloRadiusEMCal = m_config.emcalRadius;  // Use function set_emcal_radius_user(float set_variable) to set this in your Fun4All macro

  // Create variables and containers etc.
  bool is_match;
  int index;
  RawCluster *cluster;

  // EMCAL******************************************************
  // Get state
  SvtxTrackState *thisState = track->get_state(caloRadiusEMCal);

  // Print out additional track states with high verbosity:
  if(m_config.track2CaloVerbosity > 10 && thisState != nullptr)
  {
    thisState->identify();
    std::cout << "size states " << track->size_states() << std::endl;
    for (auto state_iter = track->begin_states();
          state_iter != track->end_states();
          ++state_iter)
      {
        SvtxTrackState* tstate = state_iter->second;
        tstate->identify();
        std::cout<< "radius " << std::sqrt((tstate->get_x())*(tstate->get_x()) + (tstate->get_y())*(tstate->get_y())) << std::endl;
      }
    assert(thisState);
    // Useful for debugging:
    // clustersEM->identify();
  }

  // Declare variables:
  // Track
  float _track_phi_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_eta_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_x_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_y_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_z_emc = std::numeric_limits<float>::quiet_NaN();

  // Projected track state info
  float projection_x_emc = std::numeric_limits<float>::quiet_NaN();
  float projection_y_emc = std::numeric_limits<float>::quiet_NaN();
  float projection_z_emc = std::numeric_limits<float>::quiet_NaN();
  float projection_phi_emc = std::numeric_limits<float>::quiet_NaN();
  float projection_eta_emc = std::numeric_limits<float>::quiet_NaN();

  // Basic EMCal vectors
  std::vector<float> v_emcal_clusE;
  std::vector<float> v_emcal_dphi;
  std::vector<float> v_emcal_deta;
  std::vector<float> v_emcal_dr;
  std::vector<float> v_emcal_dz;

  // cluster prob (from RawClusterV1, cluster template probability for EM shower
  /* as per John H. can be used as a cut to improve electron/photon signal purity in EMCAL */
  std::vector<float> v_emcal_cluster_prob;

  // Detailed Calo Info vectors
  std::vector<float> v_projection_z_emc;
  std::vector<float> v_projection_phi_emc;
  std::vector<float> v_projection_eta_emc;
  std::vector<float> v_emcal_cluster_z;
  std::vector<float> v_emcal_cluster_eta;
  std::vector<float> v_emcal_cluster_phi;
  std::vector<float> v_emcal_nTowers;
  std::vector<float> v_emcal_chi2;

  // 5x5 cell info
  RawClusterDefs::keytype _clus_key;
  std::vector<RawClusterDefs::keytype> v_clus_key;
  std::vector<unsigned int> v_detector_emcal_5x5Cell_Phi;
  std::vector<unsigned int> v_detector_emcal_5x5Cell_Eta;
  std::vector<float> v_detector_emcal_5x5Cell_E;


  // Set variables for matching
  is_match = false;
  index = -1;

  // If the state exists, start matching process
  if (thisState != nullptr){

    // Vertex x,y,z for forrecting position when matching
    float vx = vertex_in.GetX();
    float vy = vertex_in.GetY();
    float vz = vertex_in.GetZ();

    // Track state kinematics
    _track_x_emc = thisState->get_x() - vx;
    _track_y_emc = thisState->get_y() - vy;
    _track_z_emc = thisState->get_z() - vz;
    _track_phi_emc = std::atan2(_track_y_emc, _track_x_emc);
    _track_eta_emc = std::asinh(_track_z_emc / std::sqrt((_track_x_emc * _track_x_emc) + (_track_y_emc * _track_y_emc)));

    // Same quantities as above but in detector coordinates i.e. wrt (0,0,0);
    // There are to be stored in nTuple -- not used in any calculations
    projection_x_emc = thisState->get_x();
    projection_y_emc = thisState->get_y();
    projection_z_emc = thisState->get_z();
    projection_phi_emc = std::atan2(projection_y_emc, projection_x_emc);
    projection_eta_emc = std::asinh(projection_z_emc / std::sqrt((projection_x_emc * projection_x_emc) + (projection_y_emc * projection_y_emc)));


    // Create objects, containers, iterators for clusters
    cluster = nullptr;
    RawClusterContainer::Range begin_end_EMC = clustersEM->getClusters();
    RawClusterContainer::Iterator clusIter_EMC;

    // Loop over the EMCal clusters
    for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
    {
      // Get cluster with respect to PV, not (0,0,0)
      cluster = clusIter_EMC->second;
      CLHEP::Hep3Vector v_vertex;
      v_vertex.set(vx, vy, vz);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetEVec(*cluster, v_vertex);

      // CUT -- CLUSTER ENERGY
      float cluster_energy = E_vec_cluster.mag();
      if (cluster_energy < m_config.emcalEnergyLowCut)
      {
        continue;
      }

      // CUT -- DPHI
      float clPhi = E_vec_cluster.phi();
      float dphi = PiRange(_track_phi_emc - clPhi);
      if (dphi > m_config.emcalDPhiCutHigh || dphi < m_config.emcalDPhiCutLow)
      {
        continue;
      }

      // CUT -- DZ
      float clZ = cluster->get_z() - vz;
      float dz = _track_z_emc - clZ;
      if (dz > m_config.emcalDZCutHigh || dz < m_config.emcalDZCutLow)
      {
        continue;
      }

      // CUT -- DETA    // Currently not being used, but want to include functionality now
      float clEta = E_vec_cluster.pseudoRapidity();
      float deta = _track_eta_emc - clEta;
      // if (deta > m_deta_cut_high || deta < m_deta_cut_low) continue;


      // Calculate dr
      float tmparg = caloRadiusEMCal * dphi;
      float dr = std::sqrt((tmparg * tmparg) + (dz * dz));  // sqrt((R*dphi)^2 + (dz)^2
      // float dr = sqrt((dphi*dphi + deta*deta)); //previous version

      // Get detailed Calo Info
      unsigned int _emcal_nTowers = cluster->getNTowers();
      float _emcal_chi2 = cluster->get_chi2();
      _clus_key = cluster->get_id();
      float detectorcoordinates_emcal_x = cluster->get_x(); // In detector coordinates!!! Not wrt to the PV. It is wrt (0,0,0)
      float detectorcoordinates_emcal_y = cluster->get_y();
      float detectorcoordinates_emcal_z = cluster->get_z();
      float detectorcoordinates_emcal_phi = std::atan2(detectorcoordinates_emcal_y, detectorcoordinates_emcal_x);
      float arg = (detectorcoordinates_emcal_x * detectorcoordinates_emcal_x) + (detectorcoordinates_emcal_y * detectorcoordinates_emcal_y);
      float detectorcoordinates_emcal_eta = std::asinh(detectorcoordinates_emcal_z / std::sqrt(arg));


      // Add potential match's information to vectors
      // Basic info
      v_emcal_clusE.push_back(cluster_energy);
      v_emcal_dphi.push_back(dphi);
      v_emcal_dz.push_back(dz);
      v_emcal_deta.push_back(deta);
      v_emcal_dr.push_back(dr);
      v_emcal_cluster_prob.push_back( cluster->get_prob() );

      // Detailed Calo Info
      v_emcal_cluster_z.push_back(detectorcoordinates_emcal_z);
      v_emcal_cluster_eta.push_back(detectorcoordinates_emcal_eta);
      v_emcal_cluster_phi.push_back(detectorcoordinates_emcal_phi);
      v_emcal_nTowers.push_back(_emcal_nTowers);
      v_emcal_chi2.push_back(_emcal_chi2);

      // 5x5 call -- key used later in Get5x5CellInfo()
      v_clus_key.push_back(_clus_key);

      is_match = true;

      // Print outs
      if(m_config.track2CaloVerbosity > 5)
      {
        std::cout << "**********DELTA INFORMATION************" << std::endl;
        std::cout << "dphi = " << dphi << std::endl;
        std::cout << "deta = " << deta << std::endl;
        std::cout << "dz =   " << dz <<   std::endl;
        std::cout << "dr =   " << dr <<   std::endl;
        std::cout << "****************************************" << std::endl;
      }
    }

    // Find the closest match from all potential matches
    if (is_match)
    {
      float tmp = 99999;
      for (long unsigned int i = 0; i < v_emcal_dr.size(); i++)
      {
        if (v_emcal_dr[i] < tmp)
        {
          index = i;
          tmp = v_emcal_dr[i];
        }
      }
    }
  }
  // Save values to the branches!
  // If no match:
  if (index == -1)
  {
    // Basic info
    detector_emcal_deltaz[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_deltaeta[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_deltaphi[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_energy_3x3[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_energy_5x5[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_cluster_energy[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_cluster_prob[daughter_id] = std::numeric_limits<float>::quiet_NaN();

    // Detailed Calo
    if (m_config.detailedCalorimetry)
    {
      detector_emcal_projection_z[daughter_id] = std::numeric_limits<float>::quiet_NaN();
      detector_emcal_projection_eta[daughter_id] = std::numeric_limits<float>::quiet_NaN();
      detector_emcal_projection_phi[daughter_id] = std::numeric_limits<float>::quiet_NaN();
      detector_emcal_cluster_z[daughter_id] = std::numeric_limits<float>::quiet_NaN();
      detector_emcal_cluster_eta[daughter_id] = std::numeric_limits<float>::quiet_NaN();
      detector_emcal_cluster_phi[daughter_id] = std::numeric_limits<float>::quiet_NaN();
      detector_emcal_ntowers[daughter_id] = std::numeric_limits<unsigned int>::quiet_NaN();
      detector_emcal_chi2[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    }

    // 5x5 cell
    if (m_config.calo5x5CellInfo)
    {
      // Fill with nothing -- no match found, so vectors should be empty !
    }

    isTrackEMCalmatch = false;
  }

  // If there is a match:
  else
  {
    // Grab 3x3 and 5x5 energies & cell info from around the cluster center using the cluster key
    float energy_3x3 = 0;
    float energy_5x5 = 0;
    Get5x5CellInfo(
      v_clus_key[index],
      energy_3x3,
      energy_5x5,
      v_detector_emcal_5x5Cell_Phi,
      v_detector_emcal_5x5Cell_Eta,
      v_detector_emcal_5x5Cell_E
    );

    // Now fill the arrays
    // Basic info
    detector_emcal_deltaphi[daughter_id] = v_emcal_dphi[index];
    detector_emcal_deltaeta[daughter_id] = v_emcal_deta[index];
    detector_emcal_deltaz[daughter_id] = v_emcal_dz[index];
    detector_emcal_energy_3x3[daughter_id] = energy_3x3;
    detector_emcal_energy_5x5[daughter_id] = energy_5x5;
    detector_emcal_cluster_energy[daughter_id] = v_emcal_clusE[index];
    detector_emcal_cluster_prob[daughter_id] = v_emcal_cluster_prob[index];

    // Detailed Calo
    if (m_config.detailedCalorimetry)
    {
      detector_emcal_projection_z[daughter_id] = projection_z_emc;
      detector_emcal_projection_eta[daughter_id] = projection_eta_emc;
      detector_emcal_projection_phi[daughter_id] = projection_phi_emc;
      detector_emcal_cluster_z[daughter_id] = v_emcal_cluster_z[index];
      detector_emcal_cluster_eta[daughter_id] = v_emcal_cluster_eta[index];
      detector_emcal_cluster_phi[daughter_id] = v_emcal_cluster_phi[index];
      detector_emcal_ntowers[daughter_id] = v_emcal_nTowers[index];
      detector_emcal_chi2[daughter_id] = v_emcal_chi2[index];
    }

    // 5x5 cell
    if(m_config.calo5x5CellInfo)
    {
      detector_emcal_5x5Cell_Phi[daughter_id] = v_detector_emcal_5x5Cell_Phi;
      detector_emcal_5x5Cell_Eta[daughter_id] = v_detector_emcal_5x5Cell_Eta;
      detector_emcal_5x5Cell_E[daughter_id] = v_detector_emcal_5x5Cell_E;
    }


    isTrackEMCalmatch = true;
  }

  // HCAL*******************************************************

  // INNER
  /*
  std::cout << "Starting IHCAL-track matching!" << std::endl;

  // Track projection
  thisState = nullptr;
  thisState = track->get_state(caloRadiusIHCal);
  float _track_phi_ihc = std::numeric_limits<float>::quiet_NaN();
  float _track_eta_ihc = std::numeric_limits<float>::quiet_NaN();
  float _track_x_ihc = std::numeric_limits<float>::quiet_NaN();
  float _track_y_ihc = std::numeric_limits<float>::quiet_NaN();
  float _track_z_ihc = std::numeric_limits<float>::quiet_NaN();

  // IHCal variables and vectors
  float _ihcal_phi = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_eta = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_x = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_y = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_z = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_3x3 = std::numeric_limits<float>::quiet_NaN();
  //float _ihcal_5x5 = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_clusE = std::numeric_limits<float>::quiet_NaN();
  radius_scale = std::numeric_limits<float>::quiet_NaN();
  std::vector<float> v_ihcal_phi;
  std::vector<float> v_ihcal_eta;
  std::vector<float> v_ihcal_x;
  std::vector<float> v_ihcal_y;
  std::vector<float> v_ihcal_z;
  std::vector<float> v_ihcal_dphi;
  std::vector<float> v_ihcal_deta;
  std::vector<float> v_ihcal_dr;
  std::vector<float> v_ihcal_3x3;
  //std::vector<float> v_ihcal_5x5;
  std::vector<float> v_ihcal_clusE;

  if(thisState != nullptr){

  // Reset variables for matching
  is_match = false;
  index = -1;

  _track_phi_ihc = atan2(thisState->get_y(), thisState->get_x());
  _track_eta_ihc = asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y()));
  _track_x_ihc = thisState->get_x();
  _track_y_ihc = thisState->get_y();
  _track_z_ihc = thisState->get_z();


  // Create objects, containers, iterators for clusters
  cluster = nullptr;
  RawClusterContainer::Range begin_end_IHC = clustersIH->getClusters();
  RawClusterContainer::Iterator clusIter_IHC;

  // Loop over the IHCal clusters
  for (clusIter_IHC = begin_end_IHC.first; clusIter_IHC != begin_end_IHC.second; ++clusIter_IHC)
  {

  // Minimum energy cut
  cluster = clusIter_IHC->second;
  if(cluster->get_energy() < m_ihcal_e_low_cut)
  {
  continue;
  }

  std::cout << "Found a cluster about threshhold energy!" << std::endl;

  // Get cluster information
  _ihcal_phi = atan2(cluster->get_y(), cluster->get_x());
  _ihcal_eta = asinh(cluster->get_z()/sqrt(cluster->get_x()*cluster->get_x() + cluster->get_y()*cluster->get_y()));
  _ihcal_x = cluster->get_x();
  _ihcal_y = cluster->get_y();
  radius_scale = m_ihcal_radius_user / sqrt(_ihcal_x*_ihcal_x+_ihcal_y*_ihcal_y);
  _ihcal_z = radius_scale*cluster->get_z();
  // _ihcal_3x3 = get_e3x3(cluster, _towersIH, 0); //1 for ihcal
  // _ihcal_5x5 = get_e5x5(cluster, _towersIH, 0); //1 for ihcal
  _ihcal_3x3 = std::numeric_limits<float>::quiet_NaN();
  _ihcal_5x5 = std::numeric_limits<float>::quiet_NaN();
  _ihcal_clusE = cluster->get_energy();


  // Variables to determine potential matches
  float dphi = abs(PiRange(_track_phi_ihc - _ihcal_phi));
  float dz = abs(_track_z_ihc - _ihcal_z);
  float deta = abs(_ihcal_eta - _track_eta_ihc);
  float dr = sqrt((dphi*dphi + deta*deta));

  // Requirements for a possible match
  if(dphi<m_dphi_cut && dz<m_dz_cut)
  {
  //Add potential match's information to vectors
  v_ihcal_phi.push_back(_ihcal_phi);
  v_ihcal_eta.push_back(_ihcal_eta);
  v_ihcal_x.push_back(_ihcal_x);
  v_ihcal_y.push_back(_ihcal_y);
  v_ihcal_z.push_back(_ihcal_z);
  v_ihcal_dphi.push_back(dphi);
  v_ihcal_deta.push_back(deta);
  v_ihcal_dr.push_back(dr);
  v_ihcal_3x3.push_back(_ihcal_3x3);
  v_ihcal_clusE.push_back(_ihcal_clusE);

  is_match = true;
  }
  }

  // Find the closest match from all potential matches
  if (is_match == true)
  {
  float tmp = 99999;
  for(long unsigned int i = 0; i < v_ihcal_dr.size(); i++){
  if(v_ihcal_dr[i] < tmp){
  index = i;
  tmp = v_ihcal_dr[i];
  }
  }
  }
  }

  // Print out statihents
  if(index != -1){
  std::cout<<"matched tracks!!!"<<std::endl;
  std::cout<<"ihcal x = "<<v_ihcal_x[index]<<" , y = "<<v_ihcal_y[index]<<" , z = "<<v_ihcal_z[index]<<" , phi = "<<v_ihcal_phi[index]<<" , eta = "<<v_ihcal_eta[index]<<std::endl;
  std::cout<<"track projected x = "<<_track_x_ihc<<" , y = "<<_track_y_ihc<<" , z = "<<_track_z_ihc<<" , phi = "<<_track_phi_ihc<<" , eta = "<<_track_eta_ihc<<std::endl;
  std::cout<<"track px = "<<track->get_px()<<" , py = "<<track->get_py()<<" , pz = "<<track->get_pz()<<" , pt = "<<track->get_pt()<<" , p = "<<track->get_p()<<" , charge = "<<track->get_charge()<<std::endl;
  }
  */

  // Save values to the branches!
  // if(index == -1){
  detector_ihcal_deltaphi[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ihcal_deltaeta[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ihcal_energy_3x3[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ihcal_energy_5x5[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ihcal_cluster_energy[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // }
  // else{
  // detector_ihcal_deltaphi[daughter_id] = v_ihcal_dphi[index];
  // detector_ihcal_deltaeta[daughter_id] = v_ihcal_deta[index];
  // detector_ihcal_energy_3x3[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // detector_ihcal_energy_5x5[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // detector_ihcal_cluster_energy[daughter_id] = v_ihcal_clusE[index];
  // }

  // std::cout << "IHCAL CLLUSTERS MATCHED TO TRACK" << std::endl;

  // OUTER

  // std::cout << "Starting OHCAL-track matching!" << std::endl;

  /*

  // Track projection
  thisState = nullptr;
  thisState = track->get_state(caloRadiusOHCal);
  float _track_phi_ohc = std::numeric_limits<float>::quiet_NaN();
  float _track_eta_ohc = std::numeric_limits<float>::quiet_NaN();
  float _track_x_ohc = std::numeric_limits<float>::quiet_NaN();
  float _track_y_ohc = std::numeric_limits<float>::quiet_NaN();
  float _track_z_ohc = std::numeric_limits<float>::quiet_NaN();

  // OHCal variables and vectors
  float _ohcal_phi = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_eta = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_x = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_y = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_z = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_3x3 = std::numeric_limits<float>::quiet_NaN();
  //float _ohcal_5x5 = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_clusE = std::numeric_limits<float>::quiet_NaN();
  radius_scale = std::numeric_limits<float>::quiet_NaN();
  std::vector<float> v_ohcal_phi;
  std::vector<float> v_ohcal_eta;
  std::vector<float> v_ohcal_x;
  std::vector<float> v_ohcal_y;
  std::vector<float> v_ohcal_z;
  std::vector<float> v_ohcal_dphi;
  std::vector<float> v_ohcal_deta;
  std::vector<float> v_ohcal_dr;
  std::vector<float> v_ohcal_3x3;
  //std::vector<float> v_ohcal_5x5;
  std::vector<float> v_ohcal_clusE;

  // Reset variables for matching
  is_match = false;
  index = -1;

  if(thisState != nullptr)
  {
  _track_phi_ohc = atan2(thisState->get_y(), thisState->get_x());
  _track_eta_ohc = asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y()));
  _track_x_ohc = thisState->get_x();
  _track_y_ohc = thisState->get_y();
  _track_z_ohc = thisState->get_z();


  // Create objects, containers, iterators for clusters
  cluster = nullptr;
  RawClusterContainer::Range begin_end_OHC = clustersOH->getClusters();
  RawClusterContainer::Iterator clusIter_OHC;

  // Loop over the OHCal clusters
  for (clusIter_OHC = begin_end_OHC.first; clusIter_OHC != begin_end_OHC.second; ++clusIter_OHC)
  {

  // Minimum energy cut
  cluster = clusIter_OHC->second;
  if(cluster->get_energy() < m_ohcal_e_low_cut)
  {
  continue;
  }

  // Get cluster information
  _ohcal_phi = atan2(cluster->get_y(), cluster->get_x());
  _ohcal_eta = asinh(cluster->get_z()/sqrt(cluster->get_x()*cluster->get_x() + cluster->get_y()*cluster->get_y()));
  _ohcal_x = cluster->get_x();
  _ohcal_y = cluster->get_y();
  radius_scale = m_ohcal_radius_user / sqrt(_ohcal_x*_ohcal_x+_ohcal_y*_ohcal_y);
  _ohcal_z = radius_scale*cluster->get_z();
  // _ohcal_3x3 = get_e3x3(cluster, _towersOH, 0); //2 for ohcal
  // _ohcal_5x5 = get_e5x5(cluster, _towersOH, 0); //2 for ohcal
  _ohcal_3x3 = std::numeric_limits<float>::quiet_NaN();
  _ohcal_5x5 = std::numeric_limits<float>::quiet_NaN();
  _ohcal_clusE = cluster->get_energy();

  // Variables to determine potential matches
  float dphi = abs(PiRange(_track_phi_ohc - _ohcal_phi));
  float dz = abs(_track_z_ohc - _ohcal_z);
  float deta = abs(_ohcal_eta - _track_eta_ohc);
  float dr = sqrt((dphi*dphi + deta*deta));

  // Requirohents for a possible match
  if(dphi<m_dphi_cut && dz<m_dz_cut)
  {
  //Add potential match's information to vectors
  v_ohcal_phi.push_back(_ohcal_phi);
  v_ohcal_eta.push_back(_ohcal_eta);
  v_ohcal_x.push_back(_ohcal_x);
  v_ohcal_y.push_back(_ohcal_y);
  v_ohcal_z.push_back(_ohcal_z);
  v_ohcal_dphi.push_back(dphi);
  v_ohcal_deta.push_back(deta);
  v_ohcal_dr.push_back(dr);
  v_ohcal_3x3.push_back(_ohcal_3x3);
  v_ohcal_clusE.push_back(_ohcal_clusE);

  is_match = true;
  }
  }


  // Find the closest match from all potential matches
  if (is_match == true)
  {
  float tmp = 99999;
  for(long unsigned int i = 0; i < v_ohcal_dr.size(); i++){
  if(v_ohcal_dr[i] < tmp){
  index = i;
  tmp = v_ohcal_dr[i];
  }
  }
  }
  }

  std::cout << "match identified" << std::endl;

  // Print out statements
  if(index != -1){
  std::cout<<"matched tracks!!!"<<std::endl;
  std::cout<<"ohcal x = "<<v_ohcal_x[index]<<" , y = "<<v_ohcal_y[index]<<" , z = "<<v_ohcal_z[index]<<" , phi = "<<v_ohcal_phi[index]<<" , eta = "<<v_ohcal_eta[index]<<std::endl;
  std::cout<<"track projected x = "<<_track_x_ohc<<" , y = "<<_track_y_ohc<<" , z = "<<_track_z_ohc<<" , phi = "<<_track_phi_ohc<<" , eta = "<<_track_eta_ohc<<std::endl;
  std::cout<<"track px = "<<track->get_px()<<" , py = "<<track->get_py()<<" , pz = "<<track->get_pz()<<" , pt = "<<track->get_pt()<<" , p = "<<track->get_p()<<" , charge = "<<track->get_charge()<<std::endl;
  }
  */

  // Save values to the branches!
  // if(index == -1){
  detector_ohcal_deltaphi[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ohcal_deltaeta[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ohcal_energy_3x3[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ohcal_energy_5x5[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ohcal_cluster_energy[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // }
  // else{
  // detector_ohcal_deltaphi[daughter_id] = v_ohcal_dphi[index];
  // detector_ohcal_deltaeta[daughter_id] = v_ohcal_deta[index];
  // detector_ohcal_energy_3x3[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // detector_ohcal_energy_5x5[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // detector_ohcal_cluster_energy[daughter_id] = v_ohcal_clusE[index];
  // }
  // std::cout << "OHCAL CLLUSTERS MATCHED TO TRACK" << std::endl;
}

void KFParticle_caloInfo::Get5x5CellInfo(
  RawClusterDefs::keytype key_in,
  float &E3x3,
  float &E5x5,
  std::vector<unsigned int> &v_emcal_5x5Cell_Phi,
  std::vector<unsigned int> &v_emcal_5x5Cell_Eta,
  std::vector<float> &v_emcal_5x5Cell_E
  ){
  // Make sure energy values are reset to start at 0
  E5x5 = 0.0;
  E3x3 = 0.0;

  // Get cluster from cluster key
  RawCluster *clus = clustersEM->getCluster(key_in);

  // Set up variables for looping over towers
  unsigned int center_phi = std::numeric_limits<int>::quiet_NaN();
  unsigned int center_eta = std::numeric_limits<int>::quiet_NaN();
  float center_energy = -9999;
  const RawCluster::TowerConstRange begin_end_Towers = clus->get_towers();
  RawCluster::TowerConstIterator towersIter;

  // Loop over the towers in the cluster to find the cluster center
  for (towersIter = begin_end_Towers.first; towersIter != begin_end_Towers.second; ++towersIter)
  {
    if (towersIter->second > center_energy)
    {
      // Setup for phi and eta
      RawTowerGeom *tower_geom = EMCalGeo->get_tower_geometry(towersIter->first);
      // Get energy, phi, and eta of center tower
      center_energy = towersIter->second;
      center_phi = EMCalGeo->get_phibin(tower_geom->get_phi());
      center_eta = EMCalGeo->get_etabin(tower_geom->get_eta());
    }
  }  // Now we have the center tower and can grab the 5x5 grid info

  // Vectors used in next for loops
  std::vector<unsigned int> v_phi_tmp;
  std::vector<unsigned int> v_eta_tmp;
  std::vector<float> v_Tower_E;
//   bool printthis = false;

  // Get phi bins, and account for wrap-around case
  for (int i = -2; i <= 2; i++)
  {
    // Wraparound in phi (256 EMCal phi bins)
    int phibin = (static_cast<int>(center_phi) + i + 256) % 256;
    v_phi_tmp.push_back(phibin);
  }

  // Get eta bins, and account for end of acceptance
  for (int i = -2; i <= 2; i++)
  {
    int etabin = center_eta + i;
    // End of eta acceptance case
    if (etabin < 0 || etabin > 95) {
      continue;
    }
    v_eta_tmp.push_back(etabin);
  }

  // Ranges for for-loops
  int phi_range = v_phi_tmp.size();
  int eta_range = v_eta_tmp.size();

  // for loop to fill vectors which will be stored in the KFP output nTuple
  for (int i = 0; i < phi_range; i++)
  {
    // Get ith phi bin
    const unsigned int iphi = v_phi_tmp[i];

    // Add towers to 3x3 if not first or last rows in phi
    bool addto3x3 = false;
    if(i < 4 && i > 0)
    {
      addto3x3 = true;
    }

    // Loop over eta
    for (int j = 0; j < eta_range; j++)
    {
      // Get jth eta bin
      const unsigned int jeta = v_eta_tmp[j];
      // Get key at phi and eta bin
      unsigned int towerinfo_key = TowerInfoDefs::encode_emcal(jeta, iphi);
      // Get energy
      TowerInfo *tower = _towersEM->get_tower_at_key(towerinfo_key);
      float energytmp = tower->get_energy();
      // Add to 5x5 energy float
      E5x5+=energytmp;
      v_Tower_E.push_back(energytmp);
       // Store in vectors
      if(m_config.calo5x5CellInfo)
      {
        v_emcal_5x5Cell_Phi.push_back(iphi);
        v_emcal_5x5Cell_Eta.push_back(jeta);
        v_emcal_5x5Cell_E.push_back(energytmp);
      }

      if(!addto3x3)
      {
        continue; // Skip adding to E3x3 if phi is outside 3x3 grid
      }
      int etadiff = (int) jeta-center_eta;
      if(etadiff < 2 && etadiff > -2)
      {
        E3x3+=energytmp;
      }

    }
  }
}

void KFParticle_caloInfo::clear(int nTracks)
{
  for (int i = 0; i < nTracks; ++i)
  {
    // Detailed Calo
    if (m_config.calo5x5CellInfo)
    {
      detector_emcal_5x5Cell_Phi[i].clear();
      detector_emcal_5x5Cell_Eta[i].clear();
      detector_emcal_5x5Cell_E[i].clear();
    }
  }
}
