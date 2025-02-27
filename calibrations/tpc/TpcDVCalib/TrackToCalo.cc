/*!
 *  \file               TrackToCalo.cc
 *  \brief              Track To Calo matching, for TPC drift velocity calibration
 *  \author Xudong Yu <xyu3@bnl.gov>
 */
#include "TrackToCalo.h"

#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertex.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssocv1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackreco/ActsPropagator.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>

#include <CLHEP/Vector/ThreeVector.h>
#include <cmath>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

//____________________________________________________________________________..
TrackToCalo::TrackToCalo(const std::string &name, const std::string &file):
 SubsysReco(name),
 _outfilename(file),
 _outfile(nullptr),
 _tree(nullptr)
{
  std::cout << "TrackToCalo::TrackToCalo(const std::string &name, const std::string &file) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
TrackToCalo::~TrackToCalo()
{
  std::cout << "TrackToCalo::~TrackToCalo() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int TrackToCalo::Init(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  std::cout << "TrackToCalo::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  delete _outfile;
  _outfile = new TFile(_outfilename.c_str(), "RECREATE");
  delete _tree;
  _tree = new TTree("tree", "A tree with track/calo info");
  _tree->Branch("calo_z",&m_calo_z,"calo_z/F");
  _tree->Branch("calo_phi",&m_calo_phi,"calo_phi/F");
  _tree->Branch("track_z",&m_track_z,"track_z/F");
  _tree->Branch("track_phi",&m_track_phi,"track_phi/F");
  _tree->Branch("dphi",&m_dphi,"dphi/F");
  _tree->Branch("dz",&m_dz,"dz/F");
  cnt=0;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrackToCalo::process_event(PHCompositeNode *topNode)
{
  std::cout<<"TrackToCalo::process_event event "<<cnt<<std::endl;
  cnt++;
  ResetTreeVectors();

  if(!svtxTrackMap)
  {
    svtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    if(!svtxTrackMap)
    {
      std::cout << "TrackToCalo::process_event : SvtxTrackMap (SvtxTrackMap) not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if ( !rawClusterContainer )
  {
    rawClusterContainer = findNode::getClass<RawClusterContainer>(topNode, m_RawClusCont_EM_name);
    if (!rawClusterContainer)
    {
      std::cout << "TrackToCalo::process_event : RawClusterContainer (" << m_RawClusCont_EM_name << ") not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(!trkrClusterContainer)
  {
    trkrClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if(!trkrClusterContainer)
    {
      std::cout << "TrackToCalo::process_event : TrkrClusterContainer (TRKR_CLUSTER) not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(!rawTowerGeomContainer)
  {
    rawTowerGeomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if(!rawTowerGeomContainer)
    {
      std::cout << "TrackToCalo::process_event : RawTowerGeomContainer (TOWERGEOM_CEMC) not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  double caloRadiusEMCal;
  if (m_use_emcal_radius)
  {
    caloRadiusEMCal = m_emcal_radius_user;
  }
  else
  {
    caloRadiusEMCal = rawTowerGeomContainer->get_radius();
  }

  for (auto &iter : *svtxTrackMap)
  {
    track = iter.second;

    if(!track)
    {
      continue;
    }

    if(track->get_pt() < m_track_pt_low_cut)
    {
      continue;
    }

    int n_tpc_clusters = 0;
    tpc_seed = track->get_tpc_seed();
    if(tpc_seed)
    {
      for(auto key_iter = tpc_seed->begin_cluster_keys(); key_iter != tpc_seed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        trkrCluster = trkrClusterContainer->findCluster(cluster_key);
        if(!trkrCluster)
        {
          continue;
        }
        if(TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::TrkrId::tpcId)
        {
          n_tpc_clusters++;
        }
      }
    }

    if(n_tpc_clusters<m_ntpc_low_cut) 
    {
      continue;
    }

    // project to R_EMCAL
    thisState = track->get_state(caloRadiusEMCal);

    if(!thisState)
    {
      _track_phi_emc.push_back(NAN);
      _track_z_emc.push_back(NAN);
    }
    else
    {
      _track_phi_emc.push_back(std::atan2(thisState->get_y(), thisState->get_x()));
      _track_z_emc.push_back(thisState->get_z());
    }
  }

  RawClusterContainer::Range begin_end_EMC = rawClusterContainer->getClusters();
  RawClusterContainer::Iterator clusIter_EMC;

  /// Loop over the EMCal clusters
  for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
  {
    cluster = clusIter_EMC->second;

    if(cluster->get_energy() < m_emcal_e_low_cut)
    {
      continue;
    }

    _emcal_phi.push_back( std::atan2(cluster->get_y(), cluster->get_x()) );
    radius_scale = caloRadiusEMCal / sqrt( pow(cluster->get_x(),2) + pow(cluster->get_y(),2) );
    _emcal_z.push_back( radius_scale*cluster->get_z() );
  }

  for(unsigned int itrack = 0; itrack < _track_phi_emc.size(); itrack++)
  {
    if(std::isnan(_track_phi_emc.at(itrack)) || std::isnan(_track_z_emc.at(itrack)))
    {
      continue;
    }

    for(unsigned int iem = 0; iem < _emcal_z.size(); iem++)
    {
      m_calo_z = _emcal_z.at(iem);
      m_track_z = _track_z_emc.at(itrack);
      m_calo_phi = _emcal_phi.at(iem);
      m_track_phi = _track_phi_emc.at(itrack);
      m_dphi = PiRange(_track_phi_emc.at(itrack) - _emcal_phi.at(iem));
      m_dz = _track_z_emc.at(itrack) - _emcal_z.at(iem);
      _tree->Fill();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrackToCalo::End(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  _outfile->cd();
  _outfile->Write();
  _outfile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackToCalo::ResetTreeVectors()
{
  _track_phi_emc.clear();
  _track_z_emc.clear();
  _emcal_phi.clear();
  _emcal_z.clear();
}

float TrackToCalo::PiRange(float deltaPhi)
{
  if(deltaPhi > M_PI) 
  {
    deltaPhi -= 2*M_PI;
  }
  if(deltaPhi < -M_PI)
  {
    deltaPhi += 2*M_PI;
  }

  return deltaPhi;
}
