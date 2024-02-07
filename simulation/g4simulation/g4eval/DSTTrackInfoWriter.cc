/*!
 * \file DSTTrackInfoWriter.cc
 * \author Alex Patton <aopatton@mit.edu>
 */

#include "DSTTrackInfoWriter.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <micromegas/MicromegasDefs.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrDefs.h>
// include new cluster
#include <trackbase/TrkrClusterContainer.h>
// #include <trackbase/TrkrClusterContainerv5.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed_v1.h>
// include SvtxTrackMap_v1 to write data to it
#include <trackbase_historic/SvtxTrackInfo_v1.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/TrackInfoContainer_v1.h>
#include <trackbase_historic/TrackStateInfo_v1.h>

#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <numeric>

#include <TClonesArray.h>
#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
//_____________________________________________________________________

//_____________________________________________________________________
DSTTrackInfoWriter::DSTTrackInfoWriter(const std::string& name)
  : SubsysReco(name)
{
}

//_____________________________________________________________________
int DSTTrackInfoWriter::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "Writer Init start" << std::endl;
  }
  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "DSTTrackInfoWriter::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if (!evalNode)
  {
    // create
    std::cout << "DSTTrackInfoWriter::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode("EVAL");
    dstNode->addNode(evalNode);
  }

  // // TClonesArary
  // m_container->arrClsDST = new TClonesArray("DSTContainerv3::ClusterStruct");
  // m_container->trkrClsDST = new TClonesArray("TrkrClusterv4");

  auto newInfoNode = new PHIODataNode<PHObject>(new TrackInfoContainer_v1, "TrackInfoContainer", "PHObject");
  evalNode->addNode(newInfoNode);

  // auto newTrackNode = new PHIODataNode<PHObject>(new SvtxTrackArray_v1, "TRACK_ARRAYV1", "PHObject");
  // evalNode->addNode(newTrackNode);

  // TClonesArary
  // fcl = new TFile("dstcl.root", "recreate");
  // tcl = new TTree("tcl", "dst tree");
  // arrEvt = new TClonesArray("DSTContainerv3::EventStruct");
  // arrTrk = new TClonesArray("DSTContainerv3::TrackStruct");
  // arrCls = new TClonesArray("DSTContainerv3::ClusterStruct");
  // tcl->Branch("evt", &arrEvt);
  // tcl->Branch("trk", &arrTrk);
  // tcl->Branch("cls", &arrCls);
  if (Verbosity() > 1)
  {
    std::cout << "Writer Init end" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTTrackInfoWriter::process_event(PHCompositeNode* topNode)
{
  // make topNode run in Init
  // Init(topNode);
  //  load nodes
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << "::" << __func__ << "::" << __LINE__ << std::endl;
    std::cout << "DSTTrackInfoWriter::process_event" << std::endl;
  }
  auto res = load_nodes(topNode);
  if (res != Fun4AllReturnCodes::EVENT_OK)
  {
    return res;
  }
  if (Verbosity() > 1)
  {
    std::cout << "Return codes  start" << Fun4AllReturnCodes::EVENT_OK << std::endl;
  }
  // cleanup output container
  if (m_track_info_container) m_track_info_container->Reset();
  if (Verbosity() > 1)
  {
    std::cout << "Evalutate track info" << std::endl;
  }
  evaluate_track_info();
  if (Verbosity() > 1)
  {
    std::cout << "exiting event"
              << "\n";

    std::cout << "Return codes end" << Fun4AllReturnCodes::EVENT_OK << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTTrackInfoWriter::load_nodes(PHCompositeNode* topNode)
{
  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  m_track_info_container = findNode::getClass<TrackInfoContainer_v1>(topNode, "TrackInfoContainer");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________
void DSTTrackInfoWriter::evaluate_track_info()
{
  if (!(m_track_info_container))
  {
    return;
  }

  m_track_info_container->Reset();
  // get track into track info
  Int_t iTrk = 0;
  // long unsigned int iKey = 0;
  if (Verbosity() > 1)
  {
    std::cout << "Before loop"
              << "\n";
  }
  for (const auto& trackpair : *m_track_map)
  {
    const auto track = trackpair.second;
    // this track will have a TPC and Silicon seed

    uint64_t hitbitmap = 0;

    SvtxTrackInfo_v1* trackInfo = new SvtxTrackInfo_v1();
    if (Verbosity() > 1)
    {
      std::cout << "Before seeds"
                << "\n";
    }
    TrackSeed* TPCSeed = track->get_tpc_seed();
    TrackSeed* SiliconSeed = track->get_silicon_seed();
    if (TPCSeed)
    {
      for (auto key_iter = TPCSeed->begin_cluster_keys(); key_iter != TPCSeed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;

        // store information in track array
        // std::cout << "TPC clusterkey: " << cluster_key <<"\n";
        // std::cout << "TPC subsurfkey: " << cluster->getSubSurfKey() << std::endl;
        uint8_t layer = TrkrDefs::getLayer(cluster_key);
        if (Verbosity() > 1)
        {
          std::cout << "Layer is: " << unsigned(layer) << std::endl;
        }
        hitbitmap = hitbitmap + ((uint64_t) 1 << layer);

        // TrkrDefs::
      }
    }
    if (Verbosity() > 1)
    {
      std::cout << "Before Silicon seeds"
                << "\n";
    }
    if (!SiliconSeed)
    {
      if (Verbosity() > 1)
      {
        std::cout << "Silicon Seed does not exist" << std::endl;
      }
    }

    if (SiliconSeed)
    {
      for (auto key_iter = SiliconSeed->begin_cluster_keys(); key_iter != SiliconSeed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;

        uint8_t layer = TrkrDefs::getLayer(cluster_key);
        if (Verbosity() > 1)
        {
          std::cout << "Layer is: " << unsigned(layer) << std::endl;
        }
        hitbitmap = hitbitmap + ((uint64_t) 1 << layer);
      }
    }
    if (Verbosity() > 1)
    {
      std::cout << "After Track seeds"
                << "\n";
    }
    trackInfo->set_hitbitmap(hitbitmap);
    trackInfo->set_chisq(track->get_chisq());
    trackInfo->set_ndf(track->get_ndf());
    trackInfo->set_crossing(track->get_crossing());

    trackInfo->set_x(track->get_x());
    trackInfo->set_y(track->get_y());
    trackInfo->set_z(track->get_z());
    trackInfo->set_px(track->get_px());
    trackInfo->set_py(track->get_py());
    trackInfo->set_pz(track->get_pz());
    if (Verbosity() > 1)
    {
      std::cout << "track crossing: " << track->get_crossing() << std::endl;

      std::cout << "track.get_z(): " << track->get_z() << std::endl;
      std::cout << "trackInfo.get_z(): " << trackInfo->get_z() << std::endl;
      std::cout << "Hitbitmap: " << trackInfo->get_hitbitmap() << std::endl;
      std::cout << "crossing: " << trackInfo->get_crossing() << std::endl;
      std::cout << "chi^2: " << trackInfo->get_chisq() << std::endl;
      std::cout << "ndf: " << unsigned(trackInfo->get_ndf()) << std::endl;
    }
    int covarianceIndex = 0;
    for (int i = 0; i < 6; i++)
    {
      for (int j = i; j < 6; j++)
      {
        if (Verbosity() > 1)
        {
          std::cout << "covariance index: " << covarianceIndex << std::endl;
        }
        trackInfo->set_covariance(covarianceIndex, track->get_error(i, j));
        covarianceIndex++;
      }
    }
    if (Verbosity() > 1)
    {
      std::cout << "Right before adding track info" << iTrk << std::endl;
    }
    m_track_info_container->add_trackinfo(iTrk, trackInfo);
    if (Verbosity() > 1)
    {
      std::cout << "Right after adding track info" << std::endl;
    }
    delete trackInfo;
    ++iTrk;
  }

  // add trackinfo to trackinfocontainer
}
//_____________________________________________________________________

//_____________________________________________________________________
