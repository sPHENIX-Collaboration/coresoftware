#include "GlobalVertexReco.h"

//#include "GlobalVertex.h"     // for GlobalVertex, GlobalVe...
#include "GlobalVertexMap.h"  // for GlobalVertexMap
#include "GlobalVertexMapv1.h"
#include "GlobalVertexv2.h"
#include "MbdVertex.h"
#include "MbdVertexMap.h"
#include "CaloVertex.h"
#include "CaloVertexMap.h"
#include "SvtxVertex.h"
#include "SvtxVertexMap.h"
#include "TruthVertex.h"
#include "TruthVertexMap.h"
#include "TruthVertexMap_v1.h"
#include "TruthVertex_v1.h"

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

// truth info, for truth z vertex
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>
#include <limits>
#include <set>      // for _Rb_tree_const_iterator
#include <utility>  // for pair

GlobalVertexReco::GlobalVertexReco(const std::string &name)
  : SubsysReco(name)
{
}

int GlobalVertexReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "======================= GlobalVertexReco::InitRun() =======================" << std::endl;
    std::cout << "===========================================================================" << std::endl;
  }

  return CreateNodes(topNode);
}

int GlobalVertexReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "GlobalVertexReco::process_event -- entered" << std::endl;
  }

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  GlobalVertexMap *globalmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  
  if (!globalmap)
  {
    std::cout << PHWHERE << "::ERROR - cannot find GlobalVertexMap" << std::endl;
    exit(-1);
  }

  SvtxVertexMap *svtxmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  MbdVertexMap *mbdmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  CaloVertexMap *calomap = findNode::getClass<CaloVertexMap>(topNode, "CaloVertexMap");
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  TruthVertexMap *truthmap = findNode::getClass<TruthVertexMap>(topNode, "TruthVertexMap");
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  // we will make 4 different kinds of global vertexes
  //  (1) SVTX+MBD vertexes - we match SVTX vertex to the nearest MBD vertex within 3 sigma in zvertex
  //      the spatial point comes from the SVTX, the timing from the MBD
  //      number of SVTX+MBD vertexes <= number of SVTX vertexes

  //  (2) SVTX only vertexes - for those cases where the MBD wasn't simulated,
  //      or all the MBD vertexes are outside the 3 sigma matching requirement
  //      we pass forward the 3d point from the SVTX with the default timing info

  //  (3) MBD only vertexes - use the default x,y positions on this module and
  //      pull in the mbd z and mbd t

  //  (4) Truth vertexes - if the G4TruthInfo node is present (should be for simulation only), we will
  //      create a GlobalVertex from the primary truth vertex, and create+fill a TruthVertexMap node

  // there may be some quirks as we get to large luminosity and the MBD becomes
  // untrust worthy, I'm guessing analyzers would resort exclusively to (1) or (2)
  // in those cases

  std::set<unsigned int> used_svtx_vtxids;
  std::set<unsigned int> used_mbd_vtxids;
  std::set<unsigned int> used_calo_vtxids;

  if (svtxmap && mbdmap && useVertexType(GlobalVertex::VTXTYPE::SVTX_MBD))
  {
    if (Verbosity())
    {
      std::cout << "GlobalVertexReco::process_event - svtxmap && mbdmap" << std::endl;
    }

    for (SvtxVertexMap::ConstIter svtxiter = svtxmap->begin();
         svtxiter != svtxmap->end();
         ++svtxiter)
    {
      const SvtxVertex *svtx = svtxiter->second;

      const MbdVertex *mbd_best = nullptr;
      float min_sigma = std::numeric_limits<float>::max();
      for (MbdVertexMap::ConstIter mbditer = mbdmap->begin();
           mbditer != mbdmap->end();
           ++mbditer)
      {
        const MbdVertex *mbd = mbditer->second;

        float combined_error = sqrt(svtx->get_error(2, 2) + pow(mbd->get_z_err(), 2));
        float sigma = std::fabs(svtx->get_z() - mbd->get_z()) / combined_error;
        if (sigma < min_sigma)
        {
          min_sigma = sigma;
          mbd_best = mbd;
        }
      }

      if (min_sigma > 3.0 || !mbd_best)
      {
        continue;
      }

      // we have a matching pair
      GlobalVertex *vertex = new GlobalVertexv2();
      vertex->set_id(globalmap->size());

      vertex->clone_insert_vtx(GlobalVertex::SVTX, svtx);
      vertex->clone_insert_vtx(GlobalVertex::MBD, mbd_best);
      used_svtx_vtxids.insert(svtx->get_id());
      used_mbd_vtxids.insert(mbd_best->get_id());
      vertex->set_id(globalmap->size());

      globalmap->insert(vertex);

      //! Reset track ids to the new vertex object
      if (trackmap)
      {
        for (auto iter = svtx->begin_tracks(); iter != svtx->end_tracks();
             ++iter)
        {
          auto *track = trackmap->find(*iter)->second;
          track->set_vertex_id(vertex->get_id());
        }
      }

      if (Verbosity() > 1)
      {
        vertex->identify();
      }
    }
  }

  // okay now loop over all unused SVTX vertexes (2nd class)...
  if (svtxmap && useVertexType(GlobalVertex::VTXTYPE::SVTX))
  {
    if (Verbosity())
    {
      std::cout << "GlobalVertexReco::process_event - svtxmap " << std::endl;
    }

    for (SvtxVertexMap::ConstIter svtxiter = svtxmap->begin();
         svtxiter != svtxmap->end();
         ++svtxiter)
    {
      const SvtxVertex *svtx = svtxiter->second;

      if (used_svtx_vtxids.contains(svtx->get_id()))
      {
        continue;
      }
      if (std::isnan(svtx->get_z()))
      {
        continue;
      }

      // we have a standalone SVTX vertex
      GlobalVertex *vertex = new GlobalVertexv2();

      vertex->set_id(globalmap->size());

      vertex->clone_insert_vtx(GlobalVertex::SVTX, svtx);
      used_svtx_vtxids.insert(svtx->get_id());

      //! Reset track ids to the new vertex object
      if (trackmap)
      {
        for (auto iter = svtx->begin_tracks(); iter != svtx->end_tracks();
             ++iter)
        {
          auto *track = trackmap->find(*iter)->second;
          track->set_vertex_id(vertex->get_id());
        }
      }
      globalmap->insert(vertex);

      if (Verbosity() > 1)
      {
        vertex->identify();
      }
    }
  }

  // okay now loop over all unused MBD vertexes (3rd class)...
  if (mbdmap && useVertexType(GlobalVertex::VTXTYPE::MBD))
  {
    if (Verbosity())
    {
      std::cout << "GlobalVertexReco::process_event -  mbdmap" << std::endl;
    }

    for (MbdVertexMap::ConstIter mbditer = mbdmap->begin();
         mbditer != mbdmap->end();
         ++mbditer)
    {
      const MbdVertex *mbd = mbditer->second;

      if (used_mbd_vtxids.contains(mbd->get_id()))
      {
        continue;
      }

      if (std::isnan(mbd->get_z()))
      {
        continue;
      }

      GlobalVertex *vertex = new GlobalVertexv2();
      vertex->set_id(globalmap->size());

      vertex->clone_insert_vtx(GlobalVertex::MBD, mbd);
      used_mbd_vtxids.insert(mbd->get_id());

      globalmap->insert(vertex);

      if (Verbosity() > 1)
      {
        vertex->identify();
      }
    }
  }

  // okay now loop over all unused MBD vertexes (3rd class)...
  if (calomap && useVertexType(GlobalVertex::VTXTYPE::CALO))
  {
    if (Verbosity())
    {
      std::cout << "GlobalVertexReco::process_event -  calomap" << std::endl;
    }

    for (CaloVertexMap::ConstIter caloiter = calomap->begin();
         caloiter != calomap->end();
         ++caloiter)
      {
	const CaloVertex *calo = caloiter->second;
	
	if (used_calo_vtxids.contains(calo->get_id()))
	  {
	    continue;
	  }
	
	if (std::isnan(calo->get_z()))
	  {
	    continue;
	  }
	
	GlobalVertex *vertex = new GlobalVertexv2();
	vertex->set_id(globalmap->size());
	
	vertex->clone_insert_vtx(GlobalVertex::CALO, calo);
	used_calo_vtxids.insert(calo->get_id());
	
	globalmap->insert(vertex);
	
	if (Verbosity() > 1)
	  {
	    vertex->identify();
	  }
      }
  }

// okay now loop over all unused MBD vertexes (3rd class)...
  if (mbdmap && calomap && useVertexType(GlobalVertex::VTXTYPE::MBD_CALO))
  {
    if (Verbosity())
    {
      std::cout << "GlobalVertexReco::process_event -  calomap + mbdmap" << std::endl;
    }

    for (MbdVertexMap::ConstIter mbditer = mbdmap->begin();
         mbditer != mbdmap->end();
         ++mbditer)
    {
      const MbdVertex *mbd = mbditer->second;

      if (used_mbd_vtxids.contains(mbd->get_id()))
      {
        continue;
      }

      
      bool getcalo =  std::isnan(mbd->get_z());

      if (getcalo)
	{
	  for (CaloVertexMap::ConstIter caloiter = calomap->begin();
	       caloiter != calomap->end();
	       ++caloiter)
	    {
	      const CaloVertex *calo = caloiter->second;
	      
	      if (used_calo_vtxids.contains(calo->get_id()))
		{
		  continue;
		}
	      
	      if (std::isnan(calo->get_z()))
		{
		  continue;
		}
	      
	      GlobalVertex *vertex = new GlobalVertexv2();
	      vertex->set_id(globalmap->size());
	
	      vertex->clone_insert_vtx(GlobalVertex::CALO, calo);
	      used_calo_vtxids.insert(calo->get_id());
	      
	      globalmap->insert(vertex);
	      
	      if (Verbosity() > 1)
		{
		  vertex->identify();
		}
	    }
	  
	}
      else
	{
	  GlobalVertex *vertex = new GlobalVertexv2();
	  vertex->set_id(globalmap->size());
	  
	  vertex->clone_insert_vtx(GlobalVertex::MBD, mbd);
	  used_mbd_vtxids.insert(mbd->get_id());
	  
	  globalmap->insert(vertex);
	  
	  if (Verbosity() > 1)
	    {
	      vertex->identify();
	    }
	}
    }
  }
// okay now add the truth vertex (4th class)...

  if (truthinfo)
  {
    PHG4VtxPoint *vtxp = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
    float truth_z_vertex = std::numeric_limits<float>::quiet_NaN();
    float truth_x_vertex = std::numeric_limits<float>::quiet_NaN();
    float truth_y_vertex = std::numeric_limits<float>::quiet_NaN();
    if (vtxp)
    {
      truth_z_vertex = vtxp->get_z();
      truth_x_vertex = vtxp->get_x();
      truth_y_vertex = vtxp->get_y();
      TruthVertex *tvertex = new TruthVertex_v1();
      tvertex->set_id(globalmap->size());
      tvertex->set_z(truth_z_vertex);
      tvertex->set_z_err(0);
      tvertex->set_x(truth_x_vertex);
      tvertex->set_x_err(0);
      tvertex->set_y(truth_y_vertex);
      tvertex->set_y_err(0);

      tvertex->set_t(0);
      tvertex->set_t_err(0);  // 0.1
      GlobalVertex *vertex = new GlobalVertexv2();
      vertex->clone_insert_vtx(GlobalVertex::TRUTH, tvertex);
      globalmap->insert(vertex);
      if (truthmap)
      {
        truthmap->insert(tvertex);
      }
      if (Verbosity() > 1)
      {
        vertex->identify();
      }
    }
    else if (Verbosity())
    {
      std::cout << "PHG4TruthInfoContainer has no primary vertex" << std::endl;
    }
  }
  else if (Verbosity())
  {
    std::cout << "PHG4TruthInfoContainer is missing" << std::endl;
  }

  /// Associate any tracks that were not assigned a track-vertex
  if (trackmap)
  {
    for (const auto &[tkey, track] : *trackmap)
    {
      //! Check that the vertex hasn't already been assigned
      auto trackvtxid = track->get_vertex_id();
      if (globalmap->find(trackvtxid) != globalmap->end())
      {
        continue;
      }

      float maxdz = std::numeric_limits<float>::max();
      unsigned int vtxid = std::numeric_limits<unsigned int>::max();

      for (const auto &[vkey, vertex] : *globalmap)
      {
        float dz = track->get_z() - vertex->get_z();
        if (std::fabs(dz) < maxdz)
        {
          maxdz = dz;
          vtxid = vkey;
        }
      }

      track->set_vertex_id(vtxid);
      if (Verbosity())
      {
        std::cout << "Associated track with z " << track->get_z() << " to GlobalVertex id " << track->get_vertex_id() << std::endl;
      }
    }
  }

  if (Verbosity())
  {
    globalmap->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalVertexReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the GLOBAL stuff under a sub-node directory
  PHCompositeNode *globalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!globalNode)
  {
    globalNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(globalNode);
  }

  // create the GlobalVertexMap
  GlobalVertexMap *vertexes = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexes)
  {
    vertexes = new GlobalVertexMapv1();
    PHIODataNode<PHObject> *VertexMapNode = new PHIODataNode<PHObject>(vertexes, "GlobalVertexMap", "PHObject");
    globalNode->addNode(VertexMapNode);
  }

  TruthVertexMap *truthmap = findNode::getClass<TruthVertexMap>(topNode, "TruthVertexMap");
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  // Create the TruthVertexMap only if the PHG4TruthInfoContainer is present
  if (!truthmap && truthinfo)
  {
    if (Verbosity())
    {
      std::cout << "Creating TruthVertexMap node" << std::endl;
    }
    truthmap = new TruthVertexMap_v1();
    PHIODataNode<PHObject> *TruthVertexMapNode = new PHIODataNode<PHObject>(truthmap, "TruthVertexMap", "PHObject");
    globalNode->addNode(TruthVertexMapNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
