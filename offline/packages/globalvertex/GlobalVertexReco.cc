#include "GlobalVertexReco.h"

#include "GlobalVertex.h"     // for GlobalVertex, GlobalVe...
#include "GlobalVertexMap.h"  // for GlobalVertexMap
#include "GlobalVertexMapv1.h"
#include "GlobalVertexv2.h"
#include "SvtxVertex.h"
#include "SvtxVertexMap.h"
#include "MbdVertex.h"
#include "MbdVertexMap.h"

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

#include <cfloat>
#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>
#include <set>      // for _Rb_tree_const_iterator
#include <utility>  // for pair

using namespace std;

GlobalVertexReco::GlobalVertexReco(const string &name)
  : SubsysReco(name)
{
}

int GlobalVertexReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    cout << "======================= GlobalVertexReco::InitRun() =======================" << endl;
    cout << "===========================================================================" << endl;
  }

  return CreateNodes(topNode);
}

int GlobalVertexReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    cout << "GlobalVertexReco::process_event -- entered" << endl;
  }

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  GlobalVertexMap *globalmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!globalmap)
  {
    cout << PHWHERE << "::ERROR - cannot find GlobalVertexMap" << endl;
    exit(-1);
  }

  SvtxVertexMap *svtxmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  MbdVertexMap *mbdmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  // we will make 3 different kinds of global vertexes
  //  (1) SVTX+MBD vertexes - we match SVTX vertex to the nearest MBD vertex within 3 sigma in zvertex
  //      the spatial point comes from the SVTX, the timing from the MBD
  //      number of SVTX+MBD vertexes <= number of SVTX vertexes

  //  (2) SVTX only vertexes - for those cases where the MBD wasn't simulated,
  //      or all the MBD vertexes are outside the 3 sigma matching requirement
  //      we pass forward the 3d point from the SVTX with the default timing info

  //  (3) MBD only vertexes - use the default x,y positions on this module and
  //      pull in the mbd z and mbd t

  // there may be some quirks as we get to large luminosity and the MBD becomes
  // untrust worthy, I'm guessing analyzers would resort exclusively to (1) or (2)
  // in those cases

  std::set<unsigned int> used_svtx_vtxids;
  std::set<unsigned int> used_mbd_vtxids;
 
  if (svtxmap && mbdmap)
  {
    if (Verbosity())
    {
      cout << "GlobalVertexReco::process_event - svtxmap && mbdmap" << endl;
    }

    for (SvtxVertexMap::ConstIter svtxiter = svtxmap->begin();
         svtxiter != svtxmap->end();
         ++svtxiter)
    {
      const SvtxVertex *svtx = svtxiter->second;

      const MbdVertex *mbd_best = nullptr;
      float min_sigma = FLT_MAX;
      for (MbdVertexMap::ConstIter mbditer = mbdmap->begin();
           mbditer != mbdmap->end();
           ++mbditer)
      {
        const MbdVertex *mbd = mbditer->second;

        float combined_error = sqrt(svtx->get_error(2, 2) + pow(mbd->get_z_err(), 2));
        float sigma = fabs(svtx->get_z() - mbd->get_z()) / combined_error;
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

      vertex->insert_vtx(GlobalVertex::SVTX, svtx);
      vertex->insert_vtx(GlobalVertex::MBD, mbd_best);
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
          auto track = trackmap->find(*iter)->second;
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
  if (svtxmap)
  {
    if (Verbosity())
    {
      cout << "GlobalVertexReco::process_event - svtxmap " << endl;
    }

    for (SvtxVertexMap::ConstIter svtxiter = svtxmap->begin();
         svtxiter != svtxmap->end();
         ++svtxiter)
    {
      const SvtxVertex *svtx = svtxiter->second;

      if (used_svtx_vtxids.find(svtx->get_id()) != used_svtx_vtxids.end())
      {
        continue;
      }
      if (isnan(svtx->get_z()))
      {
        continue;
      }

      // we have a standalone SVTX vertex
      GlobalVertex *vertex = new GlobalVertexv2();
      
      vertex->set_id(globalmap->size());

      vertex->insert_vtx(GlobalVertex::SVTX, svtx);
      used_svtx_vtxids.insert(svtx->get_id());

      //! Reset track ids to the new vertex object
      if (trackmap)
      {
        for (auto iter = svtx->begin_tracks(); iter != svtx->end_tracks();
             ++iter)
        {
          auto track = trackmap->find(*iter)->second;
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
  if (mbdmap)
  {
    if (Verbosity())
    {
      cout << "GlobalVertexReco::process_event -  mbdmap" << endl;
    }

    for (MbdVertexMap::ConstIter mbditer = mbdmap->begin();
         mbditer != mbdmap->end();
         ++mbditer)
    {
      const MbdVertex *mbd = mbditer->second;

      if (used_mbd_vtxids.find(mbd->get_id()) != used_mbd_vtxids.end())
      {
        continue;
      }
      if (isnan(mbd->get_z()))
      {
        continue;
      }

      GlobalVertex *vertex = new GlobalVertexv2();
      vertex->set_id(globalmap->size());

      vertex->insert_vtx(GlobalVertex::MBD, mbd);
      used_mbd_vtxids.insert(mbd->get_id());

      globalmap->insert(vertex);

      if (Verbosity() > 1)
      {
        vertex->identify();
      }
    }
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
        if (fabs(dz) < maxdz)
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
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
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
  return Fun4AllReturnCodes::EVENT_OK;
}
