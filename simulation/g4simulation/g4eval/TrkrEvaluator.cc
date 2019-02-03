#include "TrkrEvaluator.h"

#include <trackbase_historic/SvtxCluster.h>
#include <trackbase_historic/SvtxClusterMap.h>
#include <trackbase_historic/SvtxHit.h>
#include <trackbase_historic/SvtxHitMap.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrDefs.h>
#include <tpc/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>


#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHTimeServer.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <set>

using namespace std;

TrkrEvaluator::TrkrEvaluator(const string& name, const string& filename,
                             unsigned int nlayers_maps,
                             unsigned int nlayers_intt,
                             unsigned int nlayers_tpc)
  : SubsysReco("TrkrEvaluator")
 {
_filename = filename;

}

int TrkrEvaluator::Init(PHCompositeNode* topNode)
{
_ievent = 0;

cout << "will write output to " << _filename.c_str() << endl;

  _tfile = new TFile(_filename.c_str(), "RECREATE");

_do_cluster_eval = true;

  if (_do_cluster_eval) _ntp_cluster = new TNtuple("ntp_cluster", "svtxcluster => max truth",
                                                   "event:hitID:x:y:z:r:phi:eta:ex:ey:ez:ephi:"
                                                   "e:adc:layer:size:phisize:"
                                                   "zsize:trackID:g4hitID:gx:"
                                                   "gy:gz:gr:gphi:geta:gt:gtrackID:gflavor:"
                                                   "gpx:gpy:gpz:gvx:gvy:gvz:gvt:"
                                                   "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                                                   "gembed:gprimary:efromtruth:nparticles:"
                                                   "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps");

  _timer = new PHTimer("_eval_timer");
  _timer->stop();

  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkrEvaluator::InitRun(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkrEvaluator::process_event(PHCompositeNode* topNode)
{
  //if ((Verbosity() > 0) && (_ievent % 100 == 0))
  {
    cout << "TrkrEvaluator::process_event - Event = " << _ievent << endl;
  }

  //-----------------------------------
  // print what is coming into the code
  //-----------------------------------

  //printInputInfo(topNode);

  //---------------------------
  // fill the Evaluator NTuples
  //---------------------------

  fillOutputNtuples(topNode);

  //--------------------------------------------------
  // Print out the ancestry information for this event
  //--------------------------------------------------

  //printOutputInfo(topNode);

  ++_ievent;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkrEvaluator::End(PHCompositeNode* topNode)
{
  _tfile->cd();

  if (_ntp_cluster) _ntp_cluster->Write();

  _tfile->Close();

  delete _tfile;

  //if (Verbosity() > 0)
  {
    cout << "========================= TrkrEvaluator::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

//if (Verbosity() > -1)
  {
    if ((_errors > 0) || (Verbosity() > 0))
    {
      cout << "TrkrEvaluator::End() - Error Count: " << _errors << endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

/*
void TrkrEvaluator::printInputInfo(PHCompositeNode* topNode)
{

  if (Verbosity() > 1) cout << "TrkrEvaluator::printInputInfo() entered" << endl;

  if (Verbosity() > 3)
  {
    // event information
    cout << endl;
    cout << PHWHERE << "   INPUT FOR EVENT " << _ievent << endl;

    cout << endl;
    cout << "---PHG4HITS-------------" << endl;
    _svtxevalstack->get_truth_eval()->set_strict(_strict);
    std::set<PHG4Hit*> g4hits = _svtxevalstack->get_truth_eval()->all_truth_hits();
    unsigned int ig4hit = 0;
    for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
         iter != g4hits.end();
         ++iter)
    {
      PHG4Hit* g4hit = *iter;
      cout << ig4hit << " of " << g4hits.size();
      cout << ": PHG4Hit: " << endl;
      g4hit->identify();
      ++ig4hit;
    }

    cout << "---SVTXCLUSTERS-------------" << endl;
    SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
    if (clustermap)
    {
      unsigned int icluster = 0;
      for (SvtxClusterMap::Iter iter = clustermap->begin();
           iter != clustermap->end();
           ++iter)
      {
        SvtxCluster* cluster = iter->second;
        cout << icluster << " of " << clustermap->size();
        cout << ": SvtxCluster: " << endl;
        cluster->identify();
        ++icluster;
      }
    }

    cout << "---SVXTRACKS-------------" << endl;
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());
    if (trackmap)
    {
      unsigned int itrack = 0;
      for (SvtxTrackMap::Iter iter = trackmap->begin();
           iter != trackmap->end();
           ++iter)
      {
        cout << itrack << " of " << trackmap->size();
        SvtxTrack* track = iter->second;
        cout << " : SvtxTrack:" << endl;
        track->identify();
        cout << endl;
      }
    }

    cout << "---SVXVERTEXES-------------" << endl;
    SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    if (vertexmap)
    {
      unsigned int ivertex = 0;
      for (SvtxVertexMap::Iter iter = vertexmap->begin();
           iter != vertexmap->end();
           ++iter)
      {
        cout << ivertex << " of " << vertexmap->size();
        SvtxVertex* vertex = iter->second;
        cout << " : SvtxVertex:" << endl;
        vertex->identify();
        cout << endl;
      }
    }
  }


  return;
}
*/
/*
void TrkrEvaluator::printOutputInfo(PHCompositeNode* topNode)
{

  if (Verbosity() > 1) cout << "TrkrEvaluator::printOutputInfo() entered" << endl;

  //==========================================
  // print out some useful stuff for debugging
  //==========================================

  if (Verbosity() > 0)
  {
    SvtxTrackEval* trackeval = _svtxevalstack->get_track_eval();
    SvtxClusterEval* clustereval = _svtxevalstack->get_cluster_eval();
    SvtxTruthEval* trutheval = _svtxevalstack->get_truth_eval();

    // event information
    cout << endl;
    cout << PHWHERE << "   NEW OUTPUT FOR EVENT " << _ievent << endl;
    cout << endl;

    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

    PHG4VtxPoint* gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
    float gvx = gvertex->get_x();
    float gvy = gvertex->get_y();
    float gvz = gvertex->get_z();

    float vx = NAN;
    float vy = NAN;
    float vz = NAN;

    SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    if (vertexmap)
    {
      if (!vertexmap->empty())
      {
        SvtxVertex* vertex = (vertexmap->begin()->second);

        vx = vertex->get_x();
        vy = vertex->get_y();
        vz = vertex->get_z();
      }
    }

    cout << "===Vertex Reconstruction=======================" << endl;
    cout << "vtrue = (" << gvx << "," << gvy << "," << gvz << ") => vreco = (" << vx << "," << vy << "," << vz << ")" << endl;
    cout << endl;

    cout << "===Tracking Summary============================" << endl;
    unsigned int ng4hits[100] = {0};
    std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits();
    for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
         iter != g4hits.end();
         ++iter)
    {
      PHG4Hit* g4hit = *iter;
      ++ng4hits[g4hit->get_layer()];
    }

    SvtxHitMap* hitmap = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
    unsigned int nhits[100] = {0};
    if (hitmap)
    {
      for (SvtxHitMap::Iter iter = hitmap->begin();
           iter != hitmap->end();
           ++iter)
      {
        SvtxHit* hit = iter->second;
        ++nhits[hit->get_layer()];
      }
    }

    SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
    unsigned int nclusters[100] = {0};
    if (clustermap)
    {
      for (SvtxClusterMap::Iter iter = clustermap->begin();
           iter != clustermap->end();
           ++iter)
      {
        SvtxCluster* cluster = iter->second;
        ++nclusters[cluster->get_layer()];
      }
    }
    for (unsigned int ilayer = 0; ilayer < _nlayers_maps + _nlayers_intt + _nlayers_tpc; ++ilayer)
    {
      cout << "layer " << ilayer << ": nG4hits = " << ng4hits[ilayer]
           << " => nHits = " << nhits[ilayer]
           << " => nClusters = " << nclusters[ilayer] << endl;
    }

    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());

    cout << "nGtracks = " << std::distance(truthinfo->GetPrimaryParticleRange().first, truthinfo->GetPrimaryParticleRange().second);
    cout << " => nTracks = ";
    if (trackmap)
      cout << trackmap->size() << endl;
    else
      cout << 0 << endl;

    // cluster wise information
    if (Verbosity() > 1)
    {
      for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
           iter != g4hits.end();
           ++iter)
      {
        PHG4Hit* g4hit = *iter;

        cout << endl;
        cout << "===PHG4Hit===================================" << endl;
        cout << " PHG4Hit: ";
        g4hit->identify();

        std::set<SvtxCluster*> clusters = clustereval->all_clusters_from(g4hit);

        for (std::set<SvtxCluster*>::iterator jter = clusters.begin();
             jter != clusters.end();
             ++jter)
        {
          SvtxCluster* cluster = *jter;
          cout << "===Created-SvtxCluster================" << endl;
          cout << "SvtxCluster: ";
          cluster->identify();
        }
      }

      PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
           iter != range.second;
           ++iter)
      {
        PHG4Particle* particle = iter->second;

        // track-wise information
        cout << endl;

        cout << "=== Gtrack ===================================================" << endl;
        cout << " PHG4Particle id = " << particle->get_track_id() << endl;
        particle->identify();
        cout << " ptrue = (";
        cout.width(5);
        cout << particle->get_px();
        cout << ",";
        cout.width(5);
        cout << particle->get_py();
        cout << ",";
        cout.width(5);
        cout << particle->get_pz();
        cout << ")" << endl;

        cout << " vtrue = (";
        cout.width(5);
        cout << truthinfo->GetVtx(particle->get_vtx_id())->get_x();
        cout << ",";
        cout.width(5);
        cout << truthinfo->GetVtx(particle->get_vtx_id())->get_y();
        cout << ",";
        cout.width(5);
        cout << truthinfo->GetVtx(particle->get_vtx_id())->get_z();
        cout << ")" << endl;

        cout << " pt = " << sqrt(pow(particle->get_px(), 2) + pow(particle->get_py(), 2)) << endl;
        cout << " phi = " << atan2(particle->get_py(), particle->get_px()) << endl;
        cout << " eta = " << asinh(particle->get_pz() / sqrt(pow(particle->get_px(), 2) + pow(particle->get_py(), 2))) << endl;

        cout << " embed flag = " << truthinfo->isEmbeded(particle->get_track_id()) << endl;

        cout << " ---Associated-PHG4Hits-----------------------------------------" << endl;
        std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits(particle);
        for (std::set<PHG4Hit*>::iterator jter = g4hits.begin();
             jter != g4hits.end();
             ++jter)
        {
          PHG4Hit* g4hit = *jter;

          float x = 0.5 * (g4hit->get_x(0) + g4hit->get_x(1));
          float y = 0.5 * (g4hit->get_y(0) + g4hit->get_y(1));
          float z = 0.5 * (g4hit->get_z(0) + g4hit->get_z(1));

          cout << " #" << g4hit->get_hit_id() << " xtrue = (";
          cout.width(5);
          cout << x;
          cout << ",";
          cout.width(5);
          cout << y;
          cout << ",";
          cout.width(5);
          cout << z;
          cout << ")";

          std::set<SvtxCluster*> clusters = clustereval->all_clusters_from(g4hit);
          for (std::set<SvtxCluster*>::iterator kter = clusters.begin();
               kter != clusters.end();
               ++kter)
          {
            SvtxCluster* cluster = *kter;

            float x = cluster->get_x();
            float y = cluster->get_y();
            float z = cluster->get_z();

            cout << " => #" << cluster->get_id();
            cout << " xreco = (";
            cout.width(5);
            cout << x;
            cout << ",";
            cout.width(5);
            cout << y;
            cout << ",";
            cout.width(5);
            cout << z;
            cout << ")";
          }

          cout << endl;
        }

        if (trackmap && clustermap)
        {
          std::set<SvtxTrack*> tracks = trackeval->all_tracks_from(particle);
          for (std::set<SvtxTrack*>::iterator jter = tracks.begin();
               jter != tracks.end();
               ++jter)
          {
            SvtxTrack* track = *jter;

            float px = track->get_px();
            float py = track->get_py();
            float pz = track->get_pz();

            cout << "===Created-SvtxTrack==========================================" << endl;
            cout << " SvtxTrack id = " << track->get_id() << endl;
            cout << " preco = (";
            cout.width(5);
            cout << px;
            cout << ",";
            cout.width(5);
            cout << py;
            cout << ",";
            cout.width(5);
            cout << pz;
            cout << ")" << endl;
            cout << " quality = " << track->get_quality() << endl;
            cout << " nfromtruth = " << trackeval->get_nclusters_contribution(track, particle) << endl;

            cout << " ---Associated-SvtxClusters-to-PHG4Hits-------------------------" << endl;

            for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
                 iter != track->end_clusters();
                 ++iter)
            {
              unsigned int cluster_id = *iter;
              SvtxCluster* cluster = clustermap->get(cluster_id);

              float x = cluster->get_x();
              float y = cluster->get_y();
              float z = cluster->get_z();

              cout << " #" << cluster->get_id() << " xreco = (";
              cout.width(5);
              cout << x;
              cout << ",";
              cout.width(5);
              cout << y;
              cout << ",";
              cout.width(5);
              cout << z;
              cout << ") =>";

              PHG4Hit* g4hit = clustereval->max_truth_hit_by_energy(cluster);
              if ((g4hit) && (g4hit->get_trkid() == particle->get_track_id()))
              {
                x = 0.5 * (g4hit->get_x(0) + g4hit->get_x(1));
                y = 0.5 * (g4hit->get_y(0) + g4hit->get_y(1));
                z = 0.5 * (g4hit->get_z(0) + g4hit->get_z(1));

                cout << " #" << g4hit->get_hit_id()
                     << " xtrue = (";
                cout.width(5);
                cout << x;
                cout << ",";
                cout.width(5);
                cout << y;
                cout << ",";
                cout.width(5);
                cout << z;
                cout << ") => Gtrack id = " << g4hit->get_trkid();
              }
              else
              {
                cout << " noise hit";
              }
            }

            cout << endl;
          }
        }
      }
    }

    cout << endl;

  }  // if Verbosity()

  return;
}
*/

void TrkrEvaluator::fillOutputNtuples(PHCompositeNode* topNode)
{
  if (Verbosity() > 0) cout << "TrkrEvaluator::fillOutputNtuples() entered" << endl;

  float nhit_tpc_all = 0;
  float nhit_tpc_in = 0;
  float nhit_tpc_mid = 0;
  float nhit_tpc_out = 0;
  float nclus_all = 0;
  float nclus_tpc = 0;
  float nclus_intt = 0;
  float nclus_maps = 0;

  //------------------------
  // fill the Cluster NTuple
  //------------------------

  if (Verbosity() > 0)
  {
    cout << "check for ntp_cluster" << endl;
    _timer->restart();
  }

  TrkrHitTruthAssoc *hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");
  if(!hittruthassoc) 
    {
      cout << PHWHERE << "Failed to find TRKR_HITTRUTHASSOC node, quit!" << endl;
      exit(1);
    }

  TrkrClusterHitAssoc *clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode,"TRKR_CLUSTERHITASSOC");
  if(!clusterhitassoc) 
    {
      cout << PHWHERE << "Failed to find TRKR_CLUSTERHITASSOC node, quit!" << endl;
      exit(1);
    }

  TrkrClusterContainer *clustercontainer =  findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if(!clustercontainer)
  {
      cout << PHWHERE << "Failed to find TRKR_CLUSTER node, quit!" << endl;
      exit(1);
    }

  PHG4HitContainer *g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  PHG4HitContainer *g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  PHG4HitContainer *g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");

  TrkrClusterContainer::ConstRange clusrange = clustercontainer->getClusters();
  for(TrkrClusterContainer::ConstIterator clusiter = clusrange.first; clusiter != clusrange.second; ++clusiter)
    {
      TrkrCluster *cluster = clusiter->second;
      TrkrDefs::cluskey cluskey = clusiter->first;
      unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);
      
      float hitID = cluskey;
      float layer = TrkrDefs::getLayer(cluskey);
      float x = cluster->getPosition(0);
      float y = cluster->getPosition(1);
      float z = cluster->getPosition(2);
      TVector3 pos(x, y, z);
      float r = pos.Perp();
      float phi = pos.Phi();
      float eta = pos.Eta();
      float ex = sqrt(cluster->getError(0, 0));
      float ey = sqrt(cluster->getError(1, 1));
      float ez = sqrt(cluster->getError(2, 2));
      float ephi = cluster->getPhiError();
      float e = 0.0;	 
      float adc = cluster->getAdc();
      float phisize = cluster->getPhiSize(); 
      float zsize = cluster->getZSize();
      float size = 0;
 
      float trackID = NAN;
      
      float g4hitID = NAN;
      float gx = NAN;
      float gy = NAN;
      float gz = NAN;
      float gr = NAN;
      float gphi = NAN;
      float geta = NAN;
      float gt = NAN;
      float gtrackID = NAN;
      float gflavor = NAN;
      float gpx = NAN;
      float gpy = NAN;
      float gpz = NAN;
      float gvx = NAN;
      float gvy = NAN;
      float gvz = NAN;
      float gvt = NAN;
      float gfpx = NAN;
      float gfpy = NAN;
      float gfpz = NAN;
      float gfx = NAN;
      float gfy = NAN;
      float gfz = NAN;
      float gembed = NAN;
      float gprimary = NAN;
      
      float efromtruth = NAN;
      
      // get the hit keys associated with this cluster
      TrkrClusterHitAssoc::ConstRange hitrange = clusterhitassoc->getHits(cluskey);
      for(TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
	{
	  cout << " Cluster key " << cluskey << " trkrid " << trkrid << " hitkey " << clushititer->second << endl;

	  TrkrDefs::hitkey hitkey = clushititer->second;
	  // TrkrHitTruthAssoc uses a map with (hitsetkey, (hitkey, g4hitkey)) - get the hitsetkey from the cluskey
	  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);	  

	  std::set<PHG4Hit*> truth_hits;

	  // returns a range of pairs of (hitsetkey, (hitkey, g4hitkey))
	  TrkrHitTruthAssoc::ConstRange hrange =  hittruthassoc->getCells(hitsetkey, hitkey);
	  for(TrkrHitTruthAssoc::ConstIterator hittruthiter = hrange.first; hittruthiter != hrange.second; ++hittruthiter)
	    {
	      // extract the g4 hit key here and add the hits to a set
	      cout << " hittruth hitsetkey " << hittruthiter->first << " hitkey " << hittruthiter->second.first << " g4hitkey " << hittruthiter->second.second << endl;
	      PHG4HitDefs::keytype g4hitkey = hittruthiter->second.second;
	      PHG4Hit * g4hit;
	      if(trkrid == TrkrDefs::tpcId)
		{
		  g4hit = g4hits_tpc->findHit(g4hitkey);
		}
	      else if(trkrid == TrkrDefs::inttId)
		{
		  g4hit = g4hits_intt->findHit(g4hitkey);
		}
	      else
		{
		  g4hit = g4hits_mvtx->findHit(g4hitkey);
		}
	      truth_hits.insert(g4hit);	      
	    
	      if(trkrid == TrkrDefs::tpcId)
		{  
		  // This calculates the truth cluster position for the TPC from all of the contributing g4hits, typically 2-4 for the TPC
		  // Complicated, since only the part of the energy that is collected within a layer contributes to the position
		  //===============================================================================
		  
		  PHG4CylinderCellGeomContainer* geom_container =
		    findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
		  if (!geom_container)
		    {
		      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
		      return;
		    }
		  
		  PHG4CylinderCellGeom* GeoLayer = geom_container->GetLayerCellGeom(layer);
		  // get layer boundaries here (for nominal layer value) for later use
		  // radii of layer boundaries
		  float rbin = GeoLayer->get_radius() - GeoLayer->get_thickness() / 2.0;
		  float rbout = GeoLayer->get_radius() + GeoLayer->get_thickness() / 2.0;
		  
		  gx = 0.0;
		  gy = 0.0;
		  gz = 0.0;
		  gt = 0.0;
		  float gwt = 0.0;
		  
		  for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
		       iter != truth_hits.end();
		       ++iter)
		    {
		      PHG4Hit* this_g4hit = *iter;
		      
		      float rbegin = sqrt(this_g4hit->get_x(0) * this_g4hit->get_x(0) + this_g4hit->get_y(0) * this_g4hit->get_y(0));
		      float rend = sqrt(this_g4hit->get_x(1) * this_g4hit->get_x(1) + this_g4hit->get_y(1) * this_g4hit->get_y(1));
		      //cout << " Eval: g4hit " << this_g4hit->get_hit_id() <<  " rbegin " << rbegin << " rend " << rend << endl;
		      
		      // make sure the entry point is at lower radius
		      float xl[2];
		      float yl[2];
		      float zl[2];
		      
		      if (rbegin < rend)
			{
			  xl[0] = this_g4hit->get_x(0);
			  yl[0] = this_g4hit->get_y(0);
			  zl[0] = this_g4hit->get_z(0);
			  xl[1] = this_g4hit->get_x(1);
			  yl[1] = this_g4hit->get_y(1);
			  zl[1] = this_g4hit->get_z(1);
			}
		      else
			{
			  xl[0] = this_g4hit->get_x(1);
			  yl[0] = this_g4hit->get_y(1);
			  zl[0] = this_g4hit->get_z(1);
			  xl[1] = this_g4hit->get_x(0);
			  yl[1] = this_g4hit->get_y(0);
			  zl[1] = this_g4hit->get_z(0);
			  swap(rbegin, rend);
			  //cout << "swapped in and out " << endl;
			}
		      
		      // check that the g4hit is not completely outside the cluster layer. Just skip this g4hit if it is
		      // this can happen because an electron moves across a layer boundary during drift and readout
		      // so the g4hit is recorded in the cell as contributing to that layer, even though it was outside the boundaries
		      if (rbegin < rbin && rend < rbin)
			continue;
		      if (rbegin > rbout && rend > rbout)
			continue;
		      
		      float xin = xl[0];
		      float yin = yl[0];
		      float zin = zl[0];
		      float xout = xl[1];
		      float yout = yl[1];
		      float zout = zl[1];
		      
		      float t = NAN;
		      
		      if (rbegin < rbin)
			{
			  // line segment begins before boundary, find where it crosses
			  t = line_circle_intersection(xl, yl, zl, rbin);
			  if (t > 0)
			    {
			      xin = xl[0] + t * (xl[1] - xl[0]);
			      yin = yl[0] + t * (yl[1] - yl[0]);
			      zin = zl[0] + t * (zl[1] - zl[0]);
			    }
			}
		      
		      if (rend > rbout)
			{
			  // line segment ends after boundary, find where it crosses
			  t = line_circle_intersection(xl, yl, zl, rbout);
			  if (t > 0)
			    {
			      xout = xl[0] + t * (xl[1] - xl[0]);
			      yout = yl[0] + t * (yl[1] - yl[0]);
			      zout = zl[0] + t * (zl[1] - zl[0]);
			    }
			}
		      
		      // we want only the fraction of edep inside the layer
		      gx += (xin + xout) * 0.5 * this_g4hit->get_edep() * (xout - xin) / (xl[1] - xl[0]);
		      gy += (yin + yout) * 0.5 * this_g4hit->get_edep() * (yout - yin) / (yl[1] - yl[0]);
		      gz += (zin + zout) * 0.5 * this_g4hit->get_edep() * (zout - zin) / (zl[1] - zl[0]);
		      gt += this_g4hit->get_avg_t() * this_g4hit->get_edep() * (zout - zin) / (zl[1] - zl[0]);
		      gwt += this_g4hit->get_edep() * (zout - zin) / (zl[1] - zl[0]);
		    }  // loop over this_g4hit
		  gx /= gwt;
		  gy /= gwt;
		  gz /= gwt;
		  gt /= gwt;
		  //cout << " weighted means: gx " << gx << " gy " << gy << " gz " << gz << endl;
		}  // if TPC
	      else
		{
		  // not TPC, one g4hit per cluster
		  gx = g4hit->get_avg_x();
		  gy = g4hit->get_avg_y();
		  gz = g4hit->get_avg_z();
		}  // not TPC
	      
	      g4hitID = g4hit->get_hit_id();  // just takes the las one in the list
	      TVector3 gpos(gx, gy, gz);
	      gr = gpos.Perp();
	      gphi = gpos.Phi();
	      geta = gpos.Eta();
	    }    //  if (g4hit) {

	  float nparticles = 0;

	  float cluster_data[] = {(float) _ievent,
				  hitID,
				  x,
				  y,
				  z,
				  r,
				  phi,
				  eta,
				  ex,
				  ey,
				  ez,
				  ephi,
				  e,
				  adc,
				  layer,
				  size,
				  phisize,
				  zsize,
				  trackID,
				  g4hitID,
				  gx,
				  gy,
				  gz,
				  gr,
				  gphi,
				  geta,
				  gt,
				  gtrackID,
				  gflavor,
				  gpx,
				  gpy,
				  gpz,
				  gvx,
				  gvy,
				  gvz,
				  gvt,
				  gfpx,
				  gfpy,
				  gfpz,
				  gfx,
				  gfy,
				  gfz,
				  gembed,
				  gprimary,
				  efromtruth,
				  nparticles,
				  nhit_tpc_all,
				  nhit_tpc_in,
				  nhit_tpc_mid,
				  nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps};
	  
	  _ntp_cluster->Fill(cluster_data);
	}
    }  // loop over trkr clusters

  if (Verbosity() >= 1)
    {
      _timer->stop();
      cout << "cluster time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
  
  if (Verbosity() >= 1)
    {
      _timer->stop();
      cout << "g4hit time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }


  return;
}

float TrkrEvaluator::line_circle_intersection(float x[], float y[], float z[], float radius)
{
  // parameterize the line in terms of t (distance along the line segment, from 0-1) as
  // x = x0 + t * (x1-x0); y=y0 + t * (y1-y0); z = z0 + t * (z1-z0)
  // parameterize the cylinder (centered at x,y = 0,0) as  x^2 + y^2 = radius^2,   then
  // (x0 + t*(x1-z0))^2 + (y0+t*(y1-y0))^2 = radius^2
  // (x0^2 + y0^2 - radius^2) + (2x0*(x1-x0) + 2y0*(y1-y0))*t +  ((x1-x0)^2 + (y1-y0)^2)*t^2 = 0 = C + B*t + A*t^2
  // quadratic with:  A = (x1-x0)^2+(y1-y0)^2 ;  B = 2x0*(x1-x0) + 2y0*(y1-y0);  C = x0^2 + y0^2 - radius^2
  // solution: t = (-B +/- sqrt(B^2 - 4*A*C)) / (2*A)

  float A = (x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]);
  float B = 2.0 * x[0] * (x[1] - x[0]) + 2.0 * y[0] * (y[1] - y[0]);
  float C = x[0] * x[0] + y[0] * y[0] - radius * radius;
  float tup = (-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
  float tdn = (-B - sqrt(B * B - 4.0 * A * C)) / (2.0 * A);

  // The limits are 0 and 1, but we allow a little for floating point precision
  float t;
  if (tdn >= -0.0e-4 && tdn <= 1.0004)
    t = tdn;
  else if (tup >= -0.0e-4 && tup <= 1.0004)
    t = tup;
  else
  {
    cout << PHWHERE << "   **** Oops! No valid solution for tup or tdn, tdn = " << tdn << " tup = " << tup << endl;
    cout << "   radius " << radius << " rbegin " << sqrt(x[0] * x[0] + y[0] * y[0]) << " rend " << sqrt(x[1] * x[1] + y[1] * y[1]) << endl;
    cout << "   x0 " << x[0] << " x1 " << x[1] << endl;
    cout << "   y0 " << y[0] << " y1 " << y[1] << endl;
    cout << "   z0 " << z[0] << " z1 " << z[1] << endl;
    cout << "   A " << A << " B " << B << " C " << C << endl;

    t = -1;
  }

  return t;
}
