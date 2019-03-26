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

 _nlayers_maps = nlayers_maps;
 _nlayers_intt = nlayers_intt;
 _nlayers_tpc = nlayers_tpc;

}

int TrkrEvaluator::Init(PHCompositeNode* topNode)
{
_ievent = 0;

cout << "will write output to " << _filename.c_str() << endl;

  _tfile = new TFile(_filename.c_str(), "RECREATE");

_do_cluster_eval = true;
_do_track_eval = true;

 _scan_for_embeded = false;
_do_track_match = true;

  if (_do_cluster_eval) _ntp_cluster = new TNtuple("ntp_cluster", "svtxcluster => max truth",
                                                   "event:hitID:x:y:z:r:phi:eta:ex:ey:ez:ephi:"
                                                   "e:adc:layer:size:phisize:"
                                                   "zsize:trackID:g4hitID:gx:"
                                                   "gy:gz:gr:gphi:geta:gt:gtrackID:gflavor:"
                                                   "gpx:gpy:gpz:gvx:gvy:gvz:gvt:"
                                                   "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                                                   "gembed:gprimary:efromtruth:nparticles:"
                                                   "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps");

 if (_do_track_eval) _ntp_track = new TNtuple("ntp_track", "svtxtrack => max truth",
                                               "event:trackID:px:py:pz:pt:eta:phi:charge:"
                                               "quality:chisq:ndf:nhits:nmaps:nintt:ntpc:nlmaps:nlintt:nltpc:layers:"
                                               "dca2d:dca2dsigma:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:"
                                               "presdphi:presdeta:prese3x3:prese:"
                                               "cemcdphi:cemcdeta:cemce3x3:cemce:"
                                               "hcalindphi:hcalindeta:hcaline3x3:hcaline:"
                                               "hcaloutdphi:hcaloutdeta:hcaloute3x3:hcaloute:"
                                               "gtrackID:gflavor:gnhits:gnmaps:gnintt:gntpc:gnlmaps:gnlintt:gnltpc:"
                                               "gpx:gpy:gpz:gpt:geta:gphi:"
                                               "gvx:gvy:gvz:gvt:"
                                               "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                                               "gembed:gprimary:nfromtruth:nwrong:ntrumaps:ntruintt:ntrutpc:layersfromtruth:"
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
  if ((Verbosity() > 0) && (_ievent % 100 == 0))
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
  if (_ntp_track) _ntp_track->Write();

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
  if(_ntp_cluster)
    {
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

      // loop over all clusters
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
	  float ephi = r * cluster->getPhiError();
	  float e = 0.0;	 
	  float adc = cluster->getAdc();
	  float phisize = cluster->getPhiSize(); 
	  float zsize = cluster->getZSize();
	  float size = 0;
	  float trackID = NAN;
	  //if(layer < 7) cout << " new eval: reco cluster layer : "<<  layer << " x " << x << " y " << y << " z " << z << " phi " << phi << endl; 
      
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
	  std::set<PHG4Hit*> truth_hits;      
	  TrkrClusterHitAssoc::ConstRange hitrange = clusterhitassoc->getHits(cluskey);  // returns range of pairs {cluster key, hit key} for this cluskey
	  for(TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
	    {
	      TrkrDefs::hitkey hitkey = clushititer->second;
	      // TrkrHitTruthAssoc uses a map with (hitsetkey, std::pair(hitkey, g4hitkey)) - get the hitsetkey from the cluskey
	      TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);	  

	      // get all of the g4hits for this hitkey
	      std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;    
	      hittruthassoc->getG4Hits(hitsetkey, hitkey, temp_map); 	  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
	      for(std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter =  temp_map.begin(); htiter != temp_map.end(); ++htiter) 
		{
		  // extract the g4 hit key here and add the hits to the set
		  PHG4HitDefs::keytype g4hitkey = htiter->second.second;
		  PHG4Hit * g4hit;
		  if(trkrid == TrkrDefs::tpcId)
		    g4hit = g4hits_tpc->findHit(g4hitkey);
		  else if(trkrid == TrkrDefs::inttId)
		    g4hit = g4hits_intt->findHit(g4hitkey);
		  else
		    g4hit = g4hits_mvtx->findHit(g4hitkey);
		  truth_hits.insert(g4hit);	      
		} // end loop over g4hits associated with hitsetkey and hitkey
	    } // end loop over hits associated with cluskey

	  // we have the g4hits associated with this cluster, get the truth centroid
	  gx = 0.0;
	  gy = 0.0;
	  gz = 0.0;
	  gt = 0.0;
	  float gwt = 0.0;
      
	  if(trkrid == TrkrDefs::tpcId)
	    {  
	      // This calculates the truth cluster position for the TPC from all of the contributing g4hits, typically 2-4 for the TPC
	      // Complicated, since only the part of the energy that is collected within a layer contributes to the position
	      // We need the TPC geometry object to get the layer boundaries
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
	  
	      for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
		   iter != truth_hits.end();
		   ++iter)
		{
		  PHG4Hit* this_g4hit = *iter;
		  g4hitID = this_g4hit->get_hit_id(); // take the last one
		  float rbegin = sqrt(this_g4hit->get_x(0) * this_g4hit->get_x(0) + this_g4hit->get_y(0) * this_g4hit->get_y(0));
		  float rend = sqrt(this_g4hit->get_x(1) * this_g4hit->get_x(1) + this_g4hit->get_y(1) * this_g4hit->get_y(1));
	      
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
	      //cout << " new eval: truth cluster averages: layer " << layer  << " gx " << gx << " gy " << gy << " gz " << gz << " gwt " << gwt << endl;
	    }  // if TPC
	  else
	    {
	      // not TPC, entire g4hit is contained in one detector element 
	      for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
		   iter != truth_hits.end();
		   ++iter)
		{
		  PHG4Hit* this_g4hit = *iter;
		  g4hitID = this_g4hit->get_hit_id(); // take the ID of the last one
		  gx += this_g4hit->get_edep() * this_g4hit->get_avg_x();
		  gy += this_g4hit->get_edep() * this_g4hit->get_avg_y();
		  gz += this_g4hit->get_edep() * this_g4hit->get_avg_z();
		  gt += this_g4hit->get_edep() * this_g4hit->get_avg_t();
		  gwt += this_g4hit->get_edep();
		}
	      gx /= gwt;
	      gy /= gwt;
	      gz /= gwt;
	      gt /= gwt;
	      //cout << " new eval: truth cluster averages: layer " << layer  << " gx " << gx << " gy " << gy << " gz " << gz << " gwt " << gwt << endl;
	    }  // not TPC

	  // This occasionally returns nan for gphi and throws an error message - the cluster is not useful, skip it
	  if(isnan(gx) || isnan(gy) || isnan(gz)) continue;     

	  TVector3 gpos(gx, gy, gz);
	  gr = gpos.Perp();
	  gphi = gpos.Phi();
	  geta = gpos.Eta();

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
													  
    
	}  // loop over trkr clusters
												  
      if (Verbosity() >= 1)
	{
	  _timer->stop();
	  cout << "cluster time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
	}
    }

  //===============
  // Fill the track ntuple
  //===============

 if (_ntp_track)
  {
    if (Verbosity() > 0)
    {
      cout << "Filling ntp_track " << endl;
      _timer->restart();
    }

    PHG4HitContainer *g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
    PHG4HitContainer *g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
    PHG4HitContainer *g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");    
    
    PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    
    // need things off of the DST...
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    
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

    PHG4CylinderCellGeomContainer* geom_container =
      findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    if (!geom_container)
      {
	std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
	return;
      }
    
    
    if (trackmap)
    {
      for (SvtxTrackMap::Iter iter = trackmap->begin();
           iter != trackmap->end();
           ++iter)
      {

        SvtxTrack* track = iter->second;

        float trackID = track->get_id();
        float charge = track->get_charge();
        float quality = track->get_quality();
        float chisq = track->get_chisq();
        float ndf = track->get_ndf();
        float nhits = track->size_clusters();
        unsigned int layers = 0x0;
        int maps[_nlayers_maps];
        int intt[_nlayers_intt];
        int tpc[_nlayers_tpc];
        if (_nlayers_maps > 0)
        {
          for (unsigned int i = 0; i < _nlayers_maps; i++) maps[i] = 0;
        }
        if (_nlayers_intt > 0)
        {
          for (unsigned int i = 0; i < _nlayers_intt; i++) intt[i] = 0;
        }
        if (_nlayers_tpc > 0)
        {
          for (unsigned int i = 0; i < _nlayers_tpc; i++) tpc[i] = 0;
        }

        float nmaps = 0;
        float nintt = 0;
        float ntpc = 0;
        float nlmaps = 0;
        float nlintt = 0;
        float nltpc = 0;

        for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
             iter != track->end_cluster_keys();
             ++iter)
        {
	  TrkrDefs::cluskey cluster_key = *iter;
          unsigned int layer = TrkrDefs::getLayer(cluster_key);

          if (_nlayers_maps > 0 && layer < _nlayers_maps)
          {
            maps[layer] = 1;
            nmaps++;
          }
          if (_nlayers_intt > 0 && layer >= _nlayers_maps && layer < _nlayers_maps + _nlayers_intt)
          {
            intt[layer - _nlayers_maps] = 1;
            nintt++;
          }
          if (_nlayers_tpc > 0 && layer >= (_nlayers_maps + _nlayers_intt) && layer < (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
          {
            tpc[layer - (_nlayers_maps + _nlayers_intt)] = 1;
            ntpc++;
          }
        }  // end loop over clusters
        if (_nlayers_maps > 0)
          for (unsigned int i = 0; i < _nlayers_maps; i++) nlmaps += maps[i];
        if (_nlayers_intt > 0)
          for (unsigned int i = 0; i < _nlayers_intt; i++) nlintt += intt[i];
        if (_nlayers_tpc > 0)
          for (unsigned int i = 0; i < _nlayers_tpc; i++) nltpc += tpc[i];
        layers = nlmaps + nlintt + nltpc;
        float dca2d = track->get_dca2d();
        float dca2dsigma = track->get_dca2d_error();
        float dca3dxy = track->get_dca3d_xy();
        float dca3dxysigma = track->get_dca3d_xy_error();
        float dca3dz = track->get_dca3d_z();
        float dca3dzsigma = track->get_dca3d_z_error();
        float px = track->get_px();
        float py = track->get_py();
        float pz = track->get_pz();
        TVector3 v(px, py, pz);
        float pt = v.Pt();
        float eta = v.Eta();
        float phi = v.Phi();
        float pcax = track->get_x();
        float pcay = track->get_y();
        float pcaz = track->get_z();

        float presdphi = track->get_cal_dphi(SvtxTrack::PRES);
        float presdeta = track->get_cal_deta(SvtxTrack::PRES);
        float prese3x3 = track->get_cal_energy_3x3(SvtxTrack::PRES);
        float prese = track->get_cal_cluster_e(SvtxTrack::PRES);

        float cemcdphi = track->get_cal_dphi(SvtxTrack::CEMC);
        float cemcdeta = track->get_cal_deta(SvtxTrack::CEMC);
        float cemce3x3 = track->get_cal_energy_3x3(SvtxTrack::CEMC);
        float cemce = track->get_cal_cluster_e(SvtxTrack::CEMC);

        float hcalindphi = track->get_cal_dphi(SvtxTrack::HCALIN);
        float hcalindeta = track->get_cal_deta(SvtxTrack::HCALIN);
        float hcaline3x3 = track->get_cal_energy_3x3(SvtxTrack::HCALIN);
        float hcaline = track->get_cal_cluster_e(SvtxTrack::HCALIN);

        float hcaloutdphi = track->get_cal_dphi(SvtxTrack::HCALOUT);
        float hcaloutdeta = track->get_cal_deta(SvtxTrack::HCALOUT);
        float hcaloute3x3 = track->get_cal_energy_3x3(SvtxTrack::HCALOUT);
        float hcaloute = track->get_cal_cluster_e(SvtxTrack::HCALOUT);

        float gtrackID = NAN;
        float gflavor = NAN;
        float ng4hits = NAN;
        unsigned int ngmaps = 0;
        unsigned int ngintt = 0;
        unsigned int ngtpc = 0;
        unsigned int nglmaps = 0;
        unsigned int nglintt = 0;
        unsigned int ngltpc = 0;
        float gpx = NAN;
        float gpy = NAN;
        float gpt = NAN;
        float geta = NAN;
        float gphi = NAN;
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

        float nfromtruth = NAN;
        float nwrong = NAN;
        float ntrumaps = NAN;
        float ntruintt = NAN;
        float ntrutpc = NAN;
        float layersfromtruth = NAN;

        if (_do_track_match)
        {
	  // to collect g4hits for this track
	  std::set<PHG4Hit*> truth_hits;      
	  // multimap of particle track ID, cluster ID
	  std::set<int> g4trkid; 
	  std::multimap<int, TrkrDefs::cluskey> particle_cluster_map;
	  std::set<pair<int, unsigned int>> particle_layer_set;

	  // get cluster keys for this track
	  for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
	       iter != track->end_cluster_keys();
	       ++iter)
	    {
	      TrkrDefs::cluskey cluskey = *iter;
	      //cout << " cluster key " << cluskey << endl;
	      unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);
	      unsigned int layer = TrkrDefs::getLayer(cluskey);

	      // get hits for this cluster
	      TrkrClusterHitAssoc::ConstRange hitrange = clusterhitassoc->getHits(cluskey);  // returns range of pairs {cluster key, hit key} for this cluskey
	      for(TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
		{
		  TrkrDefs::hitkey hitkey = clushititer->second;
		  // TrkrHitTruthAssoc uses a map with (hitsetkey, std::pair(hitkey, g4hitkey)) - get the hitsetkey from the cluskey
		  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);	  
		  
		  // get all of the g4hits for this hitkey
		  std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;    
		  hittruthassoc->getG4Hits(hitsetkey, hitkey, temp_map); 	  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
		  for(std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter =  temp_map.begin(); 
		      htiter != temp_map.end(); ++htiter) 
		    {
		      // extract the g4 hit key here and add the g4hit to the set
		      PHG4HitDefs::keytype g4hitkey = htiter->second.second;
		      PHG4Hit * g4hit;
		      if(trkrid == TrkrDefs::tpcId)
			g4hit = g4hits_tpc->findHit(g4hitkey);
		      else if(trkrid == TrkrDefs::inttId)
			g4hit = g4hits_intt->findHit(g4hitkey);
		      else
			g4hit = g4hits_mvtx->findHit(g4hitkey);
		      truth_hits.insert(g4hit);	      
		      // map of g4hitid vs trkid
		      g4trkid.insert(g4hit->get_trkid());
		      particle_cluster_map.insert(std::make_pair( g4hit->get_trkid(), cluskey ));
		      particle_layer_set.insert(std::make_pair(g4hit->get_trkid(), layer));
		      //cout << "           g4trkid " << g4hit->get_trkid() << " cluskey " << cluskey << " layer " << layer << endl;
		    } // end loop over g4hits associated with hitsetkey and hitkey
		} // end loop over hits associated with cluskey
	    } // end loop over clusters associated with track

	  // we need to get the dominant g4particle
	  std::set<std::pair<int,int>> trkid_counts;
	  for(std::set<int>::iterator iter = g4trkid.begin(); iter != g4trkid.end(); ++iter)
	    {
	      int trkid = *iter;
	      int clus_count = particle_cluster_map.count(trkid);
	      trkid_counts.insert(std::make_pair(trkid, clus_count));
	      //cout << "g4trkid " << *iter << " clus_count " << clus_count << endl;	      
	    }

	  // get the id of the dominant particle for this track
	  int g4trkid_primary = -1000000;
	  int counts = -1;
	  for(std::set<std::pair<int, int>>::iterator it = trkid_counts.begin(); it!=trkid_counts.end(); ++it)
	    {
	      if(it->second > counts)
		{
		  counts = it->second;
		  g4trkid_primary = it->first;
		}
	    }

	  PHG4Particle* g4particle = truthinfo->GetParticle(g4trkid_primary);

          if (g4particle)
          {

            if (_scan_for_embeded)
            {
              if (truthinfo->isEmbeded(g4particle->get_track_id()) <= 0) continue;
            }

            gtrackID = g4particle->get_track_id();
            gflavor = g4particle->get_pid();

            ng4hits = truth_hits.size();
            gpx = g4particle->get_px();
            gpy = g4particle->get_py();
            gpz = g4particle->get_pz();


            int lmaps[_nlayers_maps + 1];
            if (_nlayers_maps > 0)
              for (unsigned int i = 0; i < _nlayers_maps; i++) lmaps[i] = 0;

            int lintt[_nlayers_intt + 1];
            if (_nlayers_intt > 0)
              for (unsigned int i = 0; i < _nlayers_intt; i++) lintt[i] = 0;

            int ltpc[_nlayers_tpc + 1];
            if (_nlayers_tpc > 0)
              for (unsigned int i = 0; i < _nlayers_tpc; i++) ltpc[i] = 0;


	    // get the hits info for the layers for this track

	    for(std::set<std::pair<int, unsigned int>>::iterator it = particle_layer_set.begin(); it != particle_layer_set.end(); ++it)
	      {
		int this_g4trkid = it->first;
		if(this_g4trkid != g4particle->get_track_id())  continue; 
		
		unsigned int this_layer = it->second;
		if (_nlayers_maps > 0 && this_layer < _nlayers_maps)
		  {
		    lmaps[this_layer] = 1;
		    ngmaps++;
		  }
		if (_nlayers_intt > 0 && this_layer >= _nlayers_maps && this_layer < _nlayers_maps + _nlayers_intt)
		  {
		    lintt[this_layer - _nlayers_maps] = 1;
		    ngintt++;
		  }
		
		if (_nlayers_tpc > 0 && this_layer >= _nlayers_maps + _nlayers_intt && this_layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
		  {
		    ltpc[this_layer - (_nlayers_maps + _nlayers_intt)] = 1;
		    ngtpc++;
		  }
	      }
            if (_nlayers_maps > 0)
              for (unsigned int i = 0; i < _nlayers_maps; i++) nglmaps += lmaps[i];
            if (_nlayers_intt > 0)
              for (unsigned int i = 0; i < _nlayers_intt; i++) nglintt += lintt[i];
            if (_nlayers_tpc > 0)
              for (unsigned int i = 0; i < _nlayers_tpc; i++) ngltpc += ltpc[i];

            TVector3 gv(gpx, gpy, gpz);
            gpt = gv.Pt();
            geta = gv.Eta();
            gphi = gv.Phi();
	    gvx = truthinfo->GetVtx(g4particle->get_vtx_id())->get_x();
	    gvy = truthinfo->GetVtx(g4particle->get_vtx_id())->get_y();
	    gvz = truthinfo->GetVtx(g4particle->get_vtx_id())->get_z();
	    gvt = truthinfo->GetVtx(g4particle->get_vtx_id())->get_t();

	    // find outermost truth hit in the tpc for this track
	    double max_radius = 0.0; 
	    PHG4Hit * temp_g4hit = nullptr;
	    PHG4Hit * outer_g4hit = nullptr;
	    for (PHG4HitContainer::ConstIterator g4iter = g4hits_tpc->getHits().first;
		 g4iter != g4hits_tpc->getHits().second;
		 ++g4iter)
	      {
		temp_g4hit = g4iter->second;

		if(temp_g4hit->get_trkid() != g4particle->get_track_id())  continue;
		
		// get radius and use that to get outer layer
		double this_radius = sqrt(pow(temp_g4hit->get_avg_x(), 2) + pow(temp_g4hit->get_avg_y(), 2));
		if(this_radius > max_radius)
		  {
		    max_radius = this_radius;
		    outer_g4hit = temp_g4hit;
		    cout << "  trkid " << temp_g4hit->get_trkid() << " this_radius " << this_radius << " max_radius " << max_radius << endl;
		  }
	      }
	    if(max_radius > 0.0)
	      {
		for(unsigned int layer = _nlayers_maps + _nlayers_intt; layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc; ++layer)
		  {
		    PHG4CylinderCellGeom* GeoLayer = geom_container->GetLayerCellGeom(layer);
		    // radii of layer boundaries
		    float rbin = GeoLayer->get_radius() - GeoLayer->get_thickness() / 2.0;
		    float rbout = GeoLayer->get_radius() + GeoLayer->get_thickness() / 2.0;		

		    if(max_radius > rbin && max_radius <= rbout)
		      {
			cout << " max_radius " << max_radius << " layer " << layer << " outer_g4hit trkid " << outer_g4hit->get_trkid() << endl;
		      }
		  }	  
	      }
	    // the outermost truth hit for this track
            if (outer_g4hit)
            {
              gfpx = outer_g4hit->get_px(1);
              gfpy = outer_g4hit->get_py(1);
              gfpz = outer_g4hit->get_pz(1);
              gfx = outer_g4hit->get_x(1);
              gfy = outer_g4hit->get_y(1);
              gfz = outer_g4hit->get_z(1);
            }

            gembed = truthinfo->isEmbeded(g4particle->get_track_id());
	    gprimary = 0;
	    if(g4particle->get_parent_id() == 0)
	      gprimary = 1;

	    /*
            nfromtruth = trackeval->get_nclusters_contribution(track, g4particle);
            nwrong = trackeval->get_nwrongclusters_contribution(track, g4particle);
            if (_nlayers_maps == 0)
            {
              ntrumaps = 0;
            }
            else
            {
              ntrumaps = trackeval->get_layer_range_contribution(track, g4particle, 0, _nlayers_maps);
            }
            if (_nlayers_intt == 0)
            {
              ntruintt = 0;
            }
            else
            {
              ntruintt = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps, _nlayers_maps + _nlayers_intt);
            }
            ntrutpc = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt, _nlayers_maps + _nlayers_intt + _nlayers_tpc);
            layersfromtruth = trackeval->get_nclusters_contribution_by_layer(track, g4particle);
	    */
          }
	} // end if for _do_track_match

        float track_data[] = {(float) _ievent,
                              trackID,
                              px,
                              py,
                              pz,
                              pt,
                              eta,
                              phi,
                              charge,
                              quality,
                              chisq,
                              ndf,
                              nhits, nmaps, nintt, ntpc, nlmaps, nlintt, nltpc,
                              (float) layers,
                              dca2d,
                              dca2dsigma,
                              dca3dxy,
                              dca3dxysigma,
                              dca3dz,
                              dca3dzsigma,
                              pcax,
                              pcay,
                              pcaz,
                              presdphi,
                              presdeta,
                              prese3x3,
                              prese,
                              cemcdphi,
                              cemcdeta,
                              cemce3x3,
                              cemce,
                              hcalindphi,
                              hcalindeta,
                              hcaline3x3,
                              hcaline,
                              hcaloutdphi,
                              hcaloutdeta,
                              hcaloute3x3,
                              hcaloute,
                              gtrackID,
                              gflavor,
                              ng4hits,
                              (float) ngmaps,
                              (float) ngintt,
                              (float) ngtpc,
                              (float) nglmaps,
                              (float) nglintt,
                              (float) ngltpc,
                              gpx,
                              gpy,
                              gpz,
                              gpt,
                              geta,
                              gphi,
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
                              nfromtruth,
                              nwrong,
                              ntrumaps,
                              ntruintt,
                              ntrutpc,
                              layersfromtruth,
                              nhit_tpc_all,
                              nhit_tpc_in,
                              nhit_tpc_mid,
                              nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps};

        /*
	cout << "ievent " << _ievent
	     << " trackID " << trackID
	     << " nhits " << nhits
	     << " px " << px
	     << " py " << py
	     << " pz " << pz
	     << " gembed " << gembed
	     << " gprimary " << gprimary 
	     << endl;
	*/
        _ntp_track->Fill(track_data);
      }
    }
    if (Verbosity() >= 1)
    {
      _timer->stop();
      cout << "track time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
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
