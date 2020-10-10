#include "SvtxEvaluator.h"

#include "SvtxEvalStack.h"

#include "SvtxClusterEval.h"
#include "SvtxHitEval.h"
#include "SvtxTrackEval.h"
#include "SvtxTruthEval.h"
#include "SvtxVertexEval.h"

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrDefs.h>

#include <tpc/TpcDefs.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>                                       // for shared_ptr
#include <set>                                          // for _Rb_tree_cons...
#include <utility>
#include <vector>

using namespace std;

SvtxEvaluator::SvtxEvaluator(const string& name, const string& filename, const string& trackmapname,
                             unsigned int nlayers_maps,
                             unsigned int nlayers_intt,
                             unsigned int nlayers_tpc)
  : SubsysReco("SvtxEvaluator")
  , _ievent(0)
  , _iseed(0)
  , m_fSeed(NAN)
  , _svtxevalstack(nullptr)
  , _strict(false)
  , _use_initial_vertex(false)
  , _errors(0)
  , _do_info_eval(true)
  , _do_vertex_eval(true)
  , _do_gpoint_eval(true)
  , _do_g4hit_eval(true)
  , _do_hit_eval(true)
  , _do_cluster_eval(true)
  , _do_g4cluster_eval(true)
  , _do_gtrack_eval(true)
  , _do_track_eval(true)
  , _do_gseed_eval(false)
  , _do_track_match(true)
  , _do_eval_light(true)
  , _scan_for_embedded(false)
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _ntp_info(nullptr)
  , _ntp_vertex(nullptr)
  , _ntp_gpoint(nullptr)
  , _ntp_g4hit(nullptr)
  , _ntp_hit(nullptr)
  , _ntp_cluster(nullptr)
  , _ntp_g4cluster(nullptr)
  , _ntp_gtrack(nullptr)
  , _ntp_track(nullptr)
  , _ntp_gseed(nullptr)
  , _filename(filename)
  , _trackmapname(trackmapname)
  , _tfile(nullptr)
  , _timer(nullptr)
{
}

SvtxEvaluator::~SvtxEvaluator()
{
  delete _timer;
}

int SvtxEvaluator::Init(PHCompositeNode* topNode)
{
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");
  if (_do_info_eval) _ntp_info = new TNtuple("ntp_info", "event info",
                                                 "event:seed:"
					         "occ11:occ116:occ21:occ216:occ31:occ316:"
                                                 "gntrkall:gntrkprim:ntrk:"
                                                 "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps");

  if (_do_vertex_eval) _ntp_vertex = new TNtuple("ntp_vertex", "vertex => max truth",
                                                 "event:seed:vx:vy:vz:ntracks:"
                                                 "gvx:gvy:gvz:gvt:gembed:gntracks:gntracksmaps:"
                                                 "gnembed:nfromtruth:"
                                                 "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps");

  if (_do_gpoint_eval) _ntp_gpoint = new TNtuple("ntp_gpoint", "g4point => best vertex",
                                                 "event:seed:gvx:gvy:gvz:gvt:gntracks:gembed:"
                                                 "vx:vy:vz:ntracks:"
                                                 "nfromtruth:"
                                                 "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps");

  if (_do_g4hit_eval) _ntp_g4hit = new TNtuple("ntp_g4hit", "g4hit => best svtxcluster",
                                               "event:seed:g4hitID:gx:gy:gz:gt:gedep:geta:gphi:"
                                               "gdphi:gdz:"
                                               "glayer:gtrackID:gflavor:"
                                               "gpx:gpy:gpz:"
                                               "gvx:gvy:gvz:"
                                               "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                                               "gembed:gprimary:nclusters:"
                                               "clusID:x:y:z:eta:phi:e:adc:layer:size:"
                                               "phisize:zsize:efromtruth:dphitru:detatru:dztru:drtru:"
                                               "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps");

  if (_do_hit_eval) _ntp_hit = new TNtuple("ntp_hit", "svtxhit => max truth",
                                           "event:seed:hitID:e:adc:layer:"
                                           "cellID:ecell:phibin:zbin:phi:z:"
                                           "g4hitID:gedep:gx:gy:gz:gt:"
                                           "gtrackID:gflavor:"
                                           "gpx:gpy:gpz:gvx:gvy:gvz:gvt:"
                                           "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                                           "gembed:gprimary:efromtruth:"
                                           "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps");

  if (_do_cluster_eval) _ntp_cluster = new TNtuple("ntp_cluster", "svtxcluster => max truth",
                                                   "event:seed:hitID:x:y:z:r:phi:eta:theta:ex:ey:ez:ephi:"
                                                   "e:adc:layer:size:phisize:"
                                                   "zsize:trackID:g4hitID:gx:"
                                                   "gy:gz:gr:gphi:geta:gt:gtrackID:gflavor:"
                                                   "gpx:gpy:gpz:gvx:gvy:gvz:gvt:"
                                                   "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                                                   "gembed:gprimary:efromtruth:nparticles:"
                                                   "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps");

  if (_do_g4cluster_eval) _ntp_g4cluster = new TNtuple("ntp_g4cluster", "g4cluster => max truth",
						       "event:layer:gx:gy:gz:gt:gedep:gr:gphi:geta:ng4hits:gtrackID:gflavor:gembed:gprimary:gphisize:gzsize:nreco:x:y:z:r:phi:eta:ex:ey:ez:ephi:phisize:zsize:adc"); 
                                                       
  if (_do_gtrack_eval) _ntp_gtrack = new TNtuple("ntp_gtrack", "g4particle => best svtxtrack",
                                                 "event:seed:gntracks:gtrackID:gflavor:gnhits:gnmaps:gnintt:"
                                                 "gnintt1:gnintt2:gnintt3:gnintt4:"
                                                 "gnintt5:gnintt6:gnintt7:gnintt8:"
                                                 "gntpc:gnlmaps:gnlintt:gnltpc:"
                                                 "gpx:gpy:gpz:gpt:geta:gphi:"
                                                 "gvx:gvy:gvz:gvt:"
                                                 "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                                                 "gembed:gprimary:"
                                                 "trackID:px:py:pz:pt:eta:phi:deltapt:deltaeta:deltaphi:"
                                                 "charge:quality:chisq:ndf:nhits:layers:nmaps:nintt:ntpc:ntpc1:ntpc11:ntpc2:ntpc3:nlmaps:nlintt:nltpc:"
                                                 "dca2d:dca2dsigma:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:nfromtruth:nwrong:ntrumaps:ntruintt:ntrutpc:ntrutpc1:ntrutpc11:ntrutpc2:ntrutpc3:layersfromtruth:"
                                                 "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps");

  if (_do_track_eval) _ntp_track = new TNtuple("ntp_track", "svtxtrack => max truth",
                                               "event:seed:trackID:px:py:pz:pt:eta:phi:deltapt:deltaeta:deltaphi:charge:"
                                               "quality:chisq:ndf:nhits:nmaps:nintt:ntpc:ntpc1:ntpc11:ntpc2:ntpc3:nlmaps:nlintt:nltpc:layers:"
                                               "dca2d:dca2dsigma:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:"
                                               "presdphi:presdeta:prese3x3:prese:"
                                               "cemcdphi:cemcdeta:cemce3x3:cemce:"
                                               "hcalindphi:hcalindeta:hcaline3x3:hcaline:"
                                               "hcaloutdphi:hcaloutdeta:hcaloute3x3:hcaloute:"
                                               "gtrackID:gflavor:gnhits:gnmaps:gnintt:gntpc:gnlmaps:gnlintt:gnltpc:"
                                               "gpx:gpy:gpz:gpt:geta:gphi:"
                                               "gvx:gvy:gvz:gvt:"
                                               "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                                               "gembed:gprimary:nfromtruth:nwrong:ntrumaps:ntruintt:"
					       "ntrutpc:ntrutpc1:ntrutpc11:ntrutpc2:ntrutpc3:layersfromtruth:"
                                               "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps");

  if (_do_gseed_eval) _ntp_gseed = new TNtuple("ntp_gseed", "seeds from truth",
                                               "event:seed:ntrk:gx:gy:gz:gr:geta:gphi:"
                                               "glayer:"
                                               "gpx:gpy:gpz:gtpt:gtphi:gteta:"
                                               "gvx:gvy:gvz:"
                                               "gembed:gprimary:gflav:"
                                               "dphiprev:detaprev:"
                                               "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps");

  _timer = new PHTimer("_eval_timer");
  _timer->stop();

  return Fun4AllReturnCodes::EVENT_OK;
}

int SvtxEvaluator::InitRun(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int SvtxEvaluator::process_event(PHCompositeNode* topNode)
{
  if ((Verbosity() > 1) && (_ievent % 100 == 0))
  {
    cout << "SvtxEvaluator::process_event - Event = " << _ievent << endl;
  }
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED"))
  {
    _iseed = rc->get_IntFlag("RANDOMSEED");
    m_fSeed = _iseed;
  }
  else
  {
    _iseed = 0;
    m_fSeed = NAN;
  }

  if(Verbosity() > 1)
    {
      cout << "SvtxEvaluator::process_event - Seed = " << _iseed << endl;
    }

  if (!_svtxevalstack)
  {
    _svtxevalstack = new SvtxEvalStack(topNode);
    _svtxevalstack->set_strict(_strict);
    _svtxevalstack->set_verbosity(Verbosity());
    _svtxevalstack->set_use_initial_vertex(_use_initial_vertex);
  }
  else
  {
    _svtxevalstack->next_event(topNode);
  }

  //-----------------------------------
  // print what is coming into the code
  //-----------------------------------

  printInputInfo(topNode);

  //---------------------------
  // fill the Evaluator NTuples
  //---------------------------

  fillOutputNtuples(topNode);

  //--------------------------------------------------
  // Print out the ancestry information for this event
  //--------------------------------------------------

  printOutputInfo(topNode);

  ++_ievent;
  return Fun4AllReturnCodes::EVENT_OK;
}

int SvtxEvaluator::End(PHCompositeNode* topNode)
{
  _tfile->cd();

  if (_ntp_info) _ntp_info->Write();
  if (_ntp_vertex) _ntp_vertex->Write();
  if (_ntp_gpoint) _ntp_gpoint->Write();
  if (_ntp_g4hit) _ntp_g4hit->Write();
  if (_ntp_hit) _ntp_hit->Write();
  if (_ntp_cluster) _ntp_cluster->Write();
  if (_ntp_g4cluster) _ntp_g4cluster->Write();
  if (_ntp_gtrack) _ntp_gtrack->Write();
  if (_ntp_track) _ntp_track->Write();
  if (_ntp_gseed) _ntp_gseed->Write();

  _tfile->Close();

  delete _tfile;

  if (Verbosity() > 1)
  {
    cout << "========================= SvtxEvaluator::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  if(_svtxevalstack)
    {
      _errors += _svtxevalstack->get_errors();
      
      if (Verbosity() > 1)
	{
	  if ((_errors > 0) || (Verbosity() > 1))
	    {
	      cout << "SvtxEvaluator::End() - Error Count: " << _errors << endl;
	    }
	}

      delete _svtxevalstack;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void SvtxEvaluator::printInputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 1) cout << "SvtxEvaluator::printInputInfo() entered" << endl;

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
    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if (clustermap)
    {
      unsigned int icluster = 0;
      TrkrClusterContainer::ConstRange all_clusters = clustermap->getClusters();
      for (TrkrClusterContainer::ConstIterator iter = all_clusters.first;
           iter != all_clusters.second;
           ++iter)
      {
	TrkrDefs::cluskey cluster_key = iter->first;
        cout << icluster << " with key " << cluster_key << " of " << clustermap->size();
        cout << ": SvtxCluster: " << endl;
        iter->second->identify();
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
	++itrack;
      }
    }

    cout << "---SVXVERTEXES-------------" << endl;
    SvtxVertexMap* vertexmap = nullptr;
    if(_use_initial_vertex)
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    else
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapRefit");
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

void SvtxEvaluator::printOutputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 1) cout << "SvtxEvaluator::printOutputInfo() entered" << endl;

  //==========================================
  // print out some useful stuff for debugging
  //==========================================

  if (Verbosity() > 1)
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

    SvtxVertexMap* vertexmap = nullptr;
    if(_use_initial_vertex)
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    else
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapRefit");
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

    TrkrHitSetContainer *hitsetmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    unsigned int nhits[100] = {0};
    if (hitsetmap)
    {
      TrkrHitSetContainer::ConstRange all_hitsets = hitsetmap->getHitSets();
      for (TrkrHitSetContainer::ConstIterator hitsetiter = all_hitsets.first;
           hitsetiter != all_hitsets.second;
           ++hitsetiter)
      {
	// we have a single hitset, get the layer
	unsigned int layer = TrkrDefs::getLayer(hitsetiter->first);
	
	// count all hits in this hitset
	TrkrHitSet::ConstRange hitrangei = hitsetiter->second->getHits();
	for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	     hitr != hitrangei.second;
	     ++hitr)
	  {
	    ++nhits[layer];
	  }
      }
    }

    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    unsigned int nclusters[100] = {0};
    if (clustermap)
    {
      TrkrClusterContainer::ConstRange all_clusters = clustermap->getClusters();
      for (TrkrClusterContainer::ConstIterator iter = all_clusters.first;
           iter != all_clusters.second;
           ++iter)
      {
	TrkrDefs::cluskey cluster_key = iter->first;
	unsigned int layer = TrkrDefs::getLayer(cluster_key);
        ++nclusters[layer];
      }
    }

    PHG4CylinderCellGeomContainer* geom_container =
      findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    if (!geom_container)
      {
	std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
	return;
      }
    


    for (unsigned int ilayer = 0; ilayer < _nlayers_maps + _nlayers_intt + _nlayers_tpc; ++ilayer)
    {
      PHG4CylinderCellGeom* GeoLayer = geom_container->GetLayerCellGeom(ilayer);

      cout << "layer " << ilayer << ": nG4hits = " << ng4hits[ilayer]
           << " => nHits = " << nhits[ilayer]
           << " => nClusters = " << nclusters[ilayer] 	   << endl;
      if(ilayer>=_nlayers_maps + _nlayers_intt){
      cout << "layer " << ilayer
	   << " => nphi = " << GeoLayer->get_phibins()
	   << " => nz   = " << GeoLayer->get_zbins()
	   << " => ntot = " << GeoLayer->get_phibins() * GeoLayer->get_zbins()
	   << endl;
      }
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

        std::set<TrkrDefs::cluskey> clusters = clustereval->all_clusters_from(g4hit);

        for (std::set<TrkrDefs::cluskey>::iterator jter = clusters.begin();
             jter != clusters.end();
             ++jter)
        {
	  TrkrDefs::cluskey cluster_key = *jter;
          cout << "===Created-SvtxCluster================" << endl;
          cout << "SvtxCluster: ";
          TrkrCluster *cluster = clustermap->findCluster(cluster_key);
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

          std::set<TrkrDefs::cluskey> clusters = clustereval->all_clusters_from(g4hit);
          for (std::set<TrkrDefs::cluskey>::iterator kter = clusters.begin();
               kter != clusters.end();
               ++kter)
	    {
	      TrkrDefs::cluskey cluster_key = *kter;
	      
	      TrkrCluster *cluster = clustermap->findCluster(cluster_key);
	      float x = cluster->getX();
	      float y = cluster->getY();
	      float z = cluster->getZ();
	      
	      cout << " => #" << cluster_key;
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

            for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
                 iter != track->end_cluster_keys();
                 ++iter)
            {
	      TrkrDefs::cluskey cluster_key = *iter;
              TrkrCluster *cluster = clustermap->findCluster(cluster_key);

              float x = cluster->getX();
              float y = cluster->getY();
              float z = cluster->getZ();

              cout << " #" << cluster_key << " xreco = (";
              cout.width(5);
              cout << x;
              cout << ",";
              cout.width(5);
              cout << y;
              cout << ",";
              cout.width(5);
              cout << z;
              cout << ") =>";

              PHG4Hit* g4hit = clustereval->max_truth_hit_by_energy(cluster_key);
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

void SvtxEvaluator::fillOutputNtuples(PHCompositeNode* topNode)
{
  if (Verbosity() > 1) cout << "SvtxEvaluator::fillOutputNtuples() entered" << endl;
  
  SvtxVertexEval* vertexeval = _svtxevalstack->get_vertex_eval();
  SvtxClusterEval* clustereval = _svtxevalstack->get_cluster_eval();
  SvtxTrackEval* trackeval = _svtxevalstack->get_track_eval();
  SvtxHitEval* hiteval = _svtxevalstack->get_hit_eval();
  SvtxTruthEval* trutheval = _svtxevalstack->get_truth_eval();

  float nhit_tpc_all = 0;
  float nhit_tpc_in = 0;
  float nhit_tpc_mid = 0;
  float nhit_tpc_out = 0;
  float nclus_all = 0;
  float nclus_tpc = 0;
  float nclus_intt = 0;
  float nclus_maps = 0;
  float nhit[100];
  for(int i = 0; i<100;i++)nhit[i] = 0;
  float occ11  = 0;
  float occ116 = 0;
  float occ21  = 0;
  float occ216 = 0;
  float occ31  = 0;
  float occ316 = 0;

  TrkrHitSetContainer* hitmap_in = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (hitmap_in)
    {
      TrkrHitSetContainer::ConstRange all_hitsets = hitmap_in->getHitSets();
      for (TrkrHitSetContainer::ConstIterator hitsetiter = all_hitsets.first;
	   hitsetiter != all_hitsets.second;
	   ++hitsetiter)
	{
	  // we have a single hitset, get the layer
	  unsigned int layer = TrkrDefs::getLayer(hitsetiter->first);
	  if(layer >= _nlayers_maps + _nlayers_intt)
	    {	
	      // count all hits in this hitset
	      TrkrHitSet::ConstRange hitrangei = hitsetiter->second->getHits();
	      for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
		   hitr != hitrangei.second;
		   ++hitr)
		{
		  nhit[layer]++;
		  nhit_tpc_all++;
		  if ((float) layer == _nlayers_maps + _nlayers_intt) nhit_tpc_in++;
		  if ((float) layer == _nlayers_maps + _nlayers_intt + _nlayers_tpc - 1) nhit_tpc_out++;
		  if ((float) layer == _nlayers_maps + _nlayers_intt + _nlayers_tpc / 2 - 1) nhit_tpc_mid++;
		}
	    }
	}
    }
  /**********/
  PHG4CylinderCellGeomContainer* geom_container =
    findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return;
    }
  PHG4CylinderCellGeom* GeoLayer;
  int layer = _nlayers_maps + _nlayers_intt;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  int nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  float nhits = nhit[layer];
  occ11  = nhits/nbins;
  if(Verbosity() > 1)
    {
      cout << " occ11 = " << occ11
	   << " nbins = " << nbins
	   << " nhits = " << nhits
	   << " layer = " << layer
	   << endl;
    }
  layer = _nlayers_maps + _nlayers_intt + 15;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  occ116 = nhits/nbins;

  layer = _nlayers_maps + _nlayers_intt + 16;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  occ21  = nhits/nbins;
  layer = _nlayers_maps + _nlayers_intt + 31;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  occ216 = nhits/nbins;
  layer = _nlayers_maps + _nlayers_intt + 32;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  occ31  = nhits/nbins;
  layer = _nlayers_maps + _nlayers_intt + 47;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  occ316 = nhits/nbins;

  if(Verbosity() > 1)
    {
      cout << " occ11 = " << occ11
	   << " occ116 = " << occ116
	   << " occ21 = " << occ21
	   << " occ216 = " << occ216
	   << " occ31 = " << occ31
	   << " occ316 = " << occ316
	   << endl;
    }
  TrkrClusterContainer* clustermap_in = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  nclus_all = clustermap_in->size();
  TrkrClusterContainer::ConstRange all_clusters = clustermap_in->getClusters();
  for (TrkrClusterContainer::ConstIterator iter_cin = all_clusters.first;
       iter_cin != all_clusters.second;
       ++iter_cin)
  {
    TrkrDefs::cluskey cluster_key = iter_cin->first;
    unsigned int layer = TrkrDefs::getLayer(cluster_key);
    if (_nlayers_maps > 0)
      if (layer < _nlayers_maps) nclus_maps++;
    if (_nlayers_intt > 0)
      if (layer >= _nlayers_maps && layer < _nlayers_maps + _nlayers_intt) nclus_intt++;
    if (_nlayers_tpc > 0)
      if (layer >= _nlayers_maps + _nlayers_intt) nclus_tpc++;
  }

  //-----------------------
  // fill the info NTuple
  //-----------------------
  if (_ntp_info)
    {
      if (Verbosity() > 0)
	{
	  cout << "Filling ntp_info " << endl;
    }
      float ntrk = 0;
      SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());
      if (trackmap)
	ntrk = (float) trackmap->size();
      PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
      float info_data[] = {(float) _ievent,m_fSeed,
			   occ11,occ116,occ21,occ216,occ31,occ316,
			   (float)truthinfo->GetNumPrimaryVertexParticles(),
			   0,
			   ntrk,
			   nhit_tpc_all,
			   nhit_tpc_in,
			   nhit_tpc_mid,
			   nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps};

        _ntp_info->Fill(info_data);
      }
  //-----------------------
  // fill the Vertex NTuple
  //-----------------------
  bool doit = true;
  if (_ntp_vertex && doit)
  {
    if (Verbosity() > 1)
    {
      cout << "Filling ntp_vertex " << endl;
      cout << "start vertex time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
      _timer->restart();
    }

    SvtxVertexMap* vertexmap = nullptr;
    if(_use_initial_vertex)
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    else
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapRefit");
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (vertexmap && truthinfo)
    {
      const auto prange = truthinfo->GetPrimaryParticleRange();
      map<int, unsigned int> embedvtxid_particle_count;
      map<int, unsigned int> embedvtxid_maps_particle_count;
      map<int, unsigned int> vertex_particle_count;

      if (_do_eval_light == false)
      {
        for (auto iter = prange.first; iter != prange.second; ++iter)  // process all primary paricle
        {
          const int point_id = iter->second->get_vtx_id();
          int gembed = truthinfo->isEmbededVtx(iter->second->get_vtx_id());
          ++vertex_particle_count[point_id];
          ++embedvtxid_particle_count[gembed];
          PHG4Particle* g4particle = iter->second;

          if (_scan_for_embedded && gembed <= 0) continue;

          std::set<TrkrDefs::cluskey> g4clusters = clustereval->all_clusters_from(g4particle);
          unsigned int nglmaps = 0;

          int lmaps[_nlayers_maps + 1];
          if (_nlayers_maps > 0)
          {
            for (unsigned int i = 0; i < _nlayers_maps; i++)
            {
              lmaps[i] = 0;
            }
          }
          for (const TrkrDefs::cluskey g4cluster : g4clusters)
          {
            unsigned int layer = TrkrDefs::getLayer(g4cluster);
            //cout<<__LINE__<<": " << _ievent <<": " <<gtrackID << ": " << layer <<": " <<g4cluster->get_id() <<endl;
            if (_nlayers_maps > 0 && layer < _nlayers_maps)
            {
              lmaps[layer] = 1;
            }
          }
          if (_nlayers_maps > 0)
          {
            for (unsigned int i = 0; i < _nlayers_maps; i++)
            {
              nglmaps += lmaps[i];
            }
          }
          float gpx = g4particle->get_px();
          float gpy = g4particle->get_py();
          float gpz = g4particle->get_pz();
          float gpt = NAN;
          float geta = NAN;

          if (gpx != 0 && gpy != 0)
          {
            TVector3 gv(gpx, gpy, gpz);
            gpt = gv.Pt();
            geta = gv.Eta();
            //          gphi = gv.Phi();
          }

          if (nglmaps == 3 && fabs(geta) < 1.0 && gpt > 0.5)
            ++embedvtxid_maps_particle_count[gembed];
        }
      }

      auto vrange = truthinfo->GetPrimaryVtxRange();
      map<int, bool> embedvtxid_found;
      map<int, int> embedvtxid_vertex_id;
      map<int, PHG4VtxPoint*> embedvtxid_vertex;
      for (auto iter = vrange.first; iter != vrange.second; ++iter)  // process all primary vertexes
      {
        const int point_id = iter->first;
        int gembed = truthinfo->isEmbededVtx(point_id);
        if (_scan_for_embedded && gembed <= 0) continue;

        auto search = embedvtxid_found.find(gembed);
        if (search != embedvtxid_found.end())
        {
          embedvtxid_vertex_id[gembed] = point_id;
          embedvtxid_vertex[gembed] = iter->second;
        }
        else
        {
          if (vertex_particle_count[embedvtxid_vertex_id[gembed]] < vertex_particle_count[point_id])
          {
            embedvtxid_vertex_id[gembed] = point_id;
            embedvtxid_vertex[gembed] = iter->second;
          }
        }
        embedvtxid_found[gembed] = false;
      }

      unsigned int ngembed = 0;
      for (std::map<int, bool>::iterator iter = embedvtxid_found.begin();
           iter != embedvtxid_found.end();
           ++iter)
      {
        if (iter->first >= 0 || iter->first != iter->first) continue;
        ++ngembed;
      }

      for (SvtxVertexMap::Iter iter = vertexmap->begin();
           iter != vertexmap->end();
           ++iter)
      {
        SvtxVertex* vertex = iter->second;
        PHG4VtxPoint* point = vertexeval->max_truth_point_by_ntracks(vertex);
        float vx = vertex->get_x();
        float vy = vertex->get_y();
        float vz = vertex->get_z();
        float ntracks = vertex->size_tracks();
        float gvx = NAN;
        float gvy = NAN;
        float gvz = NAN;
        float gvt = NAN;
        float gembed = NAN;
        float gntracks = truthinfo->GetNumPrimaryVertexParticles();
        float gntracksmaps = NAN;
        float gnembed = NAN;
        float nfromtruth = NAN;
        if (point)
        {
          const int point_id = point->get_id();
          gvx = point->get_x();
          gvy = point->get_y();
          gvz = point->get_z();
          gvt = point->get_t();
          gembed = truthinfo->isEmbededVtx(point_id);
          gntracks = embedvtxid_particle_count[(int) gembed];
          if (embedvtxid_maps_particle_count[(int) gembed] > 0 && fabs(gvt) < 2000. && fabs(gvz) < 13.0)
            gntracksmaps = embedvtxid_maps_particle_count[(int) gembed];
          gnembed = (float) ngembed;
          nfromtruth = vertexeval->get_ntracks_contribution(vertex, point);
          embedvtxid_found[(int) gembed] = true;
        }
	
        float vertex_data[] = {(float) _ievent,m_fSeed,
                               vx,
                               vy,
                               vz,
                               ntracks,
                               gvx,
                               gvy,
                               gvz,
                               gvt,
                               gembed,
                               gntracks,
                               gntracksmaps,
                               gnembed,
                               nfromtruth,
                               nhit_tpc_all,
                               nhit_tpc_in,
                               nhit_tpc_mid,
                               nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps};

        _ntp_vertex->Fill(vertex_data);
      }

      if (!_scan_for_embedded)
      {
        for (std::map<int, bool>::iterator iter = embedvtxid_found.begin();
             iter != embedvtxid_found.end();
             ++iter)
        {
          if (embedvtxid_found[iter->first]) continue;

          float vx = NAN;
          float vy = NAN;
          float vz = NAN;
          float ntracks = NAN;

          float gvx = NAN;
          float gvy = NAN;
          float gvz = NAN;
          float gvt = NAN;
          float gembed = iter->first;
          float gntracks = NAN;
          float gntracksmaps = NAN;
          float gnembed = NAN;
          float nfromtruth = NAN;

          PHG4VtxPoint* point = embedvtxid_vertex[gembed];

          if (point)
          {
            const int point_id = point->get_id();
            gvx = point->get_x();
            gvy = point->get_y();
            gvz = point->get_z();
            gvt = point->get_t();
            gembed = truthinfo->isEmbededVtx(point_id);
            gntracks = embedvtxid_particle_count[(int) gembed];
            if (embedvtxid_maps_particle_count[(int) gembed] > 0 && fabs(gvt) < 2000 && fabs(gvz) < 13.0)
              gntracksmaps = embedvtxid_maps_particle_count[(int) gembed];
            gnembed = (float) ngembed;
            //        nfromtruth = vertexeval->get_ntracks_contribution(vertex,point);
          }

          float vertex_data[] = {(float) _ievent,m_fSeed,
                                 vx,
                                 vy,
                                 vz,
                                 ntracks,
                                 gvx,
                                 gvy,
                                 gvz,
                                 gvt,
                                 gembed,
                                 gntracks,
                                 gntracksmaps,
                                 gnembed,
                                 nfromtruth,
                                 nhit_tpc_all,
                                 nhit_tpc_in,
                                 nhit_tpc_mid,
                                 nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps};

          _ntp_vertex->Fill(vertex_data);
        }
      }

    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      cout << "vertex time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
  }

  //-----------------------
  // fill the gpoint NTuple
  //-----------------------

  if (_ntp_gpoint)
  {
    if (Verbosity() > 1)
    {
      cout << "Filling ntp_gpoint " << endl;
      _timer->restart();
    }
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

    if (truthinfo)
    {
      auto vrange = truthinfo->GetPrimaryVtxRange();
      const auto prange = truthinfo->GetPrimaryParticleRange();

      map<int, unsigned int> vertex_particle_count;
      for (auto iter = prange.first; iter != prange.second; ++iter)  // process all primary paricle
      {
        ++vertex_particle_count[iter->second->get_vtx_id()];
      }

      for (auto iter = vrange.first; iter != vrange.second; ++iter)  // process all primary vertexes
      {
        const int point_id = iter->first;
        PHG4VtxPoint* point = iter->second;

        //      PHG4VtxPoint* point =  truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());

        if (point)
        {
          SvtxVertex* vertex = vertexeval->best_vertex_from(point);

          float gvx = point->get_x();
          float gvy = point->get_y();
          float gvz = point->get_z();
          float gvt = point->get_t();
          float gntracks = vertex_particle_count[point_id];

          float gembed = truthinfo->isEmbededVtx(point_id);
          float vx = NAN;
          float vy = NAN;
          float vz = NAN;
          float ntracks = NAN;
          float nfromtruth = NAN;

          if (vertex)
          {
            vx = vertex->get_x();
            vy = vertex->get_y();
            vz = vertex->get_z();
            ntracks = vertex->size_tracks();
            nfromtruth = vertexeval->get_ntracks_contribution(vertex, point);
          }

          float gpoint_data[] = {(float) _ievent,m_fSeed,
                                 gvx,
                                 gvy,
                                 gvz,
                                 gvt,
                                 gntracks,
                                 gembed,
                                 vx,
                                 vy,
                                 vz,
                                 ntracks,
                                 nfromtruth,
                                 nhit_tpc_all,
                                 nhit_tpc_in,
                                 nhit_tpc_mid,
                                 nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps};

          _ntp_gpoint->Fill(gpoint_data);
        }
      }
    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      cout << "gpoint time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
  }

  //---------------------
  // fill the G4hit NTuple
  //---------------------

  if (_ntp_g4hit)
  {
    if (Verbosity() > 1)
    {
      cout << "Filling ntp_g4hit " << endl;
      _timer->restart();
    }
    std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits();
    for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
         iter != g4hits.end();
         ++iter)
    {
      PHG4Hit* g4hit = *iter;
      PHG4Particle* g4particle = trutheval->get_particle(g4hit);

      float g4hitID = g4hit->get_hit_id();
      float gx = g4hit->get_avg_x();
      float gy = g4hit->get_avg_y();
      float gz = g4hit->get_avg_z();
      TVector3 vg4(gx, gy, gz);
      float gt = g4hit->get_avg_t();
      TVector3 vin(g4hit->get_x(0), g4hit->get_y(0), g4hit->get_z(0));
      TVector3 vout(g4hit->get_x(1), g4hit->get_y(1), g4hit->get_z(1));
      float gdphi = vin.DeltaPhi(vout);
      float gdz = fabs(g4hit->get_z(1) - g4hit->get_z(0));
      float gedep = g4hit->get_edep();
      float glayer = g4hit->get_layer();

      float gtrackID = g4hit->get_trkid();

      float gflavor = NAN;
      float gpx = NAN;
      float gpy = NAN;
      float gpz = NAN;
      TVector3 vec(g4hit->get_avg_x(), g4hit->get_avg_y(), g4hit->get_avg_z());
      float geta = vec.Eta();
      float gphi = vec.Phi();
      float gvx = NAN;
      float gvy = NAN;
      float gvz = NAN;

      float gembed = NAN;
      float gprimary = NAN;

      float gfpx = 0.;
      float gfpy = 0.;
      float gfpz = 0.;
      float gfx = 0.;
      float gfy = 0.;
      float gfz = 0.;

      if (g4particle)
      {
        if (_scan_for_embedded)
        {
          if (trutheval->get_embed(g4particle) <= 0) continue;
        }

        gflavor = g4particle->get_pid();
        gpx = g4particle->get_px();
        gpy = g4particle->get_py();
        gpz = g4particle->get_pz();

        PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);

        if (vtx)
        {
          gvx = vtx->get_x();
          gvy = vtx->get_y();
          gvz = vtx->get_z();
        }
        PHG4Hit* outerhit = nullptr;
        if (_do_eval_light == false)
          outerhit = trutheval->get_outermost_truth_hit(g4particle);

        if (outerhit)
        {
          gfpx = outerhit->get_px(1);
          gfpy = outerhit->get_py(1);
          gfpz = outerhit->get_pz(1);
          gfx = outerhit->get_x(1);
          gfy = outerhit->get_y(1);
          gfz = outerhit->get_z(1);
        }

        gembed = trutheval->get_embed(g4particle);
        gprimary = trutheval->is_primary(g4particle);
      }  //       if (g4particle)

      std::set<TrkrDefs::cluskey> clusters = clustereval->all_clusters_from(g4hit);
      float nclusters = clusters.size();

      // best cluster reco'd
      TrkrDefs::cluskey cluster_key = clustereval->best_cluster_from(g4hit);
      
      float clusID = NAN;
      float x = NAN;
      float y = NAN;
      float z = NAN;
      float eta = NAN;
      float phi = NAN;
      float e = NAN;
      float adc = NAN;
      float layer = NAN;
      float size = NAN;
      float phisize = NAN;
      float zsize = NAN;
      float efromtruth = NAN;
      float dphitru = NAN;
      float detatru = NAN;
      float dztru = NAN;
      float drtru = NAN;

      TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
      TrkrCluster *cluster = clustermap->findCluster(cluster_key);

      if (cluster)
      {
        clusID = cluster_key;
        x = cluster->getX();
        y = cluster->getY();
        z = cluster->getZ();
        TVector3 vec2(x, y, z);
        eta = vec2.Eta();
        phi = vec2.Phi();
        e = cluster->getAdc();
        adc = cluster->getAdc();
        layer = (float) TrkrDefs::getLayer(cluster_key);
        size = 0.0;
	// count all hits for this cluster

	TrkrClusterHitAssoc *cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
	TrkrClusterHitAssoc::ConstRange hitrange = cluster_hit_map->getHits(cluster_key);  
	for(TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
	  {
	    ++size; 
	  }
        phisize = cluster->getPhiSize();
        zsize = cluster->getZSize();
        dphitru = vec2.DeltaPhi(vg4);
        detatru = eta - geta;
        dztru = z - gz;
        drtru = vec2.DeltaR(vg4);
        if (g4particle)
        {
          efromtruth = clustereval->get_energy_contribution(cluster_key, g4particle);
        }
      }

      float g4hit_data[] = {(float) _ievent,m_fSeed,
                            g4hitID,
                            gx,
                            gy,
                            gz,
                            gt,
                            gedep,
                            geta,
                            gphi,
                            gdphi,
                            gdz,
                            glayer,
                            gtrackID,
                            gflavor,
                            gpx,
                            gpy,
                            gpz,
                            gvx,
                            gvy,
                            gvz,
                            gfpx,
                            gfpy,
                            gfpz,
                            gfx,
                            gfy,
                            gfz,
                            gembed,
                            gprimary,
                            nclusters,
                            clusID,
                            x,
                            y,
                            z,
                            eta,
                            phi,
                            e,
                            adc,
                            layer,
                            size,
                            phisize,
                            zsize,
                            efromtruth,
                            dphitru,
                            detatru,
                            dztru,
                            drtru,
                            nhit_tpc_all,
                            nhit_tpc_in,
                            nhit_tpc_mid,
                            nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps};

      _ntp_g4hit->Fill(g4hit_data);
    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      cout << "g4hit time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
  }

  //--------------------
  // fill the Hit NTuple
  //--------------------

  if (_ntp_hit)
  {
    if (Verbosity() > 1)
    {
      cout << "Filling ntp_hit " << endl;
      _timer->restart();
    }
    // need things off of the DST...
    TrkrHitSetContainer *hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    PHG4CylinderCellGeomContainer* geom_container =
        findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    if (!geom_container)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return;
    }

    if (hitmap)
    {
      TrkrHitSetContainer::ConstRange all_hitsets = hitmap->getHitSets();
      for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;
           iter != all_hitsets.second;
           ++iter)
      {
	TrkrDefs::hitsetkey hitset_key = iter->first;
	TrkrHitSet *hitset = iter->second;

	// get all hits for this hitset
	TrkrHitSet::ConstRange hitrangei = hitset->getHits();
	for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	     hitr != hitrangei.second;
	     ++hitr)
	  {
	    TrkrDefs::hitkey hit_key = hitr->first;	    
	    TrkrHit *hit = hitr->second;
	    PHG4Hit* g4hit = hiteval->max_truth_hit_by_energy(hit_key);
	    PHG4Particle* g4particle = trutheval->get_particle(g4hit);
	    
	    float event = _ievent;
	    float hitID = hit_key;
	    float e = hit->getEnergy();
	    float adc = hit->getAdc();
	    float layer = TrkrDefs::getLayer(hitset_key);
	    float cellID = 0;
	    float ecell = hit->getAdc();
	    
	    float phibin = NAN;
	    float zbin = NAN;
	    float phi = NAN;
	    float z = NAN;
	    
	    if (layer >= _nlayers_maps + _nlayers_intt && layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc )
	      {
		PHG4CylinderCellGeom* GeoLayer = geom_container->GetLayerCellGeom(layer);
		phibin = (float) TpcDefs::getPad(hit_key);
		zbin = (float) TpcDefs::getTBin(hit_key);
		phi = GeoLayer->get_phicenter(phibin);
		z = GeoLayer->get_zcenter(zbin);
	      }

	    float g4hitID = NAN;
	    float gedep = NAN;
	    float gx = NAN;
	    float gy = NAN;
	    float gz = NAN;
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
	    
	    if (g4hit)
	      {
		g4hitID = g4hit->get_hit_id();
		gedep = g4hit->get_edep();
		gx = g4hit->get_avg_x();
		gy = g4hit->get_avg_y();
		gz = g4hit->get_avg_z();
		gt = g4hit->get_avg_t();
		
		if (g4particle)
		  {
		    if (_scan_for_embedded)
		      {
			if (trutheval->get_embed(g4particle) <= 0) continue;
		      }
		    
		    gtrackID = g4particle->get_track_id();
		    gflavor = g4particle->get_pid();
		    gpx = g4particle->get_px();
		    gpy = g4particle->get_py();
		    gpz = g4particle->get_pz();
		    
		    PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
		    
		    if (vtx)
		      {
			gvx = vtx->get_x();
			gvy = vtx->get_y();
			gvz = vtx->get_z();
			gvt = vtx->get_t();
		      }

		    PHG4Hit* outerhit = nullptr;
		    if (_do_eval_light == false)
		      outerhit = trutheval->get_outermost_truth_hit(g4particle);
		    if (outerhit)
		      {
			gfpx = outerhit->get_px(1);
			gfpy = outerhit->get_py(1);
			gfpz = outerhit->get_pz(1);
			gfx = outerhit->get_x(1);
			gfy = outerhit->get_y(1);
			gfz = outerhit->get_z(1);
		      }
		    gembed = trutheval->get_embed(g4particle);
		    gprimary = trutheval->is_primary(g4particle);
		  }  //   if (g4particle){
	      }
	    
	    if (g4particle)
	      {
		efromtruth = hiteval->get_energy_contribution(hit_key, g4particle);
	      }
	    
	    float hit_data[] = {
	      event,
              (float) _iseed,
	      hitID,
	      e,
	      adc,
	      layer,
	      cellID,
	      ecell,
	      (float) phibin,
	      (float) zbin,
	      phi,
	      z,
	      g4hitID,
	      gedep,
	      gx,
	      gy,
	      gz,
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
	      nhit_tpc_all,
	      nhit_tpc_in,
	      nhit_tpc_mid,
	      nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps};
	    
	    _ntp_hit->Fill(hit_data);
	  }
      }
    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      cout << "hit time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
  }

  //------------------------
  // fill the Cluster NTuple
  //------------------------

  if (Verbosity() > 1)
  {
    cout << "check for ntp_cluster" << endl;
    _timer->restart();
  }

  if (_ntp_cluster && !_scan_for_embedded)
  {
    if (Verbosity() > 1) cout << "Filling ntp_cluster (all of them) " << endl;
    // need things off of the DST...
    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    TrkrClusterHitAssoc* clusterhitmap = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
    if (clustermap && clusterhitmap)
    {
      TrkrClusterContainer::ConstRange all_clusters = clustermap->getClusters();
      for (TrkrClusterContainer::ConstIterator iter = all_clusters.first;
           iter != all_clusters.second;
           ++iter)
      {
	TrkrDefs::cluskey cluster_key = iter->first;
	TrkrCluster *cluster = clustermap->findCluster(cluster_key);
        SvtxTrack* track = trackeval->best_track_from(cluster_key);
	PHG4Particle* g4particle = clustereval->max_truth_particle_by_cluster_energy(cluster_key);

        float hitID = (float) cluster_key;
        float x = cluster->getX();
        float y = cluster->getY();
	float z = cluster->getZ();
        TVector3 pos(x, y, z);
        float r = pos.Perp();
        float phi = pos.Phi();
        float eta = pos.Eta();
        float theta = pos.Theta();
        float ex = sqrt(cluster->getError(0, 0));
        float ey = sqrt(cluster->getError(1, 1));
        float ez = cluster->getZError();
        float ephi = cluster->getRPhiError();

        float e = cluster->getAdc();
        float adc = cluster->getAdc();
        float layer = (float) TrkrDefs::getLayer(cluster_key);
        float size = 0;
	// count all hits for this cluster
	TrkrClusterHitAssoc::ConstRange hitrange = clusterhitmap->getHits(cluster_key);  
	for(TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
	  {
	    ++size; 
	  }
        float phisize = cluster->getPhiSize();
        float zsize = cluster->getZSize();

        float trackID = NAN;
        if (track) trackID = track->get_id();

        float g4hitID = NAN;
        float gx = NAN;
        float gy = NAN;
        float gz = NAN;
        float gr = NAN;
        float gphi = NAN;
        //float gedep = NAN;
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

	if(Verbosity() > 0)
	  {
	    TrkrDefs::cluskey reco_cluskey = cluster->getClusKey();		  
	    std::cout << PHWHERE << "  ****   reco: layer " << layer << std::endl;
	    cout << "              reco cluster key " << reco_cluskey << "  r " << r << "  x " << x << "  y " << y << "  z " << z << "  phi " << phi  << " adc " << adc << endl;
	  }

	// get best matching truth cluster from clustereval
	std::shared_ptr<TrkrCluster> truth_cluster = clustereval->max_truth_cluster_by_energy(cluster_key);
	if(truth_cluster)
	  {
	    if(Verbosity() > 0)
	      {
		TrkrDefs::cluskey truth_cluskey = truth_cluster->getClusKey();
		cout << "Found matching truth cluster with key " << truth_cluskey << " for reco cluster key " << cluster_key << " in layer " << layer << endl;
	      }

	    g4hitID = 0;
	    gx=truth_cluster->getX();
	    gy=truth_cluster->getY();
	    gz=truth_cluster->getZ();
	    efromtruth = truth_cluster->getError(0,0);

	    TVector3 gpos(gx, gy, gz);
	    gr = gpos.Perp();  // could also be just the center of the layer
	    gphi = gpos.Phi();
	    geta = gpos.Eta();

          if (g4particle)
          {
            gtrackID = g4particle->get_track_id();
            gflavor = g4particle->get_pid();
            gpx = g4particle->get_px();
            gpy = g4particle->get_py();
            gpz = g4particle->get_pz();

            PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
            if (vtx)
            {
              gvx = vtx->get_x();
              gvy = vtx->get_y();
              gvz = vtx->get_z();
              gvt = vtx->get_t();
            }

            PHG4Hit* outerhit = nullptr;
            if (_do_eval_light == false)
              outerhit = trutheval->get_outermost_truth_hit(g4particle);
            if (outerhit)
            {
              gfpx = outerhit->get_px(1);
              gfpy = outerhit->get_py(1);
              gfpz = outerhit->get_pz(1);
              gfx = outerhit->get_x(1);
              gfy = outerhit->get_y(1);
              gfz = outerhit->get_z(1);
            }

            gembed = trutheval->get_embed(g4particle);
            gprimary = trutheval->is_primary(g4particle);
          }  //   if (g4particle){
	  
	  if(Verbosity() > 0)
	    {
	      TrkrDefs::cluskey ckey = truth_cluster->getClusKey();		  
	      cout << "             truth cluster key " << ckey << " gr " << gr << " gx " << gx << " gy " << gy << " gz " << gz << " gphi " << gphi << " efromtruth " << efromtruth << endl;
	    }
	  }    //  if (truth_cluster) {

        if (g4particle)
        {

        }

        float nparticles = clustereval->all_truth_particles(cluster_key).size();
        float cluster_data[] = {(float) _ievent,
				(float) _iseed,
                                hitID,
                                x,
                                y,
                                z,
                                r,
                                phi,
                                eta,
				theta,
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
    }
  }
  else if (_ntp_cluster && _scan_for_embedded)
  {
    if (Verbosity() > 1) 
      cout << "Filling ntp_cluster (embedded only) " << endl;

    // if only scanning embedded signals, loop over all the tracks from
    // embedded particles and report all of their clusters, including those
    // from other sources (noise hits on the embedded track)

    // need things off of the DST...
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());
    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    TrkrClusterHitAssoc* clusterhitmap = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
    if (trackmap)
    {
      for (SvtxTrackMap::Iter iter = trackmap->begin();
           iter != trackmap->end();
           ++iter)
      {
        SvtxTrack* track = iter->second;

        PHG4Particle* truth = trackeval->max_truth_particle_by_nclusters(track);
        if (truth)
        {
          if (trutheval->get_embed(truth) <= 0) continue;
        }

        for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
             iter != track->end_cluster_keys();
             ++iter)
	  {
	  TrkrDefs::cluskey cluster_key = *iter;
          TrkrCluster* cluster = clustermap->findCluster(cluster_key);

          PHG4Hit* g4hit = clustereval->max_truth_hit_by_energy(cluster_key);
          PHG4Particle* g4particle = trutheval->get_particle(g4hit);

          float hitID = cluster_key;
          float x = cluster->getX();
          float y = cluster->getY();
          float z = cluster->getZ();
          TVector3 pos(x, y, z);
          float r = pos.Perp();
          float phi = pos.Phi();
          float eta = pos.Eta();
	  float theta = pos.Theta();
          float ex = sqrt(cluster->getError(0, 0));
          float ey = sqrt(cluster->getError(1, 1));
          float ez = cluster->getZError();

          float ephi = cluster->getRPhiError();

          float e = cluster->getAdc();
          float adc = cluster->getAdc();
          float layer = (float) TrkrDefs::getLayer(cluster_key);

	  // count all hits for this cluster
	  float size = 0;
	  TrkrClusterHitAssoc::ConstRange hitrange = clusterhitmap->getHits(cluster_key);  
	  for(TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
	    {
	      ++size; 
	    }	  
          float phisize = cluster->getPhiSize();
          float zsize = cluster->getZSize();

          float trackID = NAN;
          trackID = track->get_id();

	  float g4hitID = NAN;
	  float gx = NAN;
	  float gy = NAN;
	  float gz = NAN;
	  float gr = NAN;
	  float gphi = NAN;
	  //float gedep = NAN;
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

	  //cout << "Look for truth cluster to match reco cluster " << cluster_key << endl;

	// get best matching truth cluster from clustereval
	  std::shared_ptr<TrkrCluster> truth_cluster = clustereval->max_truth_cluster_by_energy(cluster_key);
	if(truth_cluster)
	  {
	    if(Verbosity() > 0)
	      {
		TrkrDefs::cluskey truth_cluskey = truth_cluster->getClusKey();
		cout << "         Found matching truth cluster with key " << truth_cluskey << " for reco cluster key " << cluster_key << " in layer " << layer << endl;
	      }

	    g4hitID = 0;
	    gx=truth_cluster->getX();
	    gy=truth_cluster->getY();
	    gz=truth_cluster->getZ();
	    efromtruth = truth_cluster->getError(0,0);

	    TVector3 gpos(gx, gy, gz);
	    gr = gpos.Perp();
	    gphi = gpos.Phi();
	    geta = gpos.Eta();

	    if (g4particle)
	      {
		gtrackID = g4particle->get_track_id();
		gflavor = g4particle->get_pid();
		gpx = g4particle->get_px();
		gpy = g4particle->get_py();
		gpz = g4particle->get_pz();
		
		PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
		if (vtx)
		  {
		    gvx = vtx->get_x();
		    gvy = vtx->get_y();
		    gvz = vtx->get_z();
		    gvt = vtx->get_t();
		  }
		PHG4Hit* outerhit = nullptr;
		if (_do_eval_light == false)
		  outerhit = trutheval->get_outermost_truth_hit(g4particle);
		if (outerhit)
		  {
		    gfpx = outerhit->get_px(1);
		    gfpy = outerhit->get_py(1);
		    gfpz = outerhit->get_pz(1);
		    gfx = outerhit->get_x(1);
		    gfy = outerhit->get_y(1);
		    gfz = outerhit->get_z(1);
		  }
		
		gembed = trutheval->get_embed(g4particle);
		gprimary = trutheval->is_primary(g4particle);
	      }  //   if (g4particle){
          }    //  if (g4hit) {
	
	float nparticles = clustereval->all_truth_particles(cluster_key).size();
	
	float cluster_data[] = {(float) _ievent,
				(float) _iseed,
                                hitID,
                                x,
                                y,
                                z,
                                r,
                                phi,
				eta,
				theta,
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
      }
    }
  }
  if (Verbosity() >= 1)
  {
    _timer->stop();
    cout << "cluster time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
  }

  // fill the truth cluster NTuple
  //-----------------------------------

  if (Verbosity() > 1)
  {
    cout << "check for ntp_g4cluster" << endl;
    _timer->restart();
  }

  if (_ntp_g4cluster)
    {
      if (Verbosity() > 0) 
	cout << "Filling ntp_g4cluster " << endl;

       PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");      
      PHG4TruthInfoContainer::ConstRange range = truthinfo->GetParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
           iter != range.second;
           ++iter)
	{
	  
	  PHG4Particle* g4particle = iter->second;
	  
	  if (_scan_for_embedded)
	    {
	      if (trutheval->get_embed(g4particle) <= 0) continue;
	    }

	  float gtrackID = g4particle->get_track_id();
	  float gflavor = g4particle->get_pid();	  
	  float gembed = trutheval->get_embed(g4particle);
	  float gprimary = trutheval->is_primary(g4particle);

	  if(Verbosity() > 0)
	    cout << PHWHERE << " PHG4Particle ID " << gtrackID << " gflavor " << gflavor << " gprimary " << gprimary << endl;

	  // Get the truth clusters from this particle
	  std::map<unsigned int, std::shared_ptr<TrkrCluster> > truth_clusters =   trutheval->all_truth_clusters(g4particle);

	  // loop over layers and add to ntuple
	  for ( auto it = truth_clusters.begin(); it != truth_clusters.end(); ++it  )
	    {
	      unsigned int layer = it->first;
	      std::shared_ptr<TrkrCluster> gclus = it->second;

	      float gx = gclus->getX();
	      float gy = gclus->getY();
	      float gz = gclus->getZ();
	      float gt = NAN;
	      float gedep = gclus->getError(0,0);
	      float ng4hits = gclus->getAdc();

	      TVector3 gpos(gx, gy, gz);
	      float gr = sqrt(gx*gx+gy*gy);
	      float gphi = gpos.Phi();
	      float geta = gpos.Eta();

	      if(Verbosity() > 0)
		{
		  TrkrDefs::cluskey ckey = gclus->getClusKey();		  
		  std::cout << PHWHERE << "  ****   truth: layer " << layer << std::endl;
		  cout << "             truth cluster key " << ckey << " hg4hits " << ng4hits << " gr " << gr << " gx " << gx << " gy " << gy << " gz " << gz 
		       << " gphi " << gphi << " gedep " << gedep << endl;
		}
	      
	      float gphisize = gclus->getSize(1,1);
	      float gzsize = gclus->getSize(2,2);

	      // Find the matching TrkrCluster, if it exists

	      float x = NAN;
	      float y = NAN;
	      float z = NAN;
	      float r = NAN;
	      float phi = NAN;
	      float eta = NAN;
	      float ex = NAN;
	      float ey = NAN;
	      float ez = NAN;
	      float ephi = NAN;
	      float phisize = NAN;
	      float zsize = NAN;
	      float adc = NAN;

	      float nreco = 0;

	      TrkrCluster *reco_cluster = clustereval->reco_cluster_from_truth_cluster(gclus);
	      if(reco_cluster)
		{
		  nreco = 1;
		  
		  x = reco_cluster->getX();
		  y = reco_cluster->getY();
		  z = reco_cluster->getZ();

		  TVector3 pos(x, y, z);
		  r = sqrt(x*x+y*y);
		  phi = pos.Phi();
		  eta = pos.Eta();
		  ex = sqrt(reco_cluster->getError(0, 0));
		  ey = sqrt(reco_cluster->getError(1, 1));
		  ez = reco_cluster->getZError();		  
		  ephi = reco_cluster->getRPhiError();

		  phisize = reco_cluster->getPhiSize();  
		  zsize = reco_cluster->getZSize();   
		  
		  adc = reco_cluster->getAdc();

		  if(Verbosity() > 0)
		    {
		      TrkrDefs::cluskey reco_cluskey = reco_cluster->getClusKey();		  
		      cout << "              reco cluster key " << reco_cluskey << "  r " << r << "  x " << x << "  y " << y << "  z " << z << "  phi " << phi  << " adc " << adc << endl;
		    }
		}
	      if(nreco == 0 && Verbosity() > 0)
		{
		  if(Verbosity() > 0)
		    cout << "   ----------- Failed to find matching reco cluster " << endl;
		}


	      // add this cluster to the ntuple

	      float g4cluster_data[] = {(float) _ievent,
					(float) layer,
					gx,
					gy,
					gz,
					gt,
					gedep,
					gr,
					gphi,
					geta,
					ng4hits,
					gtrackID,
					gflavor,
					gembed,
					gprimary,
					gphisize,
					gzsize,
					nreco,
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
					phisize,
					zsize,
					adc };
	      _ntp_g4cluster->Fill(g4cluster_data);
	    }
	}
    }
  
  //------------------------
  // fill the Gtrack NTuple
  //------------------------
  
  // need things off of the DST...
  
  //cout << "check for ntp_gtrack" << endl;
  
  //#ifdef FUCKER
  if (_ntp_gtrack)
  {
    if (Verbosity() > 1)
    {
      cout << "Filling ntp_gtrack " << endl;
      _timer->restart();
    }

    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (truthinfo)
    {
      PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
      Float_t gntracks = (Float_t) truthinfo->GetNumPrimaryVertexParticles();
      for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
           iter != range.second;
           ++iter)
      {

        PHG4Particle* g4particle = iter->second;

        if (_scan_for_embedded)
        {
          if (trutheval->get_embed(g4particle) <= 0) continue;
        }

        float gtrackID = g4particle->get_track_id();
        float gflavor = g4particle->get_pid();

        std::set<TrkrDefs::cluskey> g4clusters = clustereval->all_clusters_from(g4particle);

        float ng4hits = g4clusters.size();
        unsigned int ngmaps = 0;
        unsigned int ngintt = 0;
        unsigned int ngintt1 = 0;
        unsigned int ngintt2 = 0;
        unsigned int ngintt3 = 0;
        unsigned int ngintt4 = 0;
        unsigned int ngintt5 = 0;
        unsigned int ngintt6 = 0;
        unsigned int ngintt7 = 0;
        unsigned int ngintt8 = 0;
        unsigned int ngtpc = 0;
        unsigned int nglmaps = 0;
        unsigned int nglintt = 0;
        unsigned int ngltpc = 0;

        int lmaps[_nlayers_maps + 1];
        if (_nlayers_maps > 0)
          for (unsigned int i = 0; i < _nlayers_maps; i++) lmaps[i] = 0;

        int lintt[_nlayers_intt + 1];
        if (_nlayers_intt > 0)
          for (unsigned int i = 0; i < _nlayers_intt; i++) lintt[i] = 0;

        int ltpc[_nlayers_tpc + 1];
        if (_nlayers_tpc > 0)
          for (unsigned int i = 0; i < _nlayers_tpc; i++) ltpc[i] = 0;

        for (const TrkrDefs::cluskey g4cluster : g4clusters)
        {
          unsigned int layer = TrkrDefs::getLayer(g4cluster);
          //cout<<__LINE__<<": " << _ievent <<": " <<gtrackID << ": " << layer <<": " <<g4cluster->get_id() <<endl;
          if (_nlayers_maps > 0 && layer < _nlayers_maps)
          {
            lmaps[layer] = 1;
            ngmaps++;
          }

          if (_nlayers_intt > 0 && layer >= _nlayers_maps && layer < _nlayers_maps + _nlayers_intt)
          {
            lintt[layer - _nlayers_maps] = 1;
            ngintt++;
          }

          if (_nlayers_intt > 0 && layer == _nlayers_maps && layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt1++;
          }

          if (_nlayers_intt > 1 && layer == _nlayers_maps + 1 && layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt2++;
          }

          if (_nlayers_intt > 2 && layer == _nlayers_maps + 2 && layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt3++;
          }

          if (_nlayers_intt > 3 && layer == _nlayers_maps + 3 && layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt4++;
          }

          if (_nlayers_intt > 4 && layer == _nlayers_maps + 4 && layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt5++;
          }

          if (_nlayers_intt > 5 && layer == _nlayers_maps + 5 && layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt6++;
          }

          if (_nlayers_intt > 6 && layer == _nlayers_maps + 6 && layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt7++;
          }

          if (_nlayers_intt > 7 && layer == _nlayers_maps + 7 && layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt8++;
          }
          if (_nlayers_tpc > 0 && layer >= _nlayers_maps + _nlayers_intt && layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
          {
            ltpc[layer - (_nlayers_maps + _nlayers_intt)] = 1;
            ngtpc++;
          }
        }
        if (_nlayers_maps > 0)
          for (unsigned int i = 0; i < _nlayers_maps; i++) nglmaps += lmaps[i];
        if (_nlayers_intt > 0)
          for (unsigned int i = 0; i < _nlayers_intt; i++) nglintt += lintt[i];
        if (_nlayers_tpc > 0)
          for (unsigned int i = 0; i < _nlayers_tpc; i++) ngltpc += ltpc[i];

        float gpx = g4particle->get_px();
        float gpy = g4particle->get_py();
        float gpz = g4particle->get_pz();
        float gpt = NAN;
        float geta = NAN;
        float gphi = NAN;
        if (gpx != 0 && gpy != 0)
        {
          TVector3 gv(gpx, gpy, gpz);
          gpt = gv.Pt();
          geta = gv.Eta();
          gphi = gv.Phi();
        }
        PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
        float gvx = vtx->get_x();
        float gvy = vtx->get_y();
        float gvz = vtx->get_z();
        float gvt = vtx->get_t();

        float gfpx = 0.;
        float gfpy = 0.;
        float gfpz = 0.;
        float gfx = 0.;
        float gfy = 0.;
        float gfz = 0.;

        PHG4Hit* outerhit = nullptr;
        if (_do_eval_light == false)
          outerhit = trutheval->get_outermost_truth_hit(g4particle);

        if (outerhit)
        {
          gfpx = outerhit->get_px(1);
          gfpy = outerhit->get_py(1);
          gfpz = outerhit->get_pz(1);
          gfx = outerhit->get_x(1);
          gfy = outerhit->get_y(1);
          gfz = outerhit->get_z(1);
        }

        float gembed = trutheval->get_embed(g4particle);
        float gprimary = trutheval->is_primary(g4particle);

        float trackID = NAN;
        float charge = NAN;
        float quality = NAN;
        float chisq = NAN;
        float ndf = NAN;
        float nhits = NAN;
        float nmaps = 0;
        float nintt = 0;
        float ntpc = 0;
        float ntpc1 = 0;
        float ntpc11 = 0;
        float ntpc2 = 0;
        float ntpc3 = 0;
        float nlintt = 0;
        float nlmaps = 0;
        float nltpc = 0;
        unsigned int layers = 0x0;
        float dca2d = NAN;
        float dca2dsigma = NAN;
        float dca3dxy = NAN;
        float dca3dxysigma = NAN;
        float dca3dz = NAN;
        float dca3dzsigma = NAN;
        float px = NAN;
        float py = NAN;
        float pz = NAN;
        float pt = NAN;
        float eta = NAN;
        float phi = NAN;
	float deltapt = NAN;
	float deltaeta = NAN;
	float deltaphi = NAN;
        float pcax = NAN;
        float pcay = NAN;
        float pcaz = NAN;

        float nfromtruth = NAN;
        float nwrong = NAN;
        float ntrumaps = NAN;
        float ntruintt = NAN;
        float ntrutpc = NAN;
        float ntrutpc1 = NAN;
        float ntrutpc11 = NAN;
        float ntrutpc2 = NAN;
        float ntrutpc3 = NAN;
        float layersfromtruth = NAN;

        if (_do_track_match)
        {
          SvtxTrack* track = trackeval->best_track_from(g4particle);

          if (track)
          {
            trackID = track->get_id();
            charge =  track->get_charge();
            quality = track->get_quality();
            chisq =   track->get_chisq();
            ndf =     track->get_ndf();
            nhits =   track->size_cluster_keys();

            vector <int> maps(_nlayers_maps, 0);
            vector <int> intt(_nlayers_intt, 0);
            vector <int> tpc(_nlayers_tpc, 0);

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

            for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
                 iter != track->end_cluster_keys();
                 ++iter)
            {
	      TrkrDefs::cluskey cluster_key = *iter;
              //TrkrCluster* cluster = clustermap->findCluster(cluster_key);
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
              if (_nlayers_tpc > 0 &&
                  layer >= (_nlayers_maps + _nlayers_intt) &&
                  layer < (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
              {
                tpc[layer - (_nlayers_maps + _nlayers_intt)] = 1;
                ntpc++;
		if((layer - (_nlayers_maps + _nlayers_intt))<16){
		  ntpc1++;
		}
		if((layer - (_nlayers_maps + _nlayers_intt))<8){
		  ntpc11++;
		}
		else if((layer - (_nlayers_maps + _nlayers_intt))<32){
		  ntpc2++;
		}
		else if((layer - (_nlayers_maps + _nlayers_intt))<48){
		  ntpc3++;
		}
              }
            }
            if (_nlayers_maps > 0)
              for (unsigned int i = 0; i < _nlayers_maps; i++) nlmaps += maps[i];
            if (_nlayers_intt > 0)
              for (unsigned int i = 0; i < _nlayers_intt; i++) nlintt += intt[i];
            if (_nlayers_tpc > 0)
              for (unsigned int i = 0; i < _nlayers_tpc; i++) nltpc += tpc[i];

            layers = nlmaps + nlintt + nltpc;
	    /* cout << " layers " << layers 
		 << " nmaps " << nmaps 
		 << " nintt " << nintt 
		 << " ntpc  " << ntpc 
		 << " nlmaps "<< nlmaps  
		 << " nlintt " << nlintt 
		 << " nltpc  " << nltpc
		 << endl;
	    */
            dca2d = track->get_dca2d();
            dca2dsigma = track->get_dca2d_error();
            dca3dxy = track->get_dca3d_xy();
            dca3dxysigma = track->get_dca3d_xy_error();
            dca3dz = track->get_dca3d_z();
            dca3dzsigma = track->get_dca3d_z_error();
            px = track->get_px();
            py = track->get_py();
            pz = track->get_pz();
            TVector3 v(px, py, pz);
            pt = v.Pt();
            eta = v.Eta();
            phi = v.Phi();
	    float CVxx = track->get_error(3,3);
	    float CVxy = track->get_error(3,4);
	    float CVxz = track->get_error(3,5);
	    float CVyy = track->get_error(4,4);
	    float CVyz = track->get_error(4,5);
	    float CVzz = track->get_error(5,5);
	    deltapt = sqrt((CVxx*px*px+2*CVxy*px*py+CVyy*py*py)/(px*px+py*py));
	    deltaeta = sqrt((CVzz*(px*px+py*py)*(px*px+py*py)+pz*(-2*(CVxz*px+CVyz*py)*(px*px+py*py)+CVxx*px*px*pz+CVyy*py*py*pz+2*CVxy*px*py*pz))/((px*px+py*py)*(px*px+py*py)*(px*px+py*py+pz*pz)));
	    deltaphi = sqrt((CVyy*px*px-2*CVxy*px*py+CVxx*py*py)/((px*px+py*py)*(px*px+py*py)));
            pcax = track->get_x();
            pcay = track->get_y();
            pcaz = track->get_z();

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
            ntrutpc1 = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt, _nlayers_maps + _nlayers_intt + 16);
            ntrutpc11 = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt, _nlayers_maps + _nlayers_intt + 8);
            ntrutpc2 = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt+16, _nlayers_maps + _nlayers_intt + 32);
            ntrutpc3 = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt+32, _nlayers_maps + _nlayers_intt + _nlayers_tpc);

            layersfromtruth = trackeval->get_nclusters_contribution_by_layer(track, g4particle);
          }
        }
	float gtrack_data[] = {(float) _ievent,m_fSeed,
                               gntracks,
                               gtrackID,
                               gflavor,
                               ng4hits,
                               (float) ngmaps,
                               (float) ngintt,
                               (float) ngintt1,
                               (float) ngintt2,
                               (float) ngintt3,
                               (float) ngintt4,
                               (float) ngintt5,
                               (float) ngintt6,
                               (float) ngintt7,
                               (float) ngintt8,
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
                               trackID,
                               px,
                               py,
                               pz,
                               pt,
                               eta,
                               phi,
			       deltapt,
			       deltaeta,
			       deltaphi,
                               charge,
                               quality,
                               chisq,
                               ndf,
                               nhits,
                               (float) layers,
                               nmaps,
                               nintt,
			       ntpc,
                               ntpc1,
                               ntpc11,
                               ntpc2,
                               ntpc3,
                               nlmaps,
                               nlintt,
                               nltpc,
                               dca2d,
                               dca2dsigma,
                               dca3dxy,
                               dca3dxysigma,
                               dca3dz,
                               dca3dzsigma,
                               pcax,
                               pcay,
                               pcaz,
                               nfromtruth,
                               nwrong,
                               ntrumaps,
                               ntruintt,
                               ntrutpc,
                               ntrutpc1,
                               ntrutpc11,
                               ntrutpc2,
                               ntrutpc3,
                               layersfromtruth,
                               nhit_tpc_all,
                               nhit_tpc_in,
                               nhit_tpc_mid,
                               nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps};

        /*
	cout << " ievent " << _ievent
	     << " gtrackID " << gtrackID
	     << " gflavor " << gflavor
	     << " ng4hits " << ng4hits
	     << endl;
	*/

        _ntp_gtrack->Fill(gtrack_data);

      }

    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      cout << "gtrack time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
  }

  //------------------------
  // fill the Track NTuple
  //------------------------

  if (_ntp_track)
  {
    if (Verbosity() > 1)
    {
      cout << "Filling ntp_track " << endl;
      _timer->restart();
    }

    // need things off of the DST...
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());
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
        float nhits = track->size_cluster_keys();
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
        float ntpc1 = 0;
        float ntpc11 = 0;
        float ntpc2 = 0;
        float ntpc3 = 0;
        float nlmaps = 0;
        float nlintt = 0;
        float nltpc = 0;

        for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
             iter != track->end_cluster_keys();
             ++iter)
        {
	  TrkrDefs::cluskey cluster_key = *iter;
          //TrkrCluster* cluster = clustermap->findCluster(cluster_key);
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
	    if((layer - (_nlayers_maps + _nlayers_intt))<16){
	      ntpc1++;
	    }
	    if((layer - (_nlayers_maps + _nlayers_intt))<8){
	      ntpc11++;
	    }
	    else if((layer - (_nlayers_maps + _nlayers_intt))<32){
	      ntpc2++;
	    }
	    else if((layer - (_nlayers_maps + _nlayers_intt))<48){
	      ntpc3++;
	    }
          }
        }
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
	float CVxx = track->get_error(3,3);
	float CVxy = track->get_error(3,4);
	float CVxz = track->get_error(3,5);
	float CVyy = track->get_error(4,4);
	float CVyz = track->get_error(4,5);
	float CVzz = track->get_error(5,5);
	float deltapt = sqrt((CVxx*px*px+2*CVxy*px*py+CVyy*py*py)/(px*px+py*py));
        float deltaeta = sqrt((CVzz*(px*px+py*py)*(px*px+py*py)+pz*(-2*(CVxz*px+CVyz*py)*(px*px+py*py)+CVxx*px*px*pz+CVyy*py*py*pz+2*CVxy*px*py*pz))/((px*px+py*py)*(px*px+py*py)*(px*px+py*py+pz*pz)));
	float deltaphi = sqrt((CVyy*px*px-2*CVxy*px*py+CVxx*py*py)/((px*px+py*py)*(px*px+py*py)));
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
        float ntrutpc1 = NAN;
        float ntrutpc11 = NAN;
        float ntrutpc2 = NAN;
        float ntrutpc3 = NAN;
        float layersfromtruth = NAN;

        if (_do_track_match)
        {
          PHG4Particle* g4particle = trackeval->max_truth_particle_by_nclusters(track);
          if (g4particle)
          {
            if (_scan_for_embedded)
            {
              if (trutheval->get_embed(g4particle) <= 0) continue;
            }

            gtrackID = g4particle->get_track_id();
            gflavor = g4particle->get_pid();

            std::set<TrkrDefs::cluskey> g4clusters = clustereval->all_clusters_from(g4particle);
            ng4hits = g4clusters.size();
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

            for (const TrkrDefs::cluskey g4cluster : g4clusters)
            {
              unsigned int layer = TrkrDefs::getLayer(g4cluster);
              if (_nlayers_maps > 0 && layer < _nlayers_maps)
              {
                lmaps[layer] = 1;
                ngmaps++;
              }

              if (_nlayers_intt > 0 && layer >= _nlayers_maps && layer < _nlayers_maps + _nlayers_intt)
              {
                lintt[layer - _nlayers_maps] = 1;
                ngintt++;
              }

              if (_nlayers_tpc > 0 && layer >= _nlayers_maps + _nlayers_intt && layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
              {
                ltpc[layer - (_nlayers_maps + _nlayers_intt)] = 1;
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
            PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
            gvx = vtx->get_x();
            gvy = vtx->get_y();
            gvz = vtx->get_z();
            gvt = vtx->get_t();

            PHG4Hit* outerhit = nullptr;
            if (_do_eval_light == false)
              outerhit = trutheval->get_outermost_truth_hit(g4particle);
            if (outerhit)
            {
              gfpx = outerhit->get_px(1);
              gfpy = outerhit->get_py(1);
              gfpz = outerhit->get_pz(1);
              gfx = outerhit->get_x(1);
              gfy = outerhit->get_y(1);
              gfz = outerhit->get_z(1);
            }
            gembed = trutheval->get_embed(g4particle);
            gprimary = trutheval->is_primary(g4particle);

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
            ntrutpc1 = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt, _nlayers_maps + _nlayers_intt + 16);
            ntrutpc11 = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt, _nlayers_maps + _nlayers_intt + 8);
            ntrutpc2 = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt+16, _nlayers_maps + _nlayers_intt + 32);
            ntrutpc3 = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt+32, _nlayers_maps + _nlayers_intt + _nlayers_tpc);
            layersfromtruth = trackeval->get_nclusters_contribution_by_layer(track, g4particle);
          }
        }

        float track_data[] = {(float) _ievent,m_fSeed,
                              trackID,
                              px,
                              py,
                              pz,
                              pt,
                              eta,
                              phi,
			      deltapt,
			      deltaeta,
			      deltaphi,
                              charge,
                              quality,
                              chisq,
                              ndf,
                              nhits, nmaps, nintt, ntpc,
			      ntpc1,ntpc11,ntpc2,ntpc3,
			      nlmaps, nlintt, nltpc,
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
                              ntrutpc1,
                              ntrutpc11,
                              ntrutpc2,
                              ntrutpc3,
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
    if (Verbosity() > 1)
    {
      _timer->stop();
      cout << "track time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
  }

  //---------------------
  // fill the Gseed NTuple
  //---------------------
  
  if (_ntp_gseed)
    {
      if (Verbosity() > 1)
	{
	  cout << "Filling ntp_gseed " << endl;
	  _timer->restart();
	}
      
      PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
      
      float gx = NAN;
      float gy = NAN;
      float gz = NAN;
      float gr = NAN;
      float geta = NAN;
      float gphi = NAN;
      float glayer = NAN;
      float gpx = NAN;
      float gpy = NAN;
      float gpz = NAN;
      float gtpt = NAN;
      float gtphi = NAN;
      float gteta = NAN;
      float gvx = NAN;
      float gvy = NAN;
      float gvz = NAN;
      float gembed = NAN;
      float gprimary = NAN;
      float gflav = NAN;
      float dphiprev = NAN;
      float detaprev = NAN;
      
      float xval[_nlayers_maps + _nlayers_intt + _nlayers_tpc];
      float yval[_nlayers_maps + _nlayers_intt + _nlayers_tpc];
      float zval[_nlayers_maps + _nlayers_intt + _nlayers_tpc];
      if (truthinfo)
	{
	  int ntrk = 0;
	  PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
	  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
	       iter != range.second;
	       ++iter)
	    {
	      ntrk++;
	      PHG4Particle* g4particle = iter->second;
	      for (unsigned int i = 0; i < _nlayers_maps + _nlayers_intt + _nlayers_tpc; i++)
		{
		  xval[i] = 0;
		  yval[i] = 0;
		  zval[i] = 0;
		}
	      std::set<PHG4Hit*> truth_hits = trutheval->all_truth_hits(g4particle);
	      for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
		   iter != truth_hits.end();
		   ++iter)
		{
		  PHG4Hit* g4hit = *iter;
		  unsigned int layer = g4hit->get_layer();
		  //cout << "  g4hit " << g4hit->get_hit_id() << " layer = " << layer << endl;
		  if (layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc)
		    {
		      //cout << PHWHERE << " skipping out of bounds detector id " << layer << endl;
		      continue;
		    }
		  xval[layer] = g4hit->get_avg_x();
		  yval[layer] = g4hit->get_avg_y();
		  zval[layer] = g4hit->get_avg_z();
		}
	      
	      for (unsigned int i = 0; i < _nlayers_maps + _nlayers_intt + _nlayers_tpc; i++)
		{
		  gx = xval[i];
		  gy = yval[i];
		  gz = zval[i];
		  if (gx == 0 && gy == 0) continue;
		  
		  TVector3 vg4(gx, gy, gz);
		  glayer = i;
		  gr = vg4.Perp();
		  geta = vg4.Eta();
		  gphi = vg4.Phi();
		  gpx = g4particle->get_px();
		  gpy = g4particle->get_py();
		  gpz = g4particle->get_pz();
		  TVector3 vg4p(gpx, gpy, gpz);
		  
		  gtpt = vg4p.Perp();
		  gtphi = vg4p.Phi();
		  gteta = vg4p.Eta();
		  
		  PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
		  
		  if (vtx)
		    {
		      gvx = vtx->get_x();
		      gvy = vtx->get_y();
		      gvz = vtx->get_z();
		    }
		  
		  gembed = trutheval->get_embed(g4particle);
		  gprimary = trutheval->is_primary(g4particle);
		  gflav = g4particle->get_pid();
		  if (i >= 1)
		    {
		      if (xval[i - 1] != 0 && yval[i - 1] != 0)
			{
			  TVector3 vg4prev(xval[i - 1], yval[i - 1], zval[i - 1]);
			  dphiprev = vg4.DeltaPhi(vg4prev);
			  detaprev = geta - vg4prev.Eta();
			}
		    }
		  
		  float ntrk_f = ntrk;
		  float _ievent_f = _ievent;
		  float gseed_data[] = {_ievent_f,
					ntrk_f,
					gx,
					gy,
					gz,
					gr,
					geta,
					gphi,
					glayer,
					gpx,
					gpy,
					gpz,
					gtpt,
					gtphi,
					gteta,
					gvx,
					gvy,
					gvz,
					gembed,
					gprimary,
					gflav,
					dphiprev,
					detaprev,
					nhit_tpc_all,
					nhit_tpc_in,
					nhit_tpc_mid,
					nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps};
		  
		  _ntp_gseed->Fill(gseed_data);
		}
	    }
	}
      
      if (Verbosity() > 1)
	{
	  _timer->stop();
	  cout << "g4hit time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
	}
    }
  return;

}

