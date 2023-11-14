#include "TrkrNtuplizer.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrackFitUtils.h>

#include <trackbase_historic/TrackSeedContainer_v1.h>

#include <trackbase/TpcDefs.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrHitSet.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <trackermillepedealignment/HelicalFitter.h>

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
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>  // for shared_ptr
#include <set>     // for _Rb_tree_cons...
#include <utility>
#include <vector>

using namespace std;

TrkrNtuplizer::TrkrNtuplizer(const string& /*name*/, const string& filename, const string& trackmapname,
                             unsigned int nlayers_maps,
                             unsigned int nlayers_intt,
                             unsigned int nlayers_tpc,
                             unsigned int nlayers_mms)
  : SubsysReco("TrkrNtuplizer")
  , _ievent(0)
  , _iseed(0)
  , m_fSeed(NAN)
  , _do_info_eval(true)
  , _do_vertex_eval(true)
  , _do_hit_eval(true)
  , _do_cluster_eval(true)
  , _do_clus_trk_eval(true)
  , _do_track_eval(true)
  , _do_tpcseed_eval(false)
  , _do_siseed_eval(false)
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _nlayers_mms(nlayers_mms)
  , _ntp_info(nullptr)
  , _ntp_vertex(nullptr)
  , _ntp_hit(nullptr)
  , _ntp_cluster(nullptr)
  , _ntp_clus_trk(nullptr)
  , _ntp_track(nullptr)
  , _ntp_tpcseed(nullptr)
  , _ntp_siseed(nullptr)
  , _filename(filename)
  , _trackmapname(trackmapname)
  , _tfile(nullptr)
  , _timer(nullptr)
{
}

TrkrNtuplizer::~TrkrNtuplizer()
{
  delete _timer;
}

int TrkrNtuplizer::Init(PHCompositeNode* /*topNode*/)
{
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");
  _tfile->SetCompressionLevel(7);
  if (_do_info_eval)
  {
    _ntp_info = new TNtuple("ntp_info", "event info",
                            "event:seed:"
                            "occ11:occ116:occ21:occ216:occ31:occ316:"
                            "ntrk:"
			    "nhitmvtx:nhitintt:nhittpot:"
                            "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_vertex_eval)
  {
    _ntp_vertex = new TNtuple("ntp_vertex", "vertex => max truth",
                              "event:seed:vertexID:vx:vy:vz:ntracks:chi2:ndof:"
                              "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_hit_eval)
  {
    _ntp_hit = new TNtuple("ntp_hit", "svtxhit => max truth",
                           "event:seed:hitID:e:adc:layer:phielem:zelem:"
                           "cellID:ecell:phibin:tbin:phi:x:y:z:"
                           "nhitmvtx:nhitintt:nhittpot:"
			   "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_cluster_eval)
  {
    _ntp_cluster = new TNtuple("ntp_cluster", "svtxcluster => max truth",
                               "event:seed:hitID:locx:locy:x:y:z:r:phi:eta:theta:phibin:tbin:ex:ey:ez:ephi:pez:pephi:"
                               "e:adc:maxadc:layer:phielem:zelem:size:phisize:zsize:"
                               "trackID:niter:"
                               "nhitmvtx:nhitintt:nhittpot:"
			       "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }
  if (_do_clus_trk_eval)
  {
    _ntp_clus_trk = new TNtuple("ntp_clus_trk", "cluster on track",
				"event:seed:locx:locy:x:y:z:r:phi:eta:theta:phibin:tbin:ex:ey:ez:ephi:pez:pephi:"
				"e:adc:maxadc:layer:phielem:zelem:size:phisize:zsize:"
				"resphi:resz:"
				"trackID:niter:pt:eta:phi:X0:Y0:Z0:charge:nhits:"
				"vertexID:vx:vy:vz:dca2d:dca2dsigma:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:"
				"nhitmvtx:nhitintt:nhittpot:"
				"nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }
  
  if (_do_track_eval)
  {
    _ntp_track = new TNtuple("ntp_track", "svtxtrack => max truth",
                             "event:seed:trackID:crossing:px:py:pz:pt:eta:phi:deltapt:deltaeta:deltaphi:charge:"
                             "quality:chisq:ndf:nhits:nmaps:nintt:ntpc:nmms:ntpc1:ntpc11:ntpc2:ntpc3:nlmaps:nlintt:nltpc:nlmms:layers:"
                             "vertexID:vx:vy:vz:dca2d:dca2dsigma:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:"
                             "nhitmvtx:nhitintt:nhittpot:"
			     "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_tpcseed_eval)
  {
    _ntp_tpcseed = new TNtuple("ntp_tpcseed", "seeds from truth",
                               "event:seed:ntrk:gx:gy:gz:gr:geta:gphi:"
                               "glayer:"
                               "gpx:gpy:gpz:gtpt:gtphi:gteta:"
                               "gvx:gvy:gvz:"
                               "gembed:gprimary:gflav:"
                               "dphiprev:detaprev:"
                               "nhitmvtx:nhitintt:nhittpot:"
			       "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }
  if (_do_siseed_eval)
  {
    _ntp_tpcseed = new TNtuple("ntp_siseed", "seeds from truth",
                               "event:seed:ntrk:gx:gy:gz:gr:geta:gphi:"
                               "glayer:"
                               "gpx:gpy:gpz:gtpt:gtphi:gteta:"
                               "gvx:gvy:gvz:"
                               "gembed:gprimary:gflav:"
                               "dphiprev:detaprev:"
                               "nhitmvtx:nhitintt:nhittpot:"
			       "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }
  _timer = new PHTimer("_eval_timer");
  _timer->stop();

  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkrNtuplizer::InitRun(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkrNtuplizer::process_event(PHCompositeNode* topNode)
{
  if ((Verbosity() > 1) && (_ievent % 100 == 0))
  {
    cout << "TrkrNtuplizer::process_event - Event = " << _ievent << endl;
  }

  recoConsts* rc = recoConsts::instance();
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
  if(_trackmap == nullptr)
    _trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());

  if (Verbosity() > 1)
  {
    cout << "TrkrNtuplizer::process_event - Seed = " << _iseed << endl;
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

  // printOutputInfo(topNode);

  ++_ievent;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkrNtuplizer::End(PHCompositeNode* /*topNode*/)
{
  _tfile->cd();

  if (_ntp_info)
  {
    _ntp_info->Write();
  }
  if (_ntp_vertex)
  {
    _ntp_vertex->Write();
  }
  if (_ntp_hit)
  {
    _ntp_hit->Write();
  }
  if (_ntp_cluster)
  {
    _ntp_cluster->Write();
  }
  if (_ntp_clus_trk)
  {
    _ntp_clus_trk->Write();
  }
  if (_ntp_track)
  {
    _ntp_track->Write();
  }
  if (_ntp_tpcseed)
  {
    _ntp_tpcseed->Write();
  }
  if (_ntp_siseed)
  {
    _ntp_siseed->Write();
  }

  _tfile->Close();

  delete _tfile;

  if (Verbosity() > 1)
  {
    cout << "========================= TrkrNtuplizer::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TrkrNtuplizer::printInputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    cout << "TrkrNtuplizer::printInputInfo() entered" << endl;
  }

  if (Verbosity() > 3)
  {
    // event information
    cout << endl;
    cout << PHWHERE << "   INPUT FOR EVENT " << _ievent << endl;

    cout << endl;

    cout << "---SVTXCLUSTERS-------------" << endl;
    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
    if (!clustermap)
    {
      clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    }

    if (clustermap != nullptr)
    {
      unsigned int icluster = 0;
      for (const auto& hitsetkey : clustermap->getHitSetKeys())
      {
        auto range = clustermap->getClusters(hitsetkey);
        for (auto iter = range.first; iter != range.second; ++iter)
        {
          TrkrDefs::cluskey cluster_key = iter->first;
          cout << icluster << " with key " << cluster_key << " of " << clustermap->size();
          cout << ": SvtxCluster: " << endl;
          iter->second->identify();
          ++icluster;
        }
      }
    }

    cout << "---SVXTRACKS-------------" << endl;

    if (_trackmap)
    {
      unsigned int itrack = 0;
      for (SvtxTrackMap::Iter iter = _trackmap->begin();
           iter != _trackmap->end();
           ++iter)
      {
        cout << itrack << " of " << _trackmap->size();
        SvtxTrack* track = iter->second;
        cout << " : SvtxTrack:" << endl;
        track->identify();
        cout << endl;
        ++itrack;
      }
    }

    cout << "---SVXVERTEXES-------------" << endl;
    SvtxVertexMap* vertexmap = nullptr;
    vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");  // Acts vertices

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

void TrkrNtuplizer::printOutputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    cout << "TrkrNtuplizer::printOutputInfo() entered" << endl;
  }

  //==========================================
  // print out some useful stuff for debugging
  //==========================================

  if (Verbosity() > 100)
  {
    // event information
    cout << endl;
    cout << PHWHERE << "   NEW OUTPUT FOR EVENT " << _ievent << endl;
    cout << endl;

    float vx = NAN;
    float vy = NAN;
    float vz = NAN;

    SvtxVertexMap* vertexmap = nullptr;

    vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");  // Acts vertices

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
    cout << "vreco = (" << vx << "," << vy << "," << vz << ")" << endl;
    cout << endl;

    cout << "===Tracking Summary============================" << endl;

    TrkrHitSetContainer* hitsetmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
    if (!clustermap)
    {
      clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    }

    unsigned int nclusters[100] = {0};
    unsigned int nhits[100] = {0};

    ActsGeometry* tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

    if (!tgeometry)
    {
      std::cout << PHWHERE << "No Acts geometry on node tree. Can't  continue."
                << std::endl;
    }

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
        auto range = clustermap->getClusters(hitsetiter->first);
        for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
        {
          const auto cluskey = clusIter->first;
          //	    const auto cluster = clusIter->second;
          //	    unsigned int layer = TrkrDefs::getLayer(cluskey);
          nclusters[TrkrDefs::getLayer(cluskey)]++;
        }
      }
    }

    PHG4TpcCylinderGeomContainer* geom_container =
        findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    {
      if (!geom_container)
      {
        std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      }
      return;
    }

    for (unsigned int ilayer = 0; ilayer < _nlayers_maps + _nlayers_intt + _nlayers_tpc; ++ilayer)
    {
      PHG4TpcCylinderGeom* GeoLayer = geom_container->GetLayerCellGeom(ilayer);

      cout << "layer " << ilayer
           << " => nHits = " << nhits[ilayer]
           << " => nClusters = " << nclusters[ilayer] << endl;
      if (ilayer >= _nlayers_maps + _nlayers_intt && ilayer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
      {
        cout << "layer " << ilayer
             << " => nphi = " << GeoLayer->get_phibins()
             << " => nz   = " << GeoLayer->get_zbins()
             << " => ntot = " << GeoLayer->get_phibins() * GeoLayer->get_zbins()
             << endl;
      }
    }

    cout << " => nTracks = ";
    if (_trackmap)
    {
      cout << _trackmap->size() << endl;
    }
    else
    {
      cout << 0 << endl;
    }

    cout << endl;

  }  // if Verbosity()

  return;
}

void TrkrNtuplizer::fillOutputNtuples(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    cout << "TrkrNtuplizer::fillOutputNtuples() entered" << endl;
  }

  ActsGeometry* tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tgeometry)
  {
    std::cout << "No Acts geometry on node tree. Can't  continue."
              << std::endl;
    return;
  }

  float nhit_tpc_all = 0;
  float nhit_tpc_in = 0;
  float nhit_tpc_mid = 0;
  float nhit_tpc_out = 0;

  float nhit_intt = 0;
  float nhit_maps = 0;
  float nhit_mms  = 0;

  float nclus_all = 0;
  float nclus_tpc = 0;
  float nclus_intt = 0;
  float nclus_maps = 0;
  float nclus_mms = 0;
  float nhit[100];
  for (float& i : nhit) i = 0;
  float occ11 = 0;
  float occ116 = 0;
  float occ21 = 0;
  float occ216 = 0;
  float occ31 = 0;
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

      TrkrHitSet::ConstRange hitrangei = hitsetiter->second->getHits();
      for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	   hitr != hitrangei.second;
	   ++hitr)
        {
          nhit[layer]++;
          
	  if(layer<_nlayers_maps){
	    nhit_maps++;
	  }
	  if(layer>=_nlayers_maps&&layer<_nlayers_maps + _nlayers_intt){
	    nhit_intt++;
	  }
	  if ((float) layer >= _nlayers_maps + _nlayers_intt &&
	      (float) layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
	    nhit_tpc_all++;
	  if((float) layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc)
	    nhit_mms++;
          if ((float) layer == _nlayers_maps + _nlayers_intt)
	    {
	      nhit_tpc_in++;
	    }
          if ((float) layer == _nlayers_maps + _nlayers_intt + _nlayers_tpc - 1)
	    {
	      nhit_tpc_out++;
	    }
          if ((float) layer == _nlayers_maps + _nlayers_intt + _nlayers_tpc / 2 - 1)
	    {
	      nhit_tpc_mid++;
	    }
        }
    }
  }
  /**********/
  
  PHG4TpcCylinderGeomContainer* geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return;
  }
  float nhits = 0;
  /*
  PHG4TpcCylinderGeom* GeoLayer;
  int layer = _nlayers_maps + _nlayers_intt;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  int nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  float nhits = nhit[layer];
  occ11 = nhits / nbins;
  if (Verbosity() > 1)
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
  occ116 = nhits / nbins;

  layer = _nlayers_maps + _nlayers_intt + 16;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  occ21 = nhits / nbins;
  layer = _nlayers_maps + _nlayers_intt + 31;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  occ216 = nhits / nbins;
  layer = _nlayers_maps + _nlayers_intt + 32;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  occ31 = nhits / nbins;
  layer = _nlayers_maps + _nlayers_intt + 47;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  occ316 = nhits / nbins;

  if (Verbosity() > 1)
  {
    cout << " occ11 = " << occ11
         << " occ116 = " << occ116
         << " occ21 = " << occ21
         << " occ216 = " << occ216
         << " occ31 = " << occ31
         << " occ316 = " << occ316
         << endl;
  }
  */
  TrkrClusterContainer* clustermap_in = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if (!clustermap_in)
  {
    clustermap_in = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }

  if (clustermap_in)
  {
    nclus_all = clustermap_in->size();

    for (const auto& hitsetkey : clustermap_in->getHitSetKeys())
    {
      auto range = clustermap_in->getClusters(hitsetkey);
      for (auto iter_cin = range.first; iter_cin != range.second; ++iter_cin)
      {
        TrkrDefs::cluskey cluster_key = iter_cin->first;
        unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
        if (_nlayers_maps > 0)
        {
          if (layer_local < _nlayers_maps)
          {
            nclus_maps++;
          }
        }
        if (_nlayers_intt > 0)
        {
          if ((layer_local >= _nlayers_maps) && (layer_local < (_nlayers_maps + _nlayers_intt)))
          {
            nclus_intt++;
          }
        }
        if (_nlayers_tpc > 0)
        {
          if (layer_local >= (_nlayers_maps + _nlayers_intt) && layer_local < (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
          {
            nclus_tpc++;
          }
        }
        if (_nlayers_mms > 0)
        {
          if (layer_local >= (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
          {
            nclus_mms++;
          }
        }
      }
    }
  }
  //-----------------------
  // fill the info NTuple
  //-----------------------
  if (_ntp_info)
  {
    if (Verbosity() > 1)
    {
      cout << "Filling ntp_info " << endl;
    }
    float ntrk = 0;

    if (_trackmap)
    {
      ntrk = (float) _trackmap->size();
    }
    if (Verbosity() > 0)
    {
      cout << "EVENTINFO SEED: " << m_fSeed << endl;
      cout << "EVENTINFO NHIT: " << setprecision(9) << nhit_tpc_all << endl;
      cout << "EVENTINFO CLUSTPC: " << nclus_tpc << endl;
      cout << "EVENTINFO NTRKREC: " << ntrk << endl;
    }
    float info_data[] = {(float) _ievent, m_fSeed,
                         occ11, occ116, occ21, occ216, occ31, occ316,
                         ntrk,
			 nhit_maps,
			 nhit_intt,
			 nhit_mms,
                         nhit_tpc_all,
                         nhit_tpc_in,
                         nhit_tpc_mid,
                         nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

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

    //    SvtxVertexMap* vertexmap = nullptr;

    //    vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");  // Acts vertices

    float vx = NAN;
    float vy = NAN;
    float vz = NAN;
    float ntracks = NAN;

    if (Verbosity() > 1)
    {
      std::cout << " adding vertex data " << std::endl;
    }

    float vertex_data[] = {(float) _ievent, m_fSeed,
                           vx,
                           vy,
                           vz,
                           ntracks,
			   nhit_maps,
			   nhit_intt,
			   nhit_mms,
                           nhit_tpc_all,
                           nhit_tpc_in,
                           nhit_tpc_mid,
                           nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

    _ntp_vertex->Fill(vertex_data);
  }
  _timer->stop();
  cout << "vertex time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;

  //--------------------
  // fill the Hit NTuple
  //--------------------

  if (_ntp_hit)
  {
    auto m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

    if (Verbosity() >= 1)
    {
      cout << "Filling ntp_hit " << endl;
      _timer->restart();
    }
    // need things off of the DST...
    TrkrHitSetContainer* hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
 
    if (hitmap)
    {
      TrkrHitSetContainer::ConstRange all_hitsets = hitmap->getHitSets();
      for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;
           iter != all_hitsets.second;
           ++iter)
      {
        const TrkrDefs::hitsetkey hitset_key = iter->first;
        TrkrHitSet* hitset = iter->second;

        // get all hits for this hitset
        TrkrHitSet::ConstRange hitrangei = hitset->getHits();
        for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
             hitr != hitrangei.second;
             ++hitr)
        {
          TrkrDefs::hitkey hit_key = hitr->first;
          TrkrHit* hit = hitr->second;
          float event = _ievent;
          float hitID = hit_key;
          float e = hit->getEnergy();
          float adc = hit->getAdc();
          float layer_local = TrkrDefs::getLayer(hitset_key);
          float sector = TpcDefs::getSectorId(hitset_key);
          float side = TpcDefs::getSide(hitset_key);
          float cellID = 0;
          float ecell = hit->getAdc();

          float phibin = NAN;
          float tbin = NAN;
          float phi = NAN;
	  float phi_center = NAN;
          float x = NAN;
          float y = NAN;
          float z = NAN;
        
          if (TrkrDefs::getTrkrId(hitset_key) == TrkrDefs::TrkrId::tpcId)
          {
            PHG4TpcCylinderGeom* GeoLayer_local = geom_container->GetLayerCellGeom(layer_local);
 	    double radius = GeoLayer_local->get_radius();
            phibin = (float) TpcDefs::getPad(hit_key);
            tbin = (float) TpcDefs::getTBin(hit_key);
            phi = GeoLayer_local->get_phicenter(phibin);

            double zdriftlength = tbin * m_tGeometry->get_drift_velocity() * AdcClockPeriod;
            // convert z drift length to z position in the TPC
            //		cout << " tbin: " << tbin << " vdrift " <<m_tGeometry->get_drift_velocity() << " l drift: " << zdriftlength  <<endl;
            unsigned short NTBins = (unsigned short) GeoLayer_local->get_zbins();
            double m_tdriftmax = AdcClockPeriod * NTBins / 2.0;
            double clusz = (m_tdriftmax * m_tGeometry->get_drift_velocity()) - zdriftlength;
            if (side == 0)
            {
              clusz = -clusz;
            }
            z = clusz;
   phi_center = GeoLayer_local->get_phicenter(phibin);
	    x = radius * cos(phi_center);
	    y = radius * sin(phi_center);
          }

          float hit_data[] = {
              event,
              (float) _iseed,
              hitID,
              e,
              adc,
              layer_local,
              sector,
              side,
              cellID,
              ecell,
              (float) phibin,
              (float) tbin,
              phi,
	      x,
	      y,
              z,
	      nhit_maps,
	      nhit_intt,
	      nhit_mms,
              nhit_tpc_all,
              nhit_tpc_in,
              nhit_tpc_mid,
              nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

          _ntp_hit->Fill(hit_data);
        }
      }
    }
    if (Verbosity() >= 1)
    {
      _timer->stop();
      cout << "hit time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
  }

  //------------------------
  // fill the Cluster NTuple
  //------------------------

  if (Verbosity() >= 1)
  {
    cout << "check for ntp_cluster" << endl;
    _timer->restart();
  }

  if (_ntp_cluster)
  {
    if (Verbosity() > 1)
    {
      cout << "Filling ntp_cluster (all of them) " << endl;
    }
    // need things off of the DST...
    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
    if (!clustermap)
    {
      clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    }

    TrkrHitSetContainer* hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    //    TrkrClusterIterationMapv1* _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    ClusterErrorPara ClusErrPara;

    if (Verbosity() > 1)
    {
      if (clustermap != nullptr)
      {
        cout << "got clustermap" << endl;
      }
      else
      {
        cout << "no clustermap" << endl;
      }

      if (hitsets != nullptr)
      {
        cout << "got hitsets" << endl;
      }
      else
      {
        cout << "no hitsets" << endl;
      }
    }

    if (clustermap && hitsets)
    {
      for (const auto& hitsetkey : clustermap->getHitSetKeys())
      {
        int hitsetlayer = TrkrDefs::getLayer(hitsetkey);
        auto range = clustermap->getClusters(hitsetkey);
        for (auto iter = range.first; iter != range.second; ++iter)
        {
          TrkrDefs::cluskey cluster_key = iter->first;
          TrkrCluster* cluster = clustermap->findCluster(cluster_key);
	  SvtxTrack* track = nullptr;//best_track_from(cluster_key);
          float niter = 0;
         
          float hitID = (float) cluster_key;
          // auto trkrid = TrkrDefs::getTrkrId(cluster_key);
          Acts::Vector3 cglob;
          cglob = tgeometry->getGlobalPosition(cluster_key, cluster);
          float x = cglob(0);
          float y = cglob(1);
          float z = cglob(2);
	  float locx = cluster->getLocalX();
	  float locy = cluster->getLocalY(); 
          TVector3 pos(x, y, z);
          float r = pos.Perp();
          float phi = pos.Phi();
          float eta = pos.Eta();
          float theta = pos.Theta();

          float ex = 0;
          float ey = 0;
          float ez = 0;
          float ephi = 0;
          float pez = 0;
          float pephi = 0;
          float size = 0;
          float phisize = 0;
          float zsize = 0;
          float maxadc = -999;
         
	  auto para_errors = ClusErrPara.get_clusterv5_modified_error(cluster,r ,cluster_key);

	  phisize = cluster->getPhiSize();
	  zsize = cluster->getZSize();
	  // double clusRadius = r;
	  ez = sqrt(para_errors.second);
	  ephi = sqrt(para_errors.first);
	  maxadc = cluster->getMaxAdc();
          
          float e = cluster->getAdc();
          float adc = cluster->getAdc();
          float layer_local = (float) TrkrDefs::getLayer(cluster_key);
          float sector = TpcDefs::getSectorId(cluster_key);
          float side = TpcDefs::getSide(cluster_key);
	  phisize = cluster->getPhiSize();
	  zsize = cluster->getZSize();
	  size = cluster->getSize();
	  // double clusRadius = r;
	  ez = sqrt(para_errors.second);
	  ephi = sqrt(para_errors.first);
	  maxadc = cluster->getMaxAdc();
          // count all hits for this cluster
          const TrkrDefs::hitsetkey hitsetkey_local = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);
          int hitsetlayer2 = TrkrDefs::getLayer(hitsetkey_local);
          if (hitsetlayer != layer_local)
          {
            cout << "WARNING hitset layer " << hitsetlayer << "| " << hitsetlayer2 << " layer " << layer_local << endl;
          }
          /*else{
            cout << "Good    hitset layer " << hitsetlayer << "| " << hitsetlayer2 << " layer " << layer << endl;
          }
          */

          float trackID = NAN;
          if (track != nullptr)
          {
            trackID = track->get_id();
          }
	  float phibin = NAN;
	  float tbin  = NAN;
          if (TrkrDefs::getTrkrId(hitsetkey_local) == TrkrDefs::TrkrId::tpcId)
          {
	    PHG4TpcCylinderGeom* GeoLayer_local = geom_container->GetLayerCellGeom(layer_local);
	    phibin = GeoLayer_local->get_pad_float(phi,side);
	    tbin   = GeoLayer_local->get_tbin_float(locy-39.6);
	  }
	  else{
	    phibin = locx;
	    tbin = locy;
	  }
          float cluster_data[] = {(float) _ievent,
                                  (float) _iseed,
                                  hitID,
				  locx,
				  locy,
                                  x,
                                  y,
                                  z,
                                  r,
                                  phi,
                                  eta,
                                  theta,
				  phibin,
				  tbin,
                                  ex,
                                  ey,
                                  ez,
                                  ephi,
                                  pez,
                                  pephi,
                                  e,
                                  adc,
                                  maxadc,
                                  layer_local,
                                  sector,
                                  side,
                                  size,
                                  phisize,
                                  zsize,
                                  trackID,
                                  niter,
				  nhit_maps,
				  nhit_intt,
				  nhit_mms,
                                  nhit_tpc_all,
                                  nhit_tpc_in,
                                  nhit_tpc_mid,
                                  nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

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
  if (_ntp_clus_trk)
  {
    // if (Verbosity() > 1)
    {
      cout << "Filling ntp_clus_trk " << endl;
    }
    // need things off of the DST...
    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
    if (!clustermap)
    {
      clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    }
 

    ClusterErrorPara ClusErrPara;

    if (Verbosity() > 1)
    {
      if (clustermap != nullptr) cout << "got clustermap" << endl;
      else cout << "no clustermap" << endl;
    }

    TrackSeedContainer *_tpc_seeds = findNode::getClass<TrackSeedContainer>(topNode, _clustrackseedcontainer);
    TrackSeedContainer *_tpcseeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
    TrackSeedContainer *_silseeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
    if (!_tpc_seeds)
      {
	cerr << PHWHERE << " ERROR: Can't find " << _clustrackseedcontainer << endl;
	return ;
      }
   

    for (unsigned int phtrk_iter = 0; phtrk_iter < _tpc_seeds->size(); ++phtrk_iter){

      TrackSeed *tpcseed = _tpc_seeds->get(phtrk_iter);
      if(!tpcseed) { 
	continue; 
      }  
      TrackSeed *silseed = nullptr;
      if(_clustrackseedcontainer.find("SvtxTrackSeed") != std::string::npos)
	{
	  unsigned int tpcindex = tpcseed->get_tpc_seed_index();
	  unsigned int silindex = tpcseed->get_silicon_seed_index();
	  tpcseed = _tpcseeds->get(tpcindex);
	  silseed = _silseeds->get(silindex);
	}

    
      std::vector<Acts::Vector3> clusterPositions;
      std::vector<TrkrDefs::cluskey> clusterKeys;
      clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
			 tpcseed->end_cluster_keys());
      if(silseed)
	{
	  clusterKeys.insert(clusterKeys.end(), silseed->begin_cluster_keys(),
			     silseed->end_cluster_keys());
	}
      TrackFitUtils::getTrackletClusters(tgeometry, clustermap, 
					 clusterPositions, clusterKeys);
      std::vector<float> fitparams = TrackFitUtils::fitClusters(clusterPositions, clusterKeys);
 
      float charge = NAN;
      if(tpcseed->get_qOverR()>0)
	{ charge = 1; }
      else
	{ charge = -1; }
      //	      "pt:eta:phi:X0:Y0:charge:nhits:"
      float tpt = tpcseed->get_pt();
      float teta = tpcseed->get_eta();
      float tphi = tpcseed->get_phi(clustermap,tgeometry);
      float tX0 = tpcseed->get_X0();
      float tY0 = tpcseed->get_Y0();
      float tZ0 = tpcseed->get_Z0();
      
      float nhits_local = 0;
      nhits_local += tpcseed->size_cluster_keys();
      if(silseed)
	{
	  nhits_local += silseed->size_cluster_keys();
	}

      for (size_t i = 0; i < clusterPositions.size(); i++)
	{/*
	   for (SvtxTrack::ConstClusterKeyIter iter_local = tpcseed->begin_cluster_keys();
	   iter_local != tpcseed->end_cluster_keys();
	   ++iter_local){
	 */

	  const TrkrDefs::cluskey cluster_key = clusterKeys[i];//*iter_local;
	  // TrkrCluster* cluster = clustermap->findCluster(cluster_key);
	  // unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
	   Acts::Vector3 position = clusterPositions[i];
	   Acts::Vector3 pca = TrackFitUtils::get_helix_pca(fitparams,position);
	   float cluster_phi = atan2(position(1), position(0));
	   float pca_phi = atan2(pca(1), pca(0));
	   float dphi = cluster_phi - pca_phi;
	   if (dphi > M_PI)
	     {
	       dphi = 2 * M_PI - dphi;
	     }
	   if (dphi < -M_PI)
	     {
	       dphi = 2 * M_PI + dphi;
	     }
	   float dz = position(2) - pca(2);
	   
	   TrkrCluster* cluster = clustermap->findCluster(cluster_key);
	   Acts::Vector3 cglob;
	   cglob = tgeometry->getGlobalPosition(cluster_key, cluster);
	   float x = cglob(0);
	   float y = cglob(1);
	   float z = cglob(2);
	   float locx = cluster->getLocalX();
	   float locy = cluster->getLocalY();
	   TVector3 pos(x, y, z);
	   float r = pos.Perp();
	   float phi = pos.Phi();
	   float eta = pos.Eta();
	   float theta = pos.Theta();
	   
	   float ex = 0;
	   float ey = 0;
	   float ez = 0;
	   float ephi = 0;
	   float pez = 0;
	   float pephi = 0;
	   float size = 0;
	   float phisize = 0;
	   float zsize = 0;
	   float maxadc = -999;
	   
	   auto para_errors = ClusErrPara.get_clusterv5_modified_error(cluster,r ,cluster_key);
	   
	   phisize = cluster->getPhiSize();
	   zsize = cluster->getZSize();
	   // double clusRadius = r;
	   ez = sqrt(para_errors.second);
	   ephi = sqrt(para_errors.first);
	   maxadc = cluster->getMaxAdc();
	   
	   float e = cluster->getAdc();
	   float adc = cluster->getAdc();
	   float layer_local = (float) TrkrDefs::getLayer(cluster_key);
	   float sector = TpcDefs::getSectorId(cluster_key);
	   float side = TpcDefs::getSide(cluster_key);
	   
	   float trackID = phtrk_iter;
	   float phibin = NAN;
	   float tbin  = NAN;
	   const TrkrDefs::hitsetkey hsk = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);
	   if (TrkrDefs::getTrkrId(hsk) == TrkrDefs::TrkrId::tpcId)
	     {
	       PHG4TpcCylinderGeom* GeoLayer_local = geom_container->GetLayerCellGeom(layer_local);
	       phibin = GeoLayer_local->get_pad_float(phi,side);
	       tbin   = GeoLayer_local->get_tbin_float(locy-39.6);
	     }
	   else{
	     phibin = locx;
	     tbin = locy;
	   }
	   /*
	     "pt:eta:phi:X0:Y0:charge:nhits:"
	   */
	  float clus_trk_data[] = {(float) _ievent,
				 (float) _iseed,
				 locx,
				 locy,
				 x,
				 y,
				 z,
				 r,
				 phi,
				 eta,
				 theta,
				 phibin,
				 tbin,
				 ex,
				 ey,
				 ez,
				 ephi,
				 pez,
				 pephi,
				 e,
				 adc,
				 maxadc,
				 layer_local,
				 sector,
				 side,
				 size,
				 phisize,
				 zsize,
				   dphi,
				   dz,
				   trackID,
				   0,
				   tpt,
				   teta,
				   tphi,
				   tX0,
				   tY0,
				   tZ0,
				   charge,
				   nhit_maps,
				   nhit_intt,
				   nhit_mms,
				   nhits_local,
				   nhit_tpc_all,
				   nhit_tpc_in,
				   nhit_tpc_mid,
				   nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};
	  
	  _ntp_clus_trk->Fill(clus_trk_data);
	  
	}

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

    GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

    if (_trackmap)
    {
      for (auto& iter : *_trackmap)
      {
        SvtxTrack* track = iter.second;
	if(!track) continue;
        float trackID = track->get_id();
        TrackSeed* tpcseed = track->get_tpc_seed();
        TrackSeed* silseed = track->get_silicon_seed();
        short int crossing_int = track->get_crossing();
        float crossing;
        if (crossing_int == SHRT_MAX)
        {
          crossing = NAN;
        }
        else
        {
          crossing = (float) crossing_int;
        }
        float charge = track->get_charge();
        float quality = track->get_quality();
        float chisq = track->get_chisq();
        float ndf = track->get_ndf();
        float nhits_local = 0;
        if (tpcseed)
        {
          nhits_local += tpcseed->size_cluster_keys();
        }
        if (silseed)
        {
          nhits_local += silseed->size_cluster_keys();
        }
        unsigned int layers = 0x0;
        int maps[_nlayers_maps];
        int intt[_nlayers_intt];
        int tpc[_nlayers_tpc];
        int mms[_nlayers_mms];
        if (_nlayers_maps > 0)
        {
          for (unsigned int i = 0; i < _nlayers_maps; i++)
          {
            maps[i] = 0;
          }
        }
        if (_nlayers_intt > 0)
        {
          for (unsigned int i = 0; i < _nlayers_intt; i++)
          {
            intt[i] = 0;
          }
        }
        if (_nlayers_tpc > 0)
        {
          for (unsigned int i = 0; i < _nlayers_tpc; i++)
          {
            tpc[i] = 0;
          }
        }
        if (_nlayers_mms > 0)
        {
          for (unsigned int i = 0; i < _nlayers_mms; i++)
          {
            mms[i] = 0;
          }
        }

        float nmaps = 0;
        float nintt = 0;
        float nmms = 0;
        float ntpc = 0;
        float ntpc1 = 0;
        float ntpc11 = 0;
        float ntpc2 = 0;
        float ntpc3 = 0;
        float nlmaps = 0;
        float nlintt = 0;
        float nltpc = 0;
        float nlmms = 0;

        if (tpcseed)
        {
          for (SvtxTrack::ConstClusterKeyIter iter_local = tpcseed->begin_cluster_keys();
               iter_local != tpcseed->end_cluster_keys();
               ++iter_local)
          {
            TrkrDefs::cluskey cluster_key = *iter_local;
            // TrkrCluster* cluster = clustermap->findCluster(cluster_key);
            unsigned int layer_local = TrkrDefs::getLayer(cluster_key);

            if (_nlayers_maps > 0 && layer_local < _nlayers_maps)
            {
              maps[layer_local] = 1;
              nmaps++;
            }
            if (_nlayers_intt > 0 && layer_local >= _nlayers_maps && layer_local < _nlayers_maps + _nlayers_intt)
            {
              intt[layer_local - _nlayers_maps] = 1;
              nintt++;
            }
            if (_nlayers_mms > 0 && layer_local >= _nlayers_maps + _nlayers_intt + _nlayers_tpc && layer_local < _nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms)
            {
              mms[layer_local - (_nlayers_maps + _nlayers_intt + _nlayers_tpc)] = 1;
              nmms++;
            }
            if (_nlayers_tpc > 0 && layer_local >= (_nlayers_maps + _nlayers_intt) && layer_local < (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
            {
              tpc[layer_local - (_nlayers_maps + _nlayers_intt)] = 1;
              ntpc++;
              if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 8)
              {
                ntpc11++;
              }

              if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 16)
              {
                ntpc1++;
              }
              else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 32)
              {
                ntpc2++;
              }
              else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 48)
              {
                ntpc3++;
              }
            }
          }
        }

        if (silseed)
        {
          for (SvtxTrack::ConstClusterKeyIter iter_local = silseed->begin_cluster_keys();
               iter_local != silseed->end_cluster_keys();
               ++iter_local)
          {
            TrkrDefs::cluskey cluster_key = *iter_local;
            // TrkrCluster* cluster = clustermap->findCluster(cluster_key);
            unsigned int layer_local = TrkrDefs::getLayer(cluster_key);

            if (_nlayers_maps > 0 && layer_local < _nlayers_maps)
            {
              maps[layer_local] = 1;
              nmaps++;
            }
            if (_nlayers_intt > 0 && layer_local >= _nlayers_maps && layer_local < _nlayers_maps + _nlayers_intt)
            {
              intt[layer_local - _nlayers_maps] = 1;
              nintt++;
            }
            if (_nlayers_mms > 0 && layer_local >= _nlayers_maps + _nlayers_intt + _nlayers_tpc && layer_local < _nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms)
            {
              mms[layer_local - (_nlayers_maps + _nlayers_intt + _nlayers_tpc)] = 1;
              nmms++;
            }
            if (_nlayers_tpc > 0 && layer_local >= (_nlayers_maps + _nlayers_intt) && layer_local < (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
            {
              tpc[layer_local - (_nlayers_maps + _nlayers_intt)] = 1;
              ntpc++;
              if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 8)
              {
                ntpc11++;
              }

              if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 16)
              {
                ntpc1++;
              }
              else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 32)
              {
                ntpc2++;
              }
              else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 48)
              {
                ntpc3++;
              }
            }
          }
        }

        if (_nlayers_maps > 0)
        {
          for (unsigned int i = 0; i < _nlayers_maps; i++)
          {
            nlmaps += maps[i];
          }
        }
        if (_nlayers_intt > 0)
        {
          for (unsigned int i = 0; i < _nlayers_intt; i++)
          {
            nlintt += intt[i];
          }
        }
        if (_nlayers_tpc > 0)
        {
          for (unsigned int i = 0; i < _nlayers_tpc; i++)
          {
            nltpc += tpc[i];
          }
        }
        if (_nlayers_mms > 0)
        {
          for (unsigned int i = 0; i < _nlayers_mms; i++)
          {
            nlmms += mms[i];
          }
        }
        layers = nlmaps + nlintt + nltpc + nlmms;

        float dca3dxy = NAN, dca3dz = NAN,
              dca3dxysigma = NAN, dca3dzsigma = NAN;
        float dca2d = NAN, dca2dsigma = NAN;

        int vertexID = track->get_vertex_id();
	float vx = NAN;
        float vy = NAN;
        float vz = NAN;
	if(vertexmap) 
	  {
        GlobalVertex* vertex = vertexmap->get(vertexID);
  
        if (vertex)
        {
          vx = vertex->get_x();
          vy = vertex->get_y();
          vz = vertex->get_z();
	  Acts::Vector3 vert(vx,vy,vz);
          auto dcapair = TrackAnalysisUtils::get_dca(track, vert);
	  dca3dxy = dcapair.first.first;
	  dca3dxysigma = dcapair.first.second;
	  dca3dz = dcapair.second.first;
	  dca3dzsigma = dcapair.second.second;
        }
	  }
        float px = track->get_px();
        float py = track->get_py();
        float pz = track->get_pz();
        TVector3 v(px, py, pz);
        float pt = v.Pt();
        float eta = v.Eta();
        float phi = v.Phi();
        float CVxx = track->get_error(3, 3);
        float CVxy = track->get_error(3, 4);
        float CVxz = track->get_error(3, 5);
        float CVyy = track->get_error(4, 4);
        float CVyz = track->get_error(4, 5);
        float CVzz = track->get_error(5, 5);
        float deltapt = sqrt((CVxx * px * px + 2 * CVxy * px * py + CVyy * py * py) / (px * px + py * py));
        float deltaeta = sqrt((CVzz * (px * px + py * py) * (px * px + py * py) + pz * (-2 * (CVxz * px + CVyz * py) * (px * px + py * py) + CVxx * px * px * pz + CVyy * py * py * pz + 2 * CVxy * px * py * pz)) / ((px * px + py * py) * (px * px + py * py) * (px * px + py * py + pz * pz)));
        float deltaphi = sqrt((CVyy * px * px - 2 * CVxy * px * py + CVxx * py * py) / ((px * px + py * py) * (px * px + py * py)));
        float pcax = track->get_x();
        float pcay = track->get_y();
        float pcaz = track->get_z();

        float track_data[] = {(float) _ievent, m_fSeed,
                              trackID,
                              crossing,
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
                              nhits_local, nmaps, nintt, ntpc, nmms,
                              ntpc1, ntpc11, ntpc2, ntpc3,
                              nlmaps, nlintt, nltpc, nlmms,
                              (float) layers,
                              (float) vertexID,
                              vx,
                              vy,
                              vz,
                              dca2d,
                              dca2dsigma,
                              dca3dxy,
                              dca3dxysigma,
                              dca3dz,
                              dca3dzsigma,
                              pcax,
                              pcay,
                              pcaz,
			      nhit_maps,
			      nhit_intt,
			      nhit_mms,
                              nhit_tpc_all,
                              nhit_tpc_in,
                              nhit_tpc_mid,
                              nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

        if (Verbosity() >= 1)
        {
          cout << "ievent " << _ievent
               << " trackID " << trackID
               << " nhits " << nhits
               << " px " << px
               << " py " << py
               << " pz " << pz
               << endl;
        }

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

  if (_ntp_tpcseed)
  {
    if (Verbosity() > 1)
    {
      cout << "Filling ntp_tpcseed " << endl;
      _timer->restart();
    }

    if (Verbosity() > 1)
    {
      _timer->stop();
      cout << "tpcseed time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
  }

  if (_ntp_siseed)
  {
    if (Verbosity() > 1)
    {
      cout << "Filling ntp_tpcseed " << endl;
      _timer->restart();
    }

    if (Verbosity() > 1)
    {
      _timer->stop();
      cout << "tpcseed time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
  }
  return;
}

std::vector<TrkrDefs::cluskey> TrkrNtuplizer::get_track_ckeys(SvtxTrack* track)
{
  std::vector<TrkrDefs::cluskey> cluster_keys;
  TrackSeed *tpcseed = track->get_tpc_seed();
  TrackSeed *silseed = track->get_silicon_seed();
  if(silseed)
    {
      for(auto iter = silseed->begin_cluster_keys();
	  iter!= silseed->end_cluster_keys();
	  ++iter)
	{ cluster_keys.push_back(*iter); }
    }
  if(tpcseed)
    {
      for(auto iter = tpcseed->begin_cluster_keys();
          iter!= tpcseed->end_cluster_keys();
          ++iter)
        { cluster_keys.push_back(*iter); }
    }
  
  return cluster_keys;
}

SvtxTrack* TrkrNtuplizer::best_track_from(TrkrDefs::cluskey cluster_key)
{
  std::map<TrkrDefs::cluskey, SvtxTrack*>::iterator find_iter =
    _cache_best_track_from_cluster.find(cluster_key);
  if (find_iter != _cache_best_track_from_cluster.end())
    {
      return find_iter->second;
    }

  SvtxTrack* best_track = nullptr;
  float best_quality = FLT_MAX;

  std::set<SvtxTrack*> tracks = all_tracks_from(cluster_key);
  // loop over all SvtxTracks
  for (std::set<SvtxTrack*>::iterator iter = tracks.begin();
       iter != tracks.end();
       ++iter)
  {
    SvtxTrack* candidate = *iter;
    if (candidate->get_quality() < best_quality)
    {
      best_quality = candidate->get_quality();
      best_track = candidate;
    }
  }

  _cache_best_track_from_cluster.insert(make_pair(cluster_key, best_track));
  return best_track;
}

std::set<SvtxTrack*> TrkrNtuplizer::all_tracks_from(TrkrDefs::cluskey cluster_key)
{
  std::set<SvtxTrack*> tracks;

  if (_cache_track_from_cluster_exists == false) create_cache_track_from_cluster();
  std::map<TrkrDefs::cluskey, std::set<SvtxTrack*> >::iterator find_iter =
    _cache_all_tracks_from_cluster.find(cluster_key);
  if (find_iter != _cache_all_tracks_from_cluster.end())
    {
      return find_iter->second;
    }
  else
    {
      return tracks;
    }

  // loop over all SvtxTracks
  for (SvtxTrackMap::Iter iter = _trackmap->begin();
       iter != _trackmap->end();
       ++iter)
  {
    SvtxTrack* track = iter->second;
    std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);

    // loop over all clusters
    for (const auto& candidate : cluster_keys)
    {

      //      if (_strict)
      //      {
      //        assert(candidate);
      //      }
      //      else if (!candidate)
      //      {
      //        ++_errors;
      //        continue;
      //      }

      if (cluster_key == candidate)
      {
        tracks.insert(track);
      }
    }
  }

  _cache_all_tracks_from_cluster.insert(make_pair(cluster_key, tracks));

  return tracks;
}

void TrkrNtuplizer::create_cache_track_from_cluster()
{

  if(!_trackmap) return;

  // loop over all SvtxTracks
  for (SvtxTrackMap::Iter iter = _trackmap->begin();
       iter != _trackmap->end();
       ++iter)
  {
    SvtxTrack* track = iter->second;
    std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);

    // loop over all clusters
    for (const auto& candidate_key : cluster_keys)
    {
      //check if cluster has an entry in cache
      std::map<TrkrDefs::cluskey, std::set<SvtxTrack*> >::iterator cliter =
          _cache_all_tracks_from_cluster.find(candidate_key);
      if (cliter != _cache_all_tracks_from_cluster.end())
      {                                //got entry
        cliter->second.insert(track);  //add track to list;
      }
      else
      {
        std::set<SvtxTrack*> tracks;
        tracks.insert(track);
        _cache_all_tracks_from_cluster.insert(make_pair(candidate_key, tracks));
      }
    }
  }
  _cache_track_from_cluster_exists = true;

  return;
}

TMatrixF TrkrNtuplizer::calculateClusterError(TrkrCluster* c, float& clusphi)
{
  TMatrixF localErr(3, 3);
  localErr[0][0] = 0.;
  localErr[0][1] = 0.;
  localErr[0][2] = 0.;
  localErr[1][0] = 0.;
  localErr[1][1] = c->getActsLocalError(0, 0);
  localErr[1][2] = c->getActsLocalError(0, 1);
  localErr[2][0] = 0.;
  localErr[2][1] = c->getActsLocalError(1, 0);
  localErr[2][2] = c->getActsLocalError(1, 1);

  TMatrixF ROT(3, 3);
  ROT[0][0] = cos(clusphi);
  ROT[0][1] = -sin(clusphi);
  ROT[0][2] = 0.0;
  ROT[1][0] = sin(clusphi);
  ROT[1][1] = cos(clusphi);
  ROT[1][2] = 0.0;
  ROT[2][0] = 0.0;
  ROT[2][1] = 0.0;
  ROT[2][2] = 1.0;
  TMatrixF ROT_T(3, 3);
  ROT_T.Transpose(ROT);

  TMatrixF err(3, 3);
  err = ROT * localErr * ROT_T;
  return err;
}
