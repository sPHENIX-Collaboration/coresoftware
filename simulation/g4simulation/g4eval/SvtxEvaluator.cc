#include "SvtxEvaluator.h"

#include "SvtxEvalStack.h"

#include "SvtxClusterEval.h"
#include "SvtxHitEval.h"
#include "SvtxTrackEval.h"
#include "SvtxTruthEval.h"
#include "SvtxVertexEval.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackAnalysisUtils.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

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

SvtxEvaluator::SvtxEvaluator(const std::string& /*name*/, const std::string& filename, const std::string& trackmapname,
                             unsigned int nlayers_maps,
                             unsigned int nlayers_intt,
                             unsigned int nlayers_tpc,
                             unsigned int nlayers_mms)
  : SubsysReco("SvtxEvaluator")
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _nlayers_mms(nlayers_mms)
  , _filename(filename)
  , _trackmapname(trackmapname)
{
}

SvtxEvaluator::~SvtxEvaluator()
{
  delete _timer;
}

int SvtxEvaluator::Init(PHCompositeNode* /*topNode*/)
{
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");
  _tfile->SetCompressionLevel(0);
  if (_do_info_eval)
  {
    _ntp_info = new TNtuple("ntp_info", "event info",
                            "event:seed:"
                            "occ11:occ116:occ21:occ216:occ31:occ316:"
                            "gntrkall:gntrkprim:ntrk:"
                            "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_vertex_eval)
  {
    _ntp_vertex = new TNtuple("ntp_vertex", "vertex => max truth",
                              "event:seed:vertexID:vx:vy:vz:ntracks:chi2:ndof:"
                              "gvx:gvy:gvz:gvt:gembed:gntracks:gntracksmaps:"
                              "gnembed:nfromtruth:"
                              "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_gpoint_eval)
  {
    _ntp_gpoint = new TNtuple("ntp_gpoint", "g4point => best vertex",
                              "event:seed:gvx:gvy:gvz:gvt:gntracks:gembed:"
                              "vx:vy:vz:ntracks:"
                              "nfromtruth:"
                              "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_g4hit_eval)
  {
    _ntp_g4hit = new TNtuple("ntp_g4hit", "g4hit => best svtxcluster",
                             "event:seed:g4hitID:gx:gy:gz:gt:gpl:gedep:geta:gphi:"
                             "gdphi:gdz:"
                             "glayer:gtrackID:gflavor:"
                             "gpx:gpy:gpz:"
                             "gvx:gvy:gvz:"
                             "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                             "gembed:gprimary:nclusters:"
                             "clusID:x:y:z:eta:phi:e:adc:layer:size:"
                             "efromtruth:dphitru:detatru:dztru:drtru:"
                             "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_hit_eval)
  {
    _ntp_hit = new TNtuple("ntp_hit", "svtxhit => max truth",
                           "event:seed:hitID:e:adc:layer:phielem:zelem:"
                           "cellID:ecell:phibin:zbin:phi:z:"
                           "g4hitID:gedep:gx:gy:gz:gt:"
                           "gtrackID:gflavor:"
                           "gpx:gpy:gpz:gvx:gvy:gvz:gvt:"
                           "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                           "gembed:gprimary:efromtruth:"
                           "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_cluster_eval)
  {
    _ntp_cluster = new TNtuple("ntp_cluster", "svtxcluster => max truth",
                               "event:seed:hitID:x:y:z:r:phi:eta:theta:ex:ey:ez:ephi:pez:pephi:"
                               "e:adc:maxadc:layer:phielem:zelem:size:phisize:zsize:"
			       "pedge:redge:ovlp:"
                               "trackID:niter:g4hitID:gx:"
                               "gy:gz:gr:gphi:geta:gt:gtrackID:gflavor:"
                               "gpx:gpy:gpz:gvx:gvy:gvz:gvt:"
                               "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                               "gembed:gprimary:efromtruth:nparticles:"
                               "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_g4cluster_eval)
  {
    _ntp_g4cluster = new TNtuple("ntp_g4cluster", "g4cluster => max truth",
                                 "event:layer:gx:gy:gz:gt:gedep:gr:gphi:geta:gtrackID:gflavor:gembed:gprimary:gphisize:gzsize:gadc:nreco:x:y:z:r:phi:eta:ex:ey:ez:ephi:adc:phisize:zsize");
  }

  if (_do_gtrack_eval)
  {
    _ntp_gtrack = new TNtuple("ntp_gtrack", "g4particle => best svtxtrack",
                              "event:seed:gntracks:gnchghad:gtrackID:gflavor:gnhits:gnmaps:gnintt:gnmms:"
                              "gnintt1:gnintt2:gnintt3:gnintt4:"
                              "gnintt5:gnintt6:gnintt7:gnintt8:"
                              "gntpc:gnlmaps:gnlintt:gnltpc:gnlmms:"
                              "gpx:gpy:gpz:gpt:geta:gphi:"
                              "gvx:gvy:gvz:gvt:"
                              "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                              "gembed:gprimary:"
                              "trackID:px:py:pz:pt:eta:phi:deltapt:deltaeta:deltaphi:"
			      "siqr:siphi:sithe:six0:siy0:tpqr:tpphi:tpthe:tpx0:tpy0:"
                              "charge:quality:chisq:ndf:nhits:layers:nmaps:nintt:ntpc:nmms:ntpc1:ntpc11:ntpc2:ntpc3:nlmaps:nlintt:nltpc:nlmms:"
                              "vertexID:vx:vy:vz:dca2d:dca2dsigma:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:nfromtruth:nwrong:ntrumaps:nwrongmaps:ntruintt:nwrongintt:ntrutpc:nwrongtpc:ntrumms:nwrongmms:ntrutpc1:nwrongtpc1:ntrutpc11:nwrongtpc11:ntrutpc2:nwrongtpc2:ntrutpc3:nwrongtpc3:layersfromtruth:"
			      "npedge:nredge:nbig:novlp:merr:msize:"
                              "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_track_eval)
  {
    _ntp_track = new TNtuple("ntp_track", "svtxtrack => max truth",
                             "event:seed:trackID:crossing:px:py:pz:pt:eta:phi:deltapt:deltaeta:deltaphi:"
                             "siqr:siphi:sithe:six0:siy0:tpqr:tpphi:tpthe:tpx0:tpy0:"
			     "charge:quality:chisq:ndf:nhits:nmaps:nintt:ntpc:nmms:ntpc1:ntpc11:ntpc2:ntpc3:nlmaps:nlintt:nltpc:nlmms:layers:"
                             "vertexID:vx:vy:vz:dca2d:dca2dsigma:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:"
                             "gtrackID:singlematch:gflavor:gnhits:gnmaps:gnintt:gntpc:gnmms:gnlmaps:gnlintt:gnltpc:gnlmms:"
                             "gpx:gpy:gpz:gpt:geta:gphi:"
                             "gvx:gvy:gvz:gvt:"
                             "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
                             "gembed:gprimary:nfromtruth:nwrong:ntrumaps:nwrongmaps:ntruintt:nwrongintt:"
                             "ntrutpc:nwrongtpc:ntrumms:nwrongmms:ntrutpc1:nwrongtpc1:ntrutpc11:nwrongtpc11:ntrutpc2:nwrongtpc2:ntrutpc3:nwrongtpc3:layersfromtruth:"
			     "npedge:nredge:nbig:novlp:merr:msize:"
                             "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  if (_do_gseed_eval)
  {
    _ntp_gseed = new TNtuple("ntp_gseed", "seeds from truth",
                             "event:seed:ntrk:gx:gy:gz:gr:geta:gphi:"
                             "glayer:"
                             "gpx:gpy:gpz:gtpt:gtphi:gteta:"
                             "gvx:gvy:gvz:"
                             "gembed:gprimary:gflav:"
                             "dphiprev:detaprev:"
                             "nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms");
  }

  _timer = new PHTimer("_eval_timer");
  _timer->stop();

  return Fun4AllReturnCodes::EVENT_OK;
}

int SvtxEvaluator::InitRun(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int SvtxEvaluator::process_event(PHCompositeNode* topNode)
{
  if ((Verbosity() > 1) && (_ievent % 100 == 0))
  {
    std::cout << "SvtxEvaluator::process_event - Event = " << _ievent << std::endl;
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

  if (Verbosity() > 1)
  {
    std::cout << "SvtxEvaluator::process_event - Seed = " << _iseed << std::endl;
  }

  if (!_svtxevalstack)
  {
    _svtxevalstack = new SvtxEvalStack(topNode);
    _svtxevalstack->set_strict(_strict);
    _svtxevalstack->set_verbosity(Verbosity());
    _svtxevalstack->set_use_initial_vertex(_use_initial_vertex);
    _svtxevalstack->set_use_genfit_vertex(_use_genfit_vertex);
    _svtxevalstack->next_event(topNode);
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

  // printOutputInfo(topNode);

  ++_ievent;
  return Fun4AllReturnCodes::EVENT_OK;
}

int SvtxEvaluator::End(PHCompositeNode* /*topNode*/)
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
  if (_ntp_gpoint)
  {
    _ntp_gpoint->Write();
  }
  if (_ntp_g4hit)
  {
    _ntp_g4hit->Write();
  }
  if (_ntp_hit)
  {
    _ntp_hit->Write();
  }
  if (_ntp_cluster)
  {
    _ntp_cluster->Write();
  }
  if (_ntp_g4cluster)
  {
    _ntp_g4cluster->Write();
  }
  if (_ntp_gtrack)
  {
    _ntp_gtrack->Write();
  }
  if (_ntp_track)
  {
    _ntp_track->Write();
  }
  if (_ntp_gseed)
  {
    _ntp_gseed->Write();
  }

  _tfile->Close();

  delete _tfile;

  if (Verbosity() > 1)
  {
    std::cout << "========================= SvtxEvaluator::End() ============================" << std::endl;
    std::cout << " " << _ievent << " events of output written to: " << _filename << std::endl;
    std::cout << "===========================================================================" << std::endl;
  }

  if (_svtxevalstack)
  {
    _errors += _svtxevalstack->get_errors();

    if (Verbosity() > 1)
    {
      if ((_errors > 0) || (Verbosity() > 1))
      {
        std::cout << "SvtxEvaluator::End() - Error Count: " << _errors << std::endl;
      }
    }

    delete _svtxevalstack;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void SvtxEvaluator::printInputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "SvtxEvaluator::printInputInfo() entered" << std::endl;
  }

  if (Verbosity() > 3)
  {
    // event information
    std::cout << std::endl;
    std::cout << PHWHERE << "   INPUT FOR EVENT " << _ievent << std::endl;

    std::cout << std::endl;
    std::cout << "---PHG4HITS-------------" << std::endl;
    _svtxevalstack->get_truth_eval()->set_strict(_strict);
    std::set<PHG4Hit*> g4hits = _svtxevalstack->get_truth_eval()->all_truth_hits();
    unsigned int ig4hit = 0;
    for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
         iter != g4hits.end();
         ++iter)
    {
      PHG4Hit* g4hit = *iter;
      std::cout << ig4hit << " of " << g4hits.size();
      std::cout << ": PHG4Hit: " << std::endl;
      g4hit->identify();
      ++ig4hit;
    }

    std::cout << "---SVTXCLUSTERS-------------" << std::endl;
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
          std::cout << icluster << " with key " << cluster_key << " of " << clustermap->size();
          std::cout << ": SvtxCluster: " << std::endl;
          iter->second->identify();
          ++icluster;
        }
      }
    }

    std::cout << "---SVXTRACKS-------------" << std::endl;
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());
    if (trackmap)
    {
      unsigned int itrack = 0;
      for (SvtxTrackMap::Iter iter = trackmap->begin();
           iter != trackmap->end();
           ++iter)
      {
        std::cout << itrack << " of " << trackmap->size();
        SvtxTrack* track = iter->second;
        std::cout << " : SvtxTrack:" << std::endl;
        track->identify();
        std::cout << std::endl;
        ++itrack;
      }
    }

    std::cout << "---SVXVERTEXES-------------" << std::endl;
    SvtxVertexMap* vertexmap = nullptr;
    if (_use_initial_vertex)
    {
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    }
    else if (_use_genfit_vertex)
    {
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapRefit");
    }
    else
    {
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");  // Acts vertices
    }

    if (vertexmap)
    {
      unsigned int ivertex = 0;
      for (SvtxVertexMap::Iter iter = vertexmap->begin();
           iter != vertexmap->end();
           ++iter)
      {
        std::cout << ivertex << " of " << vertexmap->size();
        SvtxVertex* vertex = iter->second;
        std::cout << " : SvtxVertex:" << std::endl;
        vertex->identify();
        std::cout << std::endl;
      }
    }
  }

  return;
}

void SvtxEvaluator::printOutputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "SvtxEvaluator::printOutputInfo() entered" << std::endl;
  }

  //==========================================
  // print out some useful stuff for debugging
  //==========================================

  if (Verbosity() > 100)
  {
    SvtxTrackEval* trackeval = _svtxevalstack->get_track_eval();
    SvtxClusterEval* clustereval = _svtxevalstack->get_cluster_eval();
    SvtxTruthEval* trutheval = _svtxevalstack->get_truth_eval();

    // event information
    std::cout << std::endl;
    std::cout << PHWHERE << "   NEW OUTPUT FOR EVENT " << _ievent << std::endl;
    std::cout << std::endl;

    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

    PHG4VtxPoint* gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
    float gvx = gvertex->get_x();
    float gvy = gvertex->get_y();
    float gvz = gvertex->get_z();

    float vx = NAN;
    float vy = NAN;
    float vz = NAN;

    SvtxVertexMap* vertexmap = nullptr;
    if (_use_initial_vertex)
    {
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    }
    else if (_use_genfit_vertex)
    {
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapRefit");
    }
    else
    {
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");  // Acts vertices
    }

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

    std::cout << "===Vertex Reconstruction=======================" << std::endl;
    std::cout << "vtrue = (" << gvx << "," << gvy << "," << gvz << ") => vreco = (" << vx << "," << vy << "," << vz << ")" << std::endl;
    std::cout << std::endl;

    std::cout << "===Tracking Summary============================" << std::endl;
    unsigned int ng4hits[100] = {0};
    std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits();
    for (auto g4hit : g4hits)
    {
      ++ng4hits[g4hit->get_layer()];
    }

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
          unsigned int local_layer = TrkrDefs::getLayer(cluskey);
          nclusters[local_layer]++;
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

      std::cout << "layer " << ilayer << ": nG4hits = " << ng4hits[ilayer]
                << " => nHits = " << nhits[ilayer]
                << " => nClusters = " << nclusters[ilayer] << std::endl;
      if (ilayer >= _nlayers_maps + _nlayers_intt && ilayer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
      {
        std::cout << "layer " << ilayer
                  << " => nphi = " << GeoLayer->get_phibins()
                  << " => nz   = " << GeoLayer->get_zbins()
                  << " => ntot = " << GeoLayer->get_phibins() * GeoLayer->get_zbins()
                  << std::endl;
      }
    }

    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());

    std::cout << "nGtracks = " << std::distance(truthinfo->GetPrimaryParticleRange().first, truthinfo->GetPrimaryParticleRange().second);
    std::cout << " => nTracks = ";
    if (trackmap)
    {
      std::cout << trackmap->size() << std::endl;
    }
    else
    {
      std::cout << 0 << std::endl;
    }

    // cluster wise information
    if (Verbosity() > 1)
    {
      for (auto g4hit : g4hits)
      {
        std::cout << std::endl;
        std::cout << "===PHG4Hit===================================" << std::endl;
        std::cout << " PHG4Hit: ";
        g4hit->identify();

        std::set<TrkrDefs::cluskey> clusters = clustereval->all_clusters_from(g4hit);

        for (unsigned long cluster_key : clusters)
        {
          std::cout << "===Created-SvtxCluster================" << std::endl;
          std::cout << "SvtxCluster: ";
          TrkrCluster* cluster = clustermap->findCluster(cluster_key);
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
        std::cout << std::endl;

        std::cout << "=== Gtrack ===================================================" << std::endl;
        std::cout << " PHG4Particle id = " << particle->get_track_id() << std::endl;
        particle->identify();
        std::cout << " ptrue = (";
        std::cout.width(5);
        std::cout << particle->get_px();
        std::cout << ",";
        std::cout.width(5);
        std::cout << particle->get_py();
        std::cout << ",";
        std::cout.width(5);
        std::cout << particle->get_pz();
        std::cout << ")" << std::endl;

        std::cout << " vtrue = (";
        std::cout.width(5);
        std::cout << truthinfo->GetVtx(particle->get_vtx_id())->get_x();
        std::cout << ",";
        std::cout.width(5);
        std::cout << truthinfo->GetVtx(particle->get_vtx_id())->get_y();
        std::cout << ",";
        std::cout.width(5);
        std::cout << truthinfo->GetVtx(particle->get_vtx_id())->get_z();
        std::cout << ")" << std::endl;

        std::cout << " pt = " << sqrt(pow(particle->get_px(), 2) + pow(particle->get_py(), 2)) << std::endl;
        std::cout << " phi = " << atan2(particle->get_py(), particle->get_px()) << std::endl;
        std::cout << " eta = " << asinh(particle->get_pz() / sqrt(pow(particle->get_px(), 2) + pow(particle->get_py(), 2))) << std::endl;

        std::cout << " embed flag = " << truthinfo->isEmbeded(particle->get_track_id()) << std::endl;

        std::cout << " ---Associated-PHG4Hits-----------------------------------------" << std::endl;
        std::set<PHG4Hit*> local_g4hits = trutheval->all_truth_hits(particle);
        for (auto g4hit : local_g4hits)
        {
          float x = 0.5 * (g4hit->get_x(0) + g4hit->get_x(1));
          float y = 0.5 * (g4hit->get_y(0) + g4hit->get_y(1));
          float z = 0.5 * (g4hit->get_z(0) + g4hit->get_z(1));

          std::cout << " #" << g4hit->get_hit_id() << " xtrue = (";
          std::cout.width(5);
          std::cout << x;
          std::cout << ",";
          std::cout.width(5);
          std::cout << y;
          std::cout << ",";
          std::cout.width(5);
          std::cout << z;
          std::cout << ")";

          std::set<TrkrDefs::cluskey> clusters = clustereval->all_clusters_from(g4hit);
          for (unsigned long cluster_key : clusters)
          {
            TrkrCluster* cluster = clustermap->findCluster(cluster_key);
            // auto trkrid = TrkrDefs::getTrkrId(cluster_key);
            Acts::Vector3 glob;
            glob = tgeometry->getGlobalPosition(cluster_key, cluster);

            float local_x = glob(0);
            float local_y = glob(1);
            float local_z = glob(2);

            std::cout << " => #" << cluster_key;
            std::cout << " xreco = (";
            std::cout.width(5);
            std::cout << local_x;
            std::cout << ",";
            std::cout.width(5);
            std::cout << local_y;
            std::cout << ",";
            std::cout.width(5);
            std::cout << local_z;
            std::cout << ")";
          }

          std::cout << std::endl;
        }

        if (trackmap && clustermap)
        {
          std::set<SvtxTrack*> tracks = trackeval->all_tracks_from(particle);
          for (auto track : tracks)
          {
            float px = track->get_px();
            float py = track->get_py();
            float pz = track->get_pz();

            std::cout << "===Created-SvtxTrack==========================================" << std::endl;
            std::cout << " SvtxTrack id = " << track->get_id() << std::endl;
            std::cout << " preco = (";
            std::cout.width(5);
            std::cout << px;
            std::cout << ",";
            std::cout.width(5);
            std::cout << py;
            std::cout << ",";
            std::cout.width(5);
            std::cout << pz;
            std::cout << ")" << std::endl;
            std::cout << " quality = " << track->get_quality() << std::endl;
            std::cout << " nfromtruth = " << trackeval->get_nclusters_contribution(track, particle) << std::endl;

            std::cout << " ---Associated-SvtxClusters-to-PHG4Hits-------------------------" << std::endl;

            for (SvtxTrack::ConstClusterKeyIter local_iter = track->begin_cluster_keys();
                 local_iter != track->end_cluster_keys();
                 ++local_iter)
            {
              TrkrDefs::cluskey cluster_key = *local_iter;
              TrkrCluster* cluster = clustermap->findCluster(cluster_key);
              // auto trkrid = TrkrDefs::getTrkrId(cluster_key);
              Acts::Vector3 glob;
              glob = tgeometry->getGlobalPosition(cluster_key, cluster);

              float x = glob(0);
              float y = glob(1);
              float z = glob(2);

              std::cout << " #" << cluster_key << " xreco = (";
              std::cout.width(5);
              std::cout << x;
              std::cout << ",";
              std::cout.width(5);
              std::cout << y;
              std::cout << ",";
              std::cout.width(5);
              std::cout << z;
              std::cout << ") =>";

              PHG4Hit* g4hit = clustereval->max_truth_hit_by_energy(cluster_key);
              if ((g4hit) && (g4hit->get_trkid() == particle->get_track_id()))
              {
                x = 0.5 * (g4hit->get_x(0) + g4hit->get_x(1));
                y = 0.5 * (g4hit->get_y(0) + g4hit->get_y(1));
                z = 0.5 * (g4hit->get_z(0) + g4hit->get_z(1));

                std::cout << " #" << g4hit->get_hit_id()
                          << " xtrue = (";
                std::cout.width(5);
                std::cout << x;
                std::cout << ",";
                std::cout.width(5);
                std::cout << y;
                std::cout << ",";
                std::cout.width(5);
                std::cout << z;
                std::cout << ") => Gtrack id = " << g4hit->get_trkid();
              }
              else
              {
                std::cout << " noise hit";
              }
            }

            std::cout << std::endl;
          }
        }
      }
    }

    std::cout << std::endl;

  }  // if Verbosity()

  return;
}

void SvtxEvaluator::fillOutputNtuples(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "SvtxEvaluator::fillOutputNtuples() entered" << std::endl;
  }

  ActsGeometry* tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tgeometry)
  {
    std::cout << "No Acts geometry on node tree. Can't  continue."
              << std::endl;
    return;
  }

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
  float nclus_mms = 0;
  float nhit[100];
  for (float& i : nhit)
  {
    i = 0;
  }
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
      if (layer >= _nlayers_maps + _nlayers_intt && layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
      {
        // count all hits in this hitset
        TrkrHitSet::ConstRange hitrangei = hitsetiter->second->getHits();
        for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
             hitr != hitrangei.second;
             ++hitr)
        {
          nhit[layer]++;
          nhit_tpc_all++;
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
  }
  /**********/
  PHG4TpcCylinderGeomContainer* geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return;
  }
  PHG4TpcCylinderGeom* GeoLayer;
  int layer = _nlayers_maps + _nlayers_intt;
  GeoLayer = geom_container->GetLayerCellGeom(layer);
  int nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  float nhits = nhit[layer];
  occ11 = nhits / nbins;
  if (Verbosity() > 1)
  {
    std::cout << " occ11 = " << occ11
              << " nbins = " << nbins
              << " nhits = " << nhits
              << " layer = " << layer
              << std::endl;
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
    std::cout << " occ11 = " << occ11
              << " occ116 = " << occ116
              << " occ21 = " << occ21
              << " occ216 = " << occ216
              << " occ31 = " << occ31
              << " occ316 = " << occ316
              << std::endl;
  }
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
        unsigned int local_layer = TrkrDefs::getLayer(cluster_key);
        if (_nlayers_maps > 0)
        {
          if (local_layer < _nlayers_maps)
          {
            nclus_maps++;
          }
        }
        if (_nlayers_intt > 0)
        {
          if (local_layer >= _nlayers_maps && local_layer < _nlayers_maps + _nlayers_intt)
          {
            nclus_intt++;
          }
        }
        if (_nlayers_tpc > 0)
        {
          if (local_layer >= _nlayers_maps + _nlayers_intt && local_layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
          {
            nclus_tpc++;
          }
        }
        if (_nlayers_mms > 0)
        {
          if (local_layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc)
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
      std::cout << "Filling ntp_info " << std::endl;
    }
    float ntrk = 0;
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());
    if (trackmap)
    {
      ntrk = (float) trackmap->size();
    }
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    int nprim = truthinfo->GetNumPrimaryVertexParticles();
    if (Verbosity() > 0)
    {
      std::cout << "EVENTINFO SEED: " << m_fSeed << std::endl;
      std::cout << "EVENTINFO NHIT: " << std::setprecision(9) << nhit_tpc_all << std::endl;
      std::cout << "EVENTINFO NTRKGEN: " << nprim << std::endl;
      std::cout << "EVENTINFO NTRKREC: " << ntrk << std::endl;
    }
    float info_data[] = {(float) _ievent, m_fSeed,
                         occ11, occ116, occ21, occ216, occ31, occ316,
                         (float) nprim,
                         0,
                         ntrk,
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
      std::cout << "Filling ntp_vertex " << std::endl;
      std::cout << "start vertex time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
      _timer->restart();
    }

    GlobalVertexMap* gvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    SvtxVertexMap* vertexmap = nullptr;
    if (_use_initial_vertex)
    {
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    }
    else if (_use_genfit_vertex)
    {
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapRefit");
    }
    else
    {
      vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");  // Acts vertices
    }

    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

    if (gvertexmap && vertexmap && truthinfo)
    {
      const auto prange = truthinfo->GetPrimaryParticleRange();
      std::map<int, unsigned int> embedvtxid_particle_count;
      std::map<int, unsigned int> embedvtxid_maps_particle_count;
      std::map<int, unsigned int> vertex_particle_count;

      if (_do_vtx_eval_light == false)
      {
        for (auto iter = prange.first; iter != prange.second; ++iter)  // process all primary paricle
        {
          const int point_id = iter->second->get_vtx_id();
          int gembed = truthinfo->isEmbededVtx(iter->second->get_vtx_id());
          ++vertex_particle_count[point_id];
          ++embedvtxid_particle_count[gembed];
          PHG4Particle* g4particle = iter->second;

          if (_scan_for_embedded && gembed <= 0)
          {
            continue;
          }

          std::set<TrkrDefs::cluskey> g4clusters = clustereval->all_clusters_from(g4particle);
          unsigned int nglmaps = 0;

	  std::vector<int> lmaps(_nlayers_maps + 1,0);
          for (const TrkrDefs::cluskey g4cluster : g4clusters)
          {
            unsigned int local_layer = TrkrDefs::getLayer(g4cluster);
            // std::cout<<__LINE__<<": " << _ievent <<": " <<gtrackID << ": " << local_layer <<": " <<g4cluster->get_id() <<std::endl;
            if (_nlayers_maps > 0 && local_layer < _nlayers_maps)
            {
              lmaps[local_layer] = 1;
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
          {
            ++embedvtxid_maps_particle_count[gembed];
          }
        }
      }

      auto vrange = truthinfo->GetPrimaryVtxRange();
      std::map<int, bool> embedvtxid_found;
      std::map<int, int> embedvtxid_vertex_id;
      std::map<int, PHG4VtxPoint*> embedvtxid_vertex;
      for (auto iter = vrange.first; iter != vrange.second; ++iter)  // process all primary vertexes
      {
        const int point_id = iter->first;
        int gembed = truthinfo->isEmbededVtx(point_id);

        if (_scan_for_embedded && gembed <= 0)
        {
          continue;
        }

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
      for (auto& iter : embedvtxid_found)
      {
        if (iter.first >= 0)
        {
          continue;
        }
        ++ngembed;
      }

      for (auto& iter : *gvertexmap)
      {
        GlobalVertex* gvertex = iter.second;
        auto svtxviter = gvertex->find_vertexes(GlobalVertex::SVTX);
	if(svtxviter == gvertex->end_vertexes())
	  {
	    continue;
	  }

	GlobalVertex::VertexVector vertices = svtxviter->second;
	for(auto& vertex : vertices)
	  {
        PHG4VtxPoint* point = vertexeval->max_truth_point_by_ntracks(vertex);
        int vertexID = vertex->get_id();
        float vx = vertex->get_x();
        float vy = vertex->get_y();
        float vz = vertex->get_z();
        float ntracks = vertex->size_tracks();
        float chi2 = vertex->get_chisq();
        float ndof = vertex->get_ndof();
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
          {
            gntracksmaps = embedvtxid_maps_particle_count[(int) gembed];
          }
          gnembed = (float) ngembed;
          nfromtruth = vertexeval->get_ntracks_contribution(vertex, point);
          embedvtxid_found[(int) gembed] = true;
        }

        float vertex_data[] = {(float) _ievent, m_fSeed,
                               (float) vertexID,
                               vx,
                               vy,
                               vz,
                               ntracks,
                               chi2,
                               ndof,
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
                               nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

        _ntp_vertex->Fill(vertex_data);
	  }
      }

      if (!_scan_for_embedded)
      {
        for (std::map<int, bool>::iterator iter = embedvtxid_found.begin();
             iter != embedvtxid_found.end();
             ++iter)
        {
          if (embedvtxid_found[iter->first])
          {
            continue;
          }

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
            {
              gntracksmaps = embedvtxid_maps_particle_count[(int) gembed];
            }
            gnembed = (float) ngembed;
            //        nfromtruth = vertexeval->get_ntracks_contribution(vertex,point);
          }

          if (Verbosity() > 1)
          {
            std::cout << " adding vertex data " << std::endl;
          }

          float vertex_data[] = {(float) _ievent, m_fSeed,
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
                                 nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

          _ntp_vertex->Fill(vertex_data);
        }
      }
    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      std::cout << "vertex time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
    }
  }

  //-----------------------
  // fill the gpoint NTuple
  //-----------------------

  if (_ntp_gpoint)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_gpoint " << std::endl;
      _timer->restart();
    }
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

    if (truthinfo)
    {
      auto vrange = truthinfo->GetPrimaryVtxRange();
      const auto prange = truthinfo->GetPrimaryParticleRange();

      std::map<int, unsigned int> vertex_particle_count;
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
          const Vertex* vertex = vertexeval->best_vertex_from(point);

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

          float gpoint_data[] = {(float) _ievent, m_fSeed,
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
                                 nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

          _ntp_gpoint->Fill(gpoint_data);
        }
      }
    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      std::cout << "gpoint time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
    }
  }

  //---------------------
  // fill the G4hit NTuple
  //---------------------

  if (_ntp_g4hit)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_g4hit " << std::endl;
      _timer->restart();
    }
    std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits();
    for (auto g4hit : g4hits)
    {
      PHG4Particle* g4particle = trutheval->get_particle(g4hit);

      float g4hitID = g4hit->get_hit_id();
      float gx = g4hit->get_avg_x();
      float gy = g4hit->get_avg_y();
      float gz = g4hit->get_avg_z();
      TVector3 vg4(gx, gy, gz);
      float gt = g4hit->get_avg_t();
      float gpl = g4hit->get_path_length();
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
          if (trutheval->get_embed(g4particle) <= 0)
          {
            continue;
          }
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
        {
          outerhit = trutheval->get_outermost_truth_hit(g4particle);
        }

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
      float local_layer = NAN;
      float size = NAN;
      float efromtruth = NAN;
      float dphitru = NAN;
      float detatru = NAN;
      float dztru = NAN;
      float drtru = NAN;

      TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
      if (!clustermap)
      {
        clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
      }

      TrkrCluster* cluster = clustermap->findCluster(cluster_key);

      if (cluster)
      {
        clusID = cluster_key;
        // auto trkrid = TrkrDefs::getTrkrId(cluster_key);
        Acts::Vector3 global;
        global = tgeometry->getGlobalPosition(cluster_key, cluster);
        x = global(0);
        y = global(1);
        z = global(2);
        TVector3 vec2(x, y, z);
        eta = vec2.Eta();
        phi = vec2.Phi();
        e = cluster->getAdc();
        adc = cluster->getAdc();
        local_layer = (float) TrkrDefs::getLayer(cluster_key);
        size = 0.0;
        // count all hits for this cluster

        TrkrClusterHitAssoc* cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
        std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator>
            hitrange = cluster_hit_map->getHits(cluster_key);
        for (std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator
                 clushititer = hitrange.first;
             clushititer != hitrange.second; ++clushititer)
        {
          ++size;
        }

        dphitru = vec2.DeltaPhi(vg4);
        detatru = eta - geta;
        dztru = z - gz;
        drtru = vec2.DeltaR(vg4);
        if (g4particle)
        {
          efromtruth = clustereval->get_energy_contribution(cluster_key, g4particle);
        }
      }

      float g4hit_data[] = {(float) _ievent, m_fSeed,
                            g4hitID,
                            gx,
                            gy,
                            gz,
                            gt,
                            gpl,
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
                            local_layer,
                            size,
                            efromtruth,
                            dphitru,
                            detatru,
                            dztru,
                            drtru,
                            nhit_tpc_all,
                            nhit_tpc_in,
                            nhit_tpc_mid,
                            nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

      _ntp_g4hit->Fill(g4hit_data);
    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      std::cout << "g4hit time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
    }
  }

  //--------------------
  // fill the Hit NTuple
  //--------------------

  if (_ntp_hit)
  {
    if (Verbosity() >= 1)
    {
      std::cout << "Filling ntp_hit " << std::endl;
      _timer->restart();
    }
    // need things off of the DST...
    TrkrHitSetContainer* hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    PHG4TpcCylinderGeomContainer* local_geom_container =
        findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    if (!local_geom_container)
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
        TrkrHitSet* hitset = iter->second;

        // get all hits for this hitset
        TrkrHitSet::ConstRange hitrangei = hitset->getHits();
        for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
             hitr != hitrangei.second;
             ++hitr)
        {
          TrkrDefs::hitkey hit_key = hitr->first;
          TrkrHit* hit = hitr->second;
          PHG4Hit* g4hit = hiteval->max_truth_hit_by_energy(hit_key);
          PHG4Particle* g4particle = trutheval->get_particle(g4hit);
          float event = _ievent;
          float hitID = hit_key;
          float e = hit->getEnergy();
          float adc = hit->getAdc();
          float local_layer = TrkrDefs::getLayer(hitset_key);
          float sector = TpcDefs::getSectorId(hitset_key);
          float side = TpcDefs::getSide(hitset_key);
          float cellID = 0;
          float ecell = hit->getAdc();

          float phibin = NAN;
          float zbin = NAN;
          float phi = NAN;
          float z = NAN;

          if (local_layer >= _nlayers_maps + _nlayers_intt && local_layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
          {
            PHG4TpcCylinderGeom* local_GeoLayer = local_geom_container->GetLayerCellGeom(local_layer);
            phibin = (float) TpcDefs::getPad(hit_key);
            zbin = (float) TpcDefs::getTBin(hit_key);
            phi = local_GeoLayer->get_phicenter(phibin);
            z = local_GeoLayer->get_zcenter(zbin);
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
                if (trutheval->get_embed(g4particle) <= 0)
                {
                  continue;
                }
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
              {
                outerhit = trutheval->get_outermost_truth_hit(g4particle);
              }
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
              local_layer,
              sector,
              side,
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
              nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

          _ntp_hit->Fill(hit_data);
        }
      }
    }
    if (Verbosity() >= 1)
    {
      _timer->stop();
      std::cout << "hit time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
    }
  }

  //------------------------
  // fill the Cluster NTuple
  //------------------------

  if (Verbosity() >= 1)
  {
    std::cout << "check for ntp_cluster" << std::endl;
    _timer->restart();
  }

  if (_ntp_cluster && !_scan_for_embedded)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_cluster (all of them) " << std::endl;
    }
    // need things off of the DST...
    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
    if (!clustermap)
    {
      clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    }

    TrkrClusterHitAssoc* clusterhitmap = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
    TrkrHitSetContainer* hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    TrkrClusterIterationMapv1* _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    ClusterErrorPara ClusErrPara;

    if (Verbosity() > 1)
    {
      if (clustermap != nullptr)
      {
        std::cout << "got clustermap" << std::endl;
      }
      else
      {
        std::cout << "no clustermap" << std::endl;
      }
      if (clusterhitmap != nullptr)
      {
        std::cout << "got clusterhitmap" << std::endl;
      }
      else
      {
        std::cout << "no clusterhitmap" << std::endl;
      }
      if (hitsets != nullptr)
      {
        std::cout << "got hitsets" << std::endl;
      }
      else
      {
        std::cout << "no hitsets" << std::endl;
      }
    }

    if (clustermap && clusterhitmap && hitsets)
    {
      for (const auto& hitsetkey : clustermap->getHitSetKeys())
      {
        int hitsetlayer = TrkrDefs::getLayer(hitsetkey);
        auto range = clustermap->getClusters(hitsetkey);
        for (auto iter = range.first; iter != range.second; ++iter)
        {
          TrkrDefs::cluskey cluster_key = iter->first;
          TrkrCluster* cluster = clustermap->findCluster(cluster_key);
          SvtxTrack* track = trackeval->best_track_from(cluster_key);
          PHG4Particle* g4particle = clustereval->max_truth_particle_by_cluster_energy(cluster_key);
          float niter = 0;
          if (_iteration_map != nullptr)
          {
            niter = _iteration_map->getIteration(cluster_key);
          }
          float hitID = (float) cluster_key;
          // auto trkrid = TrkrDefs::getTrkrId(cluster_key);
          Acts::Vector3 cglob;
          cglob = tgeometry->getGlobalPosition(cluster_key, cluster);
          float x = cglob(0);
          float y = cglob(1);
          float z = cglob(2);
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
	  float redge = NAN;
	  float pedge = NAN;
	  float ovlp = NAN;


	  auto para_errors = ClusErrPara.get_clusterv5_modified_error(cluster, r, cluster_key);
	  phisize = cluster->getPhiSize();
	  zsize = cluster->getZSize();
	  // double clusRadius = r;
	  ez = sqrt(para_errors.second);
	  ephi = sqrt(para_errors.first);
	  maxadc = cluster->getMaxAdc();
	  pedge = cluster->getEdge();
	  ovlp = cluster->getOverlap();
	
	  if(hitsetlayer==7||hitsetlayer==22||hitsetlayer==23||hitsetlayer==38||hitsetlayer==39) redge = 1;

          float e = cluster->getAdc();
          float adc = cluster->getAdc();
          float local_layer = (float) TrkrDefs::getLayer(cluster_key);
          float sector = TpcDefs::getSectorId(cluster_key);
          float side = TpcDefs::getSide(cluster_key);

          // count all hits for this cluster
          TrkrDefs::hitsetkey local_hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);
          int hitsetlayer2 = TrkrDefs::getLayer(local_hitsetkey);
          if (hitsetlayer != local_layer)
          {
            std::cout << "WARNING hitset layer " << hitsetlayer << "| " << hitsetlayer2 << " layer " << local_layer << std::endl;
          }
          /*else{
            std::cout << "Good    hitset layer " << hitsetlayer << "| " << hitsetlayer2 << " layer " << local_layer << std::endl;
          }
          */
          float sumadc = 0;
          TrkrHitSetContainer::Iterator hitset = hitsets->findOrAddHitSet(local_hitsetkey);
          std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator>
              hitrange = clusterhitmap->getHits(cluster_key);
          for (std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator
                   clushititer = hitrange.first;
               clushititer != hitrange.second; ++clushititer)
          {
            TrkrHit* hit = hitset->second->getHit(clushititer->second);
            if (!hit)
            {
              continue;
            }

            ++size;
            sumadc += (hit->getAdc() - 70);
            if ((hit->getAdc() - 70) > maxadc)
            {
              maxadc = (hit->getAdc() - 70);
            }
          }
          e = sumadc;

          float trackID = NAN;
          if (track != nullptr)
          {
            trackID = track->get_id();
          }

          float g4hitID = NAN;
          float gx = NAN;
          float gy = NAN;
          float gz = NAN;
          float gr = NAN;
          float gphi = NAN;
          // float gedep = NAN;
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

          if (Verbosity() > 1)
          {
            std::cout << PHWHERE << "  ****   reco: layer " << local_layer << std::endl;
            std::cout << "              reco cluster key " << cluster_key << "  r " << r << "  x " << x << "  y " << y << "  z " << z << "  phi " << phi << " adc " << adc << std::endl;
          }
          float nparticles = NAN;

          // get best matching truth cluster from clustereval
          const auto [truth_ckey, truth_cluster] = clustereval->max_truth_cluster_by_energy(cluster_key);
          if (truth_cluster)
          {
            if (Verbosity() > 1)
            {
              std::cout << "Found matching truth cluster with key " << truth_ckey << " for reco cluster key " << cluster_key << " in layer " << local_layer << std::endl;
            }

            g4hitID = 0;
            gx = truth_cluster->getX();
            gy = truth_cluster->getY();
            gz = truth_cluster->getZ();
            efromtruth = truth_cluster->getError(0, 0);

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
              {
                outerhit = trutheval->get_outermost_truth_hit(g4particle);
              }
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

            if (Verbosity() > 1)
            {
              std::cout << "             truth cluster key " << truth_ckey << " gr " << gr << " gx " << gx << " gy " << gy << " gz " << gz << " gphi " << gphi << " efromtruth " << efromtruth << std::endl;
            }
          }  //  if (truth_cluster) {

          if (g4particle)
          {
          }
          nparticles = clustereval->all_truth_particles(cluster_key).size();

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
                                  pez,
                                  pephi,
                                  e,
                                  adc,
                                  maxadc,
                                  local_layer,
                                  sector,
                                  side,
                                  size,
                                  phisize,
                                  zsize,
				  pedge,
				  redge,
				  ovlp,
                                  trackID,
                                  niter,
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
                                  nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

          _ntp_cluster->Fill(cluster_data);
        }
      }
    }
  }
  else if (_ntp_cluster && _scan_for_embedded)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_cluster (embedded only) " << std::endl;
    }

    // if only scanning embedded signals, loop over all the tracks from
    // embedded particles and report all of their clusters, including those
    // from other sources (noise hits on the embedded track)

    // need things off of the DST...
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());

    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
    if (!clustermap)
    {
      clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    }
    ClusterErrorPara ClusErrPara;

    TrkrClusterHitAssoc* clusterhitmap = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
    TrkrHitSetContainer* hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    TrkrClusterIterationMapv1* _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");

    if (trackmap != nullptr && clustermap != nullptr && clusterhitmap != nullptr && hitsets != nullptr)
    {
      for (auto& iter : *trackmap)
      {
        SvtxTrack* track = iter.second;

        PHG4Particle* truth = trackeval->max_truth_particle_by_nclusters(track);
        if (truth)
        {
          if (trutheval->get_embed(truth) <= 0)
          {
            continue;
          }
        }

        std::vector<TrkrDefs::cluskey> clusters;
        auto siseed = track->get_silicon_seed();
        if (siseed)
        {
          for (auto local_iter = siseed->begin_cluster_keys();
               local_iter != siseed->end_cluster_keys();
               ++local_iter)
          {
            TrkrDefs::cluskey cluster_key = *local_iter;
            clusters.push_back(cluster_key);
          }
        }
        auto tpcseed = track->get_tpc_seed();
        if (tpcseed)
        {
          for (auto local_iter = tpcseed->begin_cluster_keys();
               local_iter != tpcseed->end_cluster_keys();
               ++local_iter)
          {
            TrkrDefs::cluskey cluster_key = *local_iter;
            clusters.push_back(cluster_key);
          }
        }

        // loop over all cluster keys and build ntuple
        for (unsigned long cluster_key : clusters)
        {
          TrkrCluster* cluster = clustermap->findCluster(cluster_key);
          if (!cluster)
          {
            continue;  // possible to be missing from corrected clusters if cluster mover fails
          }

          PHG4Hit* g4hit = clustereval->max_truth_hit_by_energy(cluster_key);
          PHG4Particle* g4particle = trutheval->get_particle(g4hit);

          float niter = 0;
          if (_iteration_map != nullptr)
          {
            niter = _iteration_map->getIteration(cluster_key);
          }
          float hitID = (float) TrkrDefs::getClusIndex(cluster_key);
          Acts::Vector3 glob = tgeometry->getGlobalPosition(cluster_key, cluster);
          float x = glob(0);
          float y = glob(1);
          float z = glob(2);
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
	  float redge = NAN;
	  float pedge = NAN;
	  float ovlp = NAN;

	  auto para_errors = ClusErrPara.get_clusterv5_modified_error(cluster, r, cluster_key);
	  phisize = cluster->getPhiSize();
	  zsize = cluster->getZSize();
	  // double clusRadius = r;
	  ez = sqrt(para_errors.second);
	  ephi = sqrt(para_errors.first);
	  maxadc = cluster->getMaxAdc();
	  pedge = cluster->getEdge();
	  ovlp = cluster->getOverlap();
          
	  float e = cluster->getAdc();
          float adc = cluster->getAdc();
          float local_layer = (float) TrkrDefs::getLayer(cluster_key);
          float sector = TpcDefs::getSectorId(cluster_key);
          float side = TpcDefs::getSide(cluster_key);
          // count all hits for this cluster
	  if(local_layer==7||local_layer==22||local_layer==23||local_layer==38||local_layer==39) redge = 1;

          // count all hits for this cluster
          TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);
          TrkrHitSetContainer::Iterator hitset = hitsets->findOrAddHitSet(hitsetkey);
          std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator>
              hitrange = clusterhitmap->getHits(cluster_key);
          for (std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator
                   clushititer = hitrange.first;
               clushititer != hitrange.second; ++clushititer)
          {
            TrkrHit* hit = hitset->second->getHit(clushititer->second);
            ++size;
            if (hit->getAdc() > maxadc)
            {
              maxadc = hit->getAdc();
            }
          }

          float trackID = NAN;
          trackID = track->get_id();

          float g4hitID = NAN;
          float gx = NAN;
          float gy = NAN;
          float gz = NAN;
          float gr = NAN;
          float gphi = NAN;
          // float gedep = NAN;
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

          // get best matching truth cluster from clustereval
          const auto [truth_ckey, truth_cluster] = clustereval->max_truth_cluster_by_energy(cluster_key);
          if (truth_cluster)
          {
            if (Verbosity() > 1)
            {
              std::cout << "         Found matching truth cluster with key " << truth_ckey << " for reco cluster key " << cluster_key << " in layer " << local_layer << std::endl;
            }

            g4hitID = 0;
            gx = truth_cluster->getX();
            gy = truth_cluster->getY();
            gz = truth_cluster->getZ();
            efromtruth = truth_cluster->getError(0, 0);

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
              {
                outerhit = trutheval->get_outermost_truth_hit(g4particle);
              }
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
                                  pez,
                                  pephi,
                                  e,
                                  adc,
                                  maxadc,
                                  local_layer,
                                  sector,
                                  side,
                                  size,
                                  phisize,
                                  zsize,
				  pedge,
				  redge,
				  ovlp,
                                  trackID,
                                  niter,
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
                                  nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

          _ntp_cluster->Fill(cluster_data);
        }
      }
    }
  }
  if (Verbosity() >= 1)
  {
    _timer->stop();
    std::cout << "cluster time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
  }

  // fill the truth cluster NTuple
  //-----------------------------------

  if (Verbosity() > 1)
  {
    std::cout << "check for ntp_g4cluster" << std::endl;
    _timer->restart();
  }

  if (_ntp_g4cluster)
  {
    if (Verbosity() >= 1)
    {
      std::cout << "Filling ntp_g4cluster " << std::endl;
    }
    ClusterErrorPara ClusErrPara;
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    PHG4TruthInfoContainer::ConstRange range = truthinfo->GetParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
         iter != range.second;
         ++iter)
    {
      PHG4Particle* g4particle = iter->second;

      if (_scan_for_embedded)
      {
        if (trutheval->get_embed(g4particle) <= 0)
        {
          continue;
        }
      }

      float gtrackID = g4particle->get_track_id();
      float gflavor = g4particle->get_pid();
      float gembed = trutheval->get_embed(g4particle);
      float gprimary = trutheval->is_primary(g4particle);

      if (Verbosity() > 1)
      {
        std::cout << PHWHERE << " PHG4Particle ID " << gtrackID << " gflavor " << gflavor << " gprimary " << gprimary << std::endl;
      }

      // Get the truth clusters from this particle
      std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster> > truth_clusters = trutheval->all_truth_clusters(g4particle);

      // loop over layers and add to ntuple
      for (const auto& [ckey, gclus] : truth_clusters)
      {
        unsigned int local_layer = TrkrDefs::getLayer(ckey);

        float gx = gclus->getX();
        float gy = gclus->getY();
        float gz = gclus->getZ();
        float gt = NAN;
        float gedep = gclus->getError(0, 0);
        float gadc = (float) gclus->getAdc();

        TVector3 gpos(gx, gy, gz);
        float gr = sqrt(gx * gx + gy * gy);
        float gphi = gpos.Phi();
        float geta = gpos.Eta();

        if (Verbosity() > 1)
        {
          std::cout << PHWHERE << "  ****   truth: layer " << local_layer << std::endl;
          std::cout << "             truth cluster key " << ckey << " gr " << gr << " gx " << gx << " gy " << gy << " gz " << gz
                    << " gphi " << gphi << " gedep " << gedep << " gadc " << gadc << std::endl;
        }

        float gphisize = gclus->getSize(1, 1);
        float gzsize = gclus->getSize(2, 2);

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
        float adc = NAN;

        float nreco = 0;
        float phisize = 0;
        float zsize = 0;
        const auto [reco_ckey, reco_cluster] = clustereval->reco_cluster_from_truth_cluster(ckey, gclus);
        if (reco_cluster)
        {
          nreco = 1;
          // auto trkrid = TrkrDefs::getTrkrId(reco_ckey);
          Acts::Vector3 glob;
          glob = tgeometry->getGlobalPosition(reco_ckey, reco_cluster);
          x = glob(0);
          y = glob(1);
          z = glob(2);

          TVector3 pos(x, y, z);
          r = sqrt(x * x + y * y);
          phi = pos.Phi();
          eta = pos.Eta();
                 
	  auto para_errors = ClusErrPara.get_clusterv5_modified_error(reco_cluster, r, ckey);
	  // std::cout << " ver v4 " <<  std::endl;
	  phisize = reco_cluster->getPhiSize();
	  zsize = reco_cluster->getZSize();
	  ez = sqrt(para_errors.second);
	  ephi = sqrt(para_errors.first);
          

          adc = reco_cluster->getAdc();

          if (Verbosity() > 1)
          {
            std::cout << "              reco cluster key " << reco_ckey << "  r " << r << "  x " << x << "  y " << y << "  z " << z << "  phi " << phi << " adc " << adc << std::endl;
          }
        }
        if (nreco == 0 && Verbosity() > 1)
        {
          if (Verbosity() > 1)
          {
            std::cout << "   ----------- Failed to find matching reco cluster " << std::endl;
          }
        }

        // add this cluster to the ntuple

        float g4cluster_data[] = {(float) _ievent,
                                  (float) local_layer,
                                  gx,
                                  gy,
                                  gz,
                                  gt,
                                  gedep,
                                  gr,
                                  gphi,
                                  geta,
                                  gtrackID,
                                  gflavor,
                                  gembed,
                                  gprimary,
                                  gphisize,
                                  gzsize,
                                  gadc,
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
                                  adc,
                                  phisize,
                                  zsize};
        _ntp_g4cluster->Fill(g4cluster_data);
      }
    }
  }

  //------------------------
  // fill the Gtrack NTuple
  //------------------------

  // need things off of the DST...

  // std::cout << "check for ntp_gtrack" << std::endl;

  if (_ntp_gtrack)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_gtrack " << std::endl;
      _timer->restart();
    }
    ClusterErrorPara ClusErrPara;
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    GlobalVertexMap* gvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
    if(!clustermap)
      clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    
    if (truthinfo)
    {
      PHG4TruthInfoContainer::ConstRange range = truthinfo->GetParticleRange();
      if (_scan_for_primaries)
      {
        range = truthinfo->GetPrimaryParticleRange();
      }

      Float_t gntracks = (Float_t) truthinfo->GetNumPrimaryVertexParticles();
      Float_t gnchghad = 0;
      for(auto iter = range.first; iter!= range.second; ++iter)
	{
	   PHG4Particle* g4particle = iter->second;

	   if (_scan_for_embedded)
	     {
	       if (trutheval->get_embed(g4particle) <= 0)
		 {
		   continue;
		 }
	     }

	   float gflavor = g4particle->get_pid();
	  if(fabs(gflavor)==211 || fabs(gflavor)==321 || fabs(gflavor)==2212)
	    {
	      gnchghad++;
	    }
	}

      for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
           iter != range.second;
           ++iter)
      {
        PHG4Particle* g4particle = iter->second;

        if (_scan_for_embedded)
        {
          if (trutheval->get_embed(g4particle) <= 0)
          {
            continue;
          }
        }

        float gtrackID = g4particle->get_track_id();
        float gflavor = g4particle->get_pid();

        std::set<TrkrDefs::cluskey> g4clusters = clustereval->all_clusters_from(g4particle);

        float ng4hits = g4clusters.size();
        unsigned int ngmaps = 0;
        unsigned int ngmms = 0;
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
        unsigned int nglmms = 0;

	std::vector<int> lmaps(_nlayers_maps + 1,0);
	std::vector<int> lintt(_nlayers_intt + 1,0);
	std::vector<int> ltpc(_nlayers_tpc + 1,0);
	std::vector<int> lmms(_nlayers_mms + 1,0);

        for (const TrkrDefs::cluskey g4cluster : g4clusters)
        {
          unsigned int local_layer = TrkrDefs::getLayer(g4cluster);
          // std::cout<<__LINE__<<": " << _ievent <<": " <<gtrackID << ": " << local_layer <<": " <<g4cluster->get_id() <<std::endl;
          if (_nlayers_maps > 0 && local_layer < _nlayers_maps)
          {
            lmaps[local_layer] = 1;
            ngmaps++;
          }
          if (_nlayers_mms > 0 && local_layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc)
          {
            lmms[local_layer - _nlayers_tpc - _nlayers_intt - _nlayers_maps] = 1;
            ngmms++;
          }
          if (_nlayers_intt > 0 && local_layer >= _nlayers_maps && local_layer < _nlayers_maps + _nlayers_intt)
          {
            lintt[local_layer - _nlayers_maps] = 1;
            ngintt++;
          }

          if (_nlayers_intt > 0 && local_layer == _nlayers_maps && local_layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt1++;
          }

          if (_nlayers_intt > 1 && local_layer == _nlayers_maps + 1 && local_layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt2++;
          }

          if (_nlayers_intt > 2 && local_layer == _nlayers_maps + 2 && local_layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt3++;
          }

          if (_nlayers_intt > 3 && local_layer == _nlayers_maps + 3 && local_layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt4++;
          }

          if (_nlayers_intt > 4 && local_layer == _nlayers_maps + 4 && local_layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt5++;
          }

          if (_nlayers_intt > 5 && local_layer == _nlayers_maps + 5 && local_layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt6++;
          }

          if (_nlayers_intt > 6 && local_layer == _nlayers_maps + 6 && local_layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt7++;
          }

          if (_nlayers_intt > 7 && local_layer == _nlayers_maps + 7 && local_layer < _nlayers_maps + _nlayers_intt)
          {
            ngintt8++;
          }
          if (_nlayers_tpc > 0 && local_layer >= _nlayers_maps + _nlayers_intt && local_layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
          {
            ltpc[local_layer - (_nlayers_maps + _nlayers_intt)] = 1;
            ngtpc++;
          }
        }
        if (_nlayers_maps > 0)
        {
          for (unsigned int i = 0; i < _nlayers_maps; i++)
          {
            nglmaps += lmaps[i];
          }
        }
        if (_nlayers_intt > 0)
        {
          for (unsigned int i = 0; i < _nlayers_intt; i++)
          {
            nglintt += lintt[i];
          }
        }
        if (_nlayers_tpc > 0)
        {
          for (unsigned int i = 0; i < _nlayers_tpc; i++)
          {
            ngltpc += ltpc[i];
          }
        }
        if (_nlayers_mms > 0)
        {
          for (unsigned int i = 0; i < _nlayers_mms; i++)
          {
            nglmms += lmms[i];
          }
        }

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
        {
          outerhit = trutheval->get_outermost_truth_hit(g4particle);
        }

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
        float local_nhits = 0;
        float nmaps = 0;
        float nintt = 0;
        float ntpc = 0;
        float nmms = 0;
        float ntpc1 = 0;
        float ntpc11 = 0;
        float ntpc2 = 0;
        float ntpc3 = 0;
        float nlintt = 0;
        float nlmaps = 0;
        float nltpc = 0;
        float nlmms = 0;
        unsigned int layers = 0x0;
        int vertexID = -1;
        float vx = NAN;
        float vy = NAN;
        float vz = NAN;
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
	float nwrongmaps = NAN;
        float ntruintt = NAN;
	float nwrongintt = NAN;
        float ntrumms = NAN;
	float nwrongmms = NAN;
        float ntrutpc = NAN;
	float nwrongtpc = NAN;
        float ntrutpc1 = NAN;
	float nwrongtpc1 = NAN;
        float ntrutpc11 = NAN;
	float nwrongtpc11 =NAN;
        float ntrutpc2 = NAN;
	float nwrongtpc2 = NAN;
        float ntrutpc3 = NAN;
	float nwrongtpc3 = NAN;
        float layersfromtruth = NAN;
	float npedge = 0;
	float nredge = 0;
	float nbig = 0;
	float novlp = 0;
	float merr = 0;
	float msize = 0;
	float siqr = NAN;
	float siphi = NAN;
	float sithe = NAN;
	float six0 = NAN;
	float siy0 = NAN;
	float tpqr = NAN;
	float tpphi = NAN;
	float tpthe = NAN;
	float tpx0 = NAN;
	float tpy0 = NAN;

        if (_do_track_match)
        {
          SvtxTrack* track = trackeval->best_track_from(g4particle);

          if (track)
          {
            trackID = track->get_id();
            charge = track->get_charge();
            quality = track->get_quality();
            chisq = track->get_chisq();
            ndf = track->get_ndf();
            TrackSeed* silseed = track->get_silicon_seed();
            TrackSeed* tpcseed = track->get_tpc_seed();
            if (tpcseed)
            {
              local_nhits += tpcseed->size_cluster_keys();
            }
            if (silseed)
            {
              local_nhits += silseed->size_cluster_keys();
            }
// +1 just in case _nlayers is zero
            std::vector<int> maps(_nlayers_maps+1, 0);
            std::vector<int> intt(_nlayers_intt+1, 0);
            std::vector<int> tpc(_nlayers_tpc+1, 0);
            std::vector<int> mms(_nlayers_mms+1, 0);


            if (tpcseed)
            {
	      tpqr = tpcseed->get_qOverR();
	      tpphi = tpcseed->get_phi();
	      tpthe = tpcseed->get_theta();
	      tpx0 = tpcseed->get_X0();
	      tpy0 = tpcseed->get_Y0();
              for (TrackSeed::ConstClusterKeyIter local_iter = tpcseed->begin_cluster_keys();
                   local_iter != tpcseed->end_cluster_keys();
                   ++local_iter)
              {
                TrkrDefs::cluskey cluster_key = *local_iter;
                TrkrCluster* cluster = clustermap->findCluster(cluster_key);
                unsigned int local_layer = TrkrDefs::getLayer(cluster_key);
                if (_nlayers_maps > 0 && local_layer < _nlayers_maps)
                {
                  maps[local_layer] = 1;
                  nmaps++;
                }
                if (_nlayers_intt > 0 && local_layer >= _nlayers_maps && local_layer < _nlayers_maps + _nlayers_intt)
                {
                  intt[local_layer - _nlayers_maps] = 1;
                  nintt++;
                }
                if (_nlayers_tpc > 0 &&
                    local_layer >= (_nlayers_maps + _nlayers_intt) &&
                    local_layer < (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
                {
                  tpc[local_layer - (_nlayers_maps + _nlayers_intt)] = 1;
                  ntpc++;

                  if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 8)
                  {
                    ntpc11++;
                  }

                  if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 16)
                  {
                    // std::cout << " tpc1: layer " << local_layer << std::endl;
                    ntpc1++;
                  }
                  else if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 32)
                  {
                    // std::cout << " tpc2: layer " << local_layer << std::endl;
                    ntpc2++;
                  }
                  else if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 48)
                  {
                    // std::cout << " tpc3: layer " << local_layer << std::endl;
                    ntpc3++;
                  }
                }

                if (_nlayers_mms > 0 && local_layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc)
                {
                  mms[local_layer - (_nlayers_maps + _nlayers_intt + _nlayers_tpc)] = 1;
                  nmms++;
                }
              
		Acts::Vector3 glob;
		glob = tgeometry->getGlobalPosition(cluster_key, cluster);
		float x = glob(0);
		float y = glob(1);
		float z = glob(2);
		
		TVector3 pos(x, y, z);
		float r = sqrt(x*x+y*y);
		phi = pos.Phi();
		eta = pos.Eta();
		float gphisize = 0;
		//float zsize = 0;
		float gphierr = 0;
		//float zerr = 0;
		float govlp = 0 ;
		float gedge = 0;
		if(cluster!=nullptr){
        
		  gphisize = cluster->getPhiSize();
		  //zsize = clusterv5->getZSize();
		  auto para_errors = ClusErrPara.get_clusterv5_modified_error(cluster,r,cluster_key);
		  //zerr = sqrt(para_errors.second);
		  gphierr = sqrt(para_errors.first);
		  govlp = cluster->getOverlap();
		  gedge = cluster->getEdge();
		 
		  if(gedge>0) npedge++;
		  if(gphisize>=4) nbig++;
		  if(govlp>=2) novlp++;
		  merr+=gphierr;
		  msize+=gphisize;
		}
		if(local_layer==7||local_layer==22||local_layer==23||local_layer==38||local_layer==39) nredge++;
	      }
	    }

            if (silseed)
            {
	      siqr = silseed->get_qOverR();
	      siphi = silseed->get_phi();
	      sithe = silseed->get_theta();
	      six0 = silseed->get_X0();
	      siy0 = silseed->get_Y0();
              for (TrackSeed::ConstClusterKeyIter local_iter = silseed->begin_cluster_keys();
                   local_iter != silseed->end_cluster_keys();
                   ++local_iter)
              {
                TrkrDefs::cluskey cluster_key = *local_iter;
                // TrkrCluster* cluster = clustermap->findCluster(cluster_key);
                unsigned int local_layer = TrkrDefs::getLayer(cluster_key);
                if (_nlayers_maps > 0 && local_layer < _nlayers_maps)
                {
                  maps[local_layer] = 1;
                  nmaps++;
                }
                if (_nlayers_intt > 0 && local_layer >= _nlayers_maps && local_layer < _nlayers_maps + _nlayers_intt)
                {
                  intt[local_layer - _nlayers_maps] = 1;
                  nintt++;
                }
                if (_nlayers_tpc > 0 &&
                    local_layer >= (_nlayers_maps + _nlayers_intt) &&
                    local_layer < (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
                {
                  tpc[local_layer - (_nlayers_maps + _nlayers_intt)] = 1;
                  ntpc++;

                  if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 8)
                  {
                    ntpc11++;
                  }

                  if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 16)
                  {
                    // std::cout << " tpc1: layer " << local_layer << std::endl;
                    ntpc1++;
                  }
                  else if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 32)
                  {
                    // std::cout << " tpc2: layer " << local_layer << std::endl;
                    ntpc2++;
                  }
                  else if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 48)
                  {
                    // std::cout << " tpc3: layer " << local_layer << std::endl;
                    ntpc3++;
                  }
                }

                if (_nlayers_mms > 0 && local_layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc)
                {
                  mms[local_layer - (_nlayers_maps + _nlayers_intt + _nlayers_tpc)] = 1;
                  nmms++;
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
            /* std::cout << " layers " << layers
                 << " nmaps " << nmaps
                 << " nintt " << nintt
                 << " ntpc  " << ntpc
                 << " nlmaps "<< nlmaps
                 << " nlintt " << nlintt
                 << " nltpc  " << nltpc
                 << std::endl;
            */

            // this is the global vertex id
            vertexID = track->get_vertex_id();

            GlobalVertex* vertex = gvertexmap->get(vertexID);
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

            px = track->get_px();
            py = track->get_py();
            pz = track->get_pz();
            TVector3 v(px, py, pz);
            pt = v.Pt();
            eta = v.Eta();
            phi = v.Phi();
            float CVxx = track->get_error(3, 3);
            float CVxy = track->get_error(3, 4);
            float CVxz = track->get_error(3, 5);
            float CVyy = track->get_error(4, 4);
            float CVyz = track->get_error(4, 5);
            float CVzz = track->get_error(5, 5);
            deltapt = sqrt((CVxx * px * px + 2 * CVxy * px * py + CVyy * py * py) / (px * px + py * py));
            deltaeta = sqrt((CVzz * (px * px + py * py) * (px * px + py * py) + pz * (-2 * (CVxz * px + CVyz * py) * (px * px + py * py) + CVxx * px * px * pz + CVyy * py * py * pz + 2 * CVxy * px * py * pz)) / ((px * px + py * py) * (px * px + py * py) * (px * px + py * py + pz * pz)));
            deltaphi = sqrt((CVyy * px * px - 2 * CVxy * px * py + CVxx * py * py) / ((px * px + py * py) * (px * px + py * py)));
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
	      auto pair = trackeval->get_layer_range_contribution(track, g4particle, 0, _nlayers_maps);
	      ntrumaps = pair.first;
	      nwrongmaps = pair.second;
            }
            if (_nlayers_intt == 0)
            {
              ntruintt = 0;
            }
            else
            {
              auto pair = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps, _nlayers_maps + _nlayers_intt);
	      ntruintt = pair.first;
	      nwrongintt = pair.second;
            }
            if (_nlayers_mms == 0)
            {
              ntrumms = 0;
            }
            else
            {
	      auto pair = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt + _nlayers_tpc, _nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms);
	      ntrumms = pair.first;
	      nwrongmms = pair.second;
            }
            auto pair = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt, _nlayers_maps + _nlayers_intt + _nlayers_tpc);
	    ntrutpc = pair.first;
	    nwrongtpc = pair.second;
            pair  = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt, _nlayers_maps + _nlayers_intt + 16);
	    ntrutpc1 = pair.first;
	    nwrongtpc1 = pair.second;
            pair = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt, _nlayers_maps + _nlayers_intt + 8);
	    ntrutpc11 = pair.first;
	    nwrongtpc11 = pair.second;
            pair = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt + 16, _nlayers_maps + _nlayers_intt + 32);
	    ntrutpc2 = pair.first;
	    nwrongtpc2 = pair.second;
            pair = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt + 32, _nlayers_maps + _nlayers_intt + _nlayers_tpc);
	    ntrutpc3 = pair.first;
	    nwrongtpc3 = pair.first;
            layersfromtruth = trackeval->get_nclusters_contribution_by_layer(track, g4particle);
          }
        }

        float gtrack_data[] = {(float) _ievent, m_fSeed,
                               gntracks,
			       gnchghad,
                               gtrackID,
                               gflavor,
                               ng4hits,
                               (float) ngmaps,
                               (float) ngintt,
                               (float) ngmms,
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
                               (float) nglmms,
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
			       siqr,
			       siphi,
			       sithe,
			       six0,
			       siy0,
			       tpqr,
			       tpphi,
			       tpthe,
			       tpx0,
			       tpy0,
                               charge,
                               quality,
                               chisq,
                               ndf,
                               local_nhits,
                               (float) layers,
                               nmaps,
                               nintt,
                               ntpc,
                               nmms,
                               ntpc1,
                               ntpc11,
                               ntpc2,
                               ntpc3,
                               nlmaps,
                               nlintt,
                               nltpc,
                               nlmms,
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
                               nfromtruth,
                               nwrong,
                               ntrumaps,
			       nwrongmaps,
                               ntruintt,
			       nwrongintt,
                               ntrutpc,
			       nwrongtpc,
                               ntrumms,
			       nwrongmms,
                               ntrutpc1,
			       nwrongtpc1,
                               ntrutpc11,
			       nwrongtpc11,
                               ntrutpc2,
			       nwrongtpc2,
                               ntrutpc3,
			       nwrongtpc3,
                               layersfromtruth,
			       npedge,
			       nredge,
			       nbig,
			       novlp,
			       merr,
			       msize,
                               nhit_tpc_all,
                               nhit_tpc_in,
                               nhit_tpc_mid,
                               nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

        /*
        std::cout << " ievent " << _ievent
             << " gtrackID " << gtrackID
             << " gflavor " << gflavor
             << " ng4hits " << ng4hits
             << std::endl;
        */

        _ntp_gtrack->Fill(gtrack_data);
      }
    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      std::cout << "gtrack time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
    }
  }

  //------------------------
  // fill the Track NTuple
  //------------------------

  if (_ntp_track)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_track " << std::endl;
      _timer->restart();
    }

    // need things off of the DST...
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());

    GlobalVertexMap* gvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
    if(!clustermap)
      clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    ClusterErrorPara ClusErrPara;
    if (trackmap)
    {
      for (auto& iter : *trackmap)
      {
        SvtxTrack* track = iter.second;
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
        float local_nhits = 0;
        if (tpcseed)
        {
          local_nhits += tpcseed->size_cluster_keys();
        }
        if (silseed)
        {
          local_nhits += silseed->size_cluster_keys();
        }
        unsigned int layers = 0x0;
// +1 just in case _nlayers is zero
	std::vector<int> maps(_nlayers_maps+1,0);
        std::vector<int> intt(_nlayers_intt+1,0);
        std::vector<int> tpc(_nlayers_tpc+1,0);
        std::vector<int> mms(_nlayers_mms+1,0);

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
	float npedge = 0;
	float nredge = 0;
	float nbig = 0;
	float novlp = 0;
	float merr = 0;
	float msize = 0;
	float siqr = NAN;
	float siphi = NAN;
	float sithe = NAN;
	float six0 = NAN;
	float siy0 = NAN;
	float tpqr = NAN;
	float tpphi = NAN;
	float tpthe = NAN;
	float tpx0 = NAN;
	float tpy0 = NAN;

        if (tpcseed)
        {
	  tpqr = tpcseed->get_qOverR();
	  tpphi = tpcseed->get_phi();
	  tpthe = tpcseed->get_theta();
	  tpx0 = tpcseed->get_X0();
	  tpy0 = tpcseed->get_Y0();
          for (SvtxTrack::ConstClusterKeyIter local_iter = tpcseed->begin_cluster_keys();
               local_iter != tpcseed->end_cluster_keys();
               ++local_iter){
            TrkrDefs::cluskey cluster_key = *local_iter;
            TrkrCluster* cluster = clustermap->findCluster(cluster_key);
            unsigned int local_layer = TrkrDefs::getLayer(cluster_key);
	    
            if (_nlayers_maps > 0 && local_layer < _nlayers_maps)
            {
              maps[local_layer] = 1;
              nmaps++;
            }
            if (_nlayers_intt > 0 && local_layer >= _nlayers_maps && local_layer < _nlayers_maps + _nlayers_intt)
            {
              intt[local_layer - _nlayers_maps] = 1;
              nintt++;
            }
            if (_nlayers_mms > 0 && local_layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc && local_layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms)
            {
              mms[local_layer - (_nlayers_maps + _nlayers_intt + _nlayers_tpc)] = 1;
              nmms++;
            }
            if (_nlayers_tpc > 0 && local_layer >= (_nlayers_maps + _nlayers_intt) && local_layer < (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
            {
              tpc[local_layer - (_nlayers_maps + _nlayers_intt)] = 1;
              ntpc++;
              if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 8)
              {
                ntpc11++;
              }

              if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 16)
              {
                ntpc1++;
              }
              else if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 32)
              {
                ntpc2++;
              }
              else if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 48)
              {
                ntpc3++;
              }
            }
          
	    Acts::Vector3 glob;
	    glob = tgeometry->getGlobalPosition(cluster_key, cluster);
	    float x = glob(0);
	    float y = glob(1);
	    float z = glob(2);
	    
	    TVector3 pos(x, y, z);
	    float r = sqrt(x*x+y*y);
	    //float phi = pos.Phi();
	    //float eta = pos.Eta();
	    float rphisize = 0;
	    //float zsize = 0;
	    float rphierr = 0;
	    //float zerr = 0;
	    float rovlp = 0 ;
	    float pedge = 0;
	    if(cluster!=nullptr){
	     
	      rphisize = cluster->getPhiSize();
	      //zsize = clusterv5->getZSize();
	      auto para_errors = ClusErrPara.get_clusterv5_modified_error(cluster,r,cluster_key);
	      //zerr = sqrt(para_errors.second);
	      rphierr = sqrt(para_errors.first);
	      rovlp = cluster->getOverlap();
	      pedge = cluster->getEdge();
	       
	      if(pedge>0) npedge++;
	      if(rphisize>=4) nbig++;
	      if(rovlp>=2) novlp++;
	      merr+=rphierr;
	      msize+=rphisize;
	    }
	    if(local_layer==7||local_layer==22||local_layer==23||local_layer==38||local_layer==39) nredge++;
	    if(Verbosity() > 2)
	      {
		std::cout << " lay: "  << local_layer
			  << " pedge " << pedge   
			  << " | " << npedge  
			  << " nredge " << nredge 
			  << " rphisize " << rphisize  
			  << " | " << nbig 
			  << " rovlp " << rovlp  
			  << "  | " << novlp  
			  << std::endl;
	      }
	  }
	}
      
        if (silseed)
        {
	  siqr = silseed->get_qOverR();
	  siphi = silseed->get_phi();
	  sithe = silseed->get_theta();
	  six0 = silseed->get_X0();
	  siy0 = silseed->get_Y0();
          for (SvtxTrack::ConstClusterKeyIter local_iter = silseed->begin_cluster_keys();
               local_iter != silseed->end_cluster_keys();
               ++local_iter)
          {
            TrkrDefs::cluskey cluster_key = *local_iter;
            // TrkrCluster* cluster = clustermap->findCluster(cluster_key);
            unsigned int local_layer = TrkrDefs::getLayer(cluster_key);

            if (_nlayers_maps > 0 && local_layer < _nlayers_maps)
            {
              maps[local_layer] = 1;
              nmaps++;
            }
            if (_nlayers_intt > 0 && local_layer >= _nlayers_maps && local_layer < _nlayers_maps + _nlayers_intt)
            {
              intt[local_layer - _nlayers_maps] = 1;
              nintt++;
            }
            if (_nlayers_mms > 0 && local_layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc && local_layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms)
            {
              mms[local_layer - (_nlayers_maps + _nlayers_intt + _nlayers_tpc)] = 1;
              nmms++;
            }
            if (_nlayers_tpc > 0 && local_layer >= (_nlayers_maps + _nlayers_intt) && local_layer < (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
            {
              tpc[local_layer - (_nlayers_maps + _nlayers_intt)] = 1;
              ntpc++;
              if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 8)
              {
                ntpc11++;
              }

              if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 16)
              {
                ntpc1++;
              }
              else if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 32)
              {
                ntpc2++;
              }
              else if ((local_layer - (_nlayers_maps + _nlayers_intt)) < 48)
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

        /// this is the global vertex
        int vertexID = track->get_vertex_id();
        GlobalVertex* vertex = gvertexmap->get(vertexID);
        float vx = NAN;
        float vy = NAN;
        float vz = NAN;
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

        float gtrackID = NAN;
        float gflavor = NAN;
        float ng4hits = NAN;
        unsigned int ngmaps = 0;
        unsigned int ngintt = 0;
        unsigned int ngmms = 0;
        unsigned int ngtpc = 0;
        unsigned int nglmaps = 0;
        unsigned int nglintt = 0;
        unsigned int nglmms = 0;
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

	int ispure = 0;
        float nfromtruth = NAN;
        float nwrong = NAN;
        float ntrumaps = NAN;
	float nwrongmaps = NAN;
        float ntruintt = NAN;
	float nwrongintt = NAN;
        float ntrumms = NAN;
	float nwrongmms = NAN;
        float ntrutpc = NAN;
	float nwrongtpc = NAN;
        float ntrutpc1 = NAN;
	float nwrongtpc1 = NAN;
        float ntrutpc11 = NAN;
	float nwrongtpc11 = NAN;
        float ntrutpc2 = NAN;
	float nwrongtpc2 = NAN;
        float ntrutpc3 = NAN;
	float nwrongtpc3 = NAN;
        float layersfromtruth = NAN;

        if (_do_track_match)
        {
          PHG4Particle* g4particle = trackeval->max_truth_particle_by_nclusters(track);
          if (g4particle)
          {
            if (_scan_for_embedded)
            {
              if (trutheval->get_embed(g4particle) <= 0)
              {
                continue;
              }
            }
	    SvtxTrack* truthrecotrk = trackeval->best_track_from(g4particle);
	    if(truthrecotrk)
	      {
		if(truthrecotrk->get_id() == track->get_id())
		  {
		    ispure = 1;
		  }
	      }
            gtrackID = g4particle->get_track_id();
            gflavor = g4particle->get_pid();

            std::set<TrkrDefs::cluskey> g4clusters = clustereval->all_clusters_from(g4particle);
            ng4hits = g4clusters.size();
            gpx = g4particle->get_px();
            gpy = g4particle->get_py();
            gpz = g4particle->get_pz();

	    std::vector<int> lmaps(_nlayers_maps + 1,0);
	    std::vector<int> lintt(_nlayers_intt + 1,0);
	    std::vector<int> ltpc(_nlayers_tpc + 1,0);
	    std::vector<int> lmms(_nlayers_mms + 1,0);

            for (const TrkrDefs::cluskey g4cluster : g4clusters)
            {
              unsigned int local_layer = TrkrDefs::getLayer(g4cluster);
              if (_nlayers_maps > 0 && local_layer < _nlayers_maps)
              {
                lmaps[local_layer] = 1;
                ngmaps++;
              }

              if (_nlayers_intt > 0 && local_layer >= _nlayers_maps && local_layer < _nlayers_maps + _nlayers_intt)
              {
                lintt[local_layer - _nlayers_maps] = 1;
                ngintt++;
              }

              if (_nlayers_tpc > 0 && local_layer >= _nlayers_maps + _nlayers_intt && local_layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
              {
                ltpc[local_layer - (_nlayers_maps + _nlayers_intt)] = 1;
                ngtpc++;
              }

              if (_nlayers_mms > 0 && local_layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc && local_layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms)
              {
                lmms[local_layer - (_nlayers_maps + _nlayers_intt + _nlayers_tpc)] = 1;
                ngmms++;
              }
            }
            if (_nlayers_maps > 0)
            {
              for (unsigned int i = 0; i < _nlayers_maps; i++)
              {
                nglmaps += lmaps[i];
              }
            }
            if (_nlayers_intt > 0)
            {
              for (unsigned int i = 0; i < _nlayers_intt; i++)
              {
                nglintt += lintt[i];
              }
            }
            if (_nlayers_tpc > 0)
            {
              for (unsigned int i = 0; i < _nlayers_tpc; i++)
              {
                ngltpc += ltpc[i];
              }
            }
            if (_nlayers_mms > 0)
            {
              for (unsigned int i = 0; i < _nlayers_mms; i++)
              {
                nglmms += lmms[i];
              }
            }

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
            {
              outerhit = trutheval->get_outermost_truth_hit(g4particle);
            }
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
	      auto pair = trackeval->get_layer_range_contribution(track, g4particle, 0, _nlayers_maps);
	      ntrumaps = pair.first;
	      nwrongmaps = pair.second;
            }
            if (_nlayers_intt == 0)
            {
              ntruintt = 0;
            }
            else
            {
              auto pair =  trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps, _nlayers_maps + _nlayers_intt);
	      ntruintt = pair.first;
	      nwrongintt = pair.second;
            }
            if (_nlayers_mms == 0)
            {
              ntrumms = 0;
            }
            else
            {
              auto pair = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt + _nlayers_tpc, _nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms);
	      ntrumms = pair.first;
	      nwrongmms = pair.second;
            }
            auto pair  = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt, _nlayers_maps + _nlayers_intt + _nlayers_tpc);
	    ntrutpc = pair.first;
	    nwrongtpc = pair.second;
            pair = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt, _nlayers_maps + _nlayers_intt + 16);
	    ntrutpc1 = pair.first;
	    nwrongtpc1 = pair.second;
            pair = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt, _nlayers_maps + _nlayers_intt + 8);
	    ntrutpc11 = pair.first;
	    nwrongtpc11 = pair.second;
            pair = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt + 16, _nlayers_maps + _nlayers_intt + 32);
	    ntrutpc2 = pair.first;
	    nwrongtpc2 = pair.second;
            pair = trackeval->get_layer_range_contribution(track, g4particle, _nlayers_maps + _nlayers_intt + 32, _nlayers_maps + _nlayers_intt + _nlayers_tpc);
	    ntrutpc3 = pair.first;
	    nwrongtpc3 = pair.second;
            layersfromtruth = trackeval->get_nclusters_contribution_by_layer(track, g4particle);
          }
        }
	if(Verbosity() > 2)
	  {
	    std::cout << " npedge "  << npedge  
		      << " nredge "  << nredge  
		      << " nbig " << nbig 
		      << " novlp "<< novlp  
		      << std::endl;
	  }
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
			      siqr,
			      siphi,
			      sithe,
			      six0,
			      siy0,
			      tpqr,
			      tpphi,
			      tpthe,
			      tpx0,
			      tpy0,
                              charge,
                              quality,
                              chisq,
                              ndf,
                              local_nhits, nmaps, nintt, ntpc, nmms,
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
                              gtrackID,
			      (float)ispure,
                              gflavor,
                              ng4hits,
                              (float) ngmaps,
                              (float) ngintt,
                              (float) ngtpc,
                              (float) ngmms,
                              (float) nglmaps,
                              (float) nglintt,
                              (float) ngltpc,
                              (float) nglmms,
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
			      nwrongmaps,
                              ntruintt,
			      nwrongintt,
                              ntrutpc,
			      nwrongtpc,
                              ntrumms,
			      nwrongmms,
                              ntrutpc1,
			      nwrongtpc1,
                              ntrutpc11,
			      nwrongtpc11,
                              ntrutpc2,
			      nwrongtpc2,
                              ntrutpc3,
			      nwrongtpc3,
                              layersfromtruth,
			      npedge,
			      nredge,
			      nbig,
			      novlp,
			      merr,
			      msize,
                              nhit_tpc_all,
                              nhit_tpc_in,
                              nhit_tpc_mid,
                              nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

        if (Verbosity() >= 1)
        {
          std::cout << "ievent " << _ievent
                    << " trackID " << trackID
                    << " nhits " << local_nhits
                    << " px " << px
                    << " py " << py
                    << " pz " << pz
                    << " gembed " << gembed
                    << " gprimary " << gprimary
                    << std::endl;
        }

        _ntp_track->Fill(track_data);
      }
    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      std::cout << "track time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
    }
  }

  //---------------------
  // fill the Gseed NTuple
  //---------------------

  if (_ntp_gseed)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_gseed " << std::endl;
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

    float *xval = new float[_nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms];
    float *yval = new float[_nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms];
    float *zval = new float[_nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms];
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
        for (unsigned int i = 0; i < _nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms; i++)
        {
          xval[i] = 0;
          yval[i] = 0;
          zval[i] = 0;
        }
        std::set<PHG4Hit*> truth_hits = trutheval->all_truth_hits(g4particle);
        for (auto g4hit : truth_hits)
        {
          unsigned int local_layer = g4hit->get_layer();
          // std::cout << "  g4hit " << g4hit->get_hit_id() << " layer = " << local_layer << std::endl;
          if (local_layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms)
          {
            // std::cout << PHWHERE << " skipping out of bounds detector id " << local_layer << std::endl;
            continue;
          }
          xval[local_layer] = g4hit->get_avg_x();
          yval[local_layer] = g4hit->get_avg_y();
          zval[local_layer] = g4hit->get_avg_z();
        }

        for (unsigned int i = 0; i < _nlayers_maps + _nlayers_intt + _nlayers_tpc + _nlayers_mms; i++)
        {
          gx = xval[i];
          gy = yval[i];
          gz = zval[i];
          if (gx == 0 && gy == 0)
          {
            continue;
          }

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
                                nhit_tpc_out, nclus_all, nclus_tpc, nclus_intt, nclus_maps, nclus_mms};

          _ntp_gseed->Fill(gseed_data);
        }
	delete [] xval;
	delete [] yval;
	delete [] zval;
      }
    }

    if (Verbosity() > 1)
    {
      _timer->stop();
      std::cout << "g4hit time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
    }
  }
  return;
}

TMatrixF SvtxEvaluator::calculateClusterError(TrkrCluster* c, float& clusphi)
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
