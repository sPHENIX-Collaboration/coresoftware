#include "TrkrNtuplizer.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrackFitUtils.h>

#include <trackbase_historic/TrackSeedContainer_v1.h>

#include <micromegas/MicromegasDefs.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrHitSet.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
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

enum n_event
{
  evnev,
  evnseed,
  evsize = evnseed + 1,
};

enum n_info
{
  infonocc11,
  infonocc116,
  infonocc21,
  infonocc216,
  infonocc31,
  infonocc316,
  infonhitmvtx,
  infontrk,
  infonhitintt,
  infonhittpot,
  infonhittpcall,
  infonhittpcin,
  infonhittpcmid,
  infonhittpcout,
  infonclusall,
  infonclustpc,
  infonclusintt,
  infonclusmvtx,
  infonclustpot,
  infosize = infonclustpot + 1
};

enum n_vertex
{
  vtxnvertexID,
  vtxnvx,
  vtxnvy,
  vtxnvz,
  vtxnntracks,
  vtxnchi2,
  vtxnndof,
  vtxsize = vtxnndof + 1
};

enum n_hit
{
  nhitID,
  nhite,
  nhitadc,
  nhitlayer,
  nhitphielem,
  nhitzelem,
  nhitcellID,
  nhitecell,
  nhitphibin,
  nhittbin,
  nhitphi,
  nhitr,
  nhitx,
  nhity,
  nhitz,
  hitsize = nhitz + 1
};

enum n_seed
{
  nseedtrackID,
  nseedniter,
  nseedpt,
  nseedeta,
  nseedphi,
  nseedsyxint,
  nseedsrzint,
  nseedsxyslope,
  nseedsrzslope,
  nseedX0,
  nseedY0,
  nseedZ0,
  nseedR0,
  nseedcharge,
  nseeddedx,
  nseedn1pix,
  nseednhits,
  seedsize = nseednhits + 1
};

enum n_residual
{
  nresalpha,
  nresbeta,
  nresphio,
  nresphi,
  nresz,
  ressize = nresz + 1
};

enum n_track
{
  ntrktrackID,
  ntrkcrossing,
  ntrkpx,
  ntrkpy,
  ntrkpz,
  ntrkpt,
  ntrketa,
  ntrkphi,
  ntrkdeltapt,
  ntrkdeltaeta,
  ntrkdeltaphi,
  ntrkcharge,
  ntrkquality,
  ntrkchisq,
  ntrkndf,
  ntrknhits,
  ntrknmaps,
  ntrknintt,
  ntrkntpc,
  ntrknmms,
  ntrkntpc1,
  ntrkntpc11,
  ntrkntpc2,
  ntrkntpc3,
  ntrkndedx,
  ntrkvertexID,
  ntrkvx,
  ntrkvy,
  ntrkvz,
  ntrkdca2d,
  ntrkdca2dsigma,
  ntrkdca3dxy,
  ntrkdca3dxysigma,
  ntrkdca3dz,
  ntrkdca3dzsigma,
  ntrkpcax,
  ntrkpcay,
  ntrkpcaz,
  ntrkhlxpt,
  ntrkhlxeta,
  ntrkhlxphi,
  ntrkhlxX0,
  ntrkhlxY0,
  ntrkhlxZ0,
  ntrkhlxcharge,
  trksize = ntrkhlxcharge + 1
};

enum n_cluster
{
  nclulocx,
  nclulocy,
  nclux,
  ncluy,
  ncluz,
  nclur,
  ncluphi,
  nclueta,
  nclutheta,
  ncluphibin,
  nclutbin,
  ncluex,
  ncluey,
  ncluez,
  ncluephi,
  nclupez,
  nclupephi,
  nclue,
  ncluadc,
  nclumaxadc,
  nclulayer,
  ncluphielem,
  ncluzelem,
  nclusize,
  ncluphisize,
  ncluzsize,
  nclupedge,
  ncluredge,
  ncluovlp,
  nclutrackID,
  ncluniter,
  clusize = ncluniter + 1
};

TrkrNtuplizer::TrkrNtuplizer(const string& /*name*/, const string& filename, const string& trackmapname,
                             unsigned int nlayers_maps,
                             unsigned int nlayers_intt,
                             unsigned int nlayers_tpc,
                             unsigned int nlayers_mms)
  : SubsysReco("TrkrNtuplizer")
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _nlayers_mms(nlayers_mms)
  , _filename(filename)
  , _trackmapname(trackmapname)
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
  string str_vertex = {"vertexID:vx:vy:vz:ntracks:chi2:ndof"};
  string str_event = {"event:seed"};
  string str_hit = {"hitID:e:adc:layer:phielem:zelem:cellID:ecell:phibin:tbin:phi:r:x:y:z"};
  string str_cluster = {"locx:locy:x:y:z:r:phi:eta:theta:phibin:tbin:ex:ey:ez:ephi:pez:pephi:e:adc:maxadc:layer:phielem:zelem:size:phisize:zsize:pedge:redge:ovlp:trackID:niter"};
  string str_seed = {"seedID:siter:spt:seta:sphi:syxint:srzint:sxyslope:srzslope:sX0:sY0:sdZ0:sR0:scharge:sdedx:sn1pix:snhits"};
  string str_residual = {"alpha:beta:resphio:resphi:resz"};
  string str_track = {"trackID:crossing:px:py:pz:pt:eta:phi:deltapt:deltaeta:deltaphi:charge:quality:chisq:ndf:nhits:nmaps:nintt:ntpc:nmms:ntpc1:ntpc11:ntpc2:ntpc3:dedx:nlmaps:nlintt:nltpc:nlmms:layers:vertexID:vx:vy:vz:dca2d:dca2dsigma:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:hlxpt:hlxeta:hlxphi:hlxX0:hlxY0:hlxZ0:hlxcharge"};
  string str_info = {"occ11:occ116:occ21:occ216:occ31:occ316:ntrk:nhitmvtx:nhitintt:nhittpot:nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms"};

  if (_do_info_eval)
  {
    string ntp_varlist_info = str_event + ":" + str_info;
    _ntp_info = new TNtuple("ntp_info", "event info", ntp_varlist_info.c_str());
  }

  if (_do_vertex_eval)
  {
    string ntp_varlist_vtx = str_event + ":" + str_vertex + ":" + str_info;
    _ntp_vertex = new TNtuple("ntp_vertex", "vertex => max truth", ntp_varlist_vtx.c_str());
  }

  if (_do_hit_eval)
  {
    string ntp_varlist_ev = str_event + ":" + str_hit + ":" + str_info;
    _ntp_hit = new TNtuple("ntp_hit", "svtxhit => max truth", ntp_varlist_ev.c_str());
  }

  if (_do_cluster_eval)
  {
    string ntp_varlist_clu = str_event + ":" + str_cluster + ":" + str_info;
    _ntp_cluster = new TNtuple("ntp_cluster", "svtxcluster => max truth", ntp_varlist_clu.c_str());
  }
  if (_do_clus_trk_eval)
  {
    string ntp_varlist_clut = str_event + ":" + str_cluster + ":" + str_residual + ":" + str_seed + ":" + str_info;
    _ntp_clus_trk = new TNtuple("ntp_clus_trk", "cluster on track", ntp_varlist_clut.c_str());
  }

  if (_do_track_eval)
  {
    string ntp_varlist_trk = str_event + ":" + str_track + ":" + str_info;
    _ntp_track = new TNtuple("ntp_track", "svtxtrack => max truth", ntp_varlist_trk.c_str());
  }

  if (_do_tpcseed_eval)
  {
    string ntp_varlist_tsee = str_event + ":" + str_seed + ":" + str_info;
    _ntp_tpcseed = new TNtuple("ntp_tpcseed", "seeds from truth", ntp_varlist_tsee.c_str());
  }
  if (_do_siseed_eval)
  {
    string ntp_varlist_ssee = str_event + ":" + str_seed + ":" + str_info;
    _ntp_siseed = new TNtuple("ntp_siseed", "seeds from truth", ntp_varlist_ssee.c_str());
  }
  _timer = new PHTimer("_eval_timer");
  _timer->stop();
  /**/
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkrNtuplizer::InitRun(PHCompositeNode* /*unused*/)
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
    m_fSeed = std::numeric_limits<float>::quiet_NaN();
  }
  if (_trackmap == nullptr)
  {
    _trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname.c_str());
  }

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if (!_cluster_map)
  {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }

  _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!_tgeometry)
  {
    std::cout << "No Acts geometry on node tree. Can't  continue."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  {
    if (!_geom_container)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

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

    if (_cluster_map != nullptr)
    {
      unsigned int icluster = 0;
      for (const auto& hitsetkey : _cluster_map->getHitSetKeys())
      {
        auto range = _cluster_map->getClusters(hitsetkey);
        for (auto iter = range.first; iter != range.second; ++iter)
        {
          TrkrDefs::cluskey cluster_key = iter->first;
          cout << icluster << " with key " << cluster_key << " of " << _cluster_map->size();
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
      for (auto& iter : *_trackmap)
      {
        cout << itrack << " of " << _trackmap->size();
        SvtxTrack* track = iter.second;
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

    float vx = std::numeric_limits<float>::quiet_NaN();
    float vy = std::numeric_limits<float>::quiet_NaN();
    float vz = std::numeric_limits<float>::quiet_NaN();

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

    unsigned int nclusters[100] = {0};
    unsigned int nhits[100] = {0};

    _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

    if (!_tgeometry)
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
        auto range = _cluster_map->getClusters(hitsetiter->first);
        for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
        {
          const auto cluskey = clusIter->first;
          nclusters[TrkrDefs::getLayer(cluskey)]++;
        }
      }
    }

    for (unsigned int ilayer = 0; ilayer < _nlayers_maps + _nlayers_intt + _nlayers_tpc; ++ilayer)
    {
      PHG4TpcCylinderGeom* GeoLayer = _geom_container->GetLayerCellGeom(ilayer);

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

  float fx_event[n_event::evsize]{(float) _ievent, (float) _iseed};
  float fx_info[n_info::infosize] = {0};

  float nhit[100];
  for (float& i : nhit)
  {
    i = 0;
  }

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

        if (layer < _nlayers_maps)
        {
          fx_info[n_info::infonhitmvtx]++;
        }
        if (layer >= _nlayers_maps && layer < _nlayers_maps + _nlayers_intt)
        {
          fx_info[n_info::infonhitintt]++;
        }
        if ((float) layer >= _nlayers_maps + _nlayers_intt &&
            (float) layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
        {
          fx_info[n_info::infonhittpcall]++;
        }
        if ((float) layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc)
        {
          fx_info[n_info::infonhittpot]++;
        }
        if ((float) layer == _nlayers_maps + _nlayers_intt)
        {
          fx_info[n_info::infonhittpcin]++;
        }
        if ((float) layer == _nlayers_maps + _nlayers_intt + _nlayers_tpc - 1)
        {
          fx_info[n_info::infonhittpcout]++;
        }
        // NOLINTNEXTLINE(bugprone-integer-division)
        if ((float) layer == _nlayers_maps + _nlayers_intt + _nlayers_tpc / 2 - 1)
        {
          fx_info[n_info::infonhittpcmid]++;
        }
      }
    }
  }

  _geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!_geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return;
  }
  PHG4TpcCylinderGeom* GeoLayer;
  int layer = _nlayers_maps + _nlayers_intt;
  GeoLayer = _geom_container->GetLayerCellGeom(layer);
  int nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  float nhits = nhit[layer];
  fx_info[n_info::infonocc11] = nhits / nbins;

  layer = _nlayers_maps + _nlayers_intt + 15;
  GeoLayer = _geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  fx_info[n_info::infonocc116] = nhits / nbins;

  layer = _nlayers_maps + _nlayers_intt + 16;
  GeoLayer = _geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  fx_info[n_info::infonocc21] = nhits / nbins;
  layer = _nlayers_maps + _nlayers_intt + 31;
  GeoLayer = _geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  fx_info[n_info::infonocc216] = nhits / nbins;
  layer = _nlayers_maps + _nlayers_intt + 32;
  GeoLayer = _geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  fx_info[n_info::infonocc31] = nhits / nbins;
  layer = _nlayers_maps + _nlayers_intt + 47;
  GeoLayer = _geom_container->GetLayerCellGeom(layer);
  nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
  nhits = nhit[layer];
  fx_info[n_info::infonocc316] = nhits / nbins;

  if (Verbosity() > 1)
  {
    cout << " occ11 = " << fx_info[n_info::infonocc11]
         << " occ116 = " << fx_info[n_info::infonocc116]
         << " occ21 = " << fx_info[n_info::infonocc21]
         << " occ216 = " << fx_info[n_info::infonocc216]
         << " occ31 = " << fx_info[n_info::infonocc31]
         << " occ316 = " << fx_info[n_info::infonocc316]
         << endl;
  }

  if (_cluster_map)
  {
    fx_info[n_info::infonclusall] = _cluster_map->size();

    for (const auto& hitsetkey : _cluster_map->getHitSetKeys())
    {
      auto range = _cluster_map->getClusters(hitsetkey);
      for (auto iter_cin = range.first; iter_cin != range.second; ++iter_cin)
      {
        TrkrDefs::cluskey cluster_key = iter_cin->first;
        unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
        if (_nlayers_maps > 0)
        {
          if (layer_local < _nlayers_maps)
          {
            fx_info[n_info::infonclusmvtx]++;
          }
        }
        if (_nlayers_intt > 0)
        {
          if ((layer_local >= _nlayers_maps) && (layer_local < (_nlayers_maps + _nlayers_intt)))
          {
            fx_info[n_info::infonclusintt]++;
          }
        }
        if (_nlayers_tpc > 0)
        {
          if (layer_local >= (_nlayers_maps + _nlayers_intt) && layer_local < (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
          {
            fx_info[n_info::infonclustpc]++;
          }
        }
        if (_nlayers_mms > 0)
        {
          if (layer_local >= (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
          {
            fx_info[n_info::infonclustpot]++;
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
      fx_info[n_info::infontrk] = (float) _trackmap->size();
    }
    if (Verbosity() > 0)
    {
      cout << "EVENTINFO SEED: " << m_fSeed << endl;
      cout << "EVENTINFO NHIT: " << setprecision(9) << fx_info[n_info::infonclusall] << endl;
      cout << "EVENTINFO CLUSTPC: " << fx_info[n_info::infonclustpc] << endl;
      cout << "EVENTINFO NTRKREC: " << ntrk << endl;
    }

    float* info_data = new float[n_info::infosize + n_event::evsize];
    std::copy(fx_event, fx_event + n_event::evsize, info_data);
    std::copy(fx_info, fx_info + n_info::infosize, info_data + n_event::evsize);
    /*
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
    */
    _ntp_info->Fill(info_data);
    //    delete info_data;
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
    float fx_vertex[n_vertex::vtxsize];
    for (float& i : fx_vertex)
    {
      i = 0;
    }

    //    SvtxVertexMap* vertexmap = nullptr;

    //    vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");  // Acts vertices

    float vx = std::numeric_limits<float>::quiet_NaN();
    float vy = std::numeric_limits<float>::quiet_NaN();
    float vz = std::numeric_limits<float>::quiet_NaN();
    float ntracks = std::numeric_limits<float>::quiet_NaN();
    fx_vertex[vtxnvx] = vx;
    fx_vertex[vtxnvy] = vy;
    fx_vertex[vtxnvz] = vz;
    fx_vertex[vtxnntracks] = ntracks;
    if (Verbosity() > 1)
    {
      std::cout << " adding vertex data " << std::endl;
    }
    float* vertex_data = new float[n_info::infosize + n_event::evsize + n_vertex::vtxsize];
    std::copy(fx_event, fx_event + n_event::evsize, vertex_data);
    std::copy(fx_vertex, fx_vertex + n_vertex::vtxsize, vertex_data + n_event::evsize);
    std::copy(fx_info, fx_info + n_info::infosize, vertex_data + n_event::evsize + n_vertex::vtxsize);
    _ntp_vertex->Fill(vertex_data);
  }
  if (Verbosity() > 1)
  {
    _timer->stop();
    cout << "vertex time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
  }
  //--------------------
  // fill the Hit NTuple
  //--------------------

  if (_ntp_hit)
  {
    float fx_hit[n_hit::hitsize] = {0};
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

          fx_hit[n_hit::nhitID] = hit_key;
          fx_hit[n_hit::nhite] = hit->getEnergy();
          fx_hit[n_hit::nhitadc] = hit->getAdc();
          unsigned int layer_local = TrkrDefs::getLayer(hitset_key);
          fx_hit[n_hit::nhitlayer] = (float) layer_local;
          fx_hit[n_hit::nhitphielem] = -666;
          fx_hit[n_hit::nhitzelem] = -666;

          if (layer_local >= 3 && layer_local < 7)
          {
            fx_hit[n_hit::nhitphielem] = InttDefs::getLadderPhiId(hitset_key);
            fx_hit[n_hit::nhitzelem] = InttDefs::getLadderZId(hitset_key);
          }
          if (layer_local >= 7 && layer_local < 55)
          {
            fx_hit[n_hit::nhitphielem] = TpcDefs::getSectorId(hitset_key);
            fx_hit[n_hit::nhitzelem] = TpcDefs::getSide(hitset_key);
          }
          /*
          if(layer_local>=55){
            if(MicromegasDefs::getSegmentationType(hitset_key)==MicromegasDefs::SEGMENTATION_Z){
              sector = 1;
              side = MicromegasDefs::getStrip(hit_key);
            }else{
              sector =MicromegasDefs::getStrip(hit_key);
              side =  1;
            }
          }
          */
          fx_hit[n_hit::nhitphielem] = TpcDefs::getSectorId(hitset_key);
          fx_hit[n_hit::nhitzelem] = TpcDefs::getSide(hitset_key);
          fx_hit[n_hit::nhitcellID] = 0;
          fx_hit[n_hit::nhitecell] = hit->getAdc();
          fx_hit[n_hit::nhitphibin] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhittbin] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhitphi] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhitr] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhitx] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhity] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhitz] = std::numeric_limits<float>::quiet_NaN();

          if (layer_local >= _nlayers_maps + _nlayers_intt && layer_local < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
          {
            PHG4TpcCylinderGeom* GeoLayer_local = _geom_container->GetLayerCellGeom(layer_local);
            double radius = GeoLayer_local->get_radius();
            fx_hit[n_hit::nhitphibin] = (float) TpcDefs::getPad(hit_key);
            fx_hit[n_hit::nhittbin] = (float) TpcDefs::getTBin(hit_key);
            fx_hit[n_hit::nhitphi] = GeoLayer_local->get_phicenter(fx_hit[n_hit::nhitphibin]);
            float phi = GeoLayer_local->get_phicenter(TpcDefs::getPad(hit_key));
            float clockperiod = GeoLayer_local->get_zstep();
            auto glob = m_tGeometry->getGlobalPositionTpc(hitset_key, hit_key, phi, radius, clockperiod);
            
            fx_hit[n_hit::nhitz] = glob.z();
            fx_hit[n_hit::nhitr] = radius;
            fx_hit[n_hit::nhitx] = glob.x();
            fx_hit[n_hit::nhity] = glob.y();
          }

          float* hit_data = new float[n_info::infosize + n_event::evsize + n_hit::hitsize];
          std::copy(fx_event, fx_event + n_event::evsize, hit_data);
          std::copy(fx_hit, fx_hit + n_hit::hitsize, hit_data + n_event::evsize);
          std::copy(fx_info, fx_info + n_info::infosize, hit_data + n_event::evsize + n_hit::hitsize);
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

    TrkrHitSetContainer* hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    //    TrkrClusterIterationMapv1* _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    ClusterErrorPara ClusErrPara;

    if (Verbosity() > 1)
    {
      if (_cluster_map != nullptr)
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

    if (_cluster_map && hitsets)
    {
      for (const auto& hitsetkey : _cluster_map->getHitSetKeys())
      {
        auto range = _cluster_map->getClusters(hitsetkey);
        for (auto iter = range.first; iter != range.second; ++iter)
        {
          TrkrDefs::cluskey cluster_key = iter->first;
          Float_t fx_cluster[n_cluster::clusize];
          FillCluster(&fx_cluster[0], cluster_key);

          float* cluster_data = new float[n_info::infosize + n_event::evsize + n_cluster::clusize];
          std::copy(fx_event, fx_event + n_event::evsize, cluster_data);
          std::copy(fx_cluster, fx_cluster + n_cluster::clusize, cluster_data + n_event::evsize);
          std::copy(fx_info, fx_info + n_info::infosize, cluster_data + n_event::evsize + n_cluster::clusize);
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
    
    float drphi[2][60];
 drphi[0][0] = 0.000000; drphi[1][0] = 0.000000;
 drphi[0][0] = 0.000000; drphi[1][0] = 0.000000;
 drphi[0][1] = 0.000000; drphi[1][1] = 0.000000;
 drphi[0][2] = 0.000000; drphi[1][2] = 0.000000;
 drphi[0][3] = 0.000000; drphi[1][3] = 0.000000;
 drphi[0][4] = 0.000000; drphi[1][4] = 0.000000;
 drphi[0][5] = 0.000000; drphi[1][5] = 0.000000;
 drphi[0][6] = 0.000000; drphi[1][6] = 0.000000;
 drphi[0][7] = -0.000089; drphi[1][7] = 0.021540;
 drphi[0][8] = 0.001072; drphi[1][8] = 0.013116;
 drphi[0][9] = 0.001267; drphi[1][9] = 0.014194;
 drphi[0][10] = 0.000303; drphi[1][10] = 0.006583;
 drphi[0][11] = 0.001176; drphi[1][11] = 0.008321;
 drphi[0][12] = 0.002913; drphi[1][12] = -0.000341;
 drphi[0][13] = 0.002765; drphi[1][13] = 0.001068;
 drphi[0][14] = 0.003330; drphi[1][14] = -0.003894;
 drphi[0][15] = 0.003533; drphi[1][15] = -0.003686;
 drphi[0][16] = 0.004609; drphi[1][16] = -0.009074;
 drphi[0][17] = 0.005643; drphi[1][17] = -0.009949;
 drphi[0][18] = 0.005311; drphi[1][18] = -0.015016;
 drphi[0][19] = 0.003438; drphi[1][19] = -0.016782;
 drphi[0][20] = -0.001007; drphi[1][20] = -0.027177;
 drphi[0][21] = -0.011399; drphi[1][21] = -0.037536;
 drphi[0][22] = -0.042440; drphi[1][22] = -0.075852;
 drphi[0][23] = 0.045325; drphi[1][23] = 0.058216;
 drphi[0][24] = 0.014062; drphi[1][24] = 0.024589;
 drphi[0][25] = 0.001659; drphi[1][25] = 0.016468;
 drphi[0][26] = -0.006583; drphi[1][26] = 0.012170;
 drphi[0][27] = -0.006997; drphi[1][27] = 0.010027;
 drphi[0][28] = -0.006659; drphi[1][28] = 0.007509;
 drphi[0][29] = -0.006740; drphi[1][29] = 0.006784;
 drphi[0][30] = -0.007644; drphi[1][30] = 0.005510;
 drphi[0][31] = -0.009417; drphi[1][31] = 0.004783;
 drphi[0][32] = -0.011589; drphi[1][32] = 0.002409;
 drphi[0][33] = -0.012136; drphi[1][33] = 0.003086;
 drphi[0][34] = -0.014818; drphi[1][34] = 0.000434;
 drphi[0][35] = -0.017750; drphi[1][35] = -0.000202;
 drphi[0][36] = -0.022697; drphi[1][36] = -0.000860;
 drphi[0][37] = -0.033155; drphi[1][37] = -0.009810;
 drphi[0][38] = -0.068933; drphi[1][38] = -0.043357;
 drphi[0][39] = 0.068807; drphi[1][39] = 0.032688;
 drphi[0][40] = 0.035268; drphi[1][40] = 0.001179;
 drphi[0][40] = 0.028104; drphi[1][40] = -0.002723;
 drphi[0][41] = 0.023703; drphi[1][41] = -0.002551;
 drphi[0][42] = 0.019526; drphi[1][42] = -0.001358;
 drphi[0][43] = 0.016406; drphi[1][43] = -0.001094;
 drphi[0][44] = 0.013309; drphi[1][44] = -0.000676;
 drphi[0][45] = 0.011975; drphi[1][45] = 0.001872;
 drphi[0][46] = 0.008119; drphi[1][46] = 0.001131;
 drphi[0][47] = 0.006167; drphi[1][47] = 0.001166;
 drphi[0][48] = 0.002595; drphi[1][48] = 0.003025;
 drphi[0][49] = -0.000863; drphi[1][49] = 0.004523;
 drphi[0][50] = -0.005000; drphi[1][50] = 0.003310;
 drphi[0][51] = -0.008824; drphi[1][51] = 0.003575;
 drphi[0][52] = -0.017929; drphi[1][52] = 0.001060;
 drphi[0][53] = -0.040639; drphi[1][53] = -0.019704;
 drphi[0][54] = 0.000000; drphi[1][54] = 0.000000;
 drphi[0][55] = 0.000000; drphi[1][55] = 0.000000;
 drphi[0][56] = 0.000000; drphi[1][56] = 0.000000;
 drphi[0][57] = 0.000000; drphi[1][57] = 0.000000;

    
    if (Verbosity() > 1)
    {
      cout << "Filling ntp_clus_trk " << endl;
      if (_cluster_map != nullptr)
      {
        cout << "got clustermap" << endl;
      }
      else
      {
        cout << "no clustermap" << endl;
      }
    }

    TrackSeedContainer* _tpc_seeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
    if (!_tpc_seeds)
    {
      cerr << PHWHERE << " ERROR: Can't find "
           << "TpcTrackSeedContainer" << endl;
      return;
    }

    if (_trackmap)
    {
      int trackID = 0;
      for (auto& iter : *_trackmap)
      {
	trackID++;
        SvtxTrack* track = iter.second;
        TrackSeed* tpcseed = track->get_tpc_seed();
        TrackSeed* siseed = track->get_silicon_seed();
        std::vector<Acts::Vector3> clusterPositions;
        std::vector<Acts::Vector3> clusterPositionsVtx;
        std::vector<Acts::Vector3> clusterPositionsVtx2;
        std::vector<TrkrDefs::cluskey> clusterKeys;
        std::vector<TrkrDefs::cluskey> clusterKeysVtx;

        TrackFitUtils::position_vector_t xypoints, rzpoints;

        clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
                           tpcseed->end_cluster_keys());

        clusterKeysVtx.insert(clusterKeysVtx.end(), tpcseed->begin_cluster_keys(),
                           tpcseed->end_cluster_keys());

	/*	TrkrDefs::cluskey vtxkey = 0;
	auto keysit = clusterKeysVtx.begin();
	clusterKeysVtx.insert(keysit, vtxkey);
	*/
        /*if(siseed!=nullptr)
          clusterKeys.insert(clusterKeys.end(), siseed->begin_cluster_keys(),
                           siseed->end_cluster_keys());
        */
        TrackFitUtils::getTrackletClusters(_tgeometry, _cluster_map,
                                           clusterPositions, clusterKeys);
        TrackFitUtils::getTrackletClusters(_tgeometry, _cluster_map,
                                           clusterPositionsVtx, clusterKeys);

	/*	Acts::Vector3 vtx(0,0,0);
	auto vtxit = clusterPositionsVtx.begin();
	clusterPositionsVtx.insert(vtxit, vtx);
	*/
        for (auto& pos : clusterPositions)
        {
          float clusr = sqrt(pos.x() * pos.x() + pos.y() * pos.y());
          if (pos.y() < 0)
          {
            clusr *= -1;
          }

          // exclude silicon and tpot clusters for now
          if (fabs(clusr) > 80 || fabs(clusr) < 30)
          {
            continue;
          }
          rzpoints.push_back(std::make_pair(pos.z(), clusr));
          xypoints.push_back(std::make_pair(pos.x(), pos.y()));
        }
	for (unsigned int n = 0; n < clusterPositionsVtx.size();n++){
	  auto vtxit2 = clusterPositionsVtx2.begin();
	  /* if(n==0){
	    Acts::Vector3 vtx2(0,0,0);
	    clusterPositionsVtx2.insert(vtxit2, vtx2);	    
	  }else{
	  */
	    Acts::Vector3 cpos = clusterPositionsVtx.at(n);
	    TrkrDefs::cluskey corrkey = clusterKeysVtx.at(n);
	    unsigned int corrlayer = TrkrDefs::getLayer(corrkey);
	    float clusr = sqrt(cpos.x() * cpos.x() + cpos.y() * cpos.y());
	    //	    float clusr_corr = (clusr + (drcorr[corrlayer]*100))/clusr;
	    TVector2 vin(cpos.x(),cpos.y());
	    TVector2 vout;
	    int corrside = 0;//TpcDefs::getSide(corrkey);
	    if(cpos.z()>0) {
	      corrside=1;
}
	    vout.SetMagPhi(clusr,vin.Phi()-drphi[corrside][corrlayer]);
	    Acts::Vector3 point(vout.X(),vout.Y(),cpos.z());
	    clusterPositionsVtx2.insert(vtxit2+n, point);
	    // }
	}

        std::vector<float> fitparams_org = TrackFitUtils::fitClusters(clusterPositions, clusterKeys);
        std::vector<float> fitparams = TrackFitUtils::fitClusters(clusterPositionsVtx2, clusterKeysVtx);
        if (fitparams.size() == 0)
        {
          cout << "fit failed bailing...." << endl;
          continue;
        }

	//	std::cout << " fit 0 org " << fitparams_org.at(0) << " new: " << fitparams.at(0) << std::endl;
	//	std::cout << " fit 1 org " << fitparams_org.at(1) << " new: " << fitparams.at(1) << std::endl;
	//	std::cout << " fit 2 org " << fitparams_org.at(2) << " new: " << fitparams.at(2) << std::endl;
	//	std::cout << " fit 3 org " << fitparams_org.at(3) << " new: " << fitparams.at(3) << std::endl;
	//	std::cout << " fit 4 org " << fitparams_org.at(4) << " new: " << fitparams.at(4) << std::endl;

        float charge = std::numeric_limits<float>::quiet_NaN();
        if (tpcseed->get_qOverR() > 0)
        {
          charge = 1;
        }
        else
        {
          charge = -1;
        }

        //	      "pt:eta:phi:X0:Y0:charge:nhits:"
        float tpt = tpcseed->get_pt();
        float teta = tpcseed->get_eta();
        float tphi = tpcseed->get_phi();
        auto xyparams = TrackFitUtils::line_fit(xypoints);
        auto rzparams = TrackFitUtils::line_fit(rzpoints);
        float xyint = std::get<1>(xyparams);
        float xyslope = std::get<0>(xyparams);
        float rzint = std::get<1>(rzparams);
        float rzslope = std::get<0>(rzparams);
        float R0 = abs(-1 * xyint) / sqrt((xyslope * xyslope) + 1);
        float tX0 = tpcseed->get_X0();
        float tY0 = tpcseed->get_Y0();
        float tZ0 = tpcseed->get_Z0();

        float nhits_local = clusterPositions.size();
        if (Verbosity() > 1)
        {
          cout << " tpc: " << tpcseed->size_cluster_keys() << endl;
          if (siseed)
          {
            cout << " si " << siseed->size_cluster_keys() << endl;
          }
          cout << "done seedsize" << endl;
        }
        //      nhits_local += tpcseed->size_cluster_keys();
        // fill the Gseed NTuple
        //---------------------
	float dedx = calc_dedx(tpcseed);
	float n1pix = get_n1pix(tpcseed);
        float fx_seed[n_seed::seedsize] = {(float) trackID, 0, tpt, teta, tphi, xyint, rzint, xyslope, rzslope, tX0, tY0, tZ0, R0, charge, dedx, n1pix, nhits_local};

        if (_ntp_tpcseed)
        {
          float* tpcseed_data = new float[n_info::infosize + n_cluster::clusize + n_residual::ressize + n_seed::seedsize + n_event::evsize];
          std::copy(fx_event, fx_event + n_event::evsize, tpcseed_data);
          std::copy(fx_seed, fx_seed + n_seed::seedsize, tpcseed_data + n_event::evsize);
          std::copy(fx_info, fx_info + n_info::infosize, tpcseed_data + n_event::evsize + n_seed::seedsize);
          _ntp_tpcseed->Fill(tpcseed_data);
        }
        for (unsigned int i = 1; i < clusterPositionsVtx2.size(); i++)
        {
          TrkrDefs::cluskey cluster_key = clusterKeysVtx.at(i);
	  Acts::Vector3 position = clusterPositionsVtx2[i];
          Acts::Vector3 position_org = clusterPositions[i-1];
          Acts::Vector3 pca_org = TrackFitUtils::get_helix_pca(fitparams_org, position_org);
          Acts::Vector3 pca = TrackFitUtils::get_helix_pca(fitparams, position);

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
          float cluster_phi_org = atan2(position_org(1), position_org(0));
          float pca_phi_org = atan2(pca_org(1), pca_org(0));
          float dphi_org = cluster_phi_org - pca_phi_org;
          if (dphi_org > M_PI)
          {
            dphi_org = 2 * M_PI - dphi_org;
          }
          if (dphi_org < -M_PI)
          {
            dphi_org = 2 * M_PI + dphi_org;
          }
	  /*
	  std::cout << " i: " 
		    << " xo " << position_org(0)  
		    << " | " <<  position_org(1)  
		    << " | " <<  position_org(2)
      
		    << " x " <<  position(0)  
		    << " | " <<  position(1)  
		    << " | " <<  position(2)
      
		    << std::endl;
	  */
          float dz = position(2) - pca(2);

	  float resr = sqrt(position(0)*position(0)+position(1)*position(1)+position(2)*position(2));
	  float seedR = TMath::Abs(1.0 / tpcseed->get_qOverR());
	  float alpha = (resr * resr) / (2 * resr * seedR);
	  float beta = TMath::Abs(atan(tpcseed->get_slope()));

          float fx_res[n_residual::ressize] = {alpha,beta,dphi_org,dphi, dz};
          // sphi:syxint:srzint:sxyslope:srzslope:sX0:sY0:sdZ0:sR0

          float fx_cluster[n_cluster::clusize];
          //
          FillCluster(&fx_cluster[0], cluster_key);

          float* clus_trk_data = new float[n_info::infosize + n_cluster::clusize + n_residual::ressize + n_seed::seedsize + n_event::evsize];
          std::copy(fx_event, fx_event + n_event::evsize, clus_trk_data);
          std::copy(fx_cluster, fx_cluster + n_cluster::clusize, clus_trk_data + n_event::evsize);
          std::copy(fx_res, fx_res + n_residual::ressize, clus_trk_data + n_event::evsize + n_cluster::clusize);
          std::copy(fx_seed, fx_seed + n_seed::seedsize, clus_trk_data + n_event::evsize + n_cluster::clusize + n_residual::ressize);
          std::copy(fx_info, fx_info + n_info::infosize, clus_trk_data + n_event::evsize + n_cluster::clusize + n_residual::ressize + n_seed::seedsize);
          _ntp_clus_trk->Fill(clus_trk_data);
        }
      }
    }
  }

  //------------------------
  // fill the Track NTuple
  //------------------------

  if (_ntp_track)
  {
    if (Verbosity() >= 1)
    {
      cout << "Filling ntp_track " << endl;
      _timer->restart();
    }
    GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (_trackmap)
    {
      for (auto& iter : *_trackmap)
      {
        SvtxTrack* track = iter.second;
        float fx_track[n_track::trksize];
        FillTrack(&fx_track[0], track, vertexmap);
        float* track_data = new float[n_info::infosize + n_track::trksize + n_event::evsize];
        std::copy(fx_event, fx_event + n_event::evsize, track_data);
        std::copy(fx_track, fx_track + n_track::trksize, track_data + n_event::evsize);
        std::copy(fx_info, fx_info + n_info::infosize, track_data + n_event::evsize + n_track::trksize);
        _ntp_track->Fill(track_data);
      }
    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      cout << "track time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
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

void TrkrNtuplizer::FillTrack(float fX[50], SvtxTrack* track, GlobalVertexMap* vertexmap)
{
  float trackID = track->get_id();
  TrackSeed* tpcseed = track->get_tpc_seed();
  TrackSeed* silseed = track->get_silicon_seed();
  short int crossing_int = track->get_crossing();
  fX[n_track::ntrktrackID] = trackID;
  fX[n_track::ntrkcrossing] = std::numeric_limits<float>::quiet_NaN();
  if (crossing_int == SHRT_MAX)
  {
    fX[n_track::ntrkcrossing] = std::numeric_limits<float>::quiet_NaN();
  }
  else
  {
    fX[n_track::ntrkcrossing] = (float) crossing_int;
  }
  fX[n_track::ntrkcharge] = track->get_charge();
  fX[n_track::ntrkquality] = track->get_quality();
  fX[n_track::ntrkchisq] = track->get_chisq();
  fX[n_track::ntrkndf] = track->get_ndf();
  fX[n_track::ntrknhits] = 0;
  if (tpcseed)
  {
    fX[n_track::ntrknhits] += tpcseed->size_cluster_keys();
  }
  if (silseed)
  {
    fX[n_track::ntrknhits] += silseed->size_cluster_keys();
  }

  fX[n_track::ntrknmaps] = 0;
  fX[n_track::ntrknintt] = 0;
  fX[n_track::ntrknmms] = 0;
  fX[n_track::ntrkntpc] = 0;
  fX[n_track::ntrkntpc1] = 0;
  fX[n_track::ntrkntpc11] = 0;
  fX[n_track::ntrkntpc2] = 0;
  fX[n_track::ntrkntpc3] = 0;
  fX[n_track::ntrkndedx] = 0;
  if (tpcseed)
  {
    fX[n_track::ntrkndedx] = calc_dedx(tpcseed);
    for (SvtxTrack::ConstClusterKeyIter iter_local = tpcseed->begin_cluster_keys();
         iter_local != tpcseed->end_cluster_keys();
         ++iter_local)
    {
      TrkrDefs::cluskey cluster_key = *iter_local;
      unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
      switch (TrkrDefs::getTrkrId(cluster_key))
      {
      case TrkrDefs::TrkrId::mvtxId:
        fX[n_track::ntrknmaps]++;
        break;
      case TrkrDefs::TrkrId::inttId:
        fX[n_track::ntrknintt]++;
        break;
      case TrkrDefs::TrkrId::tpcId:
        fX[n_track::ntrkntpc]++;
        if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 8)
        {
          fX[n_track::ntrkntpc11]++;
        }
        if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 16)
        {
          fX[n_track::ntrkntpc1]++;
        }
        else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 32)
        {
          fX[n_track::ntrkntpc2]++;
        }
        else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 48)
        {
          fX[n_track::ntrkntpc3]++;
        }
        break;
      case TrkrDefs::TrkrId::micromegasId:
        fX[n_track::ntrknmms]++;
        break;
      default:
        break;
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
      unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
      switch (TrkrDefs::getTrkrId(cluster_key))
      {
      case TrkrDefs::TrkrId::mvtxId:
        fX[n_track::ntrknmaps]++;
        break;
      case TrkrDefs::TrkrId::inttId:
        fX[n_track::ntrknintt]++;
        break;
      case TrkrDefs::TrkrId::tpcId:
        fX[n_track::ntrkntpc]++;
        if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 8)
        {
          fX[n_track::ntrkntpc11]++;
        }
        if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 16)
        {
          fX[n_track::ntrkntpc1]++;
        }
        else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 32)
        {
          fX[n_track::ntrkntpc2]++;
        }
        else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 48)
        {
          fX[n_track::ntrkntpc3]++;
        }
        break;
      case TrkrDefs::TrkrId::micromegasId:
        fX[n_track::ntrknmms]++;
        break;
      default:
        break;
      }
    }
  }
  fX[n_track::ntrkdca3dxy] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkdca3dz] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkdca3dxysigma] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkdca3dzsigma] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkdca2d] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkdca2dsigma] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkvx] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkvy] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkvz] = std::numeric_limits<float>::quiet_NaN();

  int vertexID = track->get_vertex_id();
  fX[n_track::ntrkvertexID] = vertexID;
  
  if (vertexID >= 0&& vertexmap !=nullptr)
  {
    if(vertexmap->size()>100000){
      cout << "too many vtx's" << endl;
    }
    /*
    auto vertexit = vertexmap->find(vertexID);
    if (vertexit != vertexmap->end())
      {
        auto vertex = vertexit->second;

	float vx = vertex->get_x();
	float vy = vertex->get_y();
	float vz = vertex->get_z();
	fX[n_track::ntrkvx] = vx;
	fX[n_track::ntrkvy] = vy;
	fX[n_track::ntrkvz] = vz;
	Acts::Vector3 vert(vx, vy, vz);
	Acts::Vector3 zero = Acts::Vector3::Zero();
	auto dcapair = TrackAnalysisUtils::get_dca(track, zero);
	fX[n_track::ntrkdca3dxy] = dcapair.first.first;
	fX[n_track::ntrkdca3dxysigma] = dcapair.first.second;
	fX[n_track::ntrkdca3dz] = dcapair.second.first;
	fX[n_track::ntrkdca3dzsigma] = dcapair.second.second;
	
      }
    */

  }
  
  fX[n_track::ntrkcharge] = track->get_charge();
  float px = track->get_px();
  float py = track->get_py();
  float pz = track->get_pz();
  fX[n_track::ntrkpx] = px;
  fX[n_track::ntrkpy] = py;
  fX[n_track::ntrkpz] = pz;
  TVector3 v(px, py, pz);
  fX[n_track::ntrkpt] = v.Pt();
  fX[n_track::ntrketa] = v.Eta();
  fX[n_track::ntrkphi] = v.Phi();
  float CVxx = track->get_error(3, 3);
  float CVxy = track->get_error(3, 4);
  float CVxz = track->get_error(3, 5);
  float CVyy = track->get_error(4, 4);
  float CVyz = track->get_error(4, 5);
  float CVzz = track->get_error(5, 5);
  fX[n_track::ntrkdeltapt] = sqrt((CVxx * px * px + 2 * CVxy * px * py + CVyy * py * py) / (px * px + py * py));
  fX[n_track::ntrkdeltaeta] = sqrt((CVzz * (px * px + py * py) * (px * px + py * py) + pz * (-2 * (CVxz * px + CVyz * py) * (px * px + py * py) + CVxx * px * px * pz + CVyy * py * py * pz + 2 * CVxy * px * py * pz)) / ((px * px + py * py) * (px * px + py * py) * (px * px + py * py + pz * pz)));
  fX[n_track::ntrkdeltaphi] = sqrt((CVyy * px * px - 2 * CVxy * px * py + CVxx * py * py) / ((px * px + py * py) * (px * px + py * py)));

  fX[n_track::ntrkpcax] = track->get_x();
  fX[n_track::ntrkpcay] = track->get_y();
  fX[n_track::ntrkpcaz] = track->get_z();
}

float TrkrNtuplizer::calc_dedx(TrackSeed *tpcseed){

    std::vector<TrkrDefs::cluskey> clusterKeys;
    clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
		       tpcseed->end_cluster_keys());

    std::vector<float> dedxlist;
    for (unsigned long cluster_key : clusterKeys){
      unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
      if(TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::TrkrId::tpcId){
	  continue;
      }
      TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);

      float adc = cluster->getAdc();
      PHG4TpcCylinderGeom* GeoLayer_local = _geom_container->GetLayerCellGeom(layer_local);
      float thick = GeoLayer_local->get_thickness();
      
      float r = GeoLayer_local->get_radius();
      float alpha = (r * r) / (2 * r * TMath::Abs(1.0 / tpcseed->get_qOverR()));
      float beta = atan(tpcseed->get_slope());
      float alphacorr = cos(alpha);
      if(alphacorr<0||alphacorr>4){
	alphacorr=4;
      }
      float betacorr = cos(beta);
      if(betacorr<0||betacorr>4){
	betacorr=4;
      }
      adc/=thick;
      adc*=alphacorr;
      adc*=betacorr;
      dedxlist.push_back(adc);
      sort(dedxlist.begin(), dedxlist.end());
    }
    int trunc_min = 0;
    int trunc_max = (int)dedxlist.size()*0.7;
    float sumdedx = 0;
    int ndedx = 0;
    for(int j = trunc_min; j<=trunc_max;j++){
      sumdedx+=dedxlist.at(j);
      ndedx++;
    }
    sumdedx/=ndedx;
    return sumdedx;
}

float TrkrNtuplizer::get_n1pix(TrackSeed *tpcseed){

    std::vector<TrkrDefs::cluskey> clusterKeys;
    clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
		       tpcseed->end_cluster_keys());

    float n1pix = 0;
    for (unsigned long cluster_key : clusterKeys){
      TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);
      if(cluster->getPhiSize()==1){
	n1pix++;
      }
    }
    return n1pix;
}

void TrkrNtuplizer::FillCluster(float fXcluster[49], TrkrDefs::cluskey cluster_key)
{
  unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
  TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);

  Acts::Vector3 cglob;
  cglob = _tgeometry->getGlobalPosition(cluster_key, cluster);
  float x = cglob(0);
  float y = cglob(1);
  float z = cglob(2);
  TVector3 pos(x, y, z);
  float r = pos.Perp();
  float phi = pos.Phi();
  auto para_errors = _ClusErrPara.get_clusterv5_modified_error(cluster, r, cluster_key);
  float phibin = std::numeric_limits<float>::quiet_NaN();
  float tbin = std::numeric_limits<float>::quiet_NaN();
  float locx = cluster->getLocalX();
  float locy = cluster->getLocalY();

  if (layer_local >= _nlayers_maps + _nlayers_intt && layer_local < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
  {
    int side_tpc = TpcDefs::getSide(cluster_key);
    PHG4TpcCylinderGeom* GeoLayer_local = _geom_container->GetLayerCellGeom(layer_local);
    phibin = GeoLayer_local->get_pad_float(phi, side_tpc);
    tbin = GeoLayer_local->get_tbin_float(locy - 39.6);
  }
  else
  {
    phibin = locx;
    tbin = locy;
  }
  fXcluster[n_cluster::nclulocx] = locx;
  fXcluster[n_cluster::nclulocy] = locy;
  fXcluster[n_cluster::nclux] = cglob(0);
  fXcluster[n_cluster::ncluy] = cglob(1);
  fXcluster[n_cluster::ncluz] = cglob(2);
  fXcluster[n_cluster::nclur] = pos.Perp();
  fXcluster[n_cluster::ncluphi] = pos.Phi();
  fXcluster[n_cluster::nclueta] = pos.Eta();
  fXcluster[n_cluster::nclutheta] = pos.Theta();

  fXcluster[n_cluster::ncluphibin] = phibin;
  fXcluster[n_cluster::nclutbin] = tbin;

  fXcluster[n_cluster::ncluex] = std::numeric_limits<float>::quiet_NaN();
  fXcluster[n_cluster::ncluey] = std::numeric_limits<float>::quiet_NaN();
  fXcluster[n_cluster::ncluez] = sqrt(para_errors.second);
  fXcluster[n_cluster::ncluephi] = sqrt(para_errors.first);
  fXcluster[n_cluster::nclupez] = std::numeric_limits<float>::quiet_NaN();
  fXcluster[n_cluster::nclupephi] = std::numeric_limits<float>::quiet_NaN();
  fXcluster[n_cluster::nclue] = cluster->getAdc();
  fXcluster[n_cluster::ncluadc] = cluster->getAdc();
  fXcluster[n_cluster::nclumaxadc] = cluster->getMaxAdc();
  fXcluster[n_cluster::nclulayer] = layer_local;

  if (layer_local < 3)
  {
    fXcluster[n_cluster::ncluphielem] = MvtxDefs::getStaveId(cluster_key);
    fXcluster[n_cluster::ncluzelem] = MvtxDefs::getChipId(cluster_key);
  }
  if (layer_local >= 3 && layer_local < 7)
  {
    fXcluster[n_cluster::ncluphielem] = InttDefs::getLadderPhiId(cluster_key);
    fXcluster[n_cluster::ncluzelem] = InttDefs::getLadderZId(cluster_key);
  }
  if (layer_local >= 7 && layer_local < 55)
  {
    fXcluster[n_cluster::ncluphielem] = TpcDefs::getSectorId(cluster_key);
    fXcluster[n_cluster::ncluzelem] = TpcDefs::getSide(cluster_key);
  } /*
      if(layer_local>=55){
      if(MicromegasDefs::getSegmentationType(hitsetkey)==MicromegasDefs::SEGMENTATION_Z){
      sector = 1;
      side = MicromegasDefs::getTileId(hitsetkey);
      }else{
      sector =MicromegasDefs::getTileId(hitsetkey);
      side =  1;
      }
      }
    */
  fXcluster[n_cluster::nclusize] = cluster->getSize();
  fXcluster[n_cluster::ncluphisize] = cluster->getPhiSize();
  fXcluster[n_cluster::ncluzsize] = cluster->getZSize();
  fXcluster[n_cluster::nclupedge] = cluster->getEdge();
  if (layer_local == 7 || layer_local == 22 ||
      layer_local == 23 || layer_local == 28 ||
      layer_local == 39 || layer_local == 54)
  {
    fXcluster[n_cluster::ncluredge] = 1;
  }

  fXcluster[n_cluster::ncluovlp] = 3;  // cluster->getOvlp();
  fXcluster[n_cluster::nclutrackID] = std::numeric_limits<float>::quiet_NaN();
  fXcluster[n_cluster::ncluniter] = 0;

  return;
}

std::vector<TrkrDefs::cluskey> TrkrNtuplizer::get_track_ckeys(SvtxTrack* track)
{
  std::vector<TrkrDefs::cluskey> cluster_keys;
  TrackSeed* tpcseed = track->get_tpc_seed();
  TrackSeed* silseed = track->get_silicon_seed();
  if (silseed)
  {
    for (auto iter = silseed->begin_cluster_keys();
         iter != silseed->end_cluster_keys();
         ++iter)
    {
      cluster_keys.push_back(*iter);
    }
  }
  if (tpcseed)
  {
    for (auto iter = tpcseed->begin_cluster_keys();
         iter != tpcseed->end_cluster_keys();
         ++iter)
    {
      cluster_keys.push_back(*iter);
    }
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
  for (auto candidate : tracks)
  {
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

  if (_cache_track_from_cluster_exists == false)
  {
    create_cache_track_from_cluster();
  }
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
  for (auto& iter : *_trackmap)
  {
    SvtxTrack* track = iter.second;
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
  if (!_trackmap)
  {
    return;
  }

  // loop over all SvtxTracks
  for (auto& iter : *_trackmap)
  {
    SvtxTrack* track = iter.second;
    std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);

    // loop over all clusters
    for (const auto& candidate_key : cluster_keys)
    {
      // check if cluster has an entry in cache
      std::map<TrkrDefs::cluskey, std::set<SvtxTrack*> >::iterator cliter =
          _cache_all_tracks_from_cluster.find(candidate_key);
      if (cliter != _cache_all_tracks_from_cluster.end())
      {                                // got entry
        cliter->second.insert(track);  // add track to list;
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
