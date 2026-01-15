#include "dEdxFitter.h"

#include <phool/phool.h>
#include <phool/getClass.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <g4detectors/PHG4TpcGeom.h>

#include <iostream>

//____________________________________
dEdxFitter::dEdxFitter(const std::string &name):
    SubsysReco(name)
{ 
  //initialize
  fitter = std::make_unique<GlobaldEdxFitter>();
}

//___________________________________
int dEdxFitter::InitRun(PHCompositeNode *topNode)
{
  std::cout << PHWHERE << " Opening file " << _outfile << std::endl;
  PHTFileServer::get().open( _outfile, "RECREATE");

  return 0;
}

//__________________________________
//Call user instructions for every event
int dEdxFitter::process_event(PHCompositeNode *topNode)
{
  _event++;
  if(_event%1000==0) std::cout << PHWHERE << "Events processed: " << _event << std::endl;
 
  GetNodes(topNode);

  if(Verbosity()>1)
  {
    std::cout << "--------------------------------" << std::endl;
    std::cout << "event " << _event << std::endl;
  }

  process_tracks(topNode);

  return 0;
}

//_____________________________________
void dEdxFitter::process_tracks(PHCompositeNode *topNode)
{

  for(const auto &[key, track] : *_trackmap)
  {
    if(!track) continue;

    double trackID = track->get_id();
    if(Verbosity()>1) std::cout << "track ID " << trackID << std::endl;
    if(std::isnan(track->get_x()) ||
       std::isnan(track->get_y()) ||
       std::isnan(track->get_z()) ||
       std::isnan(track->get_px()) ||
       std::isnan(track->get_py()) ||
       std::isnan(track->get_pz()))
    {
      std::cout << "malformed track:" << std::endl;
      track->identify();
      std::cout << "skipping..." << std::endl;
      continue;
    }

    // ignore TPC-only tracks
    if(!track->get_silicon_seed())
    {
      if(Verbosity()>1) std::cout << "TPC-only track, skipping..." << std::endl;
      continue;
    }

    std::tuple<int,int,int> nclus = get_nclus(track);
    int nmaps = std::get<0>(nclus);
    int nintt = std::get<1>(nclus);
    int ntpc = std::get<2>(nclus);

    if(nmaps>=nmaps_cut && nintt>=nintt_cut && ntpc>=ntpc_cut && fabs(track->get_eta())<eta_cut && get_dcaxy(track)<dcaxy_cut)
    {
      fitter->addTrack(get_dedx(track),track->get_p());
    }

    if(fitter->getNtracks() > ntracks_to_fit)
    {
      minima.push_back(fitter->get_minimum());
      fitter->reset();
    }
  }
}

std::tuple<int,int,int> dEdxFitter::get_nclus(SvtxTrack* track)
{
  int nmaps = 0;
  int nintt = 0;
  int ntpc = 0;

  for(auto it = track->get_silicon_seed()->begin_cluster_keys(); it != track->get_silicon_seed()->end_cluster_keys(); ++it)
  {
    TrkrDefs::cluskey ckey = *it;
    auto trkrid = TrkrDefs::getTrkrId(ckey);
    if(trkrid == TrkrDefs::mvtxId)
    {
      nmaps++;
    }
    else if(trkrid == TrkrDefs::inttId)
    {
      nintt++;
    }
  }
  for(auto it = track->get_tpc_seed()->begin_cluster_keys(); it != track->get_tpc_seed()->end_cluster_keys(); ++it)
  {
    ntpc++;
  }

  return std::make_tuple(nmaps,nintt,ntpc);
}

double dEdxFitter::get_dedx(SvtxTrack* track)
{
  float layerThicknesses[4] = {0.0, 0.0, 0.0, 0.0};
  // These are randomly chosen layer thicknesses for the TPC, to get the
  // correct region thicknesses in an easy to pass way to the helper fxn
  layerThicknesses[0] = _tpcgeom->GetLayerCellGeom(7)->get_thickness();
  layerThicknesses[1] = _tpcgeom->GetLayerCellGeom(8)->get_thickness();
  layerThicknesses[2] = _tpcgeom->GetLayerCellGeom(27)->get_thickness();
  layerThicknesses[3] = _tpcgeom->GetLayerCellGeom(50)->get_thickness();

  return TrackAnalysisUtils::calc_dedx(track->get_tpc_seed(), _clustermap, _geometry, layerThicknesses);
}

double dEdxFitter::get_dcaxy(SvtxTrack* track)
{
  auto vertexit = _vertexmap->find(track->get_vertex_id());
  if(vertexit != _vertexmap->end())
  {
    SvtxVertex* vtx = vertexit->second;
    Acts::Vector3 vertex(vtx->get_x(),vtx->get_y(),vtx->get_z());
    auto dcapair = TrackAnalysisUtils::get_dca(track,vertex);
    return dcapair.first.first;
  }
  else
  {
    return std::numeric_limits<float>::quiet_NaN();
  }
}

//___________________________________
void dEdxFitter::GetNodes(PHCompositeNode *topNode)
{

  _trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if(!_trackmap && _event<2)
  {
    std::cout << PHWHERE << " cannot find SvtxTrackMap" << std::endl;
  }

  _clustermap = findNode::getClass<TrkrClusterContainer>(topNode,"TRKR_CLUSTER");
  if(!_clustermap && _event<2)
  {
    std::cout << PHWHERE << " cannot find TrkrClusterContainer TRKR_CLUSTER" << std::endl;
  }

  _geometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!_geometry && _event<2)
  {
    std::cout << PHWHERE << " cannot find ActsGeometry" << std::endl;
  }

  _tpcgeom = findNode::getClass<PHG4TpcGeomContainer>(topNode,"TPCGEOMCONTAINER");
  if(!_tpcgeom && _event<2)
  {
    std::cout << PHWHERE << " cannot find PHG4TpcGeomContainer TPCGEOMCONTAINER" << std::endl;
  }

  _vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  if(!_vertexmap && _event<2)
  {
    std::cout << PHWHERE << " cannot find SvtxVertexMap" << std::endl;
  }
}

//______________________________________
int dEdxFitter::End(PHCompositeNode *topNode)
{
  if(minima.size()==0)
  {
    minima.push_back(fitter->get_minimum());
  }

  PHTFileServer::get().cd( _outfile );

  double avg_minimum = 0.;
  for(double m : minima)
  {
    avg_minimum += m;
  }
  avg_minimum /= (double)minima.size();

  TF1* pi_band = new TF1("pi_band","bethe_bloch_new_1D(fabs(x)/[1],[0])");
  pi_band->SetParameter(0,avg_minimum);
  pi_band->SetParameter(1,dedx_constants::m_pi);
  pi_band->Write();

  TF1* K_band = new TF1("K_band","bethe_bloch_new_1D(fabs(x)/[1],[0])");
  K_band->SetParameter(0,avg_minimum);
  K_band->SetParameter(1,dedx_constants::m_K);
  K_band->Write();

  TF1* p_band = new TF1("p_band","bethe_bloch_new_1D(fabs(x)/[1],[0])");
  p_band->SetParameter(0,avg_minimum);
  p_band->SetParameter(1,dedx_constants::m_p);
  p_band->Write();

  TF1* d_band = new TF1("d_band","bethe_bloch_new_1D(fabs(x)/[1],[0])");
  d_band->SetParameter(0,avg_minimum);
  d_band->SetParameter(1,dedx_constants::m_d);
  d_band->Write();

  if(Verbosity()>0) std::cout << "dEdxFitter extracted minimum: " << avg_minimum << std::endl;

  return 0;
}
