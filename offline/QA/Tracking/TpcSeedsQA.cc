#include "TpcSeedsQA.h"

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <tpc/TpcDistortionCorrectionContainer.h>
#include <tpc/TpcGlobalPositionWrapper.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1RawHit.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TH2.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include <boost/format.hpp>
#include <cmath>

//____________________________________________________________________________..
TpcSeedsQA::TpcSeedsQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int TpcSeedsQA::InitRun(PHCompositeNode *topNode)
{
  createHistos();

  clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  actsgeom = findNode::getClass<ActsGeometry>(topNode, m_actsGeomName);
  g4geom = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, m_g4GeomName);
  trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);

  if (!trackmap or !clustermap or !actsgeom or !vertexmap)
  {
    std::cout << PHWHERE << "Missing node(s), can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (!g4geom)
  {
    std::cout << PHWHERE << " unable to find DST node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // global position wrapper
  m_globalPositionWrapper.loadNodes(topNode);

  m_clusterMover.initialize_geometry(g4geom);
  m_clusterMover.set_verbosity(0);

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // tracks with TPC clusters/tracklets
  h_ntrack1d = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d").c_str()));
  h_ntrack1d_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d_pos").c_str()));
  h_ntrack1d_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d_neg").c_str()));
  h_ntrack1d_ptg1 = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d_ptg1").c_str()));
  h_ntrack1d_ptg1_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d_ptg1_pos").c_str()));
  h_ntrack1d_ptg1_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d_ptg1_neg").c_str()));
  h_pt = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "pt").c_str()));
  h_pt_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "pt_pos").c_str()));
  h_pt_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "pt_neg").c_str()));
  h_ntrack_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks_pos").c_str()));
  h_ntrack_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks_neg").c_str()));

  h_ntpc_fullpt_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_fullpt_pos").c_str()));
  h_ntpc_fullpt_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_fullpt_neg").c_str()));
  h_ntpc_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_pos").c_str()));
  h_ntpc_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_neg").c_str()));
  h_ntpc_quality_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_quality_pos").c_str()));
  h_ntpc_quality_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_quality_neg").c_str()));
  h_ntpot_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpot_pos").c_str()));
  h_ntpot_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpot_neg").c_str()));
  h_avgnclus_eta_phi_pos = dynamic_cast<TProfile2D *>(hm->getHisto(std::string(getHistoPrefix() + "avgnclus_eta_phi_pos").c_str()));
  h_avgnclus_eta_phi_neg = dynamic_cast<TProfile2D *>(hm->getHisto(std::string(getHistoPrefix() + "avgnclus_eta_phi_neg").c_str()));
  // h_trackcrossing_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "trackcrossing_pos").c_str()));
  // h_trackcrossing_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "trackcrossing_neg").c_str()));
  h_dcaxyorigin_phi_north_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcaxyorigin_phi_north_pos").c_str()));
  h_dcaxyorigin_phi_south_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcaxyorigin_phi_south_pos").c_str()));
  h_dcaxyorigin_phi_north_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcaxyorigin_phi_north_neg").c_str()));
  h_dcaxyorigin_phi_south_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcaxyorigin_phi_south_neg").c_str()));
  h_dcaxyvtx_phi_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcaxyvtx_phi_pos").c_str()));
  h_dcaxyvtx_phi_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcaxyvtx_phi_neg").c_str()));
  h_dcazorigin_phi_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcazorigin_phi_pos").c_str()));
  h_dcazorigin_phi_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcazorigin_phi_neg").c_str()));
  h_dcazvtx_phi_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcazvtx_phi_pos").c_str()));
  h_dcazvtx_phi_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcazvtx_phi_neg").c_str()));
  h_ntrack_isfromvtx_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntrack_isfromvtx_pos").c_str()));
  h_ntrack_isfromvtx_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntrack_isfromvtx_neg").c_str()));
  h_cluster_phisize1_fraction_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "cluster_phisize1_fraction_pos").c_str()));
  h_cluster_phisize1_fraction_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "cluster_phisize1_fraction_neg").c_str()));

  // vertex
  h_nvertex = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecovertices").c_str()));
  h_vx = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vx").c_str()));
  h_vy = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vy").c_str()));
  h_vx_vy = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "vx_vy").c_str()));
  h_vz = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vz").c_str()));
  h_vt = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vt").c_str()));
  // h_vcrossing = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vertexcrossing").c_str()));
  h_vchi2dof = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vertexchi2dof").c_str()));
  h_ntrackpervertex = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntrackspervertex").c_str()));
  h_dedx = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dedx").c_str()));
  h_mip_dedx = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "mip_dedx").c_str()));
  h_dedx_pcaz = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dedx_pcaz").c_str()));

  nt_sector_event_summary = dynamic_cast<TNtuple *>(hm->getHisto(std::string(getHistoPrefix() + "sector_event_summary").c_str()));

  // TPC has 3 regions, inner, mid and outer
  std::vector<int> region_layer_low = {7, 23, 39};
  std::vector<int> region_layer_high = {22, 38, 54};

  // make a layer to region multimap
  const auto range = g4geom->get_begin_end();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    m_layers.insert(iter->first);

    for (int region = 0; region < 3; ++region)
    {
      if (iter->first >= region_layer_low[region] && iter->first <= region_layer_high[region])
      {
        m_layerRegionMap.insert(std::make_pair(iter->first, region));
      }
    }
  }

  for (auto &region : {0, 1, 2})
  {
    PhiHistoList phihist;

    phihist.cphisize1pT_side0 = h_clusphisize1pt_side0[region];
    phihist.cphisize1pT_side1 = h_clusphisize1pt_side1[region];

    phihist.cphisizegeq1pT_side0 = h_clusphisizegeq1pt_side0[region];
    phihist.cphisizegeq1pT_side1 = h_clusphisizegeq1pt_side1[region];

    phihist.Clear();

    phihistos.insert(std::make_pair(region, phihist));
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

float TpcSeedsQA::calc_dedx(TrackSeed *tpcseed)
{
  std::vector<TrkrDefs::cluskey> clusterKeys;
  clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
                     tpcseed->end_cluster_keys());

  std::vector<float> dedxlist;
  for (unsigned long cluster_key : clusterKeys)
  {
    auto detid = TrkrDefs::getTrkrId(cluster_key);
    if (detid != TrkrDefs::TrkrId::tpcId)
    {
      continue;  // the micromegas clusters are added to the TPC seeds
    }
    unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
    TrkrCluster *cluster = clustermap->findCluster(cluster_key);
    float adc = cluster->getAdc();
    PHG4TpcCylinderGeom *GeoLayer_local = g4geom->GetLayerCellGeom(layer_local);
    float thick = GeoLayer_local->get_thickness();
    float r = GeoLayer_local->get_radius();
    float alpha = (r * r) / (2 * r * TMath::Abs(1.0 / tpcseed->get_qOverR()));
    float beta = std::atan(tpcseed->get_slope());
    float alphacorr = std::cos(alpha);
    if (alphacorr < 0 || alphacorr > 4)
    {
      alphacorr = 4;
    }
    float betacorr = std::cos(beta);
    if (betacorr < 0 || betacorr > 4)
    {
      betacorr = 4;
    }
    adc /= thick;
    adc *= alphacorr;
    adc *= betacorr;
    dedxlist.push_back(adc);
    sort(dedxlist.begin(), dedxlist.end());
  }
  int trunc_min = 0;
  if (dedxlist.size() < 1)
  {
    return std::numeric_limits<float>::quiet_NaN();
  }
  int trunc_max = (int) dedxlist.size() * 0.7;
  float sumdedx = 0;
  int ndedx = 0;
  for (int j = trunc_min; j <= trunc_max; j++)
  {
    sumdedx += dedxlist.at(j);
    ndedx++;
  }
  sumdedx /= ndedx;
  return sumdedx;
}

float TpcSeedsQA::cal_track_length(SvtxTrack *track)
{
  float minR = std::numeric_limits<float>::max();
  float maxR = 0;
  for (const auto& ckey : get_cluster_keys(track))
  {
    auto cluster = clustermap->findCluster(ckey);

    // Fully correct the cluster positions for the crossing and all distortions
    Acts::Vector3 global = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cluster, track->get_crossing() );

    // add the global positions to a vector to give to the cluster mover
    float R = std::sqrt(pow(global.x(),2) + pow(global.y(),2));
    if (R < minR)
    {
      minR = R;
    }
    if (R > maxR)
    {
      maxR = R;
    }
  }
  float tracklength = maxR - minR;
  return tracklength;
}

float* TpcSeedsQA::cal_dedx_cluster(SvtxTrack *track)
{
  // get the fully corrected cluster global positions
  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;
  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_moved;
  float minR = std::numeric_limits<float>::max();
  float maxR = 0;
  for (const auto& ckey : get_cluster_keys(track))
  {
    auto cluster = clustermap->findCluster(ckey);

    // Fully correct the cluster positions for the crossing and all distortions
    Acts::Vector3 global = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cluster, track->get_crossing() );

    // add the global positions to a vector to give to the cluster mover
    global_raw.emplace_back(std::make_pair(ckey, global));
    float R = std::sqrt(pow(global.x(),2) + pow(global.y(),2));
    if (R < minR)
    {
      minR = R;
    }
    if (R > maxR)
    {
      maxR = R;
    }
  }
  float tracklength = maxR - minR;
  if (collision_or_cosmics == true && tracklength < 25)
  {
    float* dedxarray = new float[10];
    for (int i = 0; i < 10; ++i)
    {
        dedxarray[i] = -1;
    }
    return dedxarray;
  }

  // move the corrected cluster positions back to the original readout surface
  global_moved = m_clusterMover.processTrack(global_raw);

  float fcorr = std::fabs(std::sin(eta_to_theta(track->get_eta())));
  Acts::Vector3 clusglob_moved(0, 0, 0);
  float adc_z0=0; int nclus_z0=0;
  float adc_z1=0; int nclus_z1=0;
  float adc_z2=0; int nclus_z2=0;
  float adc_z3=0; int nclus_z3=0;
  float adc_z4=0; int nclus_z4=0;
  float adc_z5=0; int nclus_z5=0;
  float adc_z6=0; int nclus_z6=0;
  float adc_z7=0; int nclus_z7=0;
  float adc_z8=0; int nclus_z8=0;
  float adc_z9=0; int nclus_z9=0;
  for (const auto& pair : global_moved)
  {
    auto ckey = pair.first;
    auto cluster = clustermap->findCluster(ckey);
    clusglob_moved = pair.second;

    auto detid = TrkrDefs::getTrkrId(ckey);
    if (detid != TrkrDefs::TrkrId::tpcId)
    {
      continue;  // the micromegas clusters are added to the TPC seeds
    }

    // only counts TPC R2 and R3
    auto layer = TrkrDefs::getLayer(ckey);
    if (layer<23)
    {
      continue;
    }

    float clusgz = clusglob_moved.z();
    if (clusgz>-100 && clusgz<=-80)
    {
      adc_z0+=cluster->getAdc() * fcorr;
      nclus_z0++;
    }
    else if (clusgz>-80 && clusgz<=-60)
    {
      adc_z1+=cluster->getAdc() * fcorr;
      nclus_z1++;
    }
    else if (clusgz>-60 && clusgz<=-40)
    {
      adc_z2+=cluster->getAdc() * fcorr;
      nclus_z2++;
    }
    else if (clusgz>-40 && clusgz<=-20)
    {
      adc_z3+=cluster->getAdc() * fcorr;
      nclus_z3++;
    }
    else if (clusgz>-20 && clusgz<=0)
    {
      adc_z4+=cluster->getAdc() * fcorr;
      nclus_z4++;
    }
    else if (clusgz>0 && clusgz<=20)
    {
      adc_z5+=cluster->getAdc() * fcorr;
      nclus_z5++;
    }
    else if (clusgz>20 && clusgz<=40)
    {
      adc_z6+=cluster->getAdc() * fcorr;
      nclus_z6++;
    }
    else if (clusgz>40 && clusgz<=60)
    {
      adc_z7+=cluster->getAdc() * fcorr;
      nclus_z7++;
    }
    else if (clusgz>60 && clusgz<=80)
    {
      adc_z8+=cluster->getAdc() * fcorr;
      nclus_z8++;
    }
    else if (clusgz>80 && clusgz<=100)
    {
      adc_z9+=cluster->getAdc() * fcorr;
      nclus_z9++;
    }
  }

  adc_z0 /= nclus_z0;
  adc_z1 /= nclus_z1;
  adc_z2 /= nclus_z2;
  adc_z3 /= nclus_z3;
  adc_z4 /= nclus_z4;
  adc_z5 /= nclus_z5;
  adc_z6 /= nclus_z6;
  adc_z7 /= nclus_z7;
  adc_z8 /= nclus_z8;
  adc_z9 /= nclus_z9;

  float* dedxarray = new float[10];
  dedxarray[0] = adc_z0;
  dedxarray[1] = adc_z1;
  dedxarray[2] = adc_z2;
  dedxarray[3] = adc_z3;
  dedxarray[4] = adc_z4;
  dedxarray[5] = adc_z5;
  dedxarray[6] = adc_z6;
  dedxarray[7] = adc_z7;
  dedxarray[8] = adc_z8;
  dedxarray[9] = adc_z9;

  return dedxarray;
}

//____________________________________________________________________________..
int TpcSeedsQA::process_event(PHCompositeNode *topNode)
{
  auto gl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
  if (gl1)
  {
    m_bco = gl1->get_bco();
  }
  else
  {
    Gl1Packet *gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
    if (!gl1PacketInfo)
    {
      m_bco = std::numeric_limits<uint64_t>::quiet_NaN();
    }
    if (gl1PacketInfo)
    {
      m_bco = gl1PacketInfo->getBCO();
    }
  }

  h_ntrack1d->Fill(trackmap->size());

  std::pair<int, int> ntrack_isfromvtx_pos;  // first: number of tracks not associated to a vertex, second: number of tracks associated to a vertex
  std::pair<int, int> ntrack_isfromvtx_neg;  // first: number of tracks not associated to a vertex, second: number of tracks associated to a vertex

  int ntrack1d_pos = 0;
  int ntrack1d_neg = 0;
  int ntrack1d_ptg1_pos = 0;
  int ntrack1d_ptg1_neg = 0;

  int nclus[2][3][12] = {{{0}}};
  int madc[2][3][12] = {{{0}}};
  for (const auto &[key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }

    int charge = track->get_charge();
    float quality = track->get_quality();
    float pt = track->get_pt();

    h_pt->Fill(pt);
    if (charge == 1)
    {
      ntrack1d_pos++;
      if (pt > 1)
      {
        ntrack1d_ptg1_pos++;
      }
      h_pt_pos->Fill(pt);
    }
    else if (charge == -1)
    {
      ntrack1d_neg++;
      if (pt > 1)
      {
        ntrack1d_ptg1_neg++;
      }
      h_pt_neg->Fill(pt);
    }

    auto ckeys = get_cluster_keys(track);
    std::vector<Acts::Vector3> cluspos;
    TrackFitUtils::getTrackletClusters(actsgeom, clustermap, cluspos, ckeys);
    float eta = track->get_eta();
    float phi = track->get_phi();

    // int trkcrossing = track->get_crossing();

//    int nmaps = 0;
//    int nintt = 0;
    int ntpc = 0;
    int ntpc_phisize1 = 0;
    int nmms = 0;

    for (auto &ckey : ckeys)
    {
      if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::tpcId)
      {
        TrkrCluster *cluster = clustermap->findCluster(ckey);
        if (cluster->getPhiSize() == 1)
        {
          ntpc_phisize1++;
        }
      }
      switch (TrkrDefs::getTrkrId(ckey))
      {
      // case TrkrDefs::mvtxId:
      //   nmaps++;
      //   break;
      // case TrkrDefs::inttId:
      //   nintt++;
      //   break;
      case TrkrDefs::tpcId:
        ntpc++;
        break;
      case TrkrDefs::micromegasId:
        nmms++;
        break;
      default:
	break;
      }
    }

    Acts::Vector3 zero = Acts::Vector3::Zero();
    auto dcapair_origin = TrackAnalysisUtils::get_dca(track, zero);

    auto trackvtx = vertexmap->get(track->get_vertex_id());
    if (!trackvtx)
    {
      // std::cout << "No vertex found for track " << track->get_id() << std::std::endl;
      if (charge == 1)
      {
        ntrack_isfromvtx_pos.first++;
      }
      else if (charge == -1)
      {
        ntrack_isfromvtx_neg.first++;
      }
    }
    else
    {
      Acts::Vector3 track_vtx(trackvtx->get_x(), trackvtx->get_y(), trackvtx->get_z());
      auto dcapair_vtx = TrackAnalysisUtils::get_dca(track, track_vtx);
      if (charge == 1)
      {
        ntrack_isfromvtx_pos.second++;
        h_dcaxyvtx_phi_pos->Fill(phi, dcapair_vtx.first.first);
        h_dcazvtx_phi_pos->Fill(phi, dcapair_vtx.second.first);
      }
      else if (charge == -1)
      {
        ntrack_isfromvtx_neg.second++;
        h_dcaxyvtx_phi_neg->Fill(phi, dcapair_vtx.first.first);
        h_dcazvtx_phi_neg->Fill(phi, dcapair_vtx.second.first);
      }
    }

    if (charge == 1)
    {
      h_ntpc_fullpt_pos->Fill(ntpc);
      if (dcapair_origin.second.first > 0)
      {
        h_dcaxyorigin_phi_north_pos->Fill(phi, dcapair_origin.first.first);
      }
      else if (dcapair_origin.second.first <= 0)
      {
        h_dcaxyorigin_phi_south_pos->Fill(phi, dcapair_origin.first.first);
      }
      h_dcazorigin_phi_pos->Fill(phi, dcapair_origin.second.first);
      if (pt > 1)
      {
        h_ntrack_pos->Fill(eta, phi);
        if (trackvtx)
        {
          float vz = trackvtx->get_z();
          float eta_min = cal_tpc_eta_min_max(vz).first;
          float eta_max = cal_tpc_eta_min_max(vz).second;
          if (eta > eta_min && eta < eta_max)
          {
            h_ntpc_pos->Fill(ntpc);
          }
        }
        h_ntpot_pos->Fill(nmms);
        h_ntpc_quality_pos->Fill(ntpc, quality);
        h_avgnclus_eta_phi_pos->Fill(eta, phi, ntpc);
        // h_trackcrossing_pos->Fill(trkcrossing);
        h_cluster_phisize1_fraction_pos->Fill((double) ntpc_phisize1 / (double) ntpc);
      }
    }
    else if (charge == -1)
    {
      h_ntpc_fullpt_neg->Fill(ntpc);
      if (dcapair_origin.second.first > 0)
      {
        h_dcaxyorigin_phi_north_neg->Fill(phi, dcapair_origin.first.first);
      }
      else if (dcapair_origin.second.first <= 0)
      {
        h_dcaxyorigin_phi_south_neg->Fill(phi, dcapair_origin.first.first);
      }
      h_dcazorigin_phi_neg->Fill(phi, dcapair_origin.second.first);
      if (pt > 1)
      {
        h_ntrack_neg->Fill(eta, phi);
        if (trackvtx)
        {
          float vz = trackvtx->get_z();
          float eta_min = cal_tpc_eta_min_max(vz).first;
          float eta_max = cal_tpc_eta_min_max(vz).second;
          if (eta > eta_min && eta < eta_max)
          {
            h_ntpc_neg->Fill(ntpc);
          }
        }
        h_ntpot_neg->Fill(nmms);
        h_ntpc_quality_neg->Fill(ntpc, quality);
        h_avgnclus_eta_phi_neg->Fill(eta, phi, ntpc);
        // h_trackcrossing_neg->Fill(trkcrossing);
        h_cluster_phisize1_fraction_neg->Fill((double) ntpc_phisize1 / (double) ntpc);
      }
    }
  }
  h_ntrack1d_pos->Fill(ntrack1d_pos);
  h_ntrack1d_neg->Fill(ntrack1d_neg);
  h_ntrack1d_ptg1_pos->Fill(ntrack1d_ptg1_pos);
  h_ntrack1d_ptg1_neg->Fill(ntrack1d_ptg1_neg);
  h_ntrack1d_ptg1->Fill(ntrack1d_ptg1_pos + ntrack1d_ptg1_neg);

  h_ntrack_isfromvtx_pos->SetBinContent(1, h_ntrack_isfromvtx_pos->GetBinContent(1) + ntrack_isfromvtx_pos.first);
  h_ntrack_isfromvtx_pos->SetBinContent(2, h_ntrack_isfromvtx_pos->GetBinContent(2) + ntrack_isfromvtx_pos.second);
  h_ntrack_isfromvtx_neg->SetBinContent(1, h_ntrack_isfromvtx_neg->GetBinContent(1) + ntrack_isfromvtx_neg.first);
  h_ntrack_isfromvtx_neg->SetBinContent(2, h_ntrack_isfromvtx_neg->GetBinContent(2) + ntrack_isfromvtx_neg.second);

  // vertex
  h_nvertex->Fill(vertexmap->size());
  for (const auto &[key, vertex] : *vertexmap)
  {
    if (!vertex)
    {
      continue;
    }

    float vx = vertex->get_x();
    float vy = vertex->get_y();
    float vz = vertex->get_z();
    float vt = vertex->get_t0();
    float vchi2 = vertex->get_chisq();
    int vndof = vertex->get_ndof();
    // int vcrossing = vertex->get_beam_crossing();

    // std::cout << "vertex (x,y,z,t,chi2,ndof,crossing)=(" << vx << "," << vy << "," << vz << "," << vt << "," << vchi2 << "," << vndof << "," << vcrossing << ")" << std::endl;

    h_vx->Fill(vx);
    h_vy->Fill(vy);
    h_vx_vy->Fill(vx, vy);
    h_vz->Fill(vz);
    h_vt->Fill(vt);
    h_vchi2dof->Fill(float(vchi2 / vndof));
    // h_vcrossing->Fill(vcrossing);

    h_ntrackpervertex->Fill(vertex->size_tracks());
  }

  auto fill = [](TH1 *h, float val)
  { if (h) { h->Fill(val); } };

  std::set<unsigned int> tpc_seed_ids;
  for (const auto &[key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }
    m_px = track->get_px();
    m_py = track->get_py();
    m_pz = track->get_pz();
    m_pt = std::sqrt(m_px * m_px + m_py * m_py);
    m_ptot = std::sqrt(m_px * m_px + m_py * m_py + m_pz * m_pz);
    TrackSeed *tpcseed = track->get_tpc_seed();
    m_charge = track->get_charge();
    m_dedx = calc_dedx(tpcseed);

    m_ntpc = 0;
    m_region.clear();
    m_clusgz.clear();
    m_cluslayer.clear();
    m_clusphisize.clear();
    m_cluszsize.clear();

    for (const auto &ckey : get_cluster_keys(track))
    {
      if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::tpcId)
      {
        m_ntpc++;
      }
    }

    for (const auto &ckey : get_cluster_keys(track))
    {
      TrkrCluster *cluster = clustermap->findCluster(ckey);
      const Acts::Vector3 clusglob = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cluster, track->get_crossing());

      switch (TrkrDefs::getTrkrId(ckey))
      {
      case TrkrDefs::tpcId:
        const auto it = m_layerRegionMap.find(TrkrDefs::getLayer(ckey));
        int region = it->second;
        m_region.push_back(region);
        m_clusgz.push_back(clusglob.z());
        m_cluslayer.push_back(TrkrDefs::getLayer(ckey));
        m_clusphisize.push_back(cluster->getPhiSize());
        m_cluszsize.push_back(cluster->getZSize());
        int this_sector = (int) TpcDefs::getSectorId(ckey);
        int this_side = (int) TpcDefs::getSide(ckey);
        int is_onepad = 0;
        if (cluster->getPhiSize() <= 1)
        {
          is_onepad = 1;
        }
        if (m_ntpc > 30 && cluster->getPhiSize() > 1 && m_dedx < 1500 && m_pt > 1.0 && m_pt < 50)
        {
          h_adc_sector[region]->Fill((this_sector + 1) * (2 * (this_side - 0.5)), cluster->getAdc());
        }
        if (m_ntpc > 30 && m_dedx < 1000 && m_pt > 1.0 && m_pt < 50 && m_charge < 0)
        {
          h_onepad_frac[region]->Fill((this_sector + 1) * (2 * (this_side - 0.5)), is_onepad);
        }
        nclus[this_side][region][this_sector] += 1;
        madc[this_side][region][this_sector] += cluster->getAdc();
        break;
      }
    }

    if (m_ptot > 0.2 && m_ptot < 4 && m_ntpc > 30)
    {
      h_dedx->Fill(m_charge * m_ptot, m_dedx);
      if (collision_or_cosmics == false || (collision_or_cosmics == true && cal_track_length(track) > 25))
      {
        h_dedx_pcaz->Fill(track->get_z(), m_dedx);
      }
    }

    if (m_ptot > 1.0 && m_pt < 4 && m_ntpc > 30 && m_charge < 0 && m_dedx < 1000 && m_dedx > 50)
    {
      h_mip_dedx->Fill(m_dedx);
    }

    if (m_ntpc > 30)
    {
      float *cluster_dedx = cal_dedx_cluster(track);
      for (int iz = 0; iz < 10; iz++)
      {
        h_dedx_pq_z[iz]->Fill(m_charge * m_ptot, cluster_dedx[iz]);
      }
    }

    // if (m_pt > 1)
    //{
    //   h_ntpc->Fill(m_ntpc);
    // }

    int nClus = m_cluslayer.size();
    for (int cl = 0; cl < nClus; cl++)
    {
      if (m_pt > 1 && m_ntpc > 25)
      {
        if (m_clusphisize[cl] == 1 && m_cluszsize[cl] > 1)
        {
          if (m_clusgz[cl] < 0.)
          {
            const auto hiter = phihistos.find(m_region[cl]);
            if (hiter == phihistos.end())
            {
              continue;
            }
            fill(hiter->second.cphisize1pT_side0, m_pt);
          }
          else if (m_clusgz[cl] > 0.)
          {
            const auto hiter = phihistos.find(m_region[cl]);
            if (hiter == phihistos.end())
            {
              continue;
            }
            fill(hiter->second.cphisize1pT_side1, m_pt);
          }
        }
        if (m_clusphisize[cl] >= 1 && m_cluszsize[cl] > 1)
        {
          if (m_clusgz[cl] < 0.)
          {
            const auto hiter = phihistos.find(m_region[cl]);
            if (hiter == phihistos.end())
            {
              continue;
            }
            fill(hiter->second.cphisizegeq1pT_side0, m_pt);
          }
          else if (m_clusgz[cl] > 0.)
          {
            const auto hiter = phihistos.find(m_region[cl]);
            if (hiter == phihistos.end())
            {
              continue;
            }
            fill(hiter->second.cphisizegeq1pT_side1, m_pt);
          }
        }
      }
    }

    for (auto &pair : phihistos)
    {
      pair.second.Clear();
    }

    for (int cl = 0; cl < nClus; cl++)
    {
      if (m_pt > 1 && m_ntpc > 25)
      {
        if (m_clusgz[cl] < 0.)
        {
          const auto hiter = phihistos.find(m_region[cl]);
          if (hiter == phihistos.end())
          {
            continue;
          }
          if (m_clusphisize[cl] == 1 && m_cluszsize[cl] > 1)
          {
            hiter->second.ntpc_side0_phisize1++;
          }
          if (m_clusphisize[cl] >= 1 && m_cluszsize[cl] > 1)
          {
            hiter->second.ntpc_side0++;
          }
        }
        else if (m_clusgz[cl] > 0.)
        {
          const auto hiter = phihistos.find(m_region[cl]);
          if (hiter == phihistos.end())
          {
            continue;
          }
          if (m_clusphisize[cl] == 1 && m_cluszsize[cl] > 1)
          {
            hiter->second.ntpc_side1_phisize1++;
          }
          if (m_clusphisize[cl] >= 1 && m_cluszsize[cl] > 1)
          {
            hiter->second.ntpc_side1++;
          }
        }
      }
    }

    for (auto &region : {0, 1, 2})
    {
      if (phihistos[region].ntpc_side0 > 0)
      {
        double frac_side0 = (double) phihistos[region].ntpc_side0_phisize1 / (double) phihistos[region].ntpc_side0;
        h_cluster_phisize1_fraction_side0[region]->Fill(frac_side0);
        h_cluster_phisize1_fraction_pt_side0[region]->Fill(m_pt, frac_side0);

        int index_pt_side0 = h_cluster_phisize1_fraction_mean_side0[region]->FindBin(m_pt) - 1;
        if (index_pt_side0 < h_cluster_phisize1_fraction_mean_side0[region]->GetNbinsX())
        {
          frac_side0_pt[region][index_pt_side0] += frac_side0;
          num_track_side0_pt[region][index_pt_side0]++;
        }
      }

      if (phihistos[region].ntpc_side1 > 0)
      {
        double frac_side1 = (double) phihistos[region].ntpc_side1_phisize1 / (double) phihistos[region].ntpc_side1;
        h_cluster_phisize1_fraction_side1[region]->Fill(frac_side1);
        h_cluster_phisize1_fraction_pt_side1[region]->Fill(m_pt, frac_side1);

        int index_pt_side1 = h_cluster_phisize1_fraction_mean_side1[region]->FindBin(m_pt) - 1;
        if (index_pt_side1 < h_cluster_phisize1_fraction_mean_side1[region]->GetNbinsX())
        {
          frac_side1_pt[region][index_pt_side1] += frac_side1;
          num_track_side1_pt[region][index_pt_side1]++;
        }
      }
    }
  }

  for (int iside = 0; iside < 2; iside++)
  {
    for (int iregion = 0; iregion < 3; iregion++)
    {
      for (int isector = 0; isector < 12; isector++)
      {
        //"event:segment:bco:side:region:sector:ncluster:meanadc"
        int nCluster = nclus[iside][iregion][isector];
        float meanAdc = 0;
        if (nCluster > 0)
        {
          meanAdc = madc[iside][iregion][isector] / (nCluster*1.);
        }
        nt_sector_event_summary->Fill(m_event, m_segment, m_bco, iside, iregion, isector, nCluster, meanAdc);
      }
    }
  }
  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<TrkrDefs::cluskey> TpcSeedsQA::get_cluster_keys(SvtxTrack *track)
{
  std::vector<TrkrDefs::cluskey> out;
  for (const auto &seed : {track->get_silicon_seed(), track->get_tpc_seed()})
  {
    if (seed)
    {
      std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
    }
  }
  return out;
}

//____________________________________________________________________________..
int TpcSeedsQA::EndRun(const int /*runnumber*/)
{
  for (auto &region : {0, 1, 2})
  {
    for (auto &index_pt : {0, 1, 2, 3})
    {
      if (num_track_side0_pt[region][index_pt] > 0)
      {
        h_cluster_phisize1_fraction_mean_side0[region]->SetBinContent(index_pt + 1, frac_side0_pt[region][index_pt] / num_track_side0_pt[region][index_pt]);
        h_cluster_phisize1_fraction_mean_numerator_side0[region]->SetBinContent(index_pt + 1, frac_side0_pt[region][index_pt]);
        h_cluster_phisize1_fraction_mean_denominator_side0[region]->SetBinContent(index_pt + 1, num_track_side0_pt[region][index_pt]);
      }
      if (num_track_side1_pt[region][index_pt] > 0)
      {
        h_cluster_phisize1_fraction_mean_side1[region]->SetBinContent(index_pt + 1, frac_side1_pt[region][index_pt] / num_track_side1_pt[region][index_pt]);
        h_cluster_phisize1_fraction_mean_numerator_side1[region]->SetBinContent(index_pt + 1, frac_side1_pt[region][index_pt]);
        h_cluster_phisize1_fraction_mean_denominator_side1[region]->SetBinContent(index_pt + 1, num_track_side1_pt[region][index_pt]);
      }
    }
  }

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcSeedsQA::End(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string TpcSeedsQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

void TpcSeedsQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpc_fullpt_pos").c_str(), "TPC clusters per positive track;Number of TPC clusters per positive track;Entries", 55, -0.5, 54.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpc_fullpt_neg").c_str(), "TPC clusters per negative track;Number of TPC clusters per negative track;Entries", 55, -0.5, 54.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpc_pos").c_str(), "TPC clusters per positive track (pT>1GeV,eta cut);Number of TPC clusters per positive track;Entries", 55, -0.5, 54.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpc_neg").c_str(), "TPC clusters per negative track (pT>1GeV,eta cut);Number of TPC clusters per negative track;Entries", 55, -0.5, 54.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpot_pos").c_str(), "TPOT clusters per positive track (pT>1GeV);Number of TPOT clusters per positive track;Entries", 2, -0.5, 1.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpot_neg").c_str(), "TPOT clusters per negative track (pT>1GeV);Number of TPOT clusters per negative track;Entries", 2, -0.5, 1.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "ntpc_quality_pos").c_str(), "Number of TPC clusters per positive track (pT>1GeV);Number of TPC clusters per positive track;Quality", 55, -0.5, 54.5, 100, 0, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "ntpc_quality_neg").c_str(), "Number of TPC clusters per negative track (pT>1GeV);Number of TPC clusters per negative track;Quality", 55, -0.5, 54.5, 100, 0, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d").c_str(), "Number of reconstructed tracks;Number of TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d_pos").c_str(), "Number of reconstructed positive tracks;Number of positive TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d_neg").c_str(), "Number of reconstructed negative tracks;Number of negative TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d_ptg1").c_str(), "Number of reconstructed tracks (pT>1GeV);Number of TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d_ptg1_pos").c_str(), "Number of reconstructed positive tracks (pT>1GeV);Number of positive TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d_ptg1_neg").c_str(), "Number of reconstructed negative tracks (pT>1GeV);Number of negative TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "pt").c_str(), "p_{T} distribution of reconstructed tracks;Track p_{T};Entries", 100, 0, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "pt_pos").c_str(), "p_{T} distribution of reconstructed positive tracks;Track p_{T};Entries", 100, 0, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "pt_neg").c_str(), "p_{T} distribution of reconstructed negative tracks;Track p_{T};Entries", 100, 0, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "nrecotracks_pos").c_str(), "Number of reconstructed positive tracks (pT>1GeV);#eta;#phi [rad];Entries", 100, -1.1, 1.1, 300, -3.14159, 3.1459);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "nrecotracks_neg").c_str(), "Number of reconstructed negative tracks (pT>1GeV);#eta;#phi [rad];Entries", 100, -1.1, 1.1, 300, -3.14159, 3.1459);
    hm->registerHisto(h);
  }

  {
    auto h = new TProfile2D(std::string(getHistoPrefix() + "avgnclus_eta_phi_pos").c_str(), "Average number of clusters per positive track (pT>1GeV);#eta;#phi [rad];Average number of clusters per positive track", 100, -1.1, 1.1, 300, -3.14159, 3.1459, 0, 55);
    hm->registerHisto(h);
  }

  {
    auto h = new TProfile2D(std::string(getHistoPrefix() + "avgnclus_eta_phi_neg").c_str(), "Average number of clusters per negative track (pT>1GeV);#eta;#phi [rad];Average number of clusters per negative track", 100, -1.1, 1.1, 300, -3.14159, 3.1459, 0, 55);
    hm->registerHisto(h);
  }

  //  {
  //    auto h = new TH1F(std::string(getHistoPrefix() + "trackcrossing_pos").c_str(), "Positive track beam bunch crossing (pT>1GeV);Positive track crossing;Entries", 100, -100, 300);
  //    hm->registerHisto(h);
  //  }

  //  {
  //    auto h = new TH1F(std::string(getHistoPrefix() + "trackcrossing_neg").c_str(), "Negative track beam bunch crossing (pT>1GeV);Negative track crossing;Entries", 100, -100, 300);
  //    hm->registerHisto(h);
  //  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcaxyorigin_phi_north_pos").c_str(), "DCA xy origin vs phi for positive track (dcaz>0);#phi [rad];DCA_{xy} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcaxyorigin_phi_south_pos").c_str(), "DCA xy origin vs phi for positive track (dcaz<0);#phi [rad];DCA_{xy} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcaxyorigin_phi_north_neg").c_str(), "DCA xy origin vs phi for negative track (dcaz>0);#phi [rad];DCA_{xy} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcaxyorigin_phi_south_neg").c_str(), "DCA xy origin vs phi for negative track (dcaz<0);#phi [rad];DCA_{xy} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcaxyvtx_phi_pos").c_str(), "DCA xy vertex vs phi for positive track;#phi [rad];DCA_{xy} wrt vertex [cm];Entries", 300, -3.14159, 3.1459, 90, -3, 3);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcaxyvtx_phi_neg").c_str(), "DCA xy vertex vs phi for negative track;#phi [rad];DCA_{xy} wrt vertex [cm];Entries", 300, -3.14159, 3.1459, 90, -3, 3);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcazorigin_phi_pos").c_str(), "DCA z origin vs phi for positive track;#phi [rad];DCA_{z} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcazorigin_phi_neg").c_str(), "DCA z origin vs phi for negative track;#phi [rad];DCA_{z} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcazvtx_phi_pos").c_str(), "DCA z vertex vs phi for positive track;#phi [rad];DCA_{z} wrt vertex [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcazvtx_phi_neg").c_str(), "DCA z vertex vs phi for negative track;#phi [rad];DCA_{z} wrt vertex [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntrack_isfromvtx_pos").c_str(), "Num of positive tracks associated to a vertex;Is track associated to a vertex;Entries", 2, -0.5, 1.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntrack_isfromvtx_neg").c_str(), "Num of negative tracks associated to a vertex;Is track associated to a vertex;Entries", 2, -0.5, 1.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "cluster_phisize1_fraction_pos").c_str(), "Fraction of TPC clusters per positive track with phi size of 1 (pT>1GeV);Fraction of TPC clusters phi size of 1;Entries", 100, 0, 1);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "cluster_phisize1_fraction_neg").c_str(), "Fraction of TPC clusters per negative track with phi size of 1 (pT>1GeV);Fraction of TPC clusters phi size of 1;Entries", 100, 0, 1);
    hm->registerHisto(h);
  }

  // vertex
  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecovertices").c_str(), "Num of reco vertices per event;Number of vertices;Entries", 20, 0, 20);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "vx").c_str(), "Vertex x;Vertex x [cm];Entries", 100, -2.5, 2.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "vy").c_str(), "Vertex y;Vertex y [cm];Entries", 100, -2.5, 2.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "vx_vy").c_str(), "Vertex x vs y;Vertex x [cm];Vertex y [cm];Entries", 100, -2.5, 2.5, 100, -2.5, 2.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "vz").c_str(), "Vertex z;Vertex z [cm];Entries", 50, -25, 25);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "vt").c_str(), "Vertex t;Vertex t [ns];Entries", 100, -1000, 20000);
    hm->registerHisto(h);
  }

  //  {
  //    auto h = new TH1F(std::string(getHistoPrefix() + "vertexcrossing").c_str(), "Vertex beam bunch crossing;Vertex crossing;Entries", 100, -100, 300);
  //    hm->registerHisto(h);
  //  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "vertexchi2dof").c_str(), "Vertex chi2/ndof;Vertex #chi2/ndof;Entries", 100, 0, 20);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntrackspervertex").c_str(), "Num of tracks per vertex;Number of tracks per vertex;Entries", 50, 0, 50);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dedx").c_str(),
                      "Num of tracks per vertex;Number of tracks per vertex;Entries",
                      500, -2, 2, 500, 0, 3000);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "mip_dedx").c_str(),
                      "dEdx of MIPs",
                      100, 0, 1000);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dedx_pcaz").c_str(),
                      "track dEdx vs pcaz",
                      100, -100, 100, 50, 0, 2000);
    hm->registerHisto(h);
  }

  for (int i = 0; i < 10; i++)
  {
    h_dedx_pq_z[i] = new TH2F((boost::format("%sdedx_pq_%i") % getHistoPrefix() % i).str().c_str(),
                              (boost::format("mean cluster dEdx at R2 & R3 corrected by path length according to track eta in z-region %i") % i).str().c_str(), 100, -3, 3, 100, 0, 10000);
    hm->registerHisto(h_dedx_pq_z[i]);
  }

  {
    auto nt = new TNtuple(std::string(getHistoPrefix() + "sector_event_summary").c_str(),
		      "sector_event_summary","event:segment:bco:side:region:sector:ncluster:meanadc");
    hm->registerHisto(nt);
  }

  for (auto &region : {0, 1, 2})
  {
    h_adc_sector[region] = new TH2F((boost::format("%sadc_sector_%i") % getHistoPrefix() % region).str().c_str(),
                                    (boost::format("ADC spectrum per, region_%i") % region).str().c_str(), 25, -12.5, 12.5, 50, 0, 1500);
    hm->registerHisto(h_adc_sector[region]);

    h_onepad_frac[region] = new TProfile((boost::format("%sonepad_frac_%i") % getHistoPrefix() % region).str().c_str(),
                                         (boost::format("TPC Cluster Phi Size == 1 fraction per sector, region_%i") % region).str().c_str(), 25, -12.5, 12.5);
    hm->registerHisto(h_onepad_frac[region]);

    h_clusphisize1pt_side0[region] = new TH1F((boost::format("%sclusphisize1pT_side0_%i") % getHistoPrefix() % region).str().c_str(),
                                              (boost::format("TPC Cluster Phi Size == 1, side 0, region_%i") % region).str().c_str(), 4, 1, 3.2);
    h_clusphisize1pt_side0[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hm->registerHisto(h_clusphisize1pt_side0[region]);

    h_clusphisize1pt_side1[region] = new TH1F((boost::format("%sclusphisize1pT_side1_%i") % getHistoPrefix() % region).str().c_str(),
                                              (boost::format("TPC Cluster Phi Size == 1, side 1, region_%i") % region).str().c_str(), 4, 1, 3.2);
    h_clusphisize1pt_side1[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hm->registerHisto(h_clusphisize1pt_side1[region]);

    h_clusphisizegeq1pt_side0[region] = new TH1F((boost::format("%sclusphisizegeq1pT_side0_%i") % getHistoPrefix() % region).str().c_str(),
                                                 (boost::format("TPC Cluster Phi Size >= 1, side 0, region_%i") % region).str().c_str(), 4, 1, 3.2);
    h_clusphisizegeq1pt_side0[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hm->registerHisto(h_clusphisizegeq1pt_side0[region]);

    h_clusphisizegeq1pt_side1[region] = new TH1F((boost::format("%sclusphisizegeq1pT_side1_%i") % getHistoPrefix() % region).str().c_str(),
                                                 (boost::format("TPC Cluster Phi Size >= 1, side 1, region_%i") % region).str().c_str(), 4, 1, 3.2);
    h_clusphisizegeq1pt_side1[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hm->registerHisto(h_clusphisizegeq1pt_side1[region]);

    h_cluster_phisize1_fraction_side0[region] = new TH1F((boost::format("%sclusphisize1frac_side0_%i") % getHistoPrefix() % region).str().c_str(),
                                                         (boost::format("Fraction of TPC Cluster Phi Size == 1, side 0, region_%i") % region).str().c_str(), 100, 0, 1);
    h_cluster_phisize1_fraction_side0[region]->GetXaxis()->SetTitle("Fraction");
    hm->registerHisto(h_cluster_phisize1_fraction_side0[region]);

    h_cluster_phisize1_fraction_side1[region] = new TH1F((boost::format("%sclusphisize1frac_side1_%i") % getHistoPrefix() % region).str().c_str(),
                                                         (boost::format("Fraction of TPC Cluster Phi Size == 1, side 1, region_%i") % region).str().c_str(), 100, 0, 1);
    h_cluster_phisize1_fraction_side1[region]->GetXaxis()->SetTitle("Fraction");
    hm->registerHisto(h_cluster_phisize1_fraction_side1[region]);

    h_cluster_phisize1_fraction_pt_side0[region] = new TH2F((boost::format("%sclusphisize1frac_pt_side0_%i") % getHistoPrefix() % region).str().c_str(),
                                                            (boost::format("Pt vs. Fraction of TPC Cluster Phi Size == 1, side 0, region_%i") % region).str().c_str(), 4, 1, 3.2, 100, 0, 1);
    h_cluster_phisize1_fraction_pt_side0[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_cluster_phisize1_fraction_pt_side0[region]->GetYaxis()->SetTitle("Fraction");
    hm->registerHisto(h_cluster_phisize1_fraction_pt_side0[region]);

    h_cluster_phisize1_fraction_pt_side1[region] = new TH2F((boost::format("%sclusphisize1frac_pt_side1_%i") % getHistoPrefix() % region).str().c_str(),
                                                            (boost::format("Pt vs. Fraction of TPC Cluster Phi Size == 1, side 1, region_%i") % region).str().c_str(), 4, 1, 3.2, 100, 0, 1);
    h_cluster_phisize1_fraction_pt_side1[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_cluster_phisize1_fraction_pt_side1[region]->GetYaxis()->SetTitle("Fraction");
    hm->registerHisto(h_cluster_phisize1_fraction_pt_side1[region]);

    h_cluster_phisize1_fraction_mean_side0[region] = new TH1F((boost::format("%sclusphisize1frac_mean_side0_%i") % getHistoPrefix() % region).str().c_str(),
                                                              (boost::format("Pt vs. Average fraction of TPC Cluster Phi Size == 1, side 0, region_%i") % region).str().c_str(), 4, 1, 3.2);
    h_cluster_phisize1_fraction_mean_side0[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_cluster_phisize1_fraction_mean_side0[region]->GetYaxis()->SetTitle("Fraction");
    hm->registerHisto(h_cluster_phisize1_fraction_mean_side0[region]);

    h_cluster_phisize1_fraction_mean_side1[region] = new TH1F((boost::format("%sclusphisize1frac_mean_side1_%i") % getHistoPrefix() % region).str().c_str(),
                                                              (boost::format("Pt vs. Average fraction of TPC Cluster Phi Size == 1, side 1, region_%i") % region).str().c_str(), 4, 1, 3.2);
    h_cluster_phisize1_fraction_mean_side1[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_cluster_phisize1_fraction_mean_side1[region]->GetYaxis()->SetTitle("Fraction");
    hm->registerHisto(h_cluster_phisize1_fraction_mean_side1[region]);

    h_cluster_phisize1_fraction_mean_numerator_side0[region] = new TH1F((boost::format("%sclusphisize1frac_mean_numerator_side0_%i") % getHistoPrefix() % region).str().c_str(),
                                                                        (boost::format("Pt vs. Average fraction mean_numerator (sum of fraction) of TPC Cluster Phi Size == 1, side 0, region_%i") % region).str().c_str(), 4, 1, 3.2);
    h_cluster_phisize1_fraction_mean_numerator_side0[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_cluster_phisize1_fraction_mean_numerator_side0[region]->GetYaxis()->SetTitle("Fraction");
    hm->registerHisto(h_cluster_phisize1_fraction_mean_numerator_side0[region]);

    h_cluster_phisize1_fraction_mean_numerator_side1[region] = new TH1F((boost::format("%sclusphisize1frac_mean_numerator_side1_%i") % getHistoPrefix() % region).str().c_str(),
                                                                        (boost::format("Pt vs. Average fraction mean_numerator (sum of fraction) of TPC Cluster Phi Size == 1, side 1, region_%i") % region).str().c_str(), 4, 1, 3.2);
    h_cluster_phisize1_fraction_mean_numerator_side1[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_cluster_phisize1_fraction_mean_numerator_side1[region]->GetYaxis()->SetTitle("Fraction");
    hm->registerHisto(h_cluster_phisize1_fraction_mean_numerator_side1[region]);

    h_cluster_phisize1_fraction_mean_denominator_side0[region] = new TH1F((boost::format("%sclusphisize1frac_mean_denominator_side0_%i") % getHistoPrefix() % region).str().c_str(),
                                                                          (boost::format("Pt vs. Average fraction mean_denominator (sum of track) of TPC Cluster Phi Size == 1, side 0, region_%i") % region).str().c_str(), 4, 1, 3.2);
    h_cluster_phisize1_fraction_mean_denominator_side0[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_cluster_phisize1_fraction_mean_denominator_side0[region]->GetYaxis()->SetTitle("Fraction");
    hm->registerHisto(h_cluster_phisize1_fraction_mean_denominator_side0[region]);

    h_cluster_phisize1_fraction_mean_denominator_side1[region] = new TH1F((boost::format("%sclusphisize1frac_mean_denominator_side1_%i") % getHistoPrefix() % region).str().c_str(),
                                                                          (boost::format("Pt vs. Average fraction mean_denominator (sum of track) of TPC Cluster Phi Size == 1, side 1, region_%i") % region).str().c_str(), 4, 1, 3.2);
    h_cluster_phisize1_fraction_mean_denominator_side1[region]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_cluster_phisize1_fraction_mean_denominator_side1[region]->GetYaxis()->SetTitle("Fraction");
    hm->registerHisto(h_cluster_phisize1_fraction_mean_denominator_side1[region]);
  }
}

std::pair<float, float> TpcSeedsQA::cal_tpc_eta_min_max(float vtxz)
{
  float R = 780.;
  float HalfZ = 2110. / 2.;
  float theta_max = std::atan2(R, HalfZ - vtxz);
  float theta_min = std::atan2(R, -(vtxz + HalfZ));
  float eta_max = -std::log(std::tan(theta_max / 2));
  float eta_min = -std::log(std::tan(theta_min / 2));
  std::pair<float, float> min_max = std::make_pair(eta_min, eta_max);
  return min_max;
}

float TpcSeedsQA::eta_to_theta(float eta)
{
  return 2*std::atan(std::exp(-eta));
}
