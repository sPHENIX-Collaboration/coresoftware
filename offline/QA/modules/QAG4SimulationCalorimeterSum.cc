#include "QAG4SimulationCalorimeterSum.h"
#include "QAHistManagerDef.h"

#include <g4eval/CaloEvalStack.h>
#include <g4eval/CaloRawClusterEval.h>
#include <g4eval/SvtxEvalStack.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <calobase/RawCluster.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <trackbase_historic/SvtxTrack.h>

#include <g4eval/SvtxTrackEval.h>  // for SvtxTrackEval

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TString.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>  // for reverse_iterator
#include <utility>   // for pair
#include <vector>

QAG4SimulationCalorimeterSum::QAG4SimulationCalorimeterSum(
    QAG4SimulationCalorimeterSum::enu_flags flags)
  : SubsysReco("QAG4SimulationCalorimeterSum")
  , _flags(flags)
  , _calo_name_cemc("CEMC")
  , _calo_name_hcalin("HCALIN")
  , _calo_name_hcalout("HCALOUT")
  , _truth_container(nullptr)
  , _magField(+1.4)
{
}

int QAG4SimulationCalorimeterSum::InitRun(PHCompositeNode *topNode)
{
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                                "G4TruthInfo");
  if (!_truth_container)
  {
    std::cout << "QAG4SimulationCalorimeterSum::InitRun - Fatal Error - "
              << "unable to find DST node "
              << "G4TruthInfo" << std::endl;
    assert(_truth_container);
  }

  if (flag(kProcessCluster))
  {
    if (!_caloevalstack_cemc)
    {
      _caloevalstack_cemc.reset(
          new CaloEvalStack(topNode, _calo_name_cemc));
      _caloevalstack_cemc->set_strict(true);
      _caloevalstack_cemc->set_verbosity(Verbosity() + 1);
    }
    if (!_caloevalstack_hcalin)
    {
      _caloevalstack_hcalin.reset(
          new CaloEvalStack(topNode, _calo_name_hcalin));
      _caloevalstack_hcalin->set_strict(true);
      _caloevalstack_hcalin->set_verbosity(Verbosity() + 1);
    }
    if (!_caloevalstack_hcalout)
    {
      _caloevalstack_hcalout.reset(
          new CaloEvalStack(topNode, _calo_name_hcalout));
      _caloevalstack_hcalout->set_strict(true);
      _caloevalstack_hcalout->set_verbosity(Verbosity() + 1);
    }
  }

  if (flag(kProcessTrackProj))
  {
    if (!_svtxevalstack)
    {
      _svtxevalstack.reset(new SvtxEvalStack(topNode));
      _svtxevalstack->set_strict(true);
      _svtxevalstack->set_verbosity(Verbosity() + 1);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationCalorimeterSum::Init(PHCompositeNode *topNode)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH1D *h = new TH1D(TString(get_histo_prefix()) + "Normalization",  //
                     TString(get_histo_prefix()) + " Normalization;Items;Count", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event");
  h->GetXaxis()->SetBinLabel(i++, (_calo_name_cemc + " Tower").c_str());
  h->GetXaxis()->SetBinLabel(i++, (_calo_name_hcalin + " Tower").c_str());
  h->GetXaxis()->SetBinLabel(i++, (_calo_name_hcalout + " Tower").c_str());
  h->GetXaxis()->SetBinLabel(i++, (_calo_name_cemc + " Cluster").c_str());
  h->GetXaxis()->SetBinLabel(i++, (_calo_name_hcalin + " Cluster").c_str());
  h->GetXaxis()->SetBinLabel(i++, (_calo_name_hcalout + " Cluster").c_str());
  h->GetXaxis()->SetBinLabel(i++, "Track");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  //  if (flag(kProcessTower))
  //    {
  //      if (Verbosity() >= 1)
  //        std::cout << "QAG4SimulationCalorimeterSum::Init - Process tower occupancies"
  //            << std::endl;
  //      Init_Tower(topNode);
  //    }
  if (flag(kProcessCluster))
  {
    if (Verbosity() >= 1)
      std::cout << "QAG4SimulationCalorimeterSum::Init - Process tower occupancies"
                << std::endl;
    Init_Cluster(topNode);
  }

  if (flag(kProcessTrackProj))
  {
    if (Verbosity() >= 1)
      std::cout << "QAG4SimulationCalorimeterSum::Init - Process sampling fraction"
                << std::endl;
    Init_TrackProj(topNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationCalorimeterSum::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
    std::cout << "QAG4SimulationCalorimeterSum::process_event() entered" << std::endl;

  if (_caloevalstack_cemc)
    _caloevalstack_cemc->next_event(topNode);
  if (_caloevalstack_hcalin)
    _caloevalstack_hcalin->next_event(topNode);
  if (_caloevalstack_hcalout)
    _caloevalstack_hcalout->next_event(topNode);
  if (_svtxevalstack)
    _svtxevalstack->next_event(topNode);

  //  if (flag(kProcessTower))
  //    {
  //      int ret = process_event_Tower(topNode);
  //
  //      if (ret != Fun4AllReturnCodes::EVENT_OK)
  //        return ret;
  //    }

  if (flag(kProcessCluster))
  {
    int ret = process_event_Cluster(topNode);

    if (ret != Fun4AllReturnCodes::EVENT_OK)
      return ret;
  }

  if (flag(kProcessTrackProj))
  {
    int ret = process_event_TrackProj(topNode);

    if (ret != Fun4AllReturnCodes::EVENT_OK)
      return ret;
  }

  // at the end, count success events
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto(
      get_histo_prefix() + "Normalization"));
  assert(h_norm);
  h_norm->Fill("Event", 1);

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string
QAG4SimulationCalorimeterSum::get_histo_prefix()
{
  return "h_QAG4Sim_CalorimeterSum_";
}

PHG4Particle *
QAG4SimulationCalorimeterSum::get_truth_particle()
{
  // get last primary
  assert(_truth_container);
  assert(not _truth_container->GetMap().empty());
  PHG4Particle *last_primary = _truth_container->GetMap().rbegin()->second;
  assert(last_primary);

  return last_primary;
}

int QAG4SimulationCalorimeterSum::Init_TrackProj(PHCompositeNode * /*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  hm->registerHisto(
      new TH2F(
          TString(get_histo_prefix()) + TString(_calo_name_cemc.c_str()) + "_TrackProj",  //
          TString(_calo_name_cemc.c_str()) + " Tower Energy Distr. around Track Proj.;Polar distance / Tower width;Azimuthal distance / Tower width",
          (Max_N_Tower - 1) * 10, -Max_N_Tower / 2, Max_N_Tower / 2,
          (Max_N_Tower - 1) * 10, -Max_N_Tower / 2, Max_N_Tower / 2));

  hm->registerHisto(
      new TH2F(
          TString(get_histo_prefix()) + TString(_calo_name_hcalin.c_str()) + "_TrackProj",  //
          TString(_calo_name_hcalin.c_str()) + " Tower Energy Distr. around Track Proj.;Polar distance / Tower width;Azimuthal distance / Tower width",
          (Max_N_Tower - 1) * 10, -Max_N_Tower / 2, Max_N_Tower / 2,
          (Max_N_Tower - 1) * 10, -Max_N_Tower / 2, Max_N_Tower / 2));

  hm->registerHisto(
      new TH2F(
          TString(get_histo_prefix()) + TString(_calo_name_hcalout.c_str()) + "_TrackProj",  //
          TString(_calo_name_hcalout.c_str()) + " Tower Energy Distr. around Track Proj.;Polar distance / Tower width;Azimuthal distance / Tower width",
          (Max_N_Tower - 1) * 10, -Max_N_Tower / 2, Max_N_Tower / 2,
          (Max_N_Tower - 1) * 10, -Max_N_Tower / 2, Max_N_Tower / 2));

  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "TrackProj_3x3Tower_EP",  //
               "Tower 3x3 sum /E_{Truth};#Sigma_{3x3}[E_{Tower}] / total truth energy",
               150, 0, 1.5));

  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "TrackProj_5x5Tower_EP",  //
               "Tower 5x5 sum /E_{Truth};#Sigma_{5x5}[E_{Tower}] / total truth energy",
               150, 0, 1.5));

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationCalorimeterSum::process_event_TrackProj(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
    std::cout << "QAG4SimulationCalorimeterSum::process_event_TrackProj() entered"
              << std::endl;

  PHG4Particle *primary = get_truth_particle();
  if (!primary)
    return Fun4AllReturnCodes::DISCARDEVENT;

  SvtxTrackEval *trackeval = _svtxevalstack->get_track_eval();
  assert(trackeval);
  SvtxTrack *track = trackeval->best_track_from(primary);
  if (!track)
    return Fun4AllReturnCodes::EVENT_OK;  // not through the whole event for missing track.

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto(
      get_histo_prefix() + "Normalization"));
  assert(h_norm);
  h_norm->Fill("Track", 1);

  {
    TH1F *hsum = dynamic_cast<TH1F *>(hm->getHisto(
        (get_histo_prefix()) + "TrackProj_3x3Tower_EP"));
    assert(hsum);

    hsum->Fill(
        (track->get_cal_energy_3x3(SvtxTrack::CEMC) + track->get_cal_energy_3x3(SvtxTrack::HCALIN) + track->get_cal_energy_3x3(SvtxTrack::HCALOUT)) / (primary->get_e() + 1e-9));
  }
  {
    TH1F *hsum = dynamic_cast<TH1F *>(hm->getHisto(
        (get_histo_prefix()) + "TrackProj_5x5Tower_EP"));
    assert(hsum);

    hsum->Fill(
        (track->get_cal_energy_5x5(SvtxTrack::CEMC) + track->get_cal_energy_5x5(SvtxTrack::HCALIN) + track->get_cal_energy_5x5(SvtxTrack::HCALOUT)) / (primary->get_e() + 1e-9));
  }

  eval_trk_proj(_calo_name_cemc, track, topNode);
  eval_trk_proj(_calo_name_hcalin, track, topNode);
  eval_trk_proj(_calo_name_hcalout, track, topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

bool QAG4SimulationCalorimeterSum::eval_trk_proj(const std::string &detector, SvtxTrack *track,
                                                 PHCompositeNode *topNode)
// Track projections
{
  assert(track);
  assert(topNode);

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH2F *h2_proj = dynamic_cast<TH2F *>(hm->getHisto(
      (get_histo_prefix()) + detector + "_TrackProj"));
  assert(h2_proj);

  // pull the tower geometry
  std::string towergeonodename = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeo = findNode::getClass<RawTowerGeomContainer>(
      topNode, towergeonodename);
  assert(towergeo);

  if (Verbosity() > 2)
  {
    towergeo->identify();
  }
  // pull the towers
  std::string towernodename = "TOWER_CALIB_" + detector;
  RawTowerContainer *towerList = findNode::getClass<RawTowerContainer>(topNode,
                                                                       towernodename);
  assert(towerList);

  if (Verbosity() > 3)
  {
    std::cout << __PRETTY_FUNCTION__ << " - info - handling track track p = ("
              << track->get_px() << ", " << track->get_py() << ", "
              << track->get_pz() << "), v = (" << track->get_x() << ", "
              << track->get_z() << ", " << track->get_z() << ")"
              << " size_states "
              << track->size_states() << std::endl;
  }

  // curved tracks inside mag field
  // straight projections thereafter

  std::vector<double> point;
  point.assign(3, NAN);

  //  const double radius = towergeo->get_radius() + towergeo->get_thickness() * 0.5;

  //  PHG4HoughTransform::projectToRadius(track, _magField, radius, point);

  if (std::isnan(point[0]) or std::isnan(point[1]) or std::isnan(point[2]))
  {
    // std::cout << __PRETTY_FUNCTION__ << "::" << Name()
    //      << " - Error - track extrapolation failure:";
    // track->identify();
    return false;
  }

  assert(not std::isnan(point[0]));
  assert(not std::isnan(point[1]));
  assert(not std::isnan(point[2]));

  double x = point[0];
  double y = point[1];
  double z = point[2];

  double phi = atan2(y, x);
  double eta = asinh(z / sqrt(x * x + y * y));

  // calculate 3x3 tower energy
  int binphi = towergeo->get_phibin(phi);
  int bineta = towergeo->get_etabin(eta);

  double etabin_width = towergeo->get_etabounds(bineta).second - towergeo->get_etabounds(bineta).first;
  if (bineta > 1 and bineta < towergeo->get_etabins() - 1)
    etabin_width = (towergeo->get_etacenter(bineta + 1) - towergeo->get_etacenter(bineta - 1)) / 2.;

  double phibin_width = towergeo->get_phibounds(binphi).second - towergeo->get_phibounds(binphi).first;

  assert(etabin_width > 0);
  assert(phibin_width > 0);

  const double etabin_shift = (towergeo->get_etacenter(bineta) - eta) / etabin_width;
  const double phibin_shift = (fmod(
                                   towergeo->get_phicenter(binphi) - phi + 5 * M_PI, 2 * M_PI) -
                               M_PI) /
                              phibin_width;

  if (Verbosity() > 3)
    std::cout << __PRETTY_FUNCTION__ << " - info - handling track proj. (" << x
              << ", " << y << ", " << z << ")"
              << ", eta " << eta << ", phi" << phi
              << ", bineta " << bineta << " - ["
              << towergeo->get_etabounds(bineta).first << ", "
              << towergeo->get_etabounds(bineta).second << "], etabin_shift = "
              << etabin_shift << ", binphi " << binphi << " - ["
              << towergeo->get_phibounds(binphi).first << ", "
              << towergeo->get_phibounds(binphi).second << "], phibin_shift = "
              << phibin_shift << std::endl;

  const int bin_search_range = (Max_N_Tower - 1) / 2;
  for (int iphi = binphi - bin_search_range; iphi <= binphi + bin_search_range;
       ++iphi)

    for (int ieta = bineta - bin_search_range;
         ieta <= bineta + bin_search_range; ++ieta)
    {
      // wrap around
      int wrapphi = iphi;
      if (wrapphi < 0)
      {
        wrapphi = towergeo->get_phibins() + wrapphi;
      }
      if (wrapphi >= towergeo->get_phibins())
      {
        wrapphi = wrapphi - towergeo->get_phibins();
      }

      // edges
      if (ieta < 0)
        continue;
      if (ieta >= towergeo->get_etabins())
        continue;

      if (Verbosity() > 3)
        std::cout << __PRETTY_FUNCTION__ << " - info - processing tower"
                  << ", bineta " << ieta << " - ["
                  << towergeo->get_etabounds(ieta).first << ", "
                  << towergeo->get_etabounds(ieta).second << "]"
                  << ", binphi "
                  << wrapphi << " - [" << towergeo->get_phibounds(wrapphi).first
                  << ", " << towergeo->get_phibounds(wrapphi).second << "]" << std::endl;

      RawTower *tower = towerList->getTower(ieta, wrapphi);
      double energy = 0;
      if (tower)
      {
        if (Verbosity() > 1)
          std::cout << __PRETTY_FUNCTION__ << " - info - tower " << ieta << " "
                    << wrapphi << " energy = " << tower->get_energy() << std::endl;

        energy = tower->get_energy();
      }

      h2_proj->Fill(ieta - bineta + etabin_shift,
                    iphi - binphi + phibin_shift, energy);

    }  //            for (int ieta = bineta-1; ieta < bineta+2; ++ieta) {

  return true;
}  //       // Track projections

// int
// QAG4SimulationCalorimeterSum::Init_Tower(PHCompositeNode *topNode)
//{
//
//   Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
//   assert(hm);
//
////  TH1F * h = nullptr;
//
//  return Fun4AllReturnCodes::EVENT_OK;
//}
//
// int
// QAG4SimulationCalorimeterSum::process_event_Tower(PHCompositeNode *topNode)
//{
//  if (Verbosity() > 2)
//    std::cout << "QAG4SimulationCalorimeterSum::process_event_Tower() entered"
//        << std::endl;
//
//  return Fun4AllReturnCodes::EVENT_OK;
//}

int QAG4SimulationCalorimeterSum::Init_Cluster(PHCompositeNode * /*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  hm->registerHisto(
      new TH2F(
          TString(get_histo_prefix()) + "Cluster_" + _calo_name_cemc.c_str() + "_" + _calo_name_hcalin.c_str(),  //
          TString(_calo_name_hcalin.c_str()) + " VS " + TString(_calo_name_cemc.c_str()) + ": best cluster energy;" + TString(_calo_name_cemc.c_str()) + " cluster energy (GeV);" + TString(_calo_name_hcalin.c_str()) + " cluster energy (GeV)",
          70, 0, 70, 70, 0, 70));

  hm->registerHisto(
      new TH2F(
          TString(get_histo_prefix()) + "Cluster_" + _calo_name_cemc.c_str() + "_" + _calo_name_hcalin.c_str() + "_" + _calo_name_hcalout.c_str(),  //
          TString(_calo_name_cemc.c_str()) + " + " + TString(_calo_name_hcalin.c_str()) + " VS " + TString(_calo_name_hcalout.c_str()) + ": best cluster energy;" + TString(_calo_name_cemc.c_str()) + " + " + TString(_calo_name_hcalin.c_str()) + " cluster energy (GeV);" + TString(_calo_name_hcalout.c_str()) + " cluster energy (GeV)",
          70, 0, 70, 70, 0, 70));

  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "Cluster_EP",  //
               "Total Cluster E_{Reco}/E_{Truth};Reco cluster energy sum / total truth energy",
               150, 0, 1.5));

  hm->registerHisto(
      new TH1F(
          TString(get_histo_prefix()) + "Cluster_Ratio_" + _calo_name_cemc.c_str() + "_" + _calo_name_hcalin.c_str(),  //
          "Energy ratio " + TString(_calo_name_cemc.c_str()) + " VS " + TString(_calo_name_hcalin.c_str()) + ";Best cluster " + TString(_calo_name_cemc.c_str()) + " / (" + TString(_calo_name_cemc.c_str()) + " + " + TString(_calo_name_hcalin.c_str()) + ")", 110, 0, 1.1));

  hm->registerHisto(
      new TH1F(
          TString(get_histo_prefix()) + "Cluster_Ratio_" + _calo_name_cemc.c_str() + "_" + _calo_name_hcalin.c_str() + "_" + TString(_calo_name_hcalout.c_str()),  //
          "Energy ratio " + TString(_calo_name_cemc.c_str()) + " + " + TString(_calo_name_hcalin.c_str()) + " VS " + TString(_calo_name_hcalout.c_str()) + ";Best cluster (" + TString(_calo_name_cemc.c_str()) + " + " + TString(_calo_name_hcalin.c_str()) + ") / (" + TString(_calo_name_cemc.c_str()) + " + " + TString(_calo_name_hcalin.c_str()) + " + " + TString(_calo_name_hcalout.c_str()) + ")", 110, 0, 1.1));

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationCalorimeterSum::process_event_Cluster(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 2)
    std::cout << "QAG4SimulationCalorimeterSum::process_event_Cluster() entered"
              << std::endl;

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  PHG4Particle *primary = get_truth_particle();
  if (!primary)
    return Fun4AllReturnCodes::DISCARDEVENT;

  RawCluster *cluster_cemc =
      _caloevalstack_cemc->get_rawcluster_eval()->best_cluster_from(primary);
  RawCluster *cluster_hcalin =
      _caloevalstack_hcalin->get_rawcluster_eval()->best_cluster_from(primary);
  RawCluster *cluster_hcalout =
      _caloevalstack_hcalout->get_rawcluster_eval()->best_cluster_from(primary);

  const double cluster_cemc_e = cluster_cemc ? cluster_cemc->get_energy() : 0;
  const double cluster_hcalin_e =
      cluster_hcalin ? cluster_hcalin->get_energy() : 0;
  const double cluster_hcalout_e =
      cluster_hcalout ? cluster_hcalout->get_energy() : 0;

  if (cluster_cemc_e + cluster_hcalin_e > 0)
  {
    TH2F *h2 = dynamic_cast<TH2F *>(hm->getHisto(
        (get_histo_prefix()) + "Cluster_" + _calo_name_cemc + "_" + _calo_name_hcalin));
    assert(h2);

    h2->Fill(cluster_cemc_e, cluster_hcalin_e);

    TH1F *hr = dynamic_cast<TH1F *>(hm->getHisto(
        (get_histo_prefix()) + "Cluster_Ratio_" + _calo_name_cemc + "_" + _calo_name_hcalin));
    assert(hr);

    hr->Fill(cluster_cemc_e / (cluster_cemc_e + cluster_hcalin_e));
  }

  if (cluster_cemc_e + cluster_hcalin_e + cluster_hcalout_e > 0)
  {
    if (Verbosity())
    {
      std::cout << "QAG4SimulationCalorimeterSum::process_event_Cluster - "
                << " cluster_cemc_e = " << cluster_cemc_e
                << " cluster_hcalin_e = " << cluster_hcalin_e
                << " cluster_hcalout_e = " << cluster_hcalout_e << " hr = "
                << (cluster_cemc_e + cluster_hcalin_e) / (cluster_cemc_e + cluster_hcalin_e + cluster_hcalout_e)
                << std::endl;
    }

    TH2F *h2 = dynamic_cast<TH2F *>(hm->getHisto(
        (get_histo_prefix()) + "Cluster_" + _calo_name_cemc + "_" + _calo_name_hcalin + "_" + _calo_name_hcalout));
    assert(h2);

    h2->Fill((cluster_cemc_e + cluster_hcalin_e), cluster_hcalout_e);

    TH1F *hr = dynamic_cast<TH1F *>(hm->getHisto(
        (get_histo_prefix()) + "Cluster_Ratio_" + _calo_name_cemc + "_" + _calo_name_hcalin + "_" + _calo_name_hcalout));
    assert(hr);

    hr->Fill(
        (cluster_cemc_e + cluster_hcalin_e) / (cluster_cemc_e + cluster_hcalin_e + cluster_hcalout_e));

    TH1F *hsum = dynamic_cast<TH1F *>(hm->getHisto(
        (get_histo_prefix()) + "Cluster_EP"));
    assert(hsum);

    hsum->Fill(
        (cluster_cemc_e + cluster_hcalin_e + cluster_hcalout_e) / (primary->get_e() + 1e-9));
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
