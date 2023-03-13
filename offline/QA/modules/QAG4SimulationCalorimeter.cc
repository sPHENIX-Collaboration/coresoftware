#include "QAG4SimulationCalorimeter.h"

#include "QAHistManagerDef.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <g4eval/CaloEvalStack.h>
#include <g4eval/CaloRawClusterEval.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TString.h>
#include <TVector3.h>

#include <CLHEP/Vector/ThreeVector.h>  // for Hep3Vector

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>  // for reverse_iterator
#include <map>
#include <utility>

QAG4SimulationCalorimeter::QAG4SimulationCalorimeter(const std::string &calo_name,
                                                     QAG4SimulationCalorimeter::enu_flags flags)
  : SubsysReco("QAG4SimulationCalorimeter_" + calo_name)
  , _calo_name(calo_name)
  , _flags(flags)
  , _calo_hit_container(nullptr)
  , _calo_abs_hit_container(nullptr)
  , _truth_container(nullptr)
{
}

int QAG4SimulationCalorimeter::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  assert(dstNode);

  if (flag(kProcessG4Hit))
  {
    _calo_hit_container = findNode::getClass<PHG4HitContainer>(topNode,
                                                               "G4HIT_" + _calo_name);
    if (!_calo_hit_container)
    {
      std::cout << "QAG4SimulationCalorimeter::InitRun - Fatal Error - "
                << "unable to find DST node "
                << "G4HIT_" + _calo_name << std::endl;
      assert(_calo_hit_container);
    }

    _calo_abs_hit_container = findNode::getClass<PHG4HitContainer>(topNode,
                                                                   "G4HIT_ABSORBER_" + _calo_name);
    if (!_calo_abs_hit_container)
    {
      std::cout << "QAG4SimulationCalorimeter::InitRun - Fatal Error - "
                << "unable to find DST node "
                << "G4HIT_ABSORBER_" + _calo_name
                << std::endl;
      assert(_calo_abs_hit_container);
    }
  }

  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                                "G4TruthInfo");
  if (!_truth_container)
  {
    std::cout << "QAG4SimulationCalorimeter::InitRun - Fatal Error - "
              << "unable to find DST node "
              << "G4TruthInfo" << std::endl;
    assert(_truth_container);
  }

  if (flag(kProcessCluster))
  {
    if (!_caloevalstack)
    {
      _caloevalstack.reset(new CaloEvalStack(topNode, _calo_name));
      _caloevalstack->set_strict(true);
      _caloevalstack->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstack->next_event(topNode);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationCalorimeter::Init(PHCompositeNode *topNode)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH1D *h = new TH1D(TString(get_histo_prefix()) + "_Normalization",  //
                     TString(_calo_name) + " Normalization;Items;Count", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event");
  h->GetXaxis()->SetBinLabel(i++, "G4Hit Active");
  h->GetXaxis()->SetBinLabel(i++, "G4Hit Absor.");
  h->GetXaxis()->SetBinLabel(i++, "Tower");
  h->GetXaxis()->SetBinLabel(i++, "Tower Hit");
  h->GetXaxis()->SetBinLabel(i++, "Cluster");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  if (flag(kProcessG4Hit))
  {
    if (Verbosity() >= 1)
      std::cout << "QAG4SimulationCalorimeter::Init - Process sampling fraction"
                << std::endl;
    Init_G4Hit(topNode);
  }
  if (flag(kProcessTower))
  {
    if (Verbosity() >= 1)
      std::cout << "QAG4SimulationCalorimeter::Init - Process tower occupancies"
                << std::endl;
    Init_Tower(topNode);
  }
  if (flag(kProcessCluster))
  {
    if (Verbosity() >= 1)
      std::cout << "QAG4SimulationCalorimeter::Init - Process tower occupancies"
                << std::endl;
    Init_Cluster(topNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationCalorimeter::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
    std::cout << "QAG4SimulationCalorimeter::process_event() entered" << std::endl;

  if (_caloevalstack)
    _caloevalstack->next_event(topNode);

  if (flag(kProcessG4Hit))
  {
    int ret = process_event_G4Hit(topNode);

    if (ret != Fun4AllReturnCodes::EVENT_OK)
      return ret;
  }

  if (flag(kProcessTower))
  {
    int ret = process_event_Tower(topNode);

    if (ret != Fun4AllReturnCodes::EVENT_OK)
      return ret;
  }

  if (flag(kProcessCluster))
  {
    int ret = process_event_Cluster(topNode);

    if (ret != Fun4AllReturnCodes::EVENT_OK)
      return ret;
  }

  // at the end, count success events
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto(
      get_histo_prefix() + "_Normalization"));
  assert(h_norm);
  h_norm->Fill("Event", 1);

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string
QAG4SimulationCalorimeter::get_histo_prefix()
{
  return "h_QAG4Sim_" + std::string(_calo_name);
}

int QAG4SimulationCalorimeter::Init_G4Hit(PHCompositeNode * /*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  hm->registerHisto(
      new TH2F(TString(get_histo_prefix()) + "_G4Hit_RZ",  //
               TString(_calo_name) + " RZ projection;G4 Hit Z (cm);G4 Hit R (cm)", 1200, -300, 300,
               600, -000, 300));

  hm->registerHisto(
      new TH2F(TString(get_histo_prefix()) + "_G4Hit_XY",  //
               TString(_calo_name) + " XY projection;G4 Hit X (cm);G4 Hit Y (cm)", 1200, -300, 300,
               1200, -300, 300));

  hm->registerHisto(
      new TH2F(TString(get_histo_prefix()) + "_G4Hit_LateralTruthProjection",  //
               TString(_calo_name) + " shower lateral projection (last primary);Polar direction (cm);Azimuthal direction (cm)",
               200, -30, 30, 200, -30, 30));

  hm->registerHisto(new TH1F(TString(get_histo_prefix()) + "_G4Hit_SF",  //
                             TString(_calo_name) + " sampling fraction;Sampling fraction", 1000, 0, .2));

  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "_G4Hit_VSF",  //
               TString(_calo_name) + " visible sampling fraction;Visible sampling fraction", 1000, 0,
               .2));

  TH1F *h =
      new TH1F(TString(get_histo_prefix()) + "_G4Hit_HitTime",  //
               TString(_calo_name) + " hit time (edep weighting);Hit time - T0 (ns);Geant4 energy density",
               1000, 0.5, 10000);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);

  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "_G4Hit_FractionTruthEnergy",  //
               TString(_calo_name) + " fraction truth energy ;G4 edep / particle energy",
               1000, 0, 1));

  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "_G4Hit_FractionEMVisibleEnergy",  //
               TString(_calo_name) + " fraction visible energy from EM; visible energy from e^{#pm} / total visible energy",
               100, 0, 1));

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationCalorimeter::process_event_G4Hit(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 2)
    std::cout << "QAG4SimulationCalorimeter::process_event_G4Hit() entered" << std::endl;

  TH1F *h = nullptr;

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto(
      get_histo_prefix() + "_Normalization"));
  assert(h_norm);

  // get primary
  assert(_truth_container);
  PHG4TruthInfoContainer::ConstRange primary_range =
      _truth_container->GetPrimaryParticleRange();
  double total_primary_energy = 1e-9;  // make it zero energy epsilon samll so it can be used for denominator
  for (PHG4TruthInfoContainer::ConstIterator particle_iter = primary_range.first;
       particle_iter != primary_range.second; ++particle_iter)
  {
    const PHG4Particle *particle = particle_iter->second;
    assert(particle);
    total_primary_energy += particle->get_e();
  }

  assert(not _truth_container->GetMap().empty());
  const PHG4Particle *last_primary =
      _truth_container->GetMap().rbegin()->second;
  assert(last_primary);

  if (Verbosity() > 2)
  {
    std::cout
        << "QAG4SimulationCalorimeter::process_event_G4Hit() handle this truth particle"
        << std::endl;
    last_primary->identify();
  }
  const PHG4VtxPoint *primary_vtx =  //
      _truth_container->GetPrimaryVtx(last_primary->get_vtx_id());
  assert(primary_vtx);
  if (Verbosity() > 2)
  {
    std::cout
        << "QAG4SimulationCalorimeter::process_event_G4Hit() handle this vertex"
        << std::endl;
    primary_vtx->identify();
  }

  const double t0 = primary_vtx->get_t();
  const TVector3 vertex(primary_vtx->get_x(), primary_vtx->get_y(),
                        primary_vtx->get_z());

  // projection axis
  TVector3 axis_proj(last_primary->get_px(), last_primary->get_py(),
                     last_primary->get_pz());
  if (axis_proj.Mag() == 0)
    axis_proj.SetXYZ(0, 0, 1);
  axis_proj = axis_proj.Unit();

  // azimuthal direction axis
  TVector3 axis_azimuth = axis_proj.Cross(TVector3(0, 0, 1));
  if (axis_azimuth.Mag() == 0)
    axis_azimuth.SetXYZ(1, 0, 0);
  axis_azimuth = axis_azimuth.Unit();

  // polar direction axis
  TVector3 axis_polar = axis_proj.Cross(axis_azimuth);
  assert(axis_polar.Mag() > 0);
  axis_polar = axis_polar.Unit();

  double e_calo = 0.0;      // active energy deposition
  double ev_calo = 0.0;     // visible energy
  double ea_calo = 0.0;     // absorber energy
  double ev_calo_em = 0.0;  // EM visible energy

  if (_calo_hit_container)
  {
    TH2F *hrz = dynamic_cast<TH2F *>(hm->getHisto(
        get_histo_prefix() + "_G4Hit_RZ"));
    assert(hrz);
    TH2F *hxy = dynamic_cast<TH2F *>(hm->getHisto(
        get_histo_prefix() + "_G4Hit_XY"));
    assert(hxy);
    TH1F *ht = dynamic_cast<TH1F *>(hm->getHisto(
        get_histo_prefix() + "_G4Hit_HitTime"));
    assert(ht);
    TH2F *hlat = dynamic_cast<TH2F *>(hm->getHisto(
        get_histo_prefix() + "_G4Hit_LateralTruthProjection"));
    assert(hlat);

    h_norm->Fill("G4Hit Active", _calo_hit_container->size());
    PHG4HitContainer::ConstRange calo_hit_range =
        _calo_hit_container->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = calo_hit_range.first;
         hit_iter != calo_hit_range.second; hit_iter++)
    {
      PHG4Hit *this_hit = hit_iter->second;
      assert(this_hit);

      e_calo += this_hit->get_edep();
      ev_calo += this_hit->get_light_yield();

      // EM visible energy that is only associated with electron energy deposition
      PHG4Particle *particle = _truth_container->GetParticle(
          this_hit->get_trkid());
      if (!particle)
      {
        std::cout << __PRETTY_FUNCTION__ << " - Error - this PHG4hit missing particle: ";
        this_hit->identify();
      }
      assert(particle);
      if (abs(particle->get_pid()) == 11)
        ev_calo_em += this_hit->get_light_yield();

      const TVector3 hit(this_hit->get_avg_x(), this_hit->get_avg_y(),
                         this_hit->get_avg_z());

      hrz->Fill(hit.Z(), hit.Perp(), this_hit->get_edep());
      hxy->Fill(hit.X(), hit.Y(), this_hit->get_edep());
      ht->Fill(this_hit->get_avg_t() - t0, this_hit->get_edep());

      const double hit_azimuth = axis_azimuth.Dot(hit - vertex);
      const double hit_polar = axis_polar.Dot(hit - vertex);
      hlat->Fill(hit_polar, hit_azimuth, this_hit->get_edep());
    }
  }

  if (_calo_abs_hit_container)
  {
    h_norm->Fill("G4Hit Absor.", _calo_abs_hit_container->size());

    PHG4HitContainer::ConstRange calo_abs_hit_range =
        _calo_abs_hit_container->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = calo_abs_hit_range.first;
         hit_iter != calo_abs_hit_range.second; hit_iter++)
    {
      PHG4Hit *this_hit = hit_iter->second;
      assert(this_hit);

      ea_calo += this_hit->get_edep();
    }
  }

  if (Verbosity() > 3)
    std::cout << "QAG4SimulationCalorimeter::process_event_G4Hit::" << _calo_name
              << " - SF = " << e_calo / (e_calo + ea_calo + 1e-9) << ", VSF = "
              << ev_calo / (e_calo + ea_calo + 1e-9) << std::endl;

  if (e_calo + ea_calo > 0)
  {
    h = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "_G4Hit_SF"));
    assert(h);
    h->Fill(e_calo / (e_calo + ea_calo));

    h = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "_G4Hit_VSF"));
    assert(h);
    h->Fill(ev_calo / (e_calo + ea_calo));
  }

  h = dynamic_cast<TH1F *>(hm->getHisto(
      get_histo_prefix() + "_G4Hit_FractionTruthEnergy"));
  assert(h);
  h->Fill((e_calo + ea_calo) / total_primary_energy);

  if (ev_calo > 0)
  {
    h = dynamic_cast<TH1F *>(hm->getHisto(
        get_histo_prefix() + "_G4Hit_FractionEMVisibleEnergy"));
    assert(h);
    h->Fill(ev_calo_em / (ev_calo));
  }

  if (Verbosity() > 3)
    std::cout << "QAG4SimulationCalorimeter::process_event_G4Hit::" << _calo_name
              << " - histogram " << h->GetName() << " Get Sum = " << h->GetSum()
              << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationCalorimeter::Init_Tower(PHCompositeNode * /*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1F *h = new TH1F(TString(get_histo_prefix()) + "_Tower_1x1",  //
                     TString(_calo_name) + " 1x1 tower;1x1 TOWER Energy (GeV)", 100, 9e-4, 100);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);

  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "_Tower_1x1_max",  //
               TString(_calo_name) + " 1x1 tower max per event;1x1 tower max per event (GeV)", 5000,
               0, 50));

  h = new TH1F(TString(get_histo_prefix()) + "_Tower_2x2",  //
               TString(_calo_name) + " 2x2 tower;2x2 TOWER Energy (GeV)", 100, 9e-4, 100);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "_Tower_2x2_max",  //
               TString(_calo_name) + " 2x2 tower max per event;2x2 tower max per event (GeV)", 5000,
               0, 50));

  h = new TH1F(TString(get_histo_prefix()) + "_Tower_3x3",  //
               TString(_calo_name) + " 3x3 tower;3x3 TOWER Energy (GeV)", 100, 9e-4, 100);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "_Tower_3x3_max",  //
               TString(_calo_name) + " 3x3 tower max per event;3x3 tower max per event (GeV)", 5000,
               0, 50));

  h = new TH1F(TString(get_histo_prefix()) + "_Tower_4x4",  //
               TString(_calo_name) + " 4x4 tower;4x4 TOWER Energy (GeV)", 100, 9e-4, 100);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "_Tower_4x4_max",  //
               TString(_calo_name) + " 4x4 tower max per event;4x4 tower max per event (GeV)", 5000,
               0, 50));

  h = new TH1F(TString(get_histo_prefix()) + "_Tower_5x5",  //
               TString(_calo_name) + " 5x5 tower;5x5 TOWER Energy (GeV)", 100, 9e-4, 100);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "_Tower_5x5_max",  //
               TString(_calo_name) + " 5x5 tower max per event;5x5 tower max per event (GeV)", 5000,
               0, 50));

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationCalorimeter::process_event_Tower(PHCompositeNode *topNode)
{
  const std::string detector(_calo_name);

  if (Verbosity() > 2)
    std::cout << "QAG4SimulationCalorimeter::process_event_Tower() entered" << std::endl;

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto(
      get_histo_prefix() + "_Normalization"));
  assert(h_norm);

  std::string towernodename = "TOWER_CALIB_" + detector;
  // Grab the towers
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode,
                                                                    towernodename);
  if (!towers)
  {
    std::cout << PHWHERE << ": Could not find node " << towernodename
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  std::string towergeomnodename = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(
      topNode, towergeomnodename);
  if (!towergeom)
  {
    std::cout << PHWHERE << ": Could not find node " << towergeomnodename
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  static const int max_size = 5;
  std::map<int, std::string> size_label;
  size_label[1] = "1x1";
  size_label[2] = "2x2";
  size_label[3] = "3x3";
  size_label[4] = "4x4";
  size_label[5] = "5x5";
  std::map<int, double> max_energy;
  std::map<int, TH1F *> energy_hist_list;
  std::map<int, TH1F *> max_energy_hist_list;

  for (int size = 1; size <= max_size; ++size)
  {
    max_energy[size] = 0;

    TH1F *h = dynamic_cast<TH1F *>(hm->getHisto(
        get_histo_prefix() + "_Tower_" + size_label[size]));
    assert(h);
    energy_hist_list[size] = h;
    h = dynamic_cast<TH1F *>(hm->getHisto(
        get_histo_prefix() + "_Tower_" + size_label[size] + "_max"));
    assert(h);
    max_energy_hist_list[size] = h;
  }

  h_norm->Fill("Tower", towergeom->size());  // total tower count
  h_norm->Fill("Tower Hit", towers->size());

  for (int binphi = 0; binphi < towergeom->get_phibins(); ++binphi)
  {
    for (int bineta = 0; bineta < towergeom->get_etabins(); ++bineta)
    {
      for (int size = 1; size <= max_size; ++size)
      {
        // for 2x2 and 4x4 use slide-2 window as implemented in DAQ
        if ((size == 2 or size == 4) and ((binphi % 2 != 0) and (bineta % 2 != 0)))
          continue;

        double energy = 0;

        // sliding window made from 2x2 sums
        for (int iphi = binphi; iphi < binphi + size; ++iphi)
        {
          for (int ieta = bineta; ieta < bineta + size; ++ieta)
          {
            if (ieta > towergeom->get_etabins())
              continue;

            // wrap around
            int wrapphi = iphi;
            assert(wrapphi >= 0);
            if (wrapphi >= towergeom->get_phibins())
            {
              wrapphi = wrapphi - towergeom->get_phibins();
            }

            RawTower *tower = towers->getTower(ieta, wrapphi);

            if (tower)
            {
              const double e_intput = tower->get_energy();

              energy += e_intput;
            }
          }
        }

        energy_hist_list[size]->Fill(energy == 0 ? 9.1e-4 : energy);  // trick to fill 0 energy tower to the first bin

        if (energy > max_energy[size])
          max_energy[size] = energy;

      }  //          for (int size = 1; size <= 4; ++size)
    }
  }

  for (int size = 1; size <= max_size; ++size)
  {
    max_energy_hist_list[size]->Fill(max_energy[size]);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationCalorimeter::Init_Cluster(PHCompositeNode * /*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  hm->registerHisto(
      new TH1F(TString(get_histo_prefix()) + "_Cluster_BestMatchERatio",  //
               TString(_calo_name) + " best matched cluster E/E_{Truth};E_{Cluster}/E_{Truth}", 150,
               0, 1.5));

  hm->registerHisto(
      new TH2F(TString(get_histo_prefix()) + "_Cluster_LateralTruthProjection",  //
               TString(_calo_name) + " best cluster lateral projection (last primary);Polar direction (cm);Azimuthal direction (cm)",
               200, -15, 15, 200, -15, 15));

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationCalorimeter::process_event_Cluster(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
  {
    std::cout << "QAG4SimulationCalorimeter::process_event_Cluster() entered"
              << std::endl;
  }

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  std::string towergeomnodename = "TOWERGEOM_" + _calo_name;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(
      topNode, towergeomnodename);
  if (!towergeom)
  {
    std::cout << PHWHERE << ": Could not find node " << towergeomnodename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get a cluster count
  TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto(get_histo_prefix() + "_Normalization"));
  assert(h_norm);

  std::string nodename = "CLUSTER_" + _calo_name;
  RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, nodename);
  assert(clusters);
  h_norm->Fill("Cluster", clusters->size());

  // get primary
  assert(_truth_container);
  assert(not _truth_container->GetMap().empty());
  PHG4Particle *last_primary = _truth_container->GetMap().rbegin()->second;
  assert(last_primary);

  if (Verbosity() > 2)
  {
    std::cout
        << "QAG4SimulationCalorimeter::process_event_Cluster() handle this truth particle"
        << std::endl;
    last_primary->identify();
  }

  assert(_caloevalstack);
  CaloRawClusterEval *clustereval = _caloevalstack->get_rawcluster_eval();
  assert(clustereval);

  TH1F *h = dynamic_cast<TH1F *>(hm->getHisto(
      get_histo_prefix() + "_Cluster_BestMatchERatio"));
  assert(h);

  RawCluster *cluster = clustereval->best_cluster_from(last_primary);
  if (cluster)
  {
    // has a cluster matched and best cluster selected

    if (Verbosity() > 3)
      std::cout << "QAG4SimulationCalorimeter::process_event_Cluster::"
                << _calo_name << " - get cluster with energy "
                << cluster->get_energy() << " VS primary energy "
                << last_primary->get_e() << std::endl;

    h->Fill(cluster->get_energy() / (last_primary->get_e() + 1e-9));  // avoids divide zero

    // now work on the projection:
    const CLHEP::Hep3Vector hit(cluster->get_position());

    const PHG4VtxPoint *primary_vtx =  //
        _truth_container->GetPrimaryVtx(last_primary->get_vtx_id());
    assert(primary_vtx);
    if (Verbosity() > 2)
    {
      std::cout
          << "QAG4SimulationCalorimeter::process_event_Cluster() handle this vertex"
          << std::endl;
      primary_vtx->identify();
    }

    const CLHEP::Hep3Vector vertex(primary_vtx->get_x(), primary_vtx->get_y(),
                                   primary_vtx->get_z());

    // projection axis
    CLHEP::Hep3Vector axis_proj(last_primary->get_px(), last_primary->get_py(),
                                last_primary->get_pz());
    if (axis_proj.mag() == 0)
      axis_proj.set(0, 0, 1);
    axis_proj = axis_proj.unit();

    // azimuthal direction axis
    CLHEP::Hep3Vector axis_azimuth = axis_proj.cross(CLHEP::Hep3Vector(0, 0, 1));
    if (axis_azimuth.mag() == 0)
      axis_azimuth.set(1, 0, 0);
    axis_azimuth = axis_azimuth.unit();

    // polar direction axis
    CLHEP::Hep3Vector axis_polar = axis_proj.cross(axis_azimuth);
    assert(axis_polar.mag() > 0);
    axis_polar = axis_polar.unit();

    TH2F *hlat = dynamic_cast<TH2F *>(hm->getHisto(
        get_histo_prefix() + "_Cluster_LateralTruthProjection"));
    assert(hlat);

    const double hit_azimuth = axis_azimuth.dot(hit - vertex);
    const double hit_polar = axis_polar.dot(hit - vertex);
    hlat->Fill(hit_polar, hit_azimuth);
  }
  else
  {
    if (Verbosity() > 3)
      std::cout << "QAG4SimulationCalorimeter::process_event_Cluster::"
                << _calo_name << " - missing cluster !";
    h->Fill(0);  // no cluster matched
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
