#include "QAG4SimulationCalorimeter.h"

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <phool/PHCompositeNode.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Particle.h>

//#include <g4hough/SvtxVertexMap.h>
//#include <g4hough/PHG4HoughTransform.h>

#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTowerGeomContainer.h>
#include <g4cemc/RawTower.h>
#include <g4cemc/RawClusterContainer.h>
#include <g4cemc/RawCluster.h>

#include <g4eval/CaloEvalStack.h>
#include <g4eval/CaloRawClusterEval.h>
#include <g4eval/CaloRawTowerEval.h>
#include <g4eval/CaloTruthEval.h>
#include <g4eval/SvtxEvalStack.h>

#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <exception>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>
#include <cmath>

#include "QAHistManagerDef.h"

using namespace std;

QAG4SimulationCalorimeter::QAG4SimulationCalorimeter(string _calo_name,
    QAG4SimulationCalorimeter::enu_flags flags) :
    SubsysReco("QAG4SimulationCalorimeter_" + _calo_name), //
    _eval_stack(NULL), _magfield(1.5), //
    _calo_name(_calo_name), _flags(flags), _ievent(0), //
    _calo_hit_container(NULL), _calo_abs_hit_container(NULL)
{

}

QAG4SimulationCalorimeter::~QAG4SimulationCalorimeter()
{
  if (_eval_stack)
    {
      delete _eval_stack;
    }
}

int
QAG4SimulationCalorimeter::InitRun(PHCompositeNode *topNode)
{
  _ievent = 0;

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));
  assert(dstNode);

  _calo_hit_container = findNode::getClass<PHG4HitContainer>(topNode,
      "G4HIT_" + _calo_name);
  if (!_calo_hit_container)
    {
      cout << "QAG4SimulationCalorimeter::InitRun - Fatal Error - "
          << "unable to find DST node " << "G4HIT_" + _calo_name << endl;
      assert(_calo_hit_container);
    }

  _calo_abs_hit_container = findNode::getClass<PHG4HitContainer>(topNode,
      "G4HIT_ABSORBER_" + _calo_name);
  if (!_calo_abs_hit_container)
    {
      cout << "QAG4SimulationCalorimeter::InitRun - Fatal Error - "
          << "unable to find DST node " << "G4HIT_ABSORBER_" + _calo_name
          << endl;
      assert(_calo_abs_hit_container);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationCalorimeter::End(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationCalorimeter::Init(PHCompositeNode *topNode)
{

  _ievent = 0;

  if (flag(kProcessSF))
    {
      cout << "QAG4SimulationCalorimeter::Init - Process sampling fraction"
          << endl;
      Init_SF(topNode);
    }
  if (flag(kProcessTower))
    {
      cout << "QAG4SimulationCalorimeter::Init - Process tower occupancies"
          << endl;
      Init_Tower(topNode);
    }
  if (flag(kProcessMCPhoton))
    {
      cout << "QAG4SimulationCalorimeter::Init - Process trakcs" << endl;
//      Init_MCPhoton(topNode);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationCalorimeter::process_event(PHCompositeNode *topNode)
{

  if (verbosity > 2)
    cout << "QAG4SimulationCalorimeter::process_event() entered" << endl;

  if (_eval_stack)
    _eval_stack->next_event(topNode);

  if (flag(kProcessSF))
    process_event_SF(topNode);
  if (flag(kProcessTower))
    process_event_Tower(topNode);
//  if (flag(kProcessMCPhoton))
//    process_event_MCPhoton(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationCalorimeter::Init_SF(PHCompositeNode *topNode)
{

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_SF", //
  "_h_" + TString(_calo_name) + "_SF", 1000, 0, .2));
  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_VSF", //
  "_h_" + TString(_calo_name) + "_VSF", 1000, 0, .2));

  hm->registerHisto(
      new TH2F("EMCalAna_h_" + TString(_calo_name) + "_RZ", //
      "EMCalAna_h_" + TString(_calo_name) + "_RZ;Z (cm);R (cm)", 1200, -300, 300, 600,
          -000, 300));

  hm->registerHisto(
      new TH2F("EMCalAna_h_" + TString(_calo_name) + "_XY", //
      "EMCalAna_h_" + TString(_calo_name) + "_XY;X (cm);Y (cm)", 1200, -300, 300, 1200,
          -300, 300));

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationCalorimeter::process_event_SF(PHCompositeNode *topNode)
{

  if (verbosity > 2)
    cout << "QAG4SimulationCalorimeter::process_event_SF() entered" << endl;

  TH1F* h = NULL;

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  double e_calo = 0.0;
  double ev_calo = 0.0;
  double ea_calo = 0.0;

  if (_calo_hit_container)
    {
      TH2F * hrz = (TH2F*) hm->getHisto("EMCalAna_h_" + _calo_name + "_RZ");
      assert(hrz);
      TH2F * hxy = (TH2F*) hm->getHisto("EMCalAna_h_" + _calo_name + "_XY");
      assert(hxy);

      PHG4HitContainer::ConstRange calo_hit_range =
          _calo_hit_container->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = calo_hit_range.first;
          hit_iter != calo_hit_range.second; hit_iter++)
        {

          PHG4Hit *this_hit = hit_iter->second;
          assert(this_hit);

          e_calo += this_hit->get_edep();
          ev_calo += this_hit->get_light_yield();

          const TVector3 hit(this_hit->get_avg_x(), this_hit->get_avg_y(),
              this_hit->get_avg_z());
          hrz->Fill(hit.Z(), hit.Perp(), this_hit->get_light_yield());
          hxy->Fill(hit.X(), hit.Y(), this_hit->get_light_yield());
        }
    }

  if (_calo_abs_hit_container)
    {
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

  h = (TH1F*) hm->getHisto("EMCalAna_h_" + _calo_name + "_SF");
  assert(h);
  h->Fill(e_calo / (e_calo + ea_calo) + 1e-9);
  h = (TH1F*) hm->getHisto("EMCalAna_h_" + _calo_name + "_VSF");
  assert(h);
  h->Fill(ev_calo / (e_calo + ea_calo) + 1e-9);

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationCalorimeter::Init_Tower(PHCompositeNode *topNode)
{

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_1x1", //
  "h_" + TString(_calo_name) + "_TOWER_1x1", 5000, 0, 50));
  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_1x1_max", //
  "EMCalAna_h_" + TString(_calo_name) + "_TOWER_1x1_max", 5000, 0, 50));
  hm->registerHisto(
      new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_1x1_max_trigger_ADC", //
      "EMCalAna_h_" + TString(_calo_name) + "_TOWER_1x1_max_trigger_ADC", 5000, 0, 50));

  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_2x2", //
  "h_" + TString(_calo_name) + "_TOWER_2x2", 5000, 0, 50));
  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_2x2_max", //
  "EMCalAna_h_" + TString(_calo_name) + "_TOWER_2x2_max", 5000, 0, 50));
  hm->registerHisto(
      new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_2x2_max_trigger_ADC", //
      "EMCalAna_h_" + TString(_calo_name) + "_TOWER_2x2_max_trigger_ADC", 5000, 0, 50));
  hm->registerHisto(
      new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_2x2_slide2_max_trigger_ADC", //
      "EMCalAna_h_" + TString(_calo_name) + "_TOWER_2x2_slide2_max_trigger_ADC", 5000, 0,
          50));

  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_3x3", //
  "h_" + TString(_calo_name) + "_TOWER_3x3", 5000, 0, 50));
  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_3x3_max", //
  "EMCalAna_h_" + TString(_calo_name) + "_TOWER_3x3_max", 5000, 0, 50));
  hm->registerHisto(
      new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_3x3_max_trigger_ADC", //
      "EMCalAna_h_" + TString(_calo_name) + "_TOWER_3x3_max_trigger_ADC", 5000, 0, 50));

  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_4x4", //
  "h_" + TString(_calo_name) + "_TOWER_4x4", 5000, 0, 50));
  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_4x4_max", //
  "EMCalAna_h_" + TString(_calo_name) + "_TOWER_4x4_max", 5000, 0, 50));
  hm->registerHisto(
      new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_4x4_max_trigger_ADC", //
      "EMCalAna_h_" + TString(_calo_name) + "_TOWER_4x4_max_trigger_ADC", 5000, 0, 50));
  hm->registerHisto(
      new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_4x4_slide2_max_trigger_ADC", //
      "EMCalAna_h_" + TString(_calo_name) + "_TOWER_4x4_slide2_max_trigger_ADC", 5000, 0,
          50));

  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_5x5", //
  "h_" + TString(_calo_name) + "_TOWER_4x4", 5000, 0, 50));
  hm->registerHisto(new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_5x5_max", //
  "EMCalAna_h_" + TString(_calo_name) + "_TOWER_5x5_max", 5000, 0, 50));
  hm->registerHisto(
      new TH1F("EMCalAna_h_" + TString(_calo_name) + "_TOWER_5x5_max_trigger_ADC", //
      "EMCalAna_h_" + TString(_calo_name) + "_TOWER_5x5_max_trigger_ADC", 5000, 0, 50));

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationCalorimeter::process_event_Tower(PHCompositeNode *topNode)
{
  const string detector(_calo_name);

  if (verbosity > 2)
    cout << "QAG4SimulationCalorimeter::process_event_SF() entered" << endl;

  string towernodename = "TOWER_CALIB_" + detector;
  // Grab the towers
  RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode,
      towernodename.c_str());
  if (!towers)
    {
      std::cout << PHWHERE << ": Could not find node " << towernodename.c_str()
          << std::endl;
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
  string towergeomnodename = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(
      topNode, towergeomnodename.c_str());
  if (!towergeom)
    {
      cout << PHWHERE << ": Could not find node " << towergeomnodename.c_str()
          << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  static const double trigger_ADC_bin = 45. / 256.; //8-bit ADC max to 45 GeV
  static const int max_size = 5;
  map<int, string> size_label;
  size_label[1] = "1x1";
  size_label[2] = "2x2";
  size_label[3] = "3x3";
  size_label[4] = "4x4";
  size_label[5] = "5x5";
  map<int, double> max_energy;
  map<int, double> max_energy_trigger_ADC;
  map<int, double> slide2_max_energy_trigger_ADC;
  map<int, TH1F*> energy_hist_list;
  map<int, TH1F*> max_energy_hist_list;
  map<int, TH1F*> max_energy_trigger_ADC_hist_list;
  map<int, TH1F*> slide2_max_energy_trigger_ADC_hist_list;

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  for (int size = 1; size <= max_size; ++size)
    {
      max_energy[size] = 0;
      max_energy_trigger_ADC[size] = 0;

      TH1F* h = NULL;

      h = (TH1F*) hm->getHisto(
          "EMCalAna_h_" + _calo_name + "_TOWER_" + size_label[size]);
      assert(h);
      energy_hist_list[size] = h;
      h = (TH1F*) hm->getHisto(
          "EMCalAna_h_" + _calo_name + "_TOWER_" + size_label[size] + "_max");
      assert(h);
      max_energy_hist_list[size] = h;

      h = (TH1F*) hm->getHisto(
          "EMCalAna_h_" + _calo_name + "_TOWER_" + size_label[size]
              + "_max_trigger_ADC");
      assert(h);
      max_energy_trigger_ADC_hist_list[size] = h;

      if (size == 2 or size == 4)
        {
          // sliding window made from 2x2 sums
          slide2_max_energy_trigger_ADC[size] = 0;

          h = (TH1F*) hm->getHisto(
              "EMCalAna_h_" + _calo_name + "_TOWER_" + size_label[size]
                  + "_slide2_max_trigger_ADC");
          assert(h);
          slide2_max_energy_trigger_ADC_hist_list[size] = h;

        }

    }

  for (int binphi = 0; binphi < towergeom->get_phibins(); ++binphi)
    {
      for (int bineta = 0; bineta < towergeom->get_etabins(); ++bineta)
        {
          for (int size = 1; size <= max_size; ++size)
            {
              double energy = 0;
              double energy_trigger_ADC = 0;
              double slide2_energy_trigger_ADC = 0;

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

                      RawTower* tower = towers->getTower(ieta, wrapphi);

                      if (tower)
                        {
                          const double e_intput = tower->get_energy();

                          const double e_trigger_ADC = round(
                              e_intput / trigger_ADC_bin) * trigger_ADC_bin;

                          energy += e_intput;
                          energy_trigger_ADC += e_trigger_ADC;

                          if ((size == 2 or size == 4) and (binphi % 2 == 0)
                              and (bineta % 2 == 0))
                            {
                              // sliding window made from 2x2 sums

                              slide2_energy_trigger_ADC += e_trigger_ADC;
                            }
                        }
                    }
                }

              energy_hist_list[size]->Fill(energy);

              if (energy > max_energy[size])
                max_energy[size] = energy;
              if (energy_trigger_ADC > max_energy_trigger_ADC[size])
                max_energy_trigger_ADC[size] = energy_trigger_ADC;

              if ((size == 2 or size == 4) and (binphi % 2 == 0)
                  and (bineta % 2 == 0))
                {
                  // sliding window made from 2x2 sums

                  if (slide2_energy_trigger_ADC
                      > slide2_max_energy_trigger_ADC[size])
                    slide2_max_energy_trigger_ADC[size] =
                        slide2_energy_trigger_ADC;
                }

            } //          for (int size = 1; size <= 4; ++size)
        }
    }

  for (int size = 1; size <= max_size; ++size)
    {
      max_energy_hist_list[size]->Fill(max_energy[size]);
      max_energy_trigger_ADC_hist_list[size]->Fill(
          max_energy_trigger_ADC[size]);

      if (size == 2 or size == 4)
        {
          // sliding window made from 2x2 sums
          slide2_max_energy_trigger_ADC_hist_list[size]->Fill(
              slide2_max_energy_trigger_ADC[size]);
        }
    }
  return Fun4AllReturnCodes::EVENT_OK;
}
