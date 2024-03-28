#include "DetermineTowerRho.h"

#include "TowerRho.h"
#include "TowerRhov1.h"

#include <jetbase/Jet.h>
#include <jetbase/JetInput.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// fastjet includes
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>

// standard includes
#include <algorithm>
#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>
#include <string>
#include <vector>

DetermineTowerRho::DetermineTowerRho(const std::string &name)
  : SubsysReco(name)
{
}

DetermineTowerRho::~DetermineTowerRho()
{
  for (auto &_input : _inputs)
  {
    delete _input;
  }
  _inputs.clear();
  _outputs.clear();
}

int DetermineTowerRho::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    print_settings(std::cout);
  }
  m_fj_algo_name = get_fj_algo_name();
  return CreateNode(topNode);
}

int DetermineTowerRho::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "DetermineTowerRho::process_event -- entered" << std::endl;
  }

  if (m_abs_jet_eta_range == 5.0)
  {
    m_abs_jet_eta_range = m_abs_tower_eta_range - m_par;
  }
  if (m_abs_ghost_eta == 5.0)
  {
    m_abs_ghost_eta = m_abs_tower_eta_range;
  }

  // get fastjet input from tower nodes
  std::vector<fastjet::PseudoJet> calo_pseudojets = get_pseudojets(topNode);

  // fastjet definitions
  fastjet::JetAlgorithm _fj_algo = get_fj_algo();  // get fastjet algorithm

  // jet definition
  fastjet::JetDefinition jet_def(_fj_algo,
                                 m_par,
                                 fastjet::E_scheme,
                                 fastjet::Best);

  // area definition
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,
                                   fastjet::GhostedAreaSpec(m_abs_ghost_eta, 1, m_ghost_area));

  fastjet::Selector jet_selector = ((!fastjet::SelectorNHardest(m_omit_nhardest)) *
                                    (fastjet::SelectorAbsEtaMax(m_abs_jet_eta_range) && fastjet::SelectorPtMin(m_jet_min_pT)));

  // Not pure ghost function
  fastjet::Selector not_pure_ghost = (!fastjet::SelectorIsPureGhost());

  for (unsigned int ipos = 0; ipos < _rho_methods.size(); ipos++)
  {
    TowerRho::Method _rho_method = _rho_methods.at(ipos);
    float rho = 0;
    float sigma = 0;

    if (_rho_method == TowerRho::Method::AREA)
    {
      fastjet::JetMedianBackgroundEstimator bge{jet_selector, jet_def, area_def};
      bge.set_particles(calo_pseudojets);
      rho = bge.rho();
      sigma = bge.sigma();
    }
    else if (_rho_method == TowerRho::Method::MULT)
    {
      // reconstruct the background jets
      fastjet::ClusterSequenceArea cs(calo_pseudojets, jet_def, area_def);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(jet_selector(cs.inclusive_jets()));

      std::vector<float> pt_over_nConstituents;
      int nfj_jets = 0;
      int n_empty_jets = 0;
      float total_constituents = 0;
      for (auto &jet : jets)
      {
        float jet_pt = jet.pt();

        std::vector<fastjet::PseudoJet> constituents = not_pure_ghost(jet.constituents());

        int jet_size = constituents.size();

        total_constituents += jet_size;
        if (jet_size == 0)
        {
          n_empty_jets++;
          continue;  // only consider jets with constituents
        }
        pt_over_nConstituents.push_back(jet_pt / jet_size);
        nfj_jets++;
      }

      int total_jets = nfj_jets + n_empty_jets;
      float mean_N = total_constituents / total_jets;
      if (mean_N < 0)
      {
        std::cout << "WARNING: mean_N = " << mean_N << " is negative, setting to 0" << std::endl;
        mean_N = 0;
      }

      float med_tmp, std_tmp;
      _median_stddev(pt_over_nConstituents, n_empty_jets, med_tmp, std_tmp);
      sigma = float(std_tmp * std::sqrt(mean_N));
      rho = med_tmp;
    }
    else
    {
      rho = 0;
      sigma = 0;
    }

    FillNode(topNode, ipos, rho, sigma);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<fastjet::PseudoJet> DetermineTowerRho::get_pseudojets(PHCompositeNode *topNode)
{
  // get fastjet input from tower nodes
  std::vector<fastjet::PseudoJet> calo_pseudojets;
  for (auto &_input : _inputs)
  {
    std::vector<Jet *> parts = _input->get_input(topNode);
    // for (unsigned int i = 0; i < parts.size(); ++i)
    // {
    for (auto &part : parts)
    {
      // auto& part = parts[i];
      float this_e = part->get_e();
      if (this_e == 0.)
      {
        continue;
      }
      float this_px = part->get_px();
      float this_py = part->get_py();
      float this_pz = part->get_pz();

      // float this_px = parts[i]->get_px();
      // float this_py = parts[i]->get_py();
      // float this_pz = parts[i]->get_pz();

      // if(m_do_tower_cut)
      // {
      // if(this_e < m_tower_threshold) continue;
      if (this_e < 0)
      {
        // make energy = +1 MeV for purposes of clustering
        float e_ratio = 0.001 / this_e;
        this_e = this_e * e_ratio;
        this_px = this_px * e_ratio;
        this_py = this_py * e_ratio;
        this_pz = this_pz * e_ratio;
      }

      // }
      fastjet::PseudoJet pseudojet(this_px, this_py, this_pz, this_e);

      if (pseudojet.perp() < m_tower_min_pT)
      {
        continue;
      }
      if (fabs(pseudojet.eta()) > m_abs_tower_eta_range)
      {
        continue;
      }

      calo_pseudojets.push_back(pseudojet);
    }
    for (auto &p : parts)
    {
      delete p;
    }
    parts.clear();
  }

  return calo_pseudojets;
}

int DetermineTowerRho::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the jet background stuff under a sub-node directory
  PHCompositeNode *bkgNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "JETBACKGROUND"));
  if (!bkgNode)
  {
    bkgNode = new PHCompositeNode("JETBACKGROUND");
    dstNode->addNode(bkgNode);
  }

  // create the TowerRho nodes
  for (auto &_output : _outputs)
  {
    TowerRho *towerrho = findNode::getClass<TowerRho>(topNode, _output);
    if (!towerrho)
    {
      towerrho = new TowerRhov1();
      PHIODataNode<PHObject> *rhoDataNode = new PHIODataNode<PHObject>(towerrho, _output, "PHObject");
      bkgNode->addNode(rhoDataNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void DetermineTowerRho::FillNode(PHCompositeNode *topNode, unsigned int ipos, const float rho, const float sigma)
{
  std::string _node_name = _outputs.at(ipos);
  TowerRho::Method rho_method = _rho_methods.at(ipos);
  TowerRho *towerbackground = findNode::getClass<TowerRho>(topNode, _node_name);
  if (!towerbackground)
  {
    std::cout << " ERROR -- can't find TowerRho node after it should have been created" << std::endl;
    return;
  }
  else
  {
    if (Verbosity() > 0)
    {
      std::cout << "DetermineTowerRho::FillNode - Filling node " << _node_name << " with rho = " << rho << " and sigma = " << sigma << std::endl;
    }
    towerbackground->set_rho(rho);
    towerbackground->set_sigma(sigma);
    towerbackground->set_method(rho_method);
  }

  return;
}

float DetermineTowerRho::_percentile(const std::vector<float> &sorted_vec,
                                     const float percentile,
                                     const float nempty) const
{
  assert(percentile >= 0. && percentile <= 1.);

  int njets = sorted_vec.size();
  if (njets == 0)
  {
    return 0;
  }

  float total_njets = njets + nempty;
  float perc_pos = (total_njets) *percentile - nempty - 0.5;

  float result;
  if (perc_pos >= 0 && njets > 1)
  {
    int pindex = int(perc_pos);
    // avoid out of range
    if (pindex + 1 > njets - 1)
    {
      pindex = njets - 2;
      perc_pos = njets - 1;
    }

    result = sorted_vec.at(pindex) * (pindex + 1 - perc_pos) + sorted_vec.at(pindex + 1) * (perc_pos - pindex);
  }
  else if (perc_pos > -0.5 && njets >= 1)
  {
    result = sorted_vec.at(0);
  }
  else
  {
    result = 0;
  }

  return result;
}

void DetermineTowerRho::_median_stddev(const std::vector<float> &vec,
                                       float n_empty_jets,
                                       float &median,
                                       float &std_dev) const
{
  if (vec.size() == 0)
  {
    median = 0;
    std_dev = 0;
    return;
  }

  std::vector<float> sorted_vec = vec;
  std::sort(sorted_vec.begin(), sorted_vec.end());

  int njets = sorted_vec.size();
  if (n_empty_jets > njets / 4.0)
  {
    std::cout << "WARNING: n_empty_jets = " << n_empty_jets << " is too large, setting to " << njets / 4.0 << std::endl;
    n_empty_jets = njets / 4.0;
  }

  float posn[2] = {0.5, (1.0 - 0.6827) / 2.0};
  float res[2];
  for (int i = 0; i < 2; i++)
  {
    res[i] = _percentile(sorted_vec, posn[i], n_empty_jets);
  }

  median = res[0];
  std_dev = res[0] - res[1];
}

fastjet::JetAlgorithm DetermineTowerRho::get_fj_algo()
{
  if (_algo == Jet::ANTIKT)
  {
    return fastjet::antikt_algorithm;
  }
  else if (_algo == Jet::KT)
  {
    return fastjet::kt_algorithm;
  }
  else if (_algo == Jet::CAMBRIDGE)
  {
    return fastjet::cambridge_algorithm;
  }
  else
  {
    std::cout << PHWHERE << " jet algorithm not recognized, using default (kt)." << std::endl;
    return fastjet::kt_algorithm;
  }
}

std::string DetermineTowerRho::get_fj_algo_name()
{
  if (_algo == Jet::ANTIKT)
  {
    return "antikt_algorithm";
  }
  else if (_algo == Jet::KT)
  {
    return "kt_algorithm";
  }
  else if (_algo == Jet::CAMBRIDGE)
  {
    return "cambridge_algorithm";
  }
  else
  {
    std::cout << PHWHERE << " jet algorithm not recognized, using default (kt)." << std::endl;
    return "kt_algorithm";
  }
}

void DetermineTowerRho::print_settings(std::ostream &os) const
{
  os << "DetermineTowerRho settings: " << std::endl;
  os << "Methods: ";
  for (auto &_rho_method : _rho_methods)
  {
    os << TowerRhov1::get_method_string(_rho_method) << ", ";
  }
  os << std::endl;
  os << "Inputs:";
  for (auto &_input : _inputs)
  {
    _input->identify();
  }
  os << "Outputs: ";
  for (auto &_output : _outputs)
  {
    os << _output << ", ";
  }
  os << std::endl;
  os << "Using jet algo: " << m_fj_algo_name << " with R = " << m_par << std::endl;
  os << "Tower eta range: " << m_abs_tower_eta_range << std::endl;
  os << "Jet eta range: " << m_abs_jet_eta_range << std::endl;
  os << "Ghost area: " << m_ghost_area << std::endl;
  os << "Omit n hardest: " << m_omit_nhardest << std::endl;
  os << "Tower min pT: " << m_tower_min_pT << std::endl;
  os << "Jet min pT: " << m_jet_min_pT << std::endl;
  os << "Ghost eta: " << m_abs_ghost_eta << std::endl;
  os << "-----------------------------------" << std::endl;
  return;
}
