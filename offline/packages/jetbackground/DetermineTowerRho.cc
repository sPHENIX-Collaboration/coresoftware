#include "DetermineTowerRho.h"

#include "TowerRhov1.h"

#include <jetbase/JetInput.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// fastjet includes
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>  // for GhostedAreaSpec
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>

// standard includes
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>  // for basic_ostringstream
#include <string>
#include <vector>

DetermineTowerRho::DetermineTowerRho(const std::string &name)
  : SubsysReco(name)
{
  // silence output from fastjet
  fastjet::ClusterSequence const clusseq;
  if (Verbosity() > 0)
  {
    fastjet::ClusterSequence::print_banner();
  }
  else
  {
    std::ostringstream nullstream;
    fastjet::ClusterSequence::set_fastjet_banner_stream(&nullstream);
    fastjet::ClusterSequence::print_banner();
    fastjet::ClusterSequence::set_fastjet_banner_stream(&std::cout);
  }
}

DetermineTowerRho::~DetermineTowerRho()
{
  // clean up memory
  for (auto &input : m_inputs)
  {
    delete input;
  }
  m_inputs.clear();
  m_output_nodes.clear();
  m_rho_methods.clear();
}

int DetermineTowerRho::InitRun(PHCompositeNode *topNode)
{
  // set algo

  // set jet_eta_range
  if (m_abs_jet_eta_range < 0)
  {
    m_abs_jet_eta_range = m_abs_input_eta_range - m_par;
  }
  if (Verbosity() > 0)
  {
    print_settings();
  }

  return CreateNodes(topNode);
}

int DetermineTowerRho::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "DetermineTowerRho::process_event -- entered" << std::endl;
  }

  std::vector<Jet *> particles{};
  for (auto &input : m_inputs)
  {
    std::vector<Jet *> const parts = input->get_input(topNode);
    for (const auto &part : parts)
    {
      particles.push_back(part);
      particles.back()->set_id(particles.size() - 1);  // unique ids ensured
    }
  }

  std::vector<fastjet::PseudoJet> calo_pseudojets{};
  for (unsigned int ipart = 0; ipart < particles.size(); ++ipart)
  {
    float this_e = particles[ipart]->get_e();
    if (this_e == 0.)
    {
      continue;
    }  // skip zero energy particles
    float this_px = particles[ipart]->get_px();
    float this_py = particles[ipart]->get_py();
    float this_pz = particles[ipart]->get_pz();

    if (this_e < 0)
    {  // make energy = +1 MeV for purposes of clustering
      float const e_ratio = 0.001 / this_e;
      this_e = this_e * e_ratio;
      this_px = this_px * e_ratio;
      this_py = this_py * e_ratio;
      this_pz = this_pz * e_ratio;
    }

    fastjet::PseudoJet pseudojet(this_px, this_py, this_pz, this_e);
    // eta cut
    if (std::abs(pseudojet.eta()) > m_abs_input_eta_range)
    {
      continue;
    }

    pseudojet.set_user_index(ipart);
    calo_pseudojets.push_back(pseudojet);

  }  // end of loop over particles

  // initialize the jet selector
  auto jet_selector = get_jet_selector();
  fastjet::JetDefinition *m_jet_def = nullptr;
  if (m_bkgd_jet_algo == Jet::ALGO::ANTIKT)
  {
    m_jet_def = new fastjet::JetDefinition(fastjet::antikt_algorithm, m_par, fastjet::E_scheme, fastjet::Best);
  }
  else if (m_bkgd_jet_algo == Jet::ALGO::KT)
  {
    m_jet_def = new fastjet::JetDefinition(fastjet::kt_algorithm, m_par, fastjet::E_scheme, fastjet::Best);
  }
  else if (m_bkgd_jet_algo == Jet::ALGO::CAMBRIDGE)
  {
    m_jet_def = new fastjet::JetDefinition(fastjet::cambridge_algorithm, m_par, fastjet::E_scheme, fastjet::Best);
  }
  else
  {
    std::cout << PHWHERE << " jet algorithm not recognized, using default (kt)." << std::endl;
    m_jet_def = new fastjet::JetDefinition(fastjet::kt_algorithm, m_par, fastjet::E_scheme, fastjet::Best);
  }
  for (unsigned int ipos = 0; ipos < m_rho_methods.size(); ipos++)
  {
    float rho = 0;
    float sigma = 0;
    auto rho_method = m_rho_methods.at(ipos);

    auto *m_eventbackground = findNode::getClass<TowerRho>(topNode, m_output_nodes.at(ipos));
    if (!m_eventbackground)
    {
      std::cout << PHWHERE << " TowerRho node " << m_output_nodes.at(ipos) << " not found" << std::endl;
      continue;
    }

    if (!m_jet_def)
    {
      std::cerr << PHWHERE << " jet definition not set" << std::endl;
      exit(1);
    }

    if (rho_method == TowerRho::Method::AREA)
    {
      fastjet::AreaDefinition const area_def(fastjet::active_area_explicit_ghosts,
                                             fastjet::GhostedAreaSpec(m_abs_input_eta_range, 1, m_ghost_area));
      auto *m_cluseq = new fastjet::ClusterSequenceArea(calo_pseudojets, *m_jet_def, area_def);
      auto fastjets = jet_selector(m_cluseq->inclusive_jets());

      std::vector<float> pT_over_X{};
      float total_X = 0;
      float njets_used = 0.0;
      float empty_X = 0;
      float const njets_total = static_cast<float>(fastjets.size());

      for (auto &fastjet : fastjets)
      {
        float const this_X = fastjet.area();
        // if (this_X <= 0 || this_X != this_X || fastjet.is_pure_ghost())
        if (this_X <= 0 || this_X != this_X)
        {
          if (Verbosity() > 2)
          {
            std::cout << PHWHERE << " ::WARNING: Discarding jet with zero area. Zero-area jets may be due to (i) too large a ghost area (ii) a jet being outside the ghost range (iii) the computation not being done using an appropriate algorithm (kt;C/A)." << std::endl;
          }
          if (!std::isnan(this_X))
          {
            empty_X += this_X;
          }
          continue;  // skip this jet
        }  // end of check on X

        // add back the truth comp pT
        float px_sum = 0;
        float py_sum = 0;
        for (auto &comp : fastjet.constituents())
        {
          if (comp.is_pure_ghost())
          {  // skip pure ghosts
            continue;
          }

          auto *p = particles[comp.user_index()];  // get original particle
          px_sum += p->get_px();
          py_sum += p->get_py();

        }  // end of loop over fastjet comps

        // calculate pT and pT/X
        float const this_pT_over_X = std::sqrt((px_sum * px_sum) + (py_sum * py_sum)) / this_X;
        pT_over_X.push_back(this_pT_over_X);
        total_X += this_X;
        njets_used += 1.0;
      }  // end of loop over fastjets

      if (empty_X != 0.0)
      {
        if (Verbosity() > 0)
        {
          std::cerr << PHWHERE << " ::WARNING: Found " << empty_X << " empty jets with zero area. This may be due to (i) too large a ghost area (ii) a jet being outside the ghost range (iii) the computation not being done using an appropriate algorithm (kt;C/A)." << std::endl;
        }
        total_X += empty_X;
      }

      float const n_empty_jets = njets_total - njets_used;
      float mean_X = (1.0 * total_X) / (njets_total);
      if (mean_X < 0)
      {
        std::cerr << PHWHERE << " mean_N < 0 , setting to 0" << std::endl;
        mean_X = 0;
      }

      float tmp_med;
      float tmp_std;
      CalcMedianStd(pT_over_X, n_empty_jets, tmp_med, tmp_std);

      sigma = std::sqrt(mean_X) * tmp_std;
      rho = tmp_med;

      // clean up
      fastjets.clear();
      pT_over_X.clear();
      delete m_cluseq;
    }
    else if (rho_method == TowerRho::Method::MULT)
    {
      // reconstruct the background jets
      auto *m_cluseq = new fastjet::ClusterSequence(calo_pseudojets, *m_jet_def);
      auto fastjets = jet_selector(m_cluseq->inclusive_jets());

      std::vector<float> pt_over_nconst{};
      int total_constituents = 0;

      for (auto &fastjet : fastjets)
      {
        auto comps = fastjet.constituents();
        if (!comps.empty())
        {
          float total_px = 0;
          float total_py = 0;
          for (auto &comp : comps)
          {
            auto *particle = particles[comp.user_index()];
            total_px += particle->get_px();
            total_py += particle->get_py();
            total_constituents++;
          }  // end of loop over constituents
          float const jet_avg_pt = (std::sqrt((total_px * total_px) + (total_py * total_py)) / (1.0 * comps.size()));
          pt_over_nconst.push_back(jet_avg_pt);

        }  // end of if comps.size() > 0
      }  // end of loop over fastjets

      // {
      //   // auto comps = fastjets[ijet].constituents();
      //   if (comps.size() > 0)
      //   {
      //     float total_px = 0;
      //     float total_py = 0;
      //     for (auto &comp : comps)
      //     {
      //       auto particle = particles[comp.user_index()];
      //       total_px += particle->get_px();
      //       total_py += particle->get_py();
      //       total_constituents++;
      //     }  // end of loop over constituents
      //     float jet_avg_pt = (std::sqrt((total_px * total_px) + (total_py * total_py)) / (1.0 * comps.size()));
      //     pt_over_nconst.push_back(jet_avg_pt);

      //   }  // end of if comps.size() > 0
      // }    // end of loop over fastjets

      float const n_empty_jets = 1.0 * (fastjets.size() - pt_over_nconst.size());
      float mean_N = (1.0 * total_constituents) / (1.0 * fastjets.size());
      if (mean_N < 0)
      {
        std::cerr << PHWHERE << " mean_N < 0 , setting to 0" << std::endl;
        mean_N = 0;
      }

      float tmp_med;
      float tmp_std;
      CalcMedianStd(pt_over_nconst, 1.0 * n_empty_jets, tmp_med, tmp_std);

      sigma = std::sqrt(mean_N) * tmp_std;
      rho = tmp_med;

      // clean up
      fastjets.clear();
      pt_over_nconst.clear();
      delete m_cluseq;
    }
    else
    {
      std::cout << PHWHERE << " rho method not recognized" << std::endl;
    }

    if (Verbosity() > 1)
    {
      std::cout << "DetermineTowerRho::process_event - Filling node "
                << m_output_nodes.at(ipos) << " with rho = " << rho
                << " and sigma = " << sigma << std::endl;
    }

    m_eventbackground->set_rho(rho);
    m_eventbackground->set_sigma(sigma);
    m_eventbackground->set_method(rho_method);

  }  // end of loop over rho methods

  delete m_jet_def;
  // clean up input vectors
  for (auto &part : particles)
  {
    delete part;
  }
  particles.clear();
  calo_pseudojets.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

void DetermineTowerRho::add_method(TowerRho::Method rho_method, std::string output_node)
{
  // get method name ( also checks if method is valid )
  std::string const method_name = TowerRhov1::get_method_string(rho_method);

  // check if method already exists
  if (std::find(m_rho_methods.begin(), m_rho_methods.end(), rho_method) != m_rho_methods.end())
  {
    std::cout << PHWHERE << " method " << method_name << " already exists, skipping" << std::endl;
    return;
  }

  m_rho_methods.push_back(rho_method);

  // if no output name is specified, use default
  if (output_node.empty())
  {
    output_node = "TowerRho_" + method_name;
  }
  m_output_nodes.push_back(output_node);
  return;
}

int DetermineTowerRho::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);  // Looking for the DST node
  auto *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  auto *bkgNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "JETBACKGROUND"));
  if (!bkgNode)
  {  // create the node if it does not exist
    bkgNode = new PHCompositeNode("JETBACKGROUND");
    dstNode->addNode(bkgNode);
  }

  // create the TowerRho nodes
  for (auto &output : m_output_nodes)
  {
    auto *rho = findNode::getClass<TowerRho>(topNode, output);
    if (!rho)
    {
      rho = new TowerRhov1();
      auto *rhoDataNode = new PHIODataNode<PHObject>(rho, output, "PHObject");
      bkgNode->addNode(rhoDataNode);
    }  // end of if TowerRho
  }  // end of loop over output nodes

  return Fun4AllReturnCodes::EVENT_OK;
}

float DetermineTowerRho::CalcPercentile(const std::vector<float> &sorted_vec, const float percentile, const float nempty)
{
  float result = 0;
  if (!sorted_vec.empty())
  {
    int const njets = sorted_vec.size();
    float const total_njets = njets + nempty;
    float perc_pos = ((total_njets) *percentile) - nempty - 0.5;

    if (perc_pos >= 0 && njets > 1)
    {
      int pindex = int(perc_pos);
      if (pindex + 1 > njets - 1)
      {  // avoid out of range
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
    }  // end of if criteria

  }  // end of if sorted_vec.size() > 0

  return result;
}

void DetermineTowerRho::CalcMedianStd(const std::vector<float> &vec, float n_empty_jets, float &median, float &std_dev)
{
  median = 0;
  std_dev = 0;
  if (!vec.empty())
  {
    // sort the vector
    std::vector<float> sorted_vec = vec;
    std::sort(sorted_vec.begin(), sorted_vec.end());

    int const njets = sorted_vec.size();
    if (n_empty_jets > njets / 4.0)
    {
      std::cout << "WARNING: n_empty_jets = " << n_empty_jets << " is too large, setting to " << njets / 4.0 << std::endl;
      n_empty_jets = njets / 4.0;
    }

    float const posn[2] = {0.5, (1.0 - 0.6827) / 2.0};
    float res[2] = {0, 0};
    for (int i = 0; i < 2; i++)
    {
      res[i] = CalcPercentile(sorted_vec, posn[i], n_empty_jets);
    }

    median = res[0];
    std_dev = res[0] - res[1];
    sorted_vec.clear();
  }  // end of if vec.size() > 0

  return;
}

fastjet::Selector DetermineTowerRho::get_jet_selector() const
{
  if (m_jet_min_pT != VOID_CUT)
  {
    return (!fastjet::SelectorNHardest(m_omit_nhardest)) * (fastjet::SelectorAbsEtaMax(m_abs_jet_eta_range)) && (fastjet::SelectorPtMin(m_jet_min_pT));
  }
  // default is no min jet pT
  return (!fastjet::SelectorNHardest(m_omit_nhardest)) * (fastjet::SelectorAbsEtaMax(m_abs_jet_eta_range));
}

void DetermineTowerRho::print_settings(std::ostream &os)
{
  os << PHWHERE << "-----------------------------------" << std::endl;
  os << "Methods: ";
  for (auto rho_method : m_rho_methods)
  {
    os << TowerRhov1::get_method_string(rho_method) << ", ";
  }
  os << std::endl;

  os << "Inputs:";
  for (auto &input : m_inputs)
  {
    input->identify();
  }

  os << "Outputs: ";
  for (const auto &output : m_output_nodes)
  {
    os << output << ", ";
  }
  os << std::endl;

  os << "Using jet algo: ";
  if (m_bkgd_jet_algo == Jet::ALGO::ANTIKT)
  {
    os << "ANTIKT r=" << m_par;
  }
  else if (m_bkgd_jet_algo == Jet::ALGO::KT)
  {
    os << "KT r=" << m_par;
  }
  else if (m_bkgd_jet_algo == Jet::ALGO::CAMBRIDGE)
  {
    os << "CAMBRIDGE r=" << m_par;
  }
  os << std::endl;

  os << "Tower eta range: " << m_abs_input_eta_range << std::endl;
  os << "Jet eta range: " << m_abs_jet_eta_range << std::endl;
  os << "Ghost area: " << m_ghost_area << std::endl;
  os << "Ghost eta: " << m_abs_input_eta_range << std::endl;
  os << "Omit n hardest: " << m_omit_nhardest << std::endl;
  if (m_jet_min_pT != VOID_CUT)
  {
    os << "Jet min pT: " << m_jet_min_pT << std::endl;
  }
  os << "-----------------------------------" << std::endl;
  return;
}
