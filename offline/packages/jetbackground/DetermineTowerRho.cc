#include "DetermineTowerRho.h"

#include "TowerRhov1.h"

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetContainerv1.h>
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
#include <limits>
#include <memory>
#include <sstream>  // for basic_ostringstream
#include <string>
#include <vector>

DetermineTowerRho::DetermineTowerRho( const std::string &name )
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
  m_jet_output_nodes.clear();
  m_rho_methods.clear();
}

int DetermineTowerRho::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << "-----------------------------------" << std::endl;
    std::cout << "Methods: ";
    for (auto rho_method : m_rho_methods)
    {
      std::cout << TowerRhov1::get_method_string(rho_method) << ", ";
    }
    std::cout << std::endl;

    std::cout << "Inputs:";
    for (auto &input : m_inputs)
    {
      input->identify(std::cout);
    }

    std::cout << "Outputs: ";
    for (const auto &output : m_output_nodes)
    {
      std::cout << output << ", ";
    }
    std::cout << std::endl;

    std::cout << "Background jet nodes: ";
    for (const auto &jet_node : m_jet_output_nodes)
    {
      std::cout << jet_node << ", ";
    }
    std::cout << std::endl;

    std::cout << "Using jet algo: ";
    if (m_bkgd_jet_algo == Jet::ALGO::ANTIKT)
    {
      std::cout << "ANTIKT r=" << m_par;
    }
    else if (m_bkgd_jet_algo == Jet::ALGO::KT)
    {
      std::cout << "KT r=" << m_par;
    }
    else if (m_bkgd_jet_algo == Jet::ALGO::CAMBRIDGE)
    {
      std::cout << "CAMBRIDGE r=" << m_par;
    }
    std::cout << std::endl;

    std::cout << "Estimator: signed scalar sum of constituent pT (prop_signedSumeT)" << std::endl;

    std::cout << "Ghosts (AREA): event tower eta range, ghost area " << m_ghost_area << std::endl;
    std::cout << "Jet eta acceptance: ";
    if (m_abs_jet_eta_range < 0)
    {
      std::cout << "event tower eta range shrunk by R";
    }
    else
    {
      std::cout << "|eta| < " << m_abs_jet_eta_range;
    }
    std::cout << std::endl;
    std::cout << "Omit n hardest: " << m_omit_nhardest << std::endl;
    std::cout << "-----------------------------------" << std::endl;
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
    if (this_e == 0.0)
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
    pseudojet.set_user_index(ipart);
    calo_pseudojets.push_back(pseudojet);

  }  // end of loop over particles

  // eta range spanned by the towers in this event: the extent of the ghosts for the
  // AREA method and the base of the default jet acceptance. Flipping a negative
  // tower to +1 MeV keeps its direction, so its eta is unchanged.
  double eta_min = std::numeric_limits<double>::max();
  double eta_max = std::numeric_limits<double>::lowest();
  for (const auto &pseudojet : calo_pseudojets)
  {
    eta_min = std::min(eta_min, pseudojet.eta());
    eta_max = std::max(eta_max, pseudojet.eta());
  }
  const bool have_towers = !calo_pseudojets.empty() && eta_max > eta_min;

  // initialize the jet selector
  auto jet_selector = get_jet_selector(eta_min, eta_max);

  // assign the jet definition based on the selected algorithm
  fastjet::JetDefinition *m_jet_def = nullptr;
  if (m_bkgd_jet_algo == Jet::ALGO::ANTIKT)
  {
    m_jet_def = new fastjet::JetDefinition(
      fastjet::antikt_algorithm,
      m_par,
      fastjet::E_scheme,
      fastjet::Best
    );
  }
  else if (m_bkgd_jet_algo == Jet::ALGO::KT)
  {
    m_jet_def = new fastjet::JetDefinition(
      fastjet::kt_algorithm,
      m_par,
      fastjet::E_scheme,
      fastjet::Best
    );
  }
  else if (m_bkgd_jet_algo == Jet::ALGO::CAMBRIDGE)
  {
    m_jet_def = new fastjet::JetDefinition(
      fastjet::cambridge_algorithm,
      m_par,
      fastjet::E_scheme,
      fastjet::Best
    );
  }
  else
  {
    std::cout << PHWHERE << " jet algorithm not recognized, using default (kt)." << std::endl;
    m_jet_def = new fastjet::JetDefinition(
      fastjet::kt_algorithm,
      m_par,
      fastjet::E_scheme,
      fastjet::Best
    );
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

    const std::string &jet_node = m_jet_output_nodes.at(ipos);
    auto *jets = findNode::getClass<JetContainer>(topNode, jet_node);
    if (!jets)
    {
      std::cout << PHWHERE << " JetContainer node " << jet_node << " not found" << std::endl;
      continue;
    }

    if (!m_jet_def)
    {
      std::cerr << PHWHERE << " jet definition not set" << std::endl;
      exit(1);
    }

    // clear last event's jets; they stay empty (and rho 0) if there is nothing to cluster
    jets->Reset();
    jets->set_rho_median(0);

    // cluster the background jets with fastjet; ghosts (and so areas) are only
    // needed, and only available, for the area method
    std::unique_ptr<fastjet::ClusterSequence> cluseq{};
    if (!have_towers)
    {
      if (Verbosity() > 1)
      {
        std::cout << "DetermineTowerRho::process_event - no towers to cluster, rho = 0" << std::endl;
      }
    }
    else if (rho_method == TowerRho::Method::AREA)
    {
      fastjet::AreaDefinition const area_def(
        fastjet::active_area_explicit_ghosts,
        fastjet::GhostedAreaSpec(eta_min, eta_max, 1, m_ghost_area)
      );
      cluseq = std::make_unique<fastjet::ClusterSequenceArea>(calo_pseudojets, *m_jet_def, area_def);
    }
    else if (rho_method == TowerRho::Method::MULT)
    {
      cluseq = std::make_unique<fastjet::ClusterSequence>(calo_pseudojets, *m_jet_def);
    }
    else
    {
      std::cout << PHWHERE << " rho method not recognized" << std::endl;
    }

    if (cluseq)
    {
      // The eta acceptance is applied on the fastjet side: the fastjet
      // axis is well defined for every jet, whereas the sPHENIX four-vector is not
      // for pure-ghost jets (zero momentum, so eta is NaN) or for net-negative-energy
      // jets (whose eta and phi flip). This is also the axis the eta acceptance is
      // defined on.
      auto fastjets = jet_selector(cluseq->inclusive_jets());

      // Convert the selected jets back to sPHENIX Jets, with the negative-energy
      // towers reincorporated, straight into the saved container. Every jet is
      // kept there, seeds included; the hardest are omitted from rho by index.
      // ConvertJets also fills prop_signedSumeT on every jet.
      const bool with_area = (rho_method == TowerRho::Method::AREA);
      ConvertJets(jets, fastjets, particles, with_area);

      // omit the hardest jets, ranked by the same quantity the median is taken over
      const std::vector<unsigned int> keep = SelectJets(jets);

      // estimate rho from the sPHENIX jets that remain
      CalcRho(jets, keep, rho_method, rho, sigma);
      jets->set_rho_median(rho);

      // prop_SeedItr: 1 = omitted from rho as one of the hardest (a seed),
      //               0 = entered the rho estimate
      std::vector<bool> is_seed(jets->size(), true);
      for (const auto ijet : keep)
      {
        is_seed[ijet] = false;
      }

      const Jet::PROPERTY seed_idx = jets->property_index(Jet::PROPERTY::prop_SeedItr);
      for (unsigned int ijet = 0; ijet < jets->size(); ijet++)
      {
        jets->get_jet(ijet)->set_property(seed_idx, is_seed[ijet] ? 1.0F : 0.0F);
      }

      if (Verbosity() > 1)
      {
        std::cout << "DetermineTowerRho::process_event - wrote " << jets->size()
          << " background jets to " << jet_node << " (rho = " << rho << ")"
        << std::endl;
      }
    }  // end of if cluseq

    if (Verbosity() > 1)
    {
      std::cout << "DetermineTowerRho::process_event - Filling node "
        << m_output_nodes.at(ipos) << " with rho = " << rho
        << " and sigma = " << sigma
      << std::endl;
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

void DetermineTowerRho::add_method(
  const TowerRho::Method rho_method,
  const std::string &output_node,
  const std::string &jet_node
)
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

  // rho node: default TowerRho_<METHOD>
  const std::string rho_node = output_node.empty() ? "TowerRho_" + method_name : output_node;
  m_output_nodes.push_back(rho_node);

  // jet node: default <rho node>_SeedJets, so that every rho node has its own
  // jet node (e.g. TowerRho_CEMC_MULT -> TowerRho_CEMC_MULT_SeedJets)
  m_jet_output_nodes.push_back(jet_node.empty() ? rho_node + "_SeedJets" : jet_node);

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

  // create a JetContainer per method, holding the background jets that method used
  for (unsigned int ipos = 0; ipos < m_jet_output_nodes.size(); ipos++)
  {
    const std::string &jet_node = m_jet_output_nodes.at(ipos);
    auto *jetcont = findNode::getClass<JetContainer>(topNode, jet_node);
    if (!jetcont)
    {
      jetcont = new JetContainerv1();
      auto *jetDataNode = new PHIODataNode<PHObject>(jetcont, jet_node, "PHObject");
      bkgNode->addNode(jetDataNode);
    }
    else
    {
      std::cout << PHWHERE << " JetContainer node " << jet_node
                << " pre-exists, will be overwritten each event" << std::endl;
    }

    jetcont->set_algo(m_bkgd_jet_algo);
    jetcont->set_jetpar_R(m_par);
    jetcont->add_property(Jet::PROPERTY::prop_SeedItr);
    jetcont->add_property(Jet::PROPERTY::prop_signedSumeT);
    if (m_rho_methods.at(ipos) == TowerRho::Method::AREA)
    {
      jetcont->add_property(Jet::PROPERTY::prop_area);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void DetermineTowerRho::ConvertJets(JetContainer *jets,
                                    const std::vector<fastjet::PseudoJet> &fastjets,
                                    const std::vector<Jet *> &particles,
                                    const bool with_area) const
{
  jets->Reset();
  jets->set_algo(m_bkgd_jet_algo);
  jets->set_jetpar_R(m_par);
  for (const auto &input : m_inputs)
  {
    jets->insert_src(input->get_src());
  }

  const Jet::PROPERTY signed_idx = jets->property_index(Jet::PROPERTY::prop_signedSumeT);
  Jet::PROPERTY area_idx = Jet::PROPERTY::no_property;
  if (with_area)
  {
    area_idx = jets->property_index(Jet::PROPERTY::prop_area);
  }

  for (unsigned int ijet = 0; ijet < fastjets.size(); ijet++)
  {
    const auto &fj = fastjets.at(ijet);
    auto *jet = jets->add_jet();

    // sum the ORIGINAL, unflipped constituent momenta -- the negative-energy
    // rescaling in process_event is a clustering device only, so summing the
    // original towers is what reincorporates them (this mirrors FastJetAlgoSub)
    float total_px = 0;
    float total_py = 0;
    float total_pz = 0;
    float total_e = 0;
    // the signed scalar sum: a negative-energy tower carries a negative pT, so it
    // subtracts here instead of being folded up as in the jet pT magnitude
    float total_signed_pt = 0;
    for (auto &comp : fj.constituents())
    {
      // is_pure_ghost() requires area information, so only ask when we have it
      if (with_area && comp.is_pure_ghost())
      {
        continue;
      }
      auto *particle = particles.at(comp.user_index());
      total_px += particle->get_px();
      total_py += particle->get_py();
      total_pz += particle->get_pz();
      total_e += particle->get_e();
      // Jet::get_pt() is a magnitude, so take the sign from the energy: the inputs
      // build the momentum as pT = E / cosh(eta), which flips sign with E
      total_signed_pt += (particle->get_e() < 0) ? -particle->get_pt() : particle->get_pt();
      jet->insert_comp(particle->get_comp_vec(), true);
    }
    jet->set_comp_sort_flag();
    jet->set_property(signed_idx, total_signed_pt);

    jet->set_px(total_px);
    jet->set_py(total_py);
    jet->set_pz(total_pz);
    jet->set_e(total_e);
    jet->set_id(ijet);

    if (with_area)
    {
      jet->set_property(area_idx, static_cast<float>(fj.area()));
    }
  }
}

float DetermineTowerRho::SignedSumeT(JetContainer *jets, const unsigned int ijet)
{
  return jets->get_jet(ijet)->get_property(jets->property_index(Jet::PROPERTY::prop_signedSumeT));
}

std::vector<unsigned int> DetermineTowerRho::SelectJets(JetContainer *jets) const
{
  const unsigned int njets = jets->size();

  // rank the jets by prop_signedSumeT, hardest first; ties keep their original order
  std::vector<float> pts(njets);
  std::vector<unsigned int> order(njets);
  for (unsigned int ijet = 0; ijet < njets; ijet++)
  {
    pts[ijet] = SignedSumeT(jets, ijet);
    order[ijet] = ijet;
  }
  const unsigned int nomit = std::min(m_omit_nhardest, njets);
  std::partial_sort(order.begin(), order.begin() + nomit, order.end(),
    [&pts](const unsigned int a, const unsigned int b)
    { return pts[a] > pts[b] || (pts[a] == pts[b] && a < b); });

  // omit the nomit hardest; the rest keep their original order
  std::vector<bool> omit(njets, false);
  for (unsigned int iomit = 0; iomit < nomit; iomit++)
  {
    omit[order[iomit]] = true;
  }

  std::vector<unsigned int> keep{};
  keep.reserve(njets - nomit);
  for (unsigned int ijet = 0; ijet < njets; ijet++)
  {
    if (!omit[ijet])
    {
      keep.push_back(ijet);
    }
  }
  return keep;
}

void DetermineTowerRho::CalcRho(JetContainer *jets,
                                const std::vector<unsigned int> &keep,
                                const TowerRho::Method rho_method, float &rho, float &sigma) const
{
  rho = 0;
  sigma = 0;

  // no jets left after the acceptance and the omitted seeds: rho = sigma = 0
  if (keep.empty())
  {
    return;
  }

  if (rho_method == TowerRho::Method::AREA)
  {
    const Jet::PROPERTY area_idx = jets->property_index(Jet::PROPERTY::prop_area);

    std::vector<float> pT_over_X{};
    float total_X = 0;
    float njets_used = 0.0;
    float empty_X = 0;
    float const njets_total = static_cast<float>(keep.size());

    for (const auto ijet : keep)
    {
      auto *jet = jets->get_jet(ijet);

      float const this_X = jet->get_property(area_idx);
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

      // signed sum pT per unit area; the negative-energy towers are already back in
      float const this_pT_over_X = SignedSumeT(jets, ijet) / this_X;
      pT_over_X.push_back(this_pT_over_X);
      total_X += this_X;
      njets_used += 1.0;
    }  // end of loop over jets

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
  }
  else if (rho_method == TowerRho::Method::MULT)
  {
    std::vector<float> pt_over_nconst{};
    int total_constituents = 0;

    for (const auto ijet : keep)
    {
      auto *jet = jets->get_jet(ijet);

      // one component per clustered input, so this is the number of constituents
      size_t const nconst = jet->size_comp();
      if (nconst > 0)
      {
        float const jet_avg_pt = SignedSumeT(jets, ijet) / (1.0 * nconst);
        pt_over_nconst.push_back(jet_avg_pt);
        total_constituents += static_cast<int>(nconst);
      }  // end of if nconst > 0
    }  // end of loop over jets

    float const n_empty_jets = 1.0 * (keep.size() - pt_over_nconst.size());
    float mean_N = (1.0 * total_constituents) / (1.0 * keep.size());
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
  }
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

fastjet::Selector DetermineTowerRho::get_jet_selector(const double eta_min, const double eta_max) const
{
  return (m_abs_jet_eta_range < 0)
    ? fastjet::SelectorEtaRange(eta_min + m_par, eta_max - m_par)
    : fastjet::SelectorAbsEtaMax(m_abs_jet_eta_range);
}
