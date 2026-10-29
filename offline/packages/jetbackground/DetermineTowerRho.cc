#include "DetermineTowerRho.h"

#include "MbdChargeKnee.h"
#include "TowerRhov1.h"

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/JetInput.h>

#include <calobase/TowerInfoDefs.h>

#include <cdbobjects/CDBHistos.h>

#include <globalvertex/GlobalVertexMap.h>

#include <mbd/MbdOut.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TH1.h>
#include <TH3.h>

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
    std::cout << "Tower shape map: ";
    if (m_smap_path.empty() && m_smap_tag.empty())
    {
      std::cout << "off";
    }
    else
    {
      std::cout << (m_smap_path.empty() ? "CDB tag " + m_smap_tag : m_smap_path) << " (seed only)";
    }
    std::cout << std::endl;
    std::cout << "-----------------------------------" << std::endl;
  }

  if (LoadTowerShape() != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return CreateNodes(topNode);
}

int DetermineTowerRho::End(PHCompositeNode * /*topNode*/)
{
  if (m_smap_on)
  {
    std::cout << "DetermineTowerRho::End(" << Name() << ") - tower shape map: " << m_n_smap_events
              << " events, s = 1 without a vertex: " << m_n_smap_novtx
              << ", outside the map's vertex range: " << m_n_smap_out_of_range << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
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

  // tower shape map: s of every particle (all 1 without a map)
  std::vector<float> scales(particles.size(), 1.0F);
  if (m_smap_on)
  {
    const int iclass = TowerShapeClass(topNode);
    if (iclass < -1)
    {
      for (auto &part : particles)
      {
        delete part;
      }
      return Fun4AllReturnCodes::ABORTRUN;
    }
    scales = TowerShapeScales(particles, iclass);
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

    // with a tower shape map the towers are clustered divided by their s; ConvertJets
    // rebuilds the jets from the unscaled towers
    const float scale = scales[ipart];
    if (scale != 1.0F)
    {
      this_e /= scale;
      this_px /= scale;
      this_py /= scale;
      this_pz /= scale;
    }

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
      std::vector<float> rank_pt{};
      ConvertJets(jets, fastjets, particles, scales, with_area, rank_pt);

      // omit the hardest jets, ranked by the same quantity the median is taken over
      // (with a tower shape map: by that quantity summed over the towers divided by s)
      const std::vector<unsigned int> keep = SelectJets(rank_pt);

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
                                    const std::vector<float> &scales,
                                    const bool with_area,
                                    std::vector<float> &rank_pt) const
{
  rank_pt.assign(fastjets.size(), 0.0F);
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
    // the same sum over the towers divided by their s (identical without a map), used to rank
    float total_rank_pt = 0;
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
      const float signed_pt = (particle->get_e() < 0) ? -particle->get_pt() : particle->get_pt();
      total_signed_pt += signed_pt;
      total_rank_pt += signed_pt / scales.at(comp.user_index());
      jet->insert_comp(particle->get_comp_vec(), true);
    }
    jet->set_comp_sort_flag();
    jet->set_property(signed_idx, total_signed_pt);
    rank_pt[ijet] = total_rank_pt;

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

int DetermineTowerRho::TowerShapeLayer(const int src)
{
  switch (src)
  {
  case Jet::CEMC_TOWERINFO_RETOWER:
    return 0;
  case Jet::HCALIN_TOWERINFO:
    return 1;
  case Jet::HCALOUT_TOWERINFO:
    return 2;
  default:
    return -1;
  }
}

int DetermineTowerRho::LoadTowerShape()
{
  m_smap_on = false;
  for (auto &layer : m_smap)
  {
    layer.clear();
  }
  if (m_smap_path.empty() && m_smap_tag.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // the layers of the inputs; every input needs one
  std::array<bool, 3> need{};
  for (auto *input : m_inputs)
  {
    const int layer = TowerShapeLayer(input->get_src());
    if (layer < 0)
    {
      std::cout << PHWHERE << " tower shape map: input source " << input->get_src()
                << " is not a TOWERINFO input of the retowered EMCal, IHCal or OHCal" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    need[layer] = true;
  }

  const std::string url = m_smap_path.empty() ? CDBInterface::instance()->getUrl(m_smap_tag) : m_smap_path;
  if (url.empty())
  {
    std::cout << PHWHERE << " no tower shape map found for CDB tag " << m_smap_tag << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // the payload is only needed here: copy the map into flat vectors and release it
  CDBHistos cdb(url);
  cdb.LoadCalibrations();
  TH1 *hzvtx = cdb.getHisto("zvtx_bins");
  TH1 *hmbdq = cdb.getHisto("mbdq_bins");
  TH1 *hknee = cdb.getHisto("mbdq_knee", false);
  if (!hzvtx || !hmbdq)
  {
    std::cout << PHWHERE << " tower shape map " << url << ": zvtx_bins or mbdq_bins missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  const int n_zvtx = hzvtx->GetNbinsX();
  const int n_mbdq = hmbdq->GetNbinsX();
  m_smap_zvtx_edges.assign(n_zvtx + 1, 0);
  for (int i = 0; i <= n_zvtx; ++i)
  {
    m_smap_zvtx_edges[i] = static_cast<float>(hzvtx->GetXaxis()->GetBinLowEdge(i + 1));
  }
  m_smap_mbdq_edges.assign(n_mbdq + 1, 0);
  for (int i = 0; i <= n_mbdq; ++i)
  {
    m_smap_mbdq_edges[i] = static_cast<float>(hmbdq->GetXaxis()->GetBinLowEdge(i + 1));
  }
  for (const auto *edges : {&m_smap_zvtx_edges, &m_smap_mbdq_edges})
  {
    for (size_t i = 1; i < edges->size(); ++i)
    {
      if (!((*edges)[i] > (*edges)[i - 1]))
      {
        std::cout << PHWHERE << " tower shape map " << url << ": bin edges not increasing" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
    }
  }

  // Q_knee of this run
  const int run = recoConsts::instance()->get_IntFlag("RUNNUMBER");
  std::string knee_source;
  bool knee_default = false;
  m_smap_knee = MbdChargeKnee::find(m_smap_knee_user, m_knee_tag, m_knee_path, hknee, run, knee_source, knee_default);
  if (!(m_smap_knee > 0))
  {
    std::cout << PHWHERE << " tower shape map: no MBD charge knee for run " << run << ": " << knee_source << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (knee_default)
  {
    std::cout << "DetermineTowerRho::LoadTowerShape - WARNING run " << run << " Q_knee = " << m_smap_knee << " (" << knee_source << ")" << std::endl;
  }

  // s per layer: TH3 (ieta, iphi, class = izbin * n_mbdQ_bins + imbd)
  const std::array<std::string, 3> name = {"s_CEMC", "s_HCALIN", "s_HCALOUT"};
  const int n_class = n_zvtx * n_mbdq;
  m_smap_neta = 24;
  m_smap_nphi = 64;
  for (int layer = 0; layer < 3; ++layer)
  {
    if (!need[layer])
    {
      continue;
    }
    auto *hs = dynamic_cast<TH3 *>(cdb.getHisto(name[layer]));
    if (!hs || hs->GetNbinsX() != m_smap_neta || hs->GetNbinsY() != m_smap_nphi || hs->GetNbinsZ() != n_class)
    {
      std::cout << PHWHERE << " tower shape map " << url << ": " << name[layer] << " missing or not "
                << m_smap_neta << " x " << m_smap_nphi << " x " << n_class << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    m_smap[layer].assign(static_cast<size_t>(n_class) * m_smap_neta * m_smap_nphi, 1.0F);
    for (int iclass = 0; iclass < n_class; ++iclass)
    {
      for (int ieta = 0; ieta < m_smap_neta; ++ieta)
      {
        for (int iphi = 0; iphi < m_smap_nphi; ++iphi)
        {
          const auto s = static_cast<float>(hs->GetBinContent(ieta + 1, iphi + 1, iclass + 1));
          if (std::isfinite(s) && s > 0)
          {
            m_smap[layer][(static_cast<size_t>(iclass) * m_smap_neta + ieta) * m_smap_nphi + iphi] = s;
          }
        }
      }
    }
  }
  m_smap_on = true;

  if (Verbosity() > 0)
  {
    std::cout << "DetermineTowerRho::LoadTowerShape - " << url << ": " << n_zvtx << " z x " << n_mbdq
              << " Q' bins, run " << run << " Q_knee = " << m_smap_knee << " (" << knee_source << ")" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int DetermineTowerRho::TowerShapeClass(PHCompositeNode *topNode)
{
  ++m_n_smap_events;

  // event vertex z
  auto *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  GlobalVertex *vtx = nullptr;
  if (vertexmap && !vertexmap->empty())
  {
    vtx = vertexmap->begin()->second;
    if (m_use_vertex_type)
    {
      const auto vertices = vertexmap->get_gvtxs_with_type({m_vertex_type});
      vtx = vertices.empty() ? nullptr : vertices.at(0);
    }
  }
  if (!vtx || !std::isfinite(vtx->get_z()))
  {
    ++m_n_smap_novtx;
    return -1;
  }
  const float vtxz = vtx->get_z();
  const auto zit = std::upper_bound(m_smap_zvtx_edges.begin(), m_smap_zvtx_edges.end(), vtxz);
  const int izbin = static_cast<int>(zit - m_smap_zvtx_edges.begin()) - 1;
  if (izbin < 0 || izbin >= static_cast<int>(m_smap_zvtx_edges.size()) - 1)
  {
    ++m_n_smap_out_of_range;
    return -1;
  }

  // Q' = (q_S + q_N) / Q_knee, outside the map's range in the first / last bin
  auto *mbdout = findNode::getClass<MbdOut>(topNode, m_mbd_node);
  if (!mbdout)
  {
    std::cout << PHWHERE << " MBD node " << m_mbd_node << " needed for the tower shape map not found" << std::endl;
    return -2;
  }
  const float qprime = (mbdout->get_q(0) + mbdout->get_q(1)) / m_smap_knee;
  const int n_mbdq = static_cast<int>(m_smap_mbdq_edges.size()) - 1;
  const auto qit = std::upper_bound(m_smap_mbdq_edges.begin(), m_smap_mbdq_edges.end(), qprime);
  const int imbd = std::clamp(static_cast<int>(qit - m_smap_mbdq_edges.begin()) - 1, 0, n_mbdq - 1);

  return (izbin * n_mbdq) + imbd;
}

std::vector<float> DetermineTowerRho::TowerShapeScales(const std::vector<Jet *> &particles, const int iclass) const
{
  std::vector<float> scales(particles.size(), 1.0F);
  if (iclass < 0)
  {
    return scales;
  }
  const size_t offset = static_cast<size_t>(iclass) * m_smap_neta * m_smap_nphi;
  for (size_t ipart = 0; ipart < particles.size(); ++ipart)
  {
    const auto &comps = particles[ipart]->get_comp_vec();
    if (comps.empty())
    {
      continue;
    }
    const int layer = TowerShapeLayer(comps.front().first);
    if (layer < 0 || m_smap[layer].empty())
    {
      continue;
    }
    // TOWERINFO channel -> (ieta, iphi); the retowered EMCal is a copy of the HCal container
    const unsigned int key = TowerInfoDefs::encode_hcal(comps.front().second);
    const int ieta = static_cast<int>(TowerInfoDefs::getCaloTowerEtaBin(key));
    const int iphi = static_cast<int>(TowerInfoDefs::getCaloTowerPhiBin(key));
    if (ieta < 0 || ieta >= m_smap_neta || iphi < 0 || iphi >= m_smap_nphi)
    {
      continue;
    }
    scales[ipart] = m_smap[layer][offset + (ieta * m_smap_nphi) + iphi];
  }
  return scales;
}

std::vector<unsigned int> DetermineTowerRho::SelectJets(const std::vector<float> &rank_pt) const
{
  const auto njets = static_cast<unsigned int>(rank_pt.size());

  // rank the jets by rank_pt (prop_signedSumeT without a tower shape map), hardest first;
  // ties keep their original order
  std::vector<unsigned int> order(njets);
  for (unsigned int ijet = 0; ijet < njets; ijet++)
  {
    order[ijet] = ijet;
  }
  const unsigned int nomit = std::min(m_omit_nhardest, njets);
  std::partial_sort(order.begin(), order.begin() + nomit, order.end(),
    [&rank_pt](const unsigned int a, const unsigned int b)
    { return rank_pt[a] > rank_pt[b] || (rank_pt[a] == rank_pt[b] && a < b); });

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
