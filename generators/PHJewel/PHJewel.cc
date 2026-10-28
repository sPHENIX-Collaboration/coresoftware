#include "PHJewel.h"
#include "JewelCore.h"
#include "JewelRun.h"
#include "PHJewelGenTrigger.h"

#include <phhepmc/PHGenIntegral.h>
#include <phhepmc/PHGenIntegralv1.h>
#include <phhepmc/PHHepMCGenHelper.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <HepMC/GenEvent.h>
#include <HepMC/HeavyIon.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/PdfInfo.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

bool PHJewel::s_Instance = false;

PHJewel::PHJewel(const std::string &name)
  : SubsysReco(name)
{
  PHHepMCGenHelper::set_embedding_id(1);  // same default as PHPythia8
}

PHJewel::~PHJewel() = default;

int PHJewel::Init(PHCompositeNode *topNode)
{
  if (s_Instance)
  {
    std::cout << PHWHERE << " only one PHJewel per process: JEWEL keeps its state in Fortran common blocks" << std::endl;
    exit(EXIT_FAILURE);
  }
  s_Instance = true;
  if (m_Card.empty() || m_Medium.empty())
  {
    std::cout << PHWHERE << " set_card() and set_medium() are both required" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (create_node_tree(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    // Init cannot stop the run: a return code here only drops this module
    exit(EXIT_FAILURE);
  }

  int seed = m_Seed;
  if (seed <= 0)
  {
    // PYTHIA 6 MRPY(1) range is 0 <= seed < 900000000
    seed = static_cast<int>(PHRandomSeed() % 900000000U);
    if (seed <= 0)
    {
      seed = 1;
    }
  }
  else if (seed >= 900000000)
  {
    std::cout << PHWHERE << " seed " << seed << " outside the PYTHIA 6 range 0 < seed < 900000000" << std::endl;
    exit(EXIT_FAILURE);
  }
  m_Seed = seed;
  if (m_Stem.empty())
  {
    m_Stem = "PHJewel_" + std::to_string(seed);
  }

  JewelRunConfig c;
  c.card = m_Card;
  c.medium = m_Medium;
  c.stem = m_Stem;
  c.seed = seed;
  c.coreLib = m_CoreLib;
  c.fortranHepMC = m_FortranHepMC;
  c.evinfo = m_EvInfo;
  c.allowTableGeneration = m_AllowTables;
  c.isospinBlock = m_IsospinBlock;

  m_Core = std::make_unique<JewelCore>();
  try
  {
    const JewelCard::Result r = startJewel(*m_Core, c, std::cout);
    m_Hash = r.hash;
  }
  catch (const std::exception &e)
  {
    // Fun4All drops a module whose Init fails and carries on writing empty
    // events; a generator that cannot start must kill the job.
    std::cout << PHWHERE << "\n"
              << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

  if (m_IntegralNode)
  {
    m_IntegralNode->set_Description("PHJewel JEWEL-2.6.0 medium=" + m_Medium + " card=" + m_Card +
                                    " hash=" + m_Hash + " seed=" + std::to_string(seed) +
                                    " embedding_id=" + std::to_string(PHHepMCGenHelper::get_embedding_id()));
  }
  if (Verbosity() > 0)
  {
    Print("CONFIG");
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHJewel::process_event(PHCompositeNode * /*topNode*/)
{
  HepMC::GenEvent *genevent = nullptr;
  long attempts = 0;
  while (!genevent)
  {
    try
    {
      m_Core->next();
    }
    catch (const std::exception &e)
    {
      std::cout << PHWHERE << " JEWEL failed during event generation: " << e.what() << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    const JewelCore::EventInfo &info = m_Core->info();
    ++m_NGenerated;
    m_SumWAll += info.weight;

    // JEWEL's own HepMC2 lines, read with HepMC's own reader: exactly what
    // Fun4AllHepMCInputManager would make of a JEWEL file
    std::istringstream is(m_Core->event_stream());
    HepMC::IO_GenEvent reader(is);
    HepMC::GenEvent *cand = reader.read_next_event();
    if (!cand)
    {
      std::cout << PHWHERE << " HepMC could not parse JEWEL event " << info.ngood << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    // JEWEL writes the centrality into the impact-parameter slot; sPHENIX code
    // reads b in fm, so put b there (Psi_RP = 0: JEWEL lays b along +x).
    // The hard-scattering info JEWEL leaves at zero goes into PdfInfo.
    if (cand->heavy_ion())
    {
      HepMC::HeavyIon hi(*cand->heavy_ion());
      hi.set_impact_parameter(static_cast<float>(info.b));
      hi.set_event_plane_angle(0.F);
      cand->set_heavy_ion(hi);
    }
    cand->set_signal_process_id(info.isub);
    cand->set_pdf_info(HepMC::PdfInfo(info.id1, info.id2, info.x1, info.x2, std::sqrt(std::max(0.0, info.Q2)), 0., 0.));

    // registered triggers: OR (default) or AND; no trigger = every event
    bool pass = m_RegisteredTriggers.empty() || !m_TriggersOR;
    for (PHJewelGenTrigger *trigger : m_RegisteredTriggers)
    {
      const bool r = trigger->Apply(cand);
      if (m_TriggersOR && r)
      {
        pass = true;
        break;
      }
      if (!m_TriggersOR && !r)
      {
        pass = false;
        break;
      }
    }
    if (Verbosity() > 1)
    {
      std::cout << Name() << ": JEWEL event " << info.ngood << " channel " << info.channel << " w=" << info.weight
                << " pass=" << pass << std::endl;
    }
    if (pass)
    {
      genevent = cand;
      m_SumWPass += info.weight;
    }
    else
    {
      delete cand;
      ++attempts;
      if (m_MaxAttempts > 0 && attempts >= m_MaxAttempts)
      {
        std::cout << PHWHERE << " no event passed the trigger in " << attempts << " attempts" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
    }
  }

  PHHepMCGenEvent *success = PHHepMCGenHelper::insert_event(genevent);
  if (!success)
  {
    std::cout << PHWHERE << " failed to add event to the HepMC record" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  ++m_NPassed;

  if (m_IntegralNode)
  {
    const double sigma_pb = m_Core->xsec_mb() * 1e9;
    m_IntegralNode->set_N_Generator_Accepted_Event(m_NGenerated);
    m_IntegralNode->set_N_Processed_Event(m_NPassed);
    m_IntegralNode->set_Sum_Of_Weight(m_SumWAll);
    // PHPythia8 convention: pb^-1
    m_IntegralNode->set_Integrated_Lumi(sigma_pb > 0 ? double(m_NGenerated) / sigma_pb : 0.0);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHJewel::End(PHCompositeNode * /*topNode*/)
{
  if (!m_Core)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  m_Core->finish();
  const JewelCore::Counts k = m_Core->counts();
  auto frac = [](double over, double tot)
  { return tot > 0 ? 100.0 * over / tot : 0.0; };
  std::cout << Name() << ": JEWEL wrote " << m_NGenerated << " events (" << m_Core->n_discarded_passes()
            << " discarded inside JEWEL), passed trigger " << m_NPassed << " (eff "
            << (m_NGenerated ? double(m_NPassed) / double(m_NGenerated) : 0.0) << ")" << std::endl
            << Name() << ": sumW all=" << m_SumWAll << " passed=" << m_SumWPass << ", sigma=" << m_Core->xsec_mb() * 1e9
            << " pb; normalise a histogram filled with weight w as sigma/sumW_all" << std::endl
            << Name() << ": isospin pp/pn/np/nn = " << m_Core->channel_events(1) << "/" << m_Core->channel_events(2)
            << "/" << m_Core->channel_events(3) << "/" << m_Core->channel_events(4) << ", PYTHIA initialisations "
            << m_Core->channel_switches() << std::endl
            << Name() << ": table extrapolations [%] splitint " << frac(k.noverspliti, k.ntotspliti) << " pdf "
            << frac(k.noverpdf, k.ntotpdf) << " xsec " << frac(k.noverxsec, k.ntotxsec) << " sudakov "
            << frac(k.noversuda, k.ntotsuda) << "  (large values: tables built for other parameters)" << std::endl
            << Name() << ": card hash " << m_Hash << ", log " << m_Stem << ".log" << std::endl;
  if (m_IntegralNode && Verbosity() > 0)
  {
    m_IntegralNode->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHJewel::create_node_tree(PHCompositeNode *topNode)
{
  PHHepMCGenHelper::create_node_tree(topNode);

  PHNodeIterator iter(topNode);
  PHCompositeNode *sumNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!sumNode)
  {
    std::cout << PHWHERE << " RUN node missing, doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (m_SaveIntegratedLuminosityFlag)
  {
    m_IntegralNode = findNode::getClass<PHGenIntegral>(sumNode, "PHGenIntegral");
    if (!m_IntegralNode)
    {
      m_IntegralNode = new PHGenIntegralv1("PHJewel with embedding ID of " + std::to_string(PHHepMCGenHelper::get_embedding_id()));
      PHIODataNode<PHObject> *newnode = new PHIODataNode<PHObject>(m_IntegralNode, "PHGenIntegral", "PHObject");
      sumNode->addNode(newnode);
    }
    else
    {
      std::cout << PHWHERE << " RUN/PHGenIntegral already exists; refusing to overwrite. "
                << "Call PHJewel::save_integrated_luminosity(false) if that is intended." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHJewel::Print(const std::string & /*what*/) const
{
  std::cout << Name() << " configuration:" << std::endl
            << "  card           : " << m_Card << std::endl
            << "  medium         : " << m_Medium << std::endl
            << "  seed           : " << m_Seed << std::endl
            << "  card hash      : " << m_Hash << std::endl
            << "  embedding id   : " << PHHepMCGenHelper::get_embedding_id() << std::endl
            << "  triggers       : " << m_RegisteredTriggers.size()
            << (m_RegisteredTriggers.empty() ? " (every event kept)" : "") << std::endl;
  if (m_RegisteredTriggers.size() > 1)
  {
    std::cout << "  combined with  : " << (m_TriggersOR ? "OR" : "AND") << std::endl;
  }
  for (const PHJewelGenTrigger *trigger : m_RegisteredTriggers)
  {
    std::cout << "    ";
    trigger->PrintConfig();
  }
}

void PHJewel::register_trigger(PHJewelGenTrigger *trigger)
{
  if (Verbosity() > 1)
  {
    std::cout << Name() << ": registering trigger " << trigger->GetName() << std::endl;
  }
  m_RegisteredTriggers.push_back(trigger);
}
