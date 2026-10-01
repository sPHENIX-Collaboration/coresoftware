#include "PHHerwig7.h"
#include "PHHerwig7Core.h"

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

// HepMC2 flavour of ThePEG's converter.  HepMCDefs.h must come first so the
// traits see HEPMC_HAS_CROSS_SECTION and carry sigma into the record.
#include <HepMC/HepMCDefs.h>
#include <HepMC/GenEvent.h>
#include <HepMC/Units.h>
#include <HepMC/WeightContainer.h>
#include <ThePEG/Config/HepMCHelper.h>
#include <ThePEG/Vectors/HepMCConverter.h>

#include <cassert>
#include <cstdlib>
#include <iostream>

PHHerwig7::PHHerwig7(const std::string &name)
  : SubsysReco(name)
{
  PHHepMCGenHelper::set_embedding_id(1);  // same default as PHPythia8
}

PHHerwig7::~PHHerwig7() = default;

int PHHerwig7::Init(PHCompositeNode *topNode)
{
  if (m_RunFile.empty())
  {
    std::cout << PHWHERE << " no Herwig .run file set (PHHerwig7::set_run_file)" << std::endl;
    exit(EXIT_FAILURE);
  }

  create_node_tree(topNode);

  int seed = m_Seed;
  if (seed <= 0)
  {
    unsigned int s = PHRandomSeed();
    // ThePEG takes an int seed through the UI; keep it positive and in range
    seed = static_cast<int>(s % 2000000000u);
    if (seed <= 0)
    {
      seed = 1;
    }
  }
  std::cout << Name() << " Herwig random seed: " << seed << std::endl;

  m_Core = std::make_unique<PHHerwig7Core>();
  m_Core->ui().m_RunFile = m_RunFile;
  m_Core->ui().m_SetupFile = m_SetupFile;
  m_Core->ui().m_Prepend = m_ReadDirs;
  m_Core->ui().m_Seed = seed;

  // Herwig's cout formatting is not restored by ThePEG
  std::ios old_state(nullptr);
  old_state.copyfmt(std::cout);
  try
  {
    m_Core->init();
  }
  catch (const std::exception &e)
  {
    std::cout.copyfmt(old_state);
    std::cout << PHWHERE << " Herwig initialisation failed: " << e.what() << std::endl;
    // Fun4All drops a module whose Init fails and carries on writing empty
    // events; a generator that cannot start must kill the job (as PHPythia8
    // does on a bad seed).
    std::cout << PHWHERE << " Check the .run file and LHAPDF_DATA_PATH (hwrun.sh sets it to Source/lhapdf)" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout.copyfmt(old_state);

  if (Verbosity() > 0)
  {
    Print("CONFIG");
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHHerwig7::process_event(PHCompositeNode * /*topNode*/)
{
  std::ios old_state(nullptr);
  old_state.copyfmt(std::cout);

  HepMC::GenEvent *genevent = nullptr;
  long attempts = 0;
  while (!genevent)
  {
    ThePEG::EventPtr ev;
    try
    {
      ev = m_Core->next();
    }
    catch (const std::exception &e)
    {
      std::cout.copyfmt(old_state);
      std::cout << PHWHERE << " Herwig threw during event generation: " << e.what() << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    if (!ev)
    {
      // shoot() returns null only when the generator has hit its own event
      // limit or given up; retrying would loop forever.
      std::cout.copyfmt(old_state);
      std::cout << PHWHERE << " Herwig returned no event after " << m_NGenerated
                << " generated (generator event limit reached?), aborting run" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    ++m_NGenerated;
    const double w = ev->weight();
    m_SumWAll += w;

    HepMC::GenEvent *cand = ThePEG::HepMCConverter<HepMC::GenEvent>::convert(
        *ev, false, ThePEG::Units::GeV, ThePEG::Units::millimeter);
    if (!cand)
    {
      std::cout << PHWHERE << " HepMC conversion failed" << std::endl;
      continue;
    }
    if (m_SaveEventWeightFlag && cand->weights().empty())
    {
      cand->weights().push_back(w);
    }

    const HepMCTruthTriggerCore::Result r = HepMCTruthTriggerCore::evaluate(cand, m_Trig);
    if (Verbosity() > 1)
    {
      std::cout << Name() << ": generated " << m_NGenerated
                << " leadJetPt=" << r.leadJetPt << " leadPhotonPt=" << r.leadPhotonPt
                << " pass=" << r.pass << std::endl;
    }
    if (r.pass)
    {
      genevent = cand;
      m_SumWPass += w;
    }
    else
    {
      delete cand;
      ++attempts;
      if (m_MaxAttempts > 0 && attempts >= m_MaxAttempts)
      {
        std::cout.copyfmt(old_state);
        std::cout << PHWHERE << " no event passed the trigger in " << attempts << " attempts" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
    }
  }

  PHHepMCGenEvent *success = PHHepMCGenHelper::insert_event(genevent);
  if (!success)
  {
    std::cout.copyfmt(old_state);
    std::cout << PHWHERE << " failed to add event to the HepMC record" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  ++m_NPassed;

  if (m_IntegralNode)
  {
    const double sigma_pb = m_Core->xsec_pb();
    m_IntegralNode->set_N_Generator_Accepted_Event(m_NGenerated);
    m_IntegralNode->set_N_Processed_Event(m_NPassed);
    m_IntegralNode->set_Sum_Of_Weight(m_SumWAll);
    // PHPythia8 convention: pb^-1
    m_IntegralNode->set_Integrated_Lumi(sigma_pb > 0 ? double(m_NGenerated) / sigma_pb : 0.0);
  }

  std::cout.copyfmt(old_state);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHHerwig7::End(PHCompositeNode * /*topNode*/)
{
  std::ios old_state(nullptr);
  old_state.copyfmt(std::cout);
  if (m_Core)
  {
    m_Core->finish();
  }
  std::cout.copyfmt(old_state);

  std::cout << Name() << ": generated " << m_NGenerated << ", passed trigger " << m_NPassed
            << " (eff " << (m_NGenerated ? double(m_NPassed) / double(m_NGenerated) : 0.0) << ")"
            << ", sumW all=" << m_SumWAll << " passed=" << m_SumWPass
            << ", sigma=" << (m_Core ? m_Core->xsec_pb() : 0.0) << " +- "
            << (m_Core ? m_Core->xsec_err_pb() : 0.0) << " pb" << std::endl;
  if (m_IntegralNode && Verbosity() > 0)
  {
    m_IntegralNode->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHHerwig7::create_node_tree(PHCompositeNode *topNode)
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
      m_IntegralNode = new PHGenIntegralv1("PHHerwig7 with embedding ID of " + std::to_string(PHHepMCGenHelper::get_embedding_id()));
      PHIODataNode<PHObject> *newnode = new PHIODataNode<PHObject>(m_IntegralNode, "PHGenIntegral", "PHObject");
      sumNode->addNode(newnode);
    }
    else
    {
      std::cout << PHWHERE << " RUN/PHGenIntegral already exists; refusing to overwrite. "
                << "Call PHHerwig7::save_integrated_luminosity(false) if that is intended." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHHerwig7::Print(const std::string & /*what*/) const
{
  std::cout << Name() << " configuration:" << std::endl
            << "  run file       : " << m_RunFile << std::endl
            << "  setup file     : " << (m_SetupFile.empty() ? "(none)" : m_SetupFile) << std::endl
            << "  embedding id   : " << PHHepMCGenHelper::get_embedding_id() << std::endl
            << "  jet trigger    : " << (m_Trig.requireJet ? "on" : "off")
            << "  R=" << m_Trig.jetR << " pT>=" << m_Trig.jetPtMin << " |eta|<" << m_Trig.jetEtaMax << std::endl
            << "  leadjet window : [" << m_Trig.windowLo << ", " << m_Trig.windowHi << ")" << std::endl
            << "  photon trigger : " << (m_Trig.requirePhoton ? "on" : "off")
            << "  pT>=" << m_Trig.photonPtMin << " |eta|<=" << m_Trig.photonEtaMax
            << " iso=" << (m_Trig.photonIso ? "on" : "off") << " prompt=" << (m_Trig.photonPrompt ? "on" : "off") << std::endl
            << "  combine        : " << (m_Trig.andMode ? "AND" : "OR") << std::endl;
}
