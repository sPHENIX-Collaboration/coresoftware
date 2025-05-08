// ----------------------------------------------------------------------------
// 'TrksInJetQA.cc'
// Derek Anderson
// 03.25.2024
//
// A "small" Fun4All module to produce QA plots for tracks,
// hits, and more.
// ----------------------------------------------------------------------------

#define TRKSINJETQA_CC

// module defintion
#include "TrksInJetQA.h"
#include <TStyle.h>
// ctor/dtor ------------------------------------------------------------------

TrksInJetQA::TrksInJetQA(const std::string& name)
  : SubsysReco(name)
  , m_moduleName(name)
{
}

TrksInJetQA::~TrksInJetQA()
{
  // print debug messages
  if (m_config.doDebug && (m_config.verbose > 4))
  {
    std::cout << "TrksInJetQA::~TrksInJetQA() Calling dtor" << std::endl;
  }

  // clean up any dangling pointers
  // deleting null ptrs is legal, setting it to null is not needed in the dtor
  delete m_outFile;
}  // end dtor

// public methods -------------------------------------------------------------

void TrksInJetQA::Configure(
    const TrksInJetQAConfig& config,
    std::optional<TrksInJetQAHist> hist)
{
  m_config = config;
  // print debug messages
  if (m_config.doDebug && (m_config.verbose > 3))
  {
    std::cout << "TrksInJetQA::~TrksInJetQA() Calling dtor" << std::endl;
  }

  if (hist.has_value())
  {
    m_hist = hist.value();
  }
  return;

}  // end 'Configure(TrksInJetQAConfig, std::optional<TrksInJetQAHist>)'

// fun4all methods ------------------------------------------------------------

int TrksInJetQA::Init(PHCompositeNode* /*topNode*/)
{
  // print debug message
  if (m_config.doDebug && (m_config.verbose > 0))
  {
    std::cout << "TrksInJetQA::Init(PHCompositeNode* /*topNode*/) Initializing" << std::endl;
  }

  // initialize output & histograms
  InitOutput();
  InitHistograms();

  // register histograms with manager if needed
  if (m_config.outMode == OutMode::QA)
  {
    RegisterHistograms();
  }

  // initialize trigger analyzer and exit
  delete m_analyzer; // make cppcheck happy
  m_analyzer = new TriggerAnalyzer();
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'Init(PHCompositeNode*)'

int TrksInJetQA::process_event(PHCompositeNode* topNode)
{
  // print debug message
  if (m_config.doDebug && (m_config.verbose > 2))
  {
    std::cout << "TrksInJetQA::process_event(PHCompositeNode* topNode) Processing Event" << std::endl;
  }

  // if needed, check if selected trigger fired
  if (m_doTrgSelect)
  {
    m_analyzer->decodeTriggers(topNode);
    bool hasTrigger = JetQADefs::DidTriggerFire(m_trgToSelect, m_analyzer);
    if (!hasTrigger)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  // run submodules
  if (m_config.doInJet)
  {
    m_inJet->Fill(topNode);
  }
  if (m_config.doInclusive)
  {
    m_inclusive->Fill(topNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'process_event(PHCompositeNode*)'

int TrksInJetQA::End(PHCompositeNode* /*topNode*/)
{
  // print debug message
  if (m_config.doDebug && (m_config.verbose > 0))
  {
    std::cout << "TrksInJetQA::End(PHCompositeNode* /*topNode*/) This is the End..." << std::endl;
  }

  // save hists to file if needed
  if (m_config.outMode == OutMode::File)
  {
    // save histograms
    if (m_config.doInJet)
    {
      m_inJet->SaveHistograms(m_outFile, "InJet");
    }
    if (m_config.doInclusive)
    {
      m_inclusive->SaveHistograms(m_outFile, "Inclusive");
    }

    // close file
    m_outFile->cd();
    m_outFile->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'End(PHCompositeNode*)'

// private methods ------------------------------------------------------------

void TrksInJetQA::InitOutput()
{
  // print debug message
  if (m_config.doDebug && (m_config.verbose > 1))
  {
    std::cout << "TrksInJetQA::InitOutput() Initializing outputs..." << std::endl;
  }

  // initialize relevent outputs
  switch (m_config.outMode)
  {
  case OutMode::File:
    m_outFile = new TFile(m_outFileName.data(), "RECREATE");
    if (!m_outFile)
    {
      std::cerr << PHWHERE << ": PANIC: couldn't create output file!" << std::endl;
      assert(m_outFile);
    }
    break;

  case OutMode::QA:
    delete m_manager;

    gStyle->SetOptTitle(0);
    m_manager = QAHistManagerDef::getHistoManager();
    if (!m_manager)
    {
      std::cerr << PHWHERE << ": PANIC: couldn't grab histogram manager!" << std::endl;
      assert(m_manager);
    }
    break;

  default:
    std::cerr << PHWHERE << ": PANIC: unknown output mode specified!\n"
              << "  Please set .outMode = OutMode::File OR OutMode::QA!"
              << std::endl;
    assert((m_config.outMode == OutMode::File) || (m_config.outMode == OutMode::QA));
    break;
  }
  return;

}  // end 'InitOutput()'

void TrksInJetQA::InitHistograms()
{
  // print debug message
  if (m_config.doDebug && (m_config.verbose > 1))
  {
    std::cout << "TrksInJetQA::InitHistograms() Initializing histograms..." << std::endl;
  }

  // histograms are always prefixed by the module name
  std::string prefix = "h_";
  prefix += m_moduleName;

  // if additional prefix provided, add it
  if (m_histPrefix.has_value())
  {
    prefix += m_histPrefix.value();
    prefix += "_";
  }

  // make suffixes
  std::string inJetSuffix = "InJet";
  std::string inclusiveSuffix = "Inclusive";
  if (m_histSuffix.has_value() && !m_histSuffix.value().empty())
  {
    inJetSuffix += "_";
    inJetSuffix += m_histSuffix.value();
    inclusiveSuffix += "_";
    inclusiveSuffix += m_histSuffix.value();
  }

  // initialize submodules, as needed
  if (m_config.doInJet)
  {
    m_inJet = std::make_unique<TrksInJetQAInJetFiller>(m_config, m_hist);
    m_inJet->MakeHistograms(prefix, inJetSuffix);
  }
  if (m_config.doInclusive)
  {
    m_inclusive = std::make_unique<TrksInJetQAInclusiveFiller>(m_config, m_hist);
    m_inclusive->MakeHistograms(prefix, inclusiveSuffix);
  }
  return;

}  // end 'InitHistograms()'

void TrksInJetQA::RegisterHistograms()
{
  // print debug message
  if (m_config.doDebug && (m_config.verbose > 1))
  {
    std::cout << "TrksInJetQA::RegisterHistograms() Registering histograms..." << std::endl;
  }

  std::vector<TH1D*> vecHist1D;
  std::vector<TH2D*> vecHist2D;
  if (m_config.doInJet)
  {
    m_inJet->GrabHistograms(vecHist1D, vecHist2D);
  }
  if (m_config.doInclusive)
  {
    m_inclusive->GrabHistograms(vecHist1D, vecHist2D);
  }

  // register w/ manager
  for (TH1D* hist1D : vecHist1D)
  {
    m_manager->registerHisto(hist1D);
  }
  for (TH2D* hist2D : vecHist2D)
  {
    m_manager->registerHisto(hist2D);
  }
  return;

}  // end 'RegisterHistograms()'

// end ------------------------------------------------------------------------
