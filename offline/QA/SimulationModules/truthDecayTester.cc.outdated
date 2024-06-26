#include "truthDecayTester.h"

#include "QAHistManagerDef.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <decayfinder/DecayFinder.h>
#include <decayfinder/DecayFinderContainerBase.h>  // for DecayFinderContainerBase::Iter
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <CLHEP/Vector/LorentzVector.h>

#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TString.h>
#include <TVector3.h>

#include <TDatabasePDG.h>

int candidateCounter = 0;

//____________________________________________________________________________..
truthDecayTester::truthDecayTester(const std::string &name)
  : SubsysReco(name)
  , m_nTracks(0)
  , m_min_pt(0.2)
  , m_min_eta(-1.1)
  , m_max_eta(1.1)
  , m_write_nTuple(true)
  , m_decay_pdg_id(0)
  , m_truth_info(nullptr)
  , m_g4particle(nullptr)
  , m_df_module_name("myFinder")
  , m_outfile_name("outputData.root")
  , m_outfile(nullptr)
  , m_tree(nullptr)
  , m_write_QAHists(true)
{
}

//____________________________________________________________________________..
truthDecayTester::~truthDecayTester() {}

//____________________________________________________________________________..
int truthDecayTester::Init(PHCompositeNode *topNode)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1 *h(nullptr);

  h = new TH1F(TString(get_histo_prefix()) + "mother_PDG_ID",  //
               ";Mother PDG ID;Entries", 1000, -500, 500);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "mother_mass",  //
               ";Mother Mass [GeV/c^{2}];Entries", 100, 0, 3);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "daughter_sum_mass",  //
               ";Daughter Sum Mass [GeV/c^{2}];Entries", 100, 0, 3);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "mother_decayLength",  //
               ";Mother Decay Length [cm];Entries", 50, 0, 1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "mother_decayTime",  //
               ";Mother Decay Time [s];Entries", 50, 0, 0.1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "mother_px",  //
               ";Mother p_{x} [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "mother_py",  //
               ";Mother p_{y} [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "mother_pz",  //
               ";Mother p_{z} [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "mother_pE",  //
               ";Mother p_{E} [GeV];Entries", 100, 0, 10);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "mother_pT",  //
               ";Mother p_{T} [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "mother_eta",  //
               ";Mother #eta;Entries", 100, -3, 3);
  hm->registerHisto(h);

  for (unsigned int i = 0; i < 4; ++i)
  {
    std::string track_number = "track_" + std::to_string(i + 1);

    h = new TH1F(TString(get_histo_prefix()) + TString(track_number) + "_PDG_ID",  //
                 ";Track PDG ID;Entries", 1000, -500, 500);
    hm->registerHisto(h);

    /*

    h = new TH1F(TString(get_histo_prefix()) + TString(track_number) + "_px",  //
                 ";Track p_{x} [GeV/c];Entries", 50, 0, 5);
    hm->registerHisto(h);
    h = new TH1F(TString(get_histo_prefix()) + TString(track_number) + "_py",  //
                 ";Track p_{y} [GeV/c];Entries", 50, 0, 5);
    hm->registerHisto(h);
    h = new TH1F(TString(get_histo_prefix()) + TString(track_number) + "_pz",  //
                 ";Track p_{z} [GeV/c];Entries", 50, 0, 5);
    hm->registerHisto(h);
    h = new TH1F(TString(get_histo_prefix()) + TString(track_number) + "_pE",  //
                 ";Track p_{E} [GeV];Entries", 50, 0, 5);
    hm->registerHisto(h);

    */

    h = new TH1F(TString(get_histo_prefix()) + TString(track_number) + "_pT",  //
                 ";Track pT [GeV/c];Entries", 50, 0, 5);
    hm->registerHisto(h);
    h = new TH1F(TString(get_histo_prefix()) + TString(track_number) + "_eta",  //
                 ";Track Eta;Entries", 100, -3, 3);
    hm->registerHisto(h);
    h = new TH1F(TString(get_histo_prefix()) + TString(track_number) + "_mass",  //
                 ";Track Mass [GeV/c^{2}];Entries", 100, 0, 3);
    hm->registerHisto(h);
  }

  h = new TH1F(TString(get_histo_prefix()) + "delta_px",  //
               ";#delta p_{x} [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "delta_py",  //
               ";#delta p_{y} [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "delta_pz",  //
               ";#delta p_{z} [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "delta_pE",  //
               ";#delta p_{E} [GeV];Entries", 100, 0, 10);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "accept_px_1percent",  //
               ";Accept p_{x} 1pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "accept_py_1percent",  //
               ";Accept p_{y} 1pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "accept_pz_1percent",  //
               ";Accept p_{z} 1pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "accept_pE_1percent",  //
               ";Accept p_{E} 1pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "accept_px_5percent",  //
               ";Accept p_{x} 5pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "accept_py_5percent",  //
               ";Accept p_{y} 5pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "accept_pz_5percent",  //
               ";Accept p_{z} 5pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "accept_pE_5percent",  //
               ";Accept p_{E} 5pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "accept_px_15percent",  //
               ";Accept p_{x} 15pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "accept_py_15percent",  //
               ";Accept p_{y} 15pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "accept_pz_15percent",  //
               ";Accept p_{z} 15pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "accept_pE_15percent",  //
               ";Accept p_{E} 15pcnt;Entries", 2, 0, 1);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "accept_pT",  //
               ";Accept p_{T};Entries", 2, 0, 1);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "accept_eta",  //
               ";Accept #eta;Entries", 2, 0, 1);
  hm->registerHisto(h);

  assert(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int truthDecayTester::process_event(PHCompositeNode *topNode)
{
  resetValues();
  CLHEP::HepLorentzVector daughterSumLV;

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1F *h_mother_PDG_ID = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "mother_PDG_ID"));
  assert(h_mother_PDG_ID);
  TH1F *h_mother_mass = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "mother_mass"));
  assert(h_mother_mass);
  TH1F *h_daughter_sum_mass = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "daughter_sum_mass"));
  assert(h_daughter_sum_mass);
  TH1F *h_mother_decayLength = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "mother_decayLength"));
  assert(h_mother_decayLength);
  TH1F *h_mother_decayTime = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "mother_decayTime"));
  assert(h_mother_decayTime);
  TH1F *h_mother_px = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "mother_px"));
  assert(h_mother_px);
  TH1F *h_mother_py = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "mother_py"));
  assert(h_mother_py);
  TH1F *h_mother_pz = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "mother_pz"));
  assert(h_mother_pz);
  TH1F *h_mother_pE = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "mother_pE"));
  assert(h_mother_pE);
  TH1F *h_mother_pT = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "mother_pT"));
  assert(h_mother_pT);
  TH1F *h_mother_eta = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "mother_eta"));
  assert(h_mother_eta);

  TH1F *h_delta_px = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "delta_px"));
  assert(h_delta_px);
  TH1F *h_delta_py = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "delta_py"));
  assert(h_delta_py);
  TH1F *h_delta_pz = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "delta_pz"));
  assert(h_delta_pz);
  TH1F *h_delta_pE = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "delta_pE"));
  assert(h_delta_pE);
  TH1F *h_accept_px_1percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_px_1percent"));
  assert(h_accept_px_1percent);
  TH1F *h_accept_py_1percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_py_1percent"));
  assert(h_accept_py_1percent);
  TH1F *h_accept_pz_1percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_pz_1percent"));
  assert(h_accept_pz_1percent);
  TH1F *h_accept_pE_1percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_pE_1percent"));
  assert(h_accept_pE_1percent);
  TH1F *h_accept_px_5percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_px_5percent"));
  assert(h_accept_px_5percent);
  TH1F *h_accept_py_5percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_py_5percent"));
  assert(h_accept_py_5percent);
  TH1F *h_accept_pz_5percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_pz_5percent"));
  assert(h_accept_pz_5percent);
  TH1F *h_accept_pE_5percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_pE_5percent"));
  assert(h_accept_pE_5percent);
  TH1F *h_accept_px_15percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_px_15percent"));
  assert(h_accept_px_15percent);
  TH1F *h_accept_py_15percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_py_15percent"));
  assert(h_accept_py_15percent);
  TH1F *h_accept_pz_15percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_pz_15percent"));
  assert(h_accept_pz_15percent);
  TH1F *h_accept_pE_15percent = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_pE_15percent"));
  assert(h_accept_pE_15percent);
  TH1F *h_accept_pT = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_pT"));
  assert(h_accept_pT);
  TH1F *h_accept_eta = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "accept_eta"));
  assert(h_accept_eta);

  TH1F *h_track_1_PDG_ID = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_1_PDG_ID"));
  assert(h_track_1_PDG_ID);
  /*
  TH1F *h_track_1_px = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_1_px"));
  assert(h_track_1_px);
  TH1F *h_track_1_py = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_1_py"));
  assert(h_track_1_py);
  TH1F *h_track_1_pz = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_1_pz"));
  assert(h_track_1_pz);
  TH1F *h_track_1_pE = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_1_pE"));
  assert(h_track_1_pE);
  */
  TH1F *h_track_1_pT = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_1_pT"));
  assert(h_track_1_pT);
  TH1F *h_track_1_eta = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_1_eta"));
  assert(h_track_1_eta);
  TH1F *h_track_1_mass = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_1_mass"));
  assert(h_track_1_mass);

  TH1F *h_track_2_PDG_ID = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_2_PDG_ID"));
  assert(h_track_2_PDG_ID);
  /*
  TH1F *h_track_2_px = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_2_px"));
  assert(h_track_2_px);
  TH1F *h_track_2_py = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_2_py"));
  assert(h_track_2_py);
  TH1F *h_track_2_pz = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_2_pz"));
  assert(h_track_2_pz);
  TH1F *h_track_2_pE = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_2_pE"));
  assert(h_track_2_pE);
  */
  TH1F *h_track_2_pT = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_2_pT"));
  assert(h_track_2_pT);
  TH1F *h_track_2_eta = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_2_eta"));
  assert(h_track_2_eta);
  TH1F *h_track_2_mass = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_2_mass"));
  assert(h_track_2_mass);

  TH1F *h_track_3_PDG_ID = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_3_PDG_ID"));
  assert(h_track_3_PDG_ID);
  /*
  TH1F *h_track_3_px = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_3_px"));
  assert(h_track_3_px);
  TH1F *h_track_3_py = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_3_py"));
  assert(h_track_3_py);
  TH1F *h_track_3_pz = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_3_pz"));
  assert(h_track_3_pz);
  TH1F *h_track_3_pE = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_3_pE"));
  assert(h_track_3_pE);
  */
  TH1F *h_track_3_pT = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_3_pT"));
  assert(h_track_3_pT);
  TH1F *h_track_3_eta = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_3_eta"));
  assert(h_track_3_eta);
  TH1F *h_track_3_mass = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_3_mass"));
  assert(h_track_3_mass);

  TH1F *h_track_4_PDG_ID = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_4_PDG_ID"));
  assert(h_track_4_PDG_ID);
  /*
  TH1F *h_track_4_px = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_4_px"));
  assert(h_track_4_px);
  TH1F *h_track_4_py = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_4_py"));
  assert(h_track_4_py);
  TH1F *h_track_4_pz = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_4_pz"));
  assert(h_track_4_pz);
  TH1F *h_track_4_pE = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_4_pE"));
  assert(h_track_4_pE);
  */
  TH1F *h_track_4_pT = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_4_pT"));
  assert(h_track_4_pT);
  TH1F *h_track_4_eta = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_4_eta"));
  assert(h_track_4_eta);
  TH1F *h_track_4_mass = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "track_4_mass"));
  assert(h_track_4_mass);

  ++m_event_number;
  if (m_decay_pdg_id == 0) getMotherPDG(topNode);
  std::vector<int> motherBarcodes = getDecayFinderMothers(topNode);

  if (motherBarcodes.size() == 1)
  {
    m_truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (!m_truth_info)
    {
      std::cout << "truthDecayTester: Missing node G4TruthInfo" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    unsigned int trackCounter = 0;
    float mother_x = 0;
    float mother_y = 0;
    float mother_z = 0;
    float daughter_x = 0;
    float daughter_y = 0;
    float daughter_z = 0;

    PHG4TruthInfoContainer::ConstRange range = m_truth_info->GetParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
         iter != range.second; ++iter)
    {
      m_g4particle = iter->second;
      if (std::find(std::begin(motherBarcodes), std::end(motherBarcodes),
                    m_g4particle->get_barcode()) != std::end(motherBarcodes) &&
          abs(m_g4particle->get_pid()) == abs(m_decay_pdg_id))
      {
        m_mother_pdg_id = m_g4particle->get_pid();
        m_mother_barcode = m_g4particle->get_barcode();
        m_mother_px = m_g4particle->get_px();
        m_mother_py = m_g4particle->get_py();
        m_mother_pz = m_g4particle->get_pz();
        m_mother_pE = m_g4particle->get_e();
        CLHEP::HepLorentzVector motherLV(m_mother_px, m_mother_py, m_mother_pz, m_mother_pE);
        m_mother_pT = motherLV.perp();
        m_mother_eta = motherLV.pseudoRapidity();
        m_mother_mass = motherLV.m();

        m_delta_px += m_mother_px;
        m_delta_py += m_mother_py;
        m_delta_pz += m_mother_pz;
        m_delta_pE += m_mother_pE;

        PHG4VtxPoint *mother_vtx = m_truth_info->GetVtx(m_g4particle->get_vtx_id());
        mother_x = mother_vtx->get_x();
        mother_y = mother_vtx->get_y();
        mother_z = mother_vtx->get_z();

        ++candidateCounter;
      }

      if (m_g4particle->get_parent_id() != 0)
      {
        PHG4Particle *mother = m_truth_info->GetParticle(m_g4particle->get_parent_id());
        if (std::find(std::begin(motherBarcodes), std::end(motherBarcodes),
                      mother->get_barcode()) != std::end(motherBarcodes) &&
            abs(mother->get_pid()) == abs(m_decay_pdg_id))
        {
          m_track_pdg_id[trackCounter] = m_g4particle->get_pid();
          m_track_mother_barcode[trackCounter] = mother->get_barcode();
          m_track_px[trackCounter] = m_g4particle->get_px();
          m_track_py[trackCounter] = m_g4particle->get_py();
          m_track_pz[trackCounter] = m_g4particle->get_pz();
          m_track_pE[trackCounter] = m_g4particle->get_e();
          CLHEP::HepLorentzVector daughterLV(m_track_px[trackCounter], m_track_py[trackCounter], m_track_pz[trackCounter], m_track_pE[trackCounter]);
          daughterSumLV += daughterLV;
          m_track_pT[trackCounter] = daughterLV.perp();
          m_track_eta[trackCounter] = daughterLV.pseudoRapidity();
          m_track_mass[trackCounter] = daughterLV.m();

          m_delta_px -= m_track_px[trackCounter];
          m_delta_py -= m_track_py[trackCounter];
          m_delta_pz -= m_track_pz[trackCounter];
          m_delta_pE -= m_track_pE[trackCounter];

          PHG4VtxPoint *daughter_vtx = m_truth_info->GetVtx(m_g4particle->get_vtx_id());
          daughter_x = daughter_vtx->get_x();
          daughter_y = daughter_vtx->get_y();
          daughter_z = daughter_vtx->get_z();

          if (m_track_pT[trackCounter] < m_min_pt) m_accept_pT = false;
          bool in_eta_range = isInRange(m_track_eta[trackCounter], m_min_eta, m_max_eta);
          if (!in_eta_range) m_accept_eta = false;

          ++trackCounter;
        }
      }
    }

    m_daughter_sum_mass = daughterSumLV.m();

    float diff_percent_px = fabs(m_delta_px / m_mother_px) * 100.;
    float diff_percent_py = fabs(m_delta_py / m_mother_py) * 100.;
    float diff_percent_pz = fabs(m_delta_pz / m_mother_pz) * 100.;
    float diff_percent_pE = fabs(m_delta_pE / m_mother_pE) * 100.;

    m_accept_px_1percent = diff_percent_px <= 1. ? 1 : 0;
    m_accept_py_1percent = diff_percent_py <= 1. ? 1 : 0;
    m_accept_pz_1percent = diff_percent_pz <= 1. ? 1 : 0;
    m_accept_pE_1percent = diff_percent_pE <= 1. ? 1 : 0;

    m_accept_px_5percent = diff_percent_px <= 5. ? 1 : 0;
    m_accept_py_5percent = diff_percent_py <= 5. ? 1 : 0;
    m_accept_pz_5percent = diff_percent_pz <= 5. ? 1 : 0;
    m_accept_pE_5percent = diff_percent_pE <= 5. ? 1 : 0;

    m_accept_px_15percent = diff_percent_px <= 15. ? 1 : 0;
    m_accept_py_15percent = diff_percent_py <= 15. ? 1 : 0;
    m_accept_pz_15percent = diff_percent_pz <= 15. ? 1 : 0;
    m_accept_pE_15percent = diff_percent_pE <= 15. ? 1 : 0;

    m_mother_decayLength = sqrt(pow(daughter_x - mother_x, 2) + pow(daughter_y - mother_y, 2) + pow(daughter_z - mother_z, 2));
    float mother_p = sqrt(pow(m_mother_px, 2) + pow(m_mother_py, 2) + pow(m_mother_pz, 2));
    m_mother_decayTime = m_mother_mass * m_mother_decayLength / mother_p;

    if (m_write_nTuple)
    {
      if (candidateCounter == 1) initializeBranches();
      m_tree->Fill();
    }
    if (m_write_QAHists)
    {
      h_mother_PDG_ID->Fill(m_mother_pdg_id);
      std::cout << m_mother_pdg_id << " is the mother PDG ID" << std::endl;
      std::cout << h_mother_PDG_ID->Integral(h_mother_PDG_ID->FindFixBin(-500.), h_mother_PDG_ID->FindFixBin(500.) - 1) << std::endl;
      h_mother_mass->Fill(m_mother_mass);
      h_daughter_sum_mass->Fill(m_daughter_sum_mass);
      h_mother_decayLength->Fill(m_mother_decayLength);
      h_mother_decayTime->Fill(m_mother_decayTime);
      h_mother_px->Fill(m_mother_px);
      h_mother_py->Fill(m_mother_py);
      h_mother_pz->Fill(m_mother_pz);
      h_mother_pE->Fill(m_mother_pE);
      h_mother_pT->Fill(m_mother_pT);
      h_mother_eta->Fill(m_mother_eta);
      h_delta_px->Fill(m_delta_px);
      h_delta_py->Fill(m_delta_py);
      h_delta_pz->Fill(m_delta_pz);
      h_delta_pE->Fill(m_delta_pE);
      h_accept_px_1percent->Fill(m_accept_px_1percent);
      h_accept_py_1percent->Fill(m_accept_py_1percent);
      h_accept_pz_1percent->Fill(m_accept_pz_1percent);
      h_accept_pE_1percent->Fill(m_accept_pE_1percent);
      h_accept_px_5percent->Fill(m_accept_px_5percent);
      h_accept_py_5percent->Fill(m_accept_py_5percent);
      h_accept_pz_5percent->Fill(m_accept_pz_5percent);
      h_accept_pE_5percent->Fill(m_accept_pE_5percent);
      h_accept_px_15percent->Fill(m_accept_px_15percent);
      h_accept_py_15percent->Fill(m_accept_py_15percent);
      h_accept_pz_15percent->Fill(m_accept_pz_15percent);
      h_accept_pE_15percent->Fill(m_accept_pE_15percent);
      h_accept_pT->Fill(m_accept_pT);
      h_accept_eta->Fill(m_accept_eta);
      if (m_nTracks >= 2)
      {
        h_track_1_PDG_ID->Fill(m_track_pdg_id[0]);
        /*
        h_track_1_px->Fill(m_track_px[0]);
        h_track_1_py->Fill(m_track_py[0]);
        h_track_1_pz->Fill(m_track_pz[0]);
        h_track_1_pE->Fill(m_track_pE[0]);
        */
        h_track_1_pT->Fill(m_track_pT[0]);
        h_track_1_eta->Fill(m_track_eta[0]);
        h_track_1_mass->Fill(m_track_mass[0]);
        h_track_2_PDG_ID->Fill(m_track_pdg_id[1]);
        /*
        h_track_2_px->Fill(m_track_px[1]);
        h_track_2_py->Fill(m_track_py[1]);
        h_track_2_pz->Fill(m_track_pz[1]);
        h_track_2_pE->Fill(m_track_pE[1]);
        */
        h_track_2_pT->Fill(m_track_pT[1]);
        h_track_2_eta->Fill(m_track_eta[1]);
        h_track_2_mass->Fill(m_track_mass[1]);
      }
      if (m_nTracks >= 3)
      {
        h_track_3_PDG_ID->Fill(m_track_pdg_id[2]);
        /*
        h_track_3_px->Fill(m_track_px[2]);
        h_track_3_py->Fill(m_track_py[2]);
        h_track_3_pz->Fill(m_track_pz[2]);
        h_track_3_pE->Fill(m_track_pE[2]);
        */
        h_track_3_pT->Fill(m_track_pT[2]);
        h_track_3_eta->Fill(m_track_eta[2]);
        h_track_3_mass->Fill(m_track_mass[2]);
      }
      if (m_nTracks == 4)
      {
        h_track_4_PDG_ID->Fill(m_track_pdg_id[3]);
        /*
        h_track_4_px->Fill(m_track_px[3]);
        h_track_4_py->Fill(m_track_py[3]);
        h_track_4_pz->Fill(m_track_pz[3]);
        h_track_4_pE->Fill(m_track_pE[3]);
        */
        h_track_4_pT->Fill(m_track_pT[3]);
        h_track_4_eta->Fill(m_track_eta[3]);
        h_track_4_mass->Fill(m_track_mass[3]);
      }
    }
  }
  else
  {
    std::cout << "You have more than one good decay in this event, this processing is not yet supported" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int truthDecayTester::End(PHCompositeNode *topNode)
{
  if (m_write_nTuple)
  {
    m_outfile->Write();
    m_outfile->Close();
    delete m_outfile;
  }

  assert(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void truthDecayTester::initializeBranches()
{
  m_outfile = new TFile(m_outfile_name.c_str(), "RECREATE");
  delete m_tree;
  m_tree = new TTree("truthDecayTester", "truthDecayTester");
  m_tree->OptimizeBaskets();
  m_tree->SetAutoSave(-5e6);  // Save the output file every 5MB

  m_tree->Branch("EventNumber", &m_event_number, "EventNumber/i");
  m_tree->Branch("mother_PDG_ID", &m_mother_pdg_id, "mother_PDG_ID/I");
  m_tree->Branch("mother_mass", &m_mother_mass, "mother_mass/F");
  m_tree->Branch("daughter_sum_mass", &m_daughter_sum_mass, "daughter_sum_mass/F");
  m_tree->Branch("mother_decayLength", &m_mother_decayLength, "mother_decayLength/F");
  m_tree->Branch("mother_decayTime", &m_mother_decayTime, "mother_decayTime/F");
  m_tree->Branch("mother_px", &m_mother_px, "mother_px/F");
  m_tree->Branch("mother_py", &m_mother_py, "mother_py/F");
  m_tree->Branch("mother_pz", &m_mother_pz, "mother_pz/F");
  m_tree->Branch("mother_pE", &m_mother_pE, "mother_pE/F");
  m_tree->Branch("mother_pT", &m_mother_pT, "mother_pT/F");
  m_tree->Branch("mother_eta", &m_mother_eta, "mother_eta/F");
  m_tree->Branch("mother_barcode", &m_mother_barcode, "mother_barcode/I");

  for (unsigned int i = 0; i < m_nTracks; ++i)
  {
    std::string track_number = "track_" + std::to_string(i + 1);

    m_tree->Branch(TString(track_number) + "_PDG_ID", &m_track_pdg_id[i], TString(track_number) + "_PDG_ID/I");
    m_tree->Branch(TString(track_number) + "_px", &m_track_px[i], TString(track_number) + "_px/F");
    m_tree->Branch(TString(track_number) + "_py", &m_track_py[i], TString(track_number) + "_py/F");
    m_tree->Branch(TString(track_number) + "_pz", &m_track_pz[i], TString(track_number) + "_pz/F");
    m_tree->Branch(TString(track_number) + "_pE", &m_track_pE[i], TString(track_number) + "_pE/F");
    m_tree->Branch(TString(track_number) + "_pT", &m_track_pT[i], TString(track_number) + "_pT/F");
    m_tree->Branch(TString(track_number) + "_eta", &m_track_eta[i], TString(track_number) + "_eta/F");
    m_tree->Branch(TString(track_number) + "_mass", &m_track_mass[i], TString(track_number) + "_mass/F");
    m_tree->Branch(TString(track_number) + "_mother_barcode", &m_track_mother_barcode[i], TString(track_number) + "_mother_barcode/I");
  }

  m_tree->Branch("delta_px", &m_delta_px, "delta_px/F");
  m_tree->Branch("delta_py", &m_delta_py, "delta_py/F");
  m_tree->Branch("delta_pz", &m_delta_pz, "delta_pz/F");
  m_tree->Branch("delta_pE", &m_delta_pE, "delta_pE/F");

  m_tree->Branch("accept_px_1percent", &m_accept_px_1percent, "accept_px_1percent/O");
  m_tree->Branch("accept_py_1percent", &m_accept_py_1percent, "accept_py_1percent/O");
  m_tree->Branch("accept_pz_1percent", &m_accept_pz_1percent, "accept_pz_1percent/O");
  m_tree->Branch("accept_pE_1percent", &m_accept_pE_1percent, "accept_pE_1percent/O");

  m_tree->Branch("accept_px_5percent", &m_accept_px_5percent, "accept_px_5percent/O");
  m_tree->Branch("accept_py_5percent", &m_accept_py_5percent, "accept_py_5percent/O");
  m_tree->Branch("accept_pz_5percent", &m_accept_pz_5percent, "accept_pz_5percent/O");
  m_tree->Branch("accept_pE_5percent", &m_accept_pE_5percent, "accept_pE_5percent/O");

  m_tree->Branch("accept_px_15percent", &m_accept_px_15percent, "accept_px_15percent/O");
  m_tree->Branch("accept_py_15percent", &m_accept_py_15percent, "accept_py_15percent/O");
  m_tree->Branch("accept_pz_15percent", &m_accept_pz_15percent, "accept_pz_15percent/O");
  m_tree->Branch("accept_pE_15percent", &m_accept_pE_15percent, "accept_pE_15percent/O");

  m_tree->Branch("accept_pT", &m_accept_pT, "accept_pT/O");
  m_tree->Branch("accept_eta", &m_accept_eta, "accept_eta/O");
}

void truthDecayTester::getMotherPDG(PHCompositeNode *topNode)
{
  PHNodeIterator nodeIter(topNode);

  std::string node_name = m_df_module_name + "_DecayMap";

  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(node_name.c_str()));
  if (findNode)
  {
    m_decayMap = findNode::getClass<DecayFinderContainer_v1>(topNode, node_name.c_str());
  }
  else
  {
  }

  std::vector<std::pair<int, int>> decay = m_decayMap->begin()->second;
  m_decay_pdg_id = decay[0].second;
}

std::vector<int> truthDecayTester::getDecayFinderMothers(PHCompositeNode *topNode)
{
  std::vector<int> m_motherBarcodes;

  PHNodeIterator nodeIter(topNode);

  std::string node_name = m_df_module_name + "_DecayMap";

  m_decayMap = findNode::getClass<DecayFinderContainer_v1>(topNode, node_name.c_str());

  for (DecayFinderContainer_v1::Iter iter = m_decayMap->begin(); iter != m_decayMap->end(); ++iter)
  {
    std::vector<std::pair<int, int>> decay = iter->second;
    m_nTracks = decay.size() - 1;
    for (unsigned int i = 0; i < decay.size(); ++i)
    {
      if (abs(decay[i].second) == abs(m_decay_pdg_id)) m_motherBarcodes.push_back(decay[i].first);
    }
  }

  return m_motherBarcodes;
}

bool truthDecayTester::isInRange(float min, float value, float max)
{
  return min <= value && value <= max;
}

void truthDecayTester::resetValues()
{
  m_delta_px = 0;
  m_delta_py = 0;
  m_delta_pz = 0;
  m_delta_pE = 0;

  m_accept_eta = true;
  m_accept_pT = true;
}

std::string truthDecayTester::get_histo_prefix()
{
  return std::string("h_") + Name() + std::string("_") + std::string("_");
}
