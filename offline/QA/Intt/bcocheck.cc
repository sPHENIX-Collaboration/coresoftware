#include "bcocheck.h"

using namespace std;

bcocheck::bcocheck(const string &name, const int run_num, const int felix_num)
  : SubsysReco(name)
{
  run_num_ = run_num;
  felix_num_ = felix_num;
  cout << "felix_num_=" << felix_num_ << "felix_num=" << felix_num << endl;
}

int bcocheck::Init(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning Init in bcocheck" << std::endl;
  }

  return 0;
}

int bcocheck::InitRun(PHCompositeNode *topNode)
{
  if (!topNode)
  {
    std::cout << "bcocheck::InitRun(PHCompositeNode* topNode)" << std::endl;
    std::cout << "\tCould not retrieve topNode; doing nothing" << std::endl;

    return 1;
  }

  cout << "felix_num_=" << felix_num_ << endl;

  // Initialize histograms
  for (int felix = 0; felix < kFelix_num_; felix++)
  {
    if (felix != felix_num_)
    {
      continue;
    }

    string name = "h2_bco_felix_" + to_string(felix);
    string title = name + Form("_Run%d", run_num_);
    h_full_bco[felix] = new TH1D(name.c_str(), title.c_str(), 128, 0, 128);
    h_full_bco[felix]->SetXTitle("BCO_FULL - BCO");
    h_full_bco[felix]->SetMinimum(0);
    tf_output_[felix] = new TFile(Form("./bco_000%d_intt%d.root", run_num_, felix), "RECREATE");
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int bcocheck::process_event(PHCompositeNode *topNode)
{
  // cout<<"1234"<<endl;

  // get raw hit
  string m_InttRawNodeName = "INTTRAWHIT";
  InttRawHitContainer *inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    cout << PHWHERE << endl;
    cout << "bcocheck::process_event(PHCompositeNode* topNode)" << endl;
    cout << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree" << endl;
    cout << "Exiting" << endl;
    gSystem->Exit(1);
    exit(1);
  }

  uint64_t longbco_full = inttcont->get_nhits() > 0
                              ? inttcont->get_hit(0)->get_bco()
                              : std::numeric_limits<uint64_t>::max();
  // uint64_t difevent_bcofull = (longbco_full &bit )-(long_prev_bcofull &bit);
  // h_interval->Fill(difevent_bcofull);
  uint64_t bco_full = longbco_full & 0x7FU;
  // cout<<"longbco_full="<<longbco_full<<endl;

  ievent_++;
  // cout<<"5678"<<endl;
  /*if(ievent_<500000){
    return Fun4AllReturnCodes::EVENT_OK;
  }*/

  if ((ievent_ % 100 == 0 && ievent_ < 1000) || ievent_ % 1000 == 0)
  {
    cout << "Process event #" << ievent_ << endl;
  }

  int nhits = inttcont->get_nhits();

  for (int i = 0; i < nhits; i++)
  {
    InttRawHit *intthit = inttcont->get_hit(i);

    int fnum = intthit->get_packetid() - 3001;  // packet id
    int bco = intthit->get_FPHX_BCO();          // FPHX bco

    int bco_diff = bco_full - bco;
    if (bco_diff < 0)
    {
      bco_diff = bco_diff + 128;
    }

    h_full_bco[fnum]->Fill(bco_diff);
  }
  // cout<<"910"<<endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int bcocheck::End(PHCompositeNode * /*topNode*/)
{
  // Mixupfraction();

  if (Verbosity() > 1)
  {
    std::cout << "Processing bcocheck done" << std::endl;
  }

  this->DrawHists();

  return Fun4AllReturnCodes::EVENT_OK;
}

void bcocheck::DrawHists()
{
  gStyle->SetOptStat(0);

  for (int felix = 0; felix < kFelix_num_; felix++)
  {
    if (felix != felix_num_)
    {
      continue;
    }

    tf_output_[felix]->cd();

    h_full_bco[felix]->Write();

    tf_output_[felix]->Close();
  }
}
