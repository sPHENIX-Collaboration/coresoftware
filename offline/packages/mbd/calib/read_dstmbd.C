//
// macro to read in MBD data after the waveforms are processed
// at this stage, we only have charge and time from each channel, not the full waveform
//
#include <mbd/MbdDefs.h>

const int NUM_PMT = MbdDefs::MBD_N_PMT;
const int NUM_ARMS = 2;

// Set up variables to read from TTree
Int_t   f_evt;
UShort_t f_clk;
UShort_t f_femclk;
Float_t f_tt[NUM_PMT];  // time from t-channels
Float_t f_tq[NUM_PMT];  // time from q-channels
Float_t f_q[NUM_PMT];   // charge

//Float_t f_bn[NUM_ARMS];   // num hit PMTs
Short_t f_bn[NUM_ARMS];   // num hit PMTs
Float_t f_bq[NUM_ARMS];   // chargesum
Float_t f_bt[NUM_ARMS];   // mean time in arm
Float_t f_bz;             // z-vertex
Float_t f_bt0;            // t-zero

TFile *tfile {nullptr};
TTree *tree {nullptr};

void read_dstmbd(const char *tfname = "beam_mbd-00009184-0000_mbd.root")
{
  cout << "tfname " << tfname << endl;

  const int NUM_PMT = 128;
  const int NUM_ARMS = 2;

  // Set up TTree
  int is_dst = 0; // whether reading from DST or private root files
  TString name;
  tfile = new TFile(tfname,"READ");
  tree = (TTree*)tfile->Get("t");
  if ( tree==nullptr )
  {
    tree = (TTree*)tfile->Get("T");
    is_dst = 1;
  }

  if ( is_dst==1 )
  {
    tree->SetBranchAddress("EvtSequence",&f_evt);
  }
  else
  {
    tree->SetBranchAddress("evt",&f_evt);
  }
  tree->SetBranchAddress("clk", &f_clk);
  tree->SetBranchAddress("femclk", &f_femclk);
  tree->SetBranchAddress("bns",&f_bn[0]);
  tree->SetBranchAddress("bnn",&f_bn[1]);
  tree->SetBranchAddress("bqs",&f_bq[0]);
  tree->SetBranchAddress("bqn",&f_bq[1]);
  tree->SetBranchAddress("bts",&f_bt[0]);
  tree->SetBranchAddress("btn",&f_bt[1]);
  tree->SetBranchAddress("bz",&f_bz);
  tree->SetBranchAddress("bt0",&f_bt0);
  if ( is_dst )
  {
    tree->SetBranchAddress("MbdPmtHits.bq",f_q);
    tree->SetBranchAddress("MbdPmtHits.btt",f_tt);
    tree->SetBranchAddress("MbdPmtHits.btq",f_tq);
  }
  else
  {
    for (int ipmt=0; ipmt<NUM_PMT; ipmt++)
    {
      name = "tt"; name += ipmt;
      tree->SetBranchAddress(name,&f_tt[ipmt]);
      name = "tq"; name += ipmt;
      tree->SetBranchAddress(name,&f_tq[ipmt]);
      name = "q"; name += ipmt;
      tree->SetBranchAddress(name,&f_q[ipmt]);
    }
  }


  // Event loop, each ientry is one triggered event
  int nentries = tree->GetEntries();
  nentries=10;
  for (int ientry=0; ientry<nentries; ientry++)
  {
    tree->GetEntry(ientry);

    if (ientry<10)
    {
      // print charge from channels 0 and 127
      cout << f_evt << "\t" << f_q[0] << "\t" << f_tt[0] << endl;
      cout <<  "\t" << f_q[127] << "\t" << f_tt[127] << endl;
    }

  }

}

