#include "InttMixupQA.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <ffarawobjects/InttRawHitContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // PHWHERE
#include <qautils/QAHistManagerDef.h>

#include <boost/format.hpp>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <TSystem.h>

#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>

using namespace std;

InttMixupQA::InttMixupQA(const string &name, const int run_num, const int felix_num)
  : SubsysReco(name)
{
  run_num_ = run_num;
  felix_num_ = felix_num;
  cout << "felix_num_=" << felix_num_ << "felix_num=" << felix_num << '\n';
}

InttMixupQA::~InttMixupQA() = default;

int InttMixupQA::Init(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning Init in InttMixupQA" << '\n';
  }

  return 0;
}

int InttMixupQA::InitRun(PHCompositeNode *topNode)
{
  if (!topNode)
  {
    std::cout << "InttMixupQA::InitRun(PHCompositeNode* topNode)" << '\n';
    std::cout << "\tCould not retrieve topNode; doing nothing" << '\n';

    return 1;
  }

  cout << "felix_num_=" << felix_num_ << '\n';

  // Initialize histograms
  for (int felix = 0; felix < kFelix_num_; felix++)
  {
    string name;

    name = (boost::format("%s_allmulti_intt%d") % getHistoPrefix() % felix).str();
    string title = name + Form("_Run%d", run_num_) + ": with clone cut";
    h_allmulti_[felix] = new TH1F(name.c_str(), title.c_str(), 200, 0, bin);
    h_allmulti_[felix]->SetXTitle("Multiplicity");
    // h_allmulti_[felix]->SetLineColor(felix+1);

    name = (boost::format("%s_allclone_intt%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_allclone_[felix] = new TH1F(name.c_str(), title.c_str(), 200, 0, 200);
    h_allclone_[felix]->SetXTitle("Clone multiplicity");
    h_allclone_[felix]->SetLineColor(felix + 1);

    // These next 4 will be plotted
    auto hm = QAHistManagerDef::getHistoManager();
    if (!hm)
    {
      std::cerr << PHWHERE << "\n"
                << "\tQAHistManagerDef::getHistoManager() returned null\n"
                << "\texit(1)" << '\n';
      gSystem->Exit(1);
      exit(1);
    }

    name = (boost::format("%s_bco_full&0x7F_prev_vs_bco_intt%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_vsprefull_bco[felix] = new TH2F(name.c_str(), title.c_str(), 128, 0, 128, 128, 0, 128);
    h_vsprefull_bco[felix]->SetXTitle("BCO");
    h_vsprefull_bco[felix]->SetYTitle("BCO_FULL previous event &0x7F");
    hm->registerHisto(h_vsprefull_bco[felix]);

    name = (boost::format("%s_bco_full&0x7F_vs_bco_intt%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_vsfull_bco[felix] = new TH2F(name.c_str(), title.c_str(), 128, 0, 128, 128, 0, 128);
    h_vsfull_bco[felix]->SetXTitle("BCO");
    h_vsfull_bco[felix]->SetYTitle("BCO_FULL &0x7F");
    hm->registerHisto(h_vsfull_bco[felix]);

    name = (boost::format("%s_bco_full&0x7F_prev_bco_intt%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_prefull_bco[felix] = new TH1F(name.c_str(), title.c_str(), 128, 0, 128);
    h_prefull_bco[felix]->SetXTitle("BCO_FULL previous event - BCO");
    h_prefull_bco[felix]->SetMinimum(0);
    hm->registerHisto(h_prefull_bco[felix]);

    name = (boost::format("%s_bco_full&0x7F_bco_intt%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_full_bco[felix] = new TH1F(name.c_str(), title.c_str(), 128, 0, 128);
    h_full_bco[felix]->SetXTitle("BCO_FULL - BCO");
    h_full_bco[felix]->SetMinimum(0);
    hm->registerHisto(h_full_bco[felix]);

    name = (boost::format("%s_bco_full&0x7F_prev_bco_intt_all%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_prefull_bco_all[felix] = new TH1F(name.c_str(), title.c_str(), 128, 0, 128);
    h_prefull_bco_all[felix]->SetXTitle("BCO_FULL previous event - BCO");
    h_prefull_bco_all[felix]->SetMinimum(0);

    name = (boost::format("%s_bco_full&0x7F_prev_vs_bco_intt_all%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_vsprefull_bco_all[felix] = new TH2F(name.c_str(), title.c_str(), 128, 0, 128, 128, 0, 128);
    h_vsprefull_bco_all[felix]->SetXTitle("BCO");
    h_vsprefull_bco_all[felix]->SetYTitle("BCO_FULL previous event &0x7F");

    string const name1 = "Mixup Multiplicity" + to_string(felix);
    string const name2 = "Mixup/all Multiplcity" + to_string(felix);
    string const title1 = name1 + Form("_Run%d", run_num_);
    string const title2 = name2 + Form("_Run%d", run_num_);
    for (int p = 0; p < divimul; p++)
    {
      h_mixupmulti[felix][p] = new TH1F((name1 + Form("_%d", p)).c_str(), (title1 + Form("_%d", p)).c_str(), 200, 0, bin);
      h_divmul[felix][p] = new TH1F((title1 + Form("_%d", p)).c_str(), (title2 + Form("_%d", p)).c_str(), 200, 0, bin);
      ;
    }

    name = (boost::format("%s_Number of mixup hit%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_mixup[felix] = new TH1F(name.c_str(), title.c_str(), 400, 0, 400);

    name = (boost::format("%s_prev_allhit vs Nmixup%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_prevsNmix[felix] = new TH2F(name.c_str(), title.c_str(), divimul + 10, 0, divimul + 10, bin, 0, bin);
    h_prevsNmix[felix]->GetXaxis()->SetNdivisions(405);
    h_prevsNmix[felix]->GetYaxis()->SetNdivisions(405);
    h_prevsNmix[felix]->SetXTitle("Number of Mixup hits");
    h_prevsNmix[felix]->SetYTitle("Number of previous event hits");

    name = (boost::format("%s_bco_full_prev-bco No copyhit%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_nocopyhit[felix] = new TH1F(name.c_str(), title.c_str(), 128, 0, 128);
    h_nocopyhit[felix]->SetXTitle("bco_full_prev-bco");

    name = (boost::format("%s_bco_full_prev-bco copyhit%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_copyhit[felix] = new TH1F(name.c_str(), title.c_str(), 128, 0, 128);
    h_copyhit[felix]->SetXTitle("bco_full_prev-bco");

    name = (boost::format("%s_Mixup & copy hit %d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_mixcopy[felix] = new TH1F(name.c_str(), title.c_str(), 50, 0, 50);
    h_mixcopy[felix]->SetXTitle("Number of Mixup copy hit");

    name = (boost::format("%s_Nmixup vs Nclone%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_mixvscopy[felix] = new TH2F(name.c_str(), title.c_str(), 100, 0, 100, 50, 0, 50);
    h_mixvscopy[felix]->SetXTitle("Nmixup");
    h_mixvscopy[felix]->SetYTitle("Ncopy");

    name = (boost::format("%s_Nmixup vs pre_Nhits%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_hitfra[felix] = new TH2F(name.c_str(), title.c_str(), bin, 0, bin, bin, 0, bin);

    name = (boost::format("%s_Nmixup vs pre_Nhits others 4bin%d") % getHistoPrefix() % felix).str();
    title = name + Form("_Run%d", run_num_);
    h_bghit[felix] = new TH2F(name.c_str(), title.c_str(), bin, 0, bin, bin, 0, bin);

    fNhit[felix].open(Form("./txtfile/NHit_%d_intt%d.txt", run_num_, felix));

    bcopeak_file[felix] = bcopeak_dir_ + Form("bco_000%d_intt%d.root", run_num_, felix);
  }

  string name = "bco_full - prev_bco_full";
  string title = name + Form("_Run%d", run_num_);
  h_interval = new TH1F(name.c_str(), title.c_str(), 200, 0, bin2);
  h_interval->SetXTitle("bco_full - prev_bco_full");

  name = "bco_full - prev_bco_full Mixup";
  title = name + Form("_Run%d", run_num_);
  h_mixinterval = new TH1F(name.c_str(), title.c_str(), 200, 0, bin2);

  name = "Interval normarize";
  title = name + Form("_Run%d", run_num_);
  h_divinter = new TH1F(name.c_str(), title.c_str(), 200, 0, bin2);
  h_divinter->SetXTitle("bco_full - prev_bco_full");

  name = "bco_full";
  title = name + Form("_Run%d", run_num_);
  h_bcofull_7 = new TH1F(name.c_str(), title.c_str(), 128, 0, 128);
  h_bcofull_7->SetXTitle("bco_full &0x7F");

  name = "bco";
  title = name + Form("_Run%d", run_num_);
  h_bco = new TH1F(name.c_str(), title.c_str(), 128, 0, 128);
  h_bco->SetXTitle("bco");

  name = "felix vs NmixupEv";
  title = name + Form("_Run%d", run_num_);
  // h_NmixEv=new TH2F(name.c_str(), title.c_str(),8,0,8,100000,0,100000);
  h_NmixEv = new TH1F(name.c_str(), title.c_str(), 8, 0, 8);

  name = "All Event";
  title = name + Form("_Run%d", run_num_);
  h_AllEv = new TH1F(name.c_str(), title.c_str(), 8, 0, 8);

  g_evfraction = new TGraph(n);

  g_cloevfraction = new TGraph(n);

  g_hitfraction = new TGraph(n);

  g_copyfraction = new TGraph(n);

  // bcopeak file set
  // bcopeak_file=bcopeak_dir_+Form("bco_000%d",run_num_)+this->GetFileSuffix()+".root";

  // hotchan file set
  // hotchan_file=hotchan_dir_+Form("InttHotDeadMap_000%d_50",run_num_)+this->GetFileSuffix()+".txt";

  tf_output_ = new TFile(output_root_.c_str(), "RECREATE");

  GetBcopeak();

  Readpeak();

  // Hotchancut();

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttMixupQA::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 5)
  {
    cout << "prev_bcofull=" << prev_bcofull << " "
         << "pre_allhit=" << pre_allhit << '\n';
  }

  // get event bco_full
  /* string node_name_intteventheader = "INTTEVENTHEADER";
   InttEventInfo *node_intteventheader_map_ = findNode::getClass<InttEventInfo>(topNode, node_name_intteventheader);

   if (!node_intteventheader_map_){
     cerr << PHWHERE << node_name_intteventheader << " node is missing." << endl;
     return Fun4AllReturnCodes::ABORTEVENT;
   }
   uint64_t longbco_full = node_intteventheader_map_->get_bco_full();
   uint64_t bco_full =longbco_full &0x7FU;
   uint64_t difevent_bcofull = (longbco_full &bit )-(long_prev_bcofull &bit);
   h_interval->Fill(difevent_bcofull);
   h_bcofull_7->Fill(bco_full);

   if(Verbosity()>5){
   cout<<"bco_full="<<bco_full<<endl;
   }*/

  // get raw hit
  string const m_InttRawNodeName = "INTTRAWHIT";
  InttRawHitContainer *inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    cout << PHWHERE << '\n';
    cout << "InttMixupQA::process_event(PHCompositeNode* topNode)" << '\n';
    cout << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree" << '\n';
    cout << "Exiting" << '\n';
    gSystem->Exit(1);
    exit(1);
  }

  uint64_t const longbco_full = (0 < inttcont->get_nhits()) ? inttcont->get_hit(0)->get_bco() : std::numeric_limits<uint64_t>::max();
  uint64_t const difevent_bcofull = (longbco_full & bit) - (long_prev_bcofull & bit);
  h_interval->Fill(difevent_bcofull);
  uint64_t const bco_full = longbco_full & 0x7FU;
  // cout<<"longbco_full="<<longbco_full<<endl;

  ievent_++;

  if (ievent_ == 1)
  {
    first_bcofull = longbco_full;
  }

  /*if(ievent_<500000){
    return Fun4AllReturnCodes::EVENT_OK;
  }*/

  if ((ievent_ % 100 == 0 && ievent_ < 1000) || ievent_ % 1000 == 0)
  {
    cout << "Process event #" << ievent_ << '\n';
  }

  int const nhits = inttcont->get_nhits();
  if (Verbosity() > 5)
  {
    cout << "Nhits = " << nhits << '\n';
  }

  map<int, int> map_hit;
  int nhit_fx[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int ncln_fx[8] = {0, 0, 0, 0, 0, 0, 0, 0};

  int Nmixup[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int Nother[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  // int Nclone[kFelix_num_] ={0,0,0,0,0,0,0,0};
  int Ncopy[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  // Loop over all INTTRAWHIT
  for (int i = 0; i < nhits; i++)
  {
    InttRawHit *intthit = inttcont->get_hit(i);

    int const fnum = intthit->get_packetid() - 3001;  // packet id
    int const fchn = intthit->get_fee();              // module
    int const adc = intthit->get_adc();               // adc
    int const chip = intthit->get_chip_id();          // chip
    int const chan = intthit->get_channel_id();       // channel
    int const bco = intthit->get_FPHX_BCO();          // FPHX bco
    int const hitID = 100000000 * (fnum + 1) + 1000000 * fchn + 10000 * chip + 10 * chan + adc;
    // int hitID_h = 10000000 * (fnum+1) + 100000 * fchn + 1000 * chip + chan;

    // clone hit cut
    auto itrHit = map_hit.find(hitID);
    if (itrHit == map_hit.end())
    {
      map_hit.insert(make_pair(hitID, 0));
      if (Verbosity() > 5)
      {
        cout << hitID << " " << fnum + 1 << " " << fchn << " " << chip << " " << chan << " " << adc << " " << bco << '\n';
      }

      // hot channel cut
      // auto itrHit_h=hotmap.find(hitID_h);
      // if (itrHit_h==hotmap.end())
      {
        nhit_fx[fnum]++;

        h_bco->Fill(bco);

        int bco_diff = bco_full - bco;
        // int bco_diff = bco - bco_full;
        if (bco_diff < 0)
        {
          bco_diff = bco_diff + 128;
        }

        int pre_bco_diff = prev_bcofull - bco;
        // int pre_bco_diff =bco - prev_bcofull ;
        if (pre_bco_diff < 0)
        {
          pre_bco_diff = pre_bco_diff + 128;
        }

        /*if(bco_diff>21 && bco_diff<80)//for Run46681
        {
          cout<<"bco_diff="<<bco_diff<<endl;
        }
        else*/
        {
          if (bcopar_[fnum].size() > 0 && bcopar_[fnum].find(bco_diff) == bcopar_[fnum].end())
          {
            if (bcopar_[fnum].find(pre_bco_diff) != bcopar_[fnum].end())
            {
              Nmixup[fnum]++;

              auto pre_itrHit = premap_hit.find(hitID);
              if (pre_itrHit != premap_hit.end())
              {
                h_copyhit[fnum]->Fill(pre_bco_diff);
                Ncopy[fnum]++;  // mixup & copy hit
              }
              else
              {
                // mixup & not copy hit
                h_nocopyhit[fnum]->Fill(pre_bco_diff);
              }
            }
            else if (otbcopar_[fnum].find(pre_bco_diff) != otbcopar_[fnum].end())
            {
              Nother[fnum]++;
            }

            h_prefull_bco[fnum]->Fill(pre_bco_diff);  // without this event collision hit
            h_vsprefull_bco[fnum]->Fill(bco, prev_bcofull);
          }
          h_prefull_bco_all[fnum]->Fill(pre_bco_diff);  // with this event collision hit
          h_vsprefull_bco_all[fnum]->Fill(bco, prev_bcofull);

          h_full_bco[fnum]->Fill(bco_diff);
          h_vsfull_bco[fnum]->Fill(bco, bco_full);
        }
      }
      /*else
      {
        cout<<"Hot channel"<<endl;
      }*/
    }
    else
    {
      ncln_fx[fnum]++;  // Number of clone hit
    }

    if (Verbosity() > 5)
    {
      cout << "bco_full=" << bco_full << " "
           << "bco=" << bco << " " << '\n';
    }
  }

  if (Verbosity() > 5)
  {
    for (int id = 0; id < 8; id++)
    {
      cout << "Felix#" << id << " "
           << "Nmixup=" << Nmixup[id] << " "
           << "Ncopy=" << Ncopy[id] << " "
           << "Nclone=" << ncln_fx[id] << '\n';
    }
  }

  for (int felix = 0; felix < kFelix_num_; felix++)
  {
    // h_allmulti_[felix]->Fill(nhit_fx[felix]);
    h_allmulti_[felix]->Fill(pre_allhit[felix]);
    h_allclone_[felix]->Fill(ncln_fx[felix]);
    h_mixup[felix]->Fill(Nmixup[felix]);
    h_mixcopy[felix]->Fill(Ncopy[felix]);
    h_hitfra[felix]->Fill(Nmixup[felix], pre_allhit[felix]);
    h_bghit[felix]->Fill(Nother[felix], pre_allhit[felix]);
    for (int p = 0; p < divimul; p++)
    {
      if (Nmixup[felix] == p + 1)
      // if (Nmixup[felix] > p * 4 && Nmixup[felix] <= (p + 1) * 4)
      {
        h_mixupmulti[felix][p]->Fill(pre_allhit[felix]);

        /*if (1500 < pre_allhit && pre_allhit < 2500)
        {
          // h_fullmixup[fnumber]->Fill(difevent_bcofull);
        }*/
      }
      else
      {
        continue;
      }

      h_divmul[felix][p]->Divide(h_mixupmulti[felix][p], h_allmulti_[felix]);
    }

    mixupfraction[felix] = (double(Nmixup[felix]) / double(pre_allhit[felix] + Nmixup[felix]));
    // copyfraction[felix] = (double(Ncopy[felix]) / double(Nmixup[felix])) * 100;
    // thisclonefraction[felix]=(double(Nthisclone)/double(map_hit.size()))*100;

    h_AllEv->Fill(felix);
    if (Nmixup[felix] > 0 && pre_allhit[felix] > 0)
    // if (Nmixup[felix] > 0 && pre_allhit[felix]>0 )
    {
      // cout<<"Nmixup="<<Nmixup[felix]<<" "<<"mixupfraction="<<mixupfraction[felix]<<endl;
      h_mixinterval->Fill(difevent_bcofull);
      h_divinter->Divide(h_mixinterval, h_interval);

      h_prevsNmix[felix]->Fill(Nmixup[felix], pre_allhit[felix]);
      NmixupEv[felix]++;
      h_NmixEv->Fill(felix);

      // mixupfraction[felix]= (double(Nmixup[felix]) / double(pre_allhit[felix] + Nmixup[felix])) * 100;
      // copyfraction[felix] = (double(Ncopy[felix]) / double(Nmixup[felix])) * 100;
      mixupfraction_sum[felix] = mixupfraction_sum[felix] + mixupfraction[felix];
      Nmixup_sum[felix] = Nmixup_sum[felix] + Nmixup[felix];
      pre_allhit_sum[felix] = pre_allhit_sum[felix] + pre_allhit[felix];
      if (Ncopy[felix] > 0)
      {
        NmixcopyEv[felix]++;
        copyfraction_sum[felix] = copyfraction_sum[felix] + copyfraction[felix];
      }
      if (Verbosity() > 5)
      {
        cout << "Felix=" << felix << " "
             << "mixupfraction=" << mixupfraction[felix] << "Nmixup=" << Nmixup[felix] << "pre_allhit=" << pre_allhit[felix] << "Ncopy=" << Ncopy[felix] << "copyfraction=" << copyfraction[felix] << '\n';
      }
    }
    // h_divinter->Divide(h_mixinterval,h_interval);
    if (felix == felix_num_)
    {
      // cout<<nhit_fx[felix]<<endl;

      fNhit[felix] << nhit_fx[felix] << '\n';
    }
    pre_allhit[felix] = nhit_fx[felix];
    // Initialize pre_allhit array

    /* pre_allhit[felix] = 0;

     pre_allhit[felix] = std::count_if(map_hit.begin(), map_hit.end(),
     [felix](const std::pair<int, int>& hit) {
     int felix_id = hit.first / 100000000 - 1;
     return felix_id == felix;
     });*/
  }

  long_prev_bcofull = longbco_full;
  prev_bcofull = bco_full;
  premap_hit = map_hit;
  //  h_allmulti_[fnum]->Fill(nhits);
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttMixupQA::End(PHCompositeNode * /*topNode*/)
{
  // Mixupfraction();

  if (Verbosity() > 1)
  {
    std::cout << "Processing InttMixupQA done" << '\n';
  }

  /*for (int felix=0;felix<8;felix++){
  cout<<"felix="<< NmixupEv[felix]<<endl;
  }*/
  this->DrawHists();

  /*for( auto& hist : h_allmulti_ ) {
    tf_output_->WriteTObject( hist, hist->GetName() );
  }

  tf_output_->Close();*/

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttMixupQA::SetOutputDir(std::string const &dir)
{
  output_dir_ = dir;

  string const run_num_str = string(8 - to_string(run_num_).size(), '0') + to_string(run_num_);
  string const fnum_str = to_string(felix_num_);

  output_root_ = output_dir_ + "root/" + output_basename_ + run_num_str + "intt" + fnum_str + ".root";
  output_pdf_ = output_dir_ + "plots/" + output_basename_ + run_num_str + "intt" + fnum_str + ".pdf";

  // output_root_ = output_dir_ + "root/" + output_basename_ + run_num_str +".root";
  // output_pdf_  = output_dir_ + "plots/" + output_basename_ + run_num_str +".pdf";

  output_txt_ = output_dir_ + "txtfile/";

  return Fun4AllReturnCodes::EVENT_OK;
}

/*int InttMixupQA::SetHistBin(string type)
{
  if (type==p)
  {
    bin=400;
    bin2=200000;
    bit=0x1FFF;//13bit
  }
  else if(type==Au)
  {
    bin=6000;
    bin2=200000;
    bit=0x7FFFFF;//21bit
  }
  else
  {
      cout<<"No set bin type"<<endl;

      return 1;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}*/

// Private functions
void InttMixupQA::GetBcopeak()
{
  // get bco finder file
  /*tf_bcopeak_=new TFile(bcopeak_file.c_str(),"READ");
  cout<<tf_bcopeak_<<endl;*/

  // get hist from INTT QA
  /*for(int i=0; i<8; i++){
   h2_bco_felix[i]=(TH2D*)gROOT->FindObject(Form("h2_bco_felix_%d",i));
   hbco[i]=h2_bco_felix[i]->ProjectionY(Form("hbco_%d",i));

   //default peakbin for can't found peak bin
   h2_bco_felix_sub[i]=(TH2D*)gROOT->FindObject(Form("h2_bco_felix_cut%d",i));
   hbco_sub[i]=h2_bco_felix_sub[i]->ProjectionY(Form("hbco_sub%d",i));
   int maxbin_sub=hbco_sub[i]->GetMaximumBin();
   DEFAULT_BCO_VALUE[i]=hbco_sub[i]->GetBinLowEdge(maxbin_sub);
   }*/

  // get hist from made by myself
  for (int i = 0; i < 8; i++)
  {
    tf_bcopeak_[i] = new TFile(bcopeak_file[i].c_str(), "READ");
    // cout<<tf_bcopeak_[i]<<endl;
    hbco[i] = dynamic_cast<TH1D *>(gROOT->FindObject(Form("h2_bco_felix_%d", i)));
    if (!hbco[i])
    {
      if (Verbosity())
      {
        std::cerr << PHWHERE << '\n';
      }
      continue;
    }

    // default peakbin for can't found peak bin
    int const maxbin_sub = hbco[i]->GetMaximumBin();
    DEFAULT_BCO_VALUE[i] = hbco[i]->GetBinLowEdge(maxbin_sub);
  }

  // check_bcohist
  int peakbco[8];
  double bg_mean[8];
  double bg_rms[8];
  /////////////////////
  // find peak and BG area
  for (int id = 0; id < 8; id++)
  {
    double const max = hbco[id]->GetMaximum();
    int const maxbin = hbco[id]->GetMaximumBin();
    peakbco[id] = hbco[id]->GetBinLowEdge(maxbin);
    cout << id << " " << max << " " << maxbin << " " << peakbco[id] << '\n';

    hbcohist[id] = new TH1F(Form("hbcohist_%d", id), Form("hbcohist_%d", id), 200, 0, max * 1.1);
    int const N = hbco[id]->GetNbinsX();
    for (int ibin = 0; ibin < N; ibin++)
    {
      double const cont = hbco[id]->GetBinContent(ibin + 1);
      hbcohist[id]->Fill(cont);
    }
  }

  for (int id = 0; id < 8; id++)
  {
    double const max = hbco[id]->GetMaximum();
    hbcohist2[id] = new TH1F(Form("hbcohist2_%d", id), Form("hbcohist2_%d", id), 200, 0, max * 1.1);

    int const N = hbco[id]->GetNbinsX();
    double const mean = hbcohist[id]->GetMean();
    double const rms = hbcohist[id]->GetRMS();

    for (int ibin = 0; ibin < N; ibin++)
    {
      double const cont = hbco[id]->GetBinContent(ibin + 1);
      if (mean + 3 * rms > cont)
      {
        hbcohist2[id]->Fill(cont);
      }
    }
    bg_mean[id] = hbcohist2[id]->GetMean();
    bg_rms[id] = hbcohist2[id]->GetRMS();
  }

  ///////
  // check peak area w/ peak +-10bin
  {
    f_felixpeak.open((output_txt_ + Form("bco_%d.txt", run_num_)).c_str());

    for (int id = 0; id < 8; id++)
    {
      f_felixpeak << id << " : ";
      // double thre   = bg_mean[id] + 6*bg_rms[id];
      double const thre = 5 * bg_mean[id];
      ;
      int const maxbin = hbco[id]->GetMaximumBin();

      int const N = hbco[id]->GetNbinsX();
      bool peak_found = false;

      for (int ibin = 0; ibin < 10; ibin++)
      {
        int binid = ibin + maxbin - 10;
        if (binid < 1)
        {
          binid += N;  // binid+=128;
        }
        double const cont = hbco[id]->GetBinContent(binid + 1);
        cout << id << " :  " << ibin << " " << binid << " " << maxbin;
        if (cont > thre)
        {
          int const peakbin = hbco[id]->GetBinLowEdge(binid + 1);
          cout << " exceed : " << peakbin;
          f_felixpeak << peakbin << " ";
          peak_found = true;
        }
        cout << '\n';
      }
      for (int ibin = 0; ibin < 10; ibin++)
      {
        int binid = ibin + maxbin;
        if (binid >= N)
        {
          binid -= N;  // binid+=128;
        }
        double const cont = hbco[id]->GetBinContent(binid + 1);
        cout << id << " :  " << ibin << " " << binid << " " << maxbin;
        if (cont > thre)
        {
          int const peakbin = hbco[id]->GetBinLowEdge(binid + 1);
          cout << " exceed : " << peakbin;
          f_felixpeak << peakbin << " ";
          peak_found = true;
        }
        cout << '\n';
      }
      if (!peak_found)
      {
        f_felixpeak << DEFAULT_BCO_VALUE[id] << " ";
      }

      f_felixpeak << '\n';
    }
    f_felixpeak.close();
  }

  // get min
  float minimum = 1000000;
  for (auto &id : hbco)
  {
    cout << id->GetMinimumBin() << '\n';
    cout << id->GetBinContent(id->GetMinimumBin()) << '\n';
    double const min = id->GetBinContent(id->GetMinimumBin());
    cout << min << '\n';
    if (min < minimum)
    {
      minimum = min;
    }
  }

  // Draw hist

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 700);
  c1->Divide(4, 4);
  for (int id = 0; id < 8; id++)
  {
    c1->cd(id + 1);
    gPad->SetLogy();
    hbco[id]->SetMinimum(minimum * 0.2);
    hbco[id]->Draw();

    TLine *la = new TLine(0, 5 * bg_mean[id], 128, 5 * bg_mean[id]);
    la->Draw("same");
    TLine *lb = new TLine(0, bg_mean[id] + 6 * bg_rms[id], 128, bg_mean[id] + 6 * bg_rms[id]);
    lb->SetLineColor(2);
    lb->Draw("same");
  }

  for (int id = 0; id < 8; id++)
  {
    c1->cd(id + 9);
    gPad->SetLogy();
    double const max = hbcohist[id]->GetMaximum();

    hbcohist[id]->Draw();
    hbcohist2[id]->SetLineColor(2);
    hbcohist2[id]->Draw("same");

    TLine *lm = new TLine(bg_mean[id], 0.001, bg_mean[id], max * 1.1);
    lm->Draw();
    TLine *l5m = new TLine(5 * bg_mean[id], 0.001, 5 * bg_mean[id], max * 1.1);
    l5m->Draw();
    TLine *l = new TLine(bg_mean[id] + 6 * bg_rms[id], 0.001, bg_mean[id] + 6 * bg_rms[id], max * 1.1);
    l->SetLineColor(2);
    l->Draw();
  }

  c1->Print((output_dir_ + "plots/" + Form("check_bco_file_%d.pdf", run_num_)).c_str());
}

void InttMixupQA::Readpeak()
{
  std::string const bcofile_ = (output_txt_ + Form("bco_%d.txt", run_num_));

  if (bcofile_.size() > 0)
  {
    cout << "BCO file : " << bcofile_.c_str() << '\n';
    ifstream file(bcofile_.c_str());
    if (!file.is_open())
    {
      std::cerr << "Failed to open file: " << bcofile_ << '\n';
      return;
    }

    string sline;

    while (getline(file, sline))
    {
      // cout<<sline<<endl;
      istringstream ss(sline);
      string sbuf;

      int felix = -1;
      int idx = 0;
      while (getline(ss, sbuf, ':'))
      {
        // 1st for felix id
        if (idx == 0)
        {
          // cout<<idx<<", id "<<sbuf<<endl;
          felix = stoi(sbuf.c_str());
          if (felix < 0 || felix > 7)
          {
            cout << "felixid out of range. id=" << felix << '\n';
            break;
          }
        }
        else
        {
          // cout<<idx<<" "<<felix<<", bco "<<sbuf<<endl;
          istringstream ss2(sbuf);
          string sbuf2;
          string sbuf3;
          string sbuf4;
          string sbuf5;
          string sbuf6;
          while (ss2 >> sbuf2)
          {
            // cout<<"   "<<sbuf2<<endl;
            bcopar_[felix].insert(stoi(sbuf2));
            int const peak = stoi(sbuf2);
            int const peak1 = peak + 1;
            int const peak2 = peak + 2;
            int const peak_1 = peak - 1;
            int const peak_2 = peak - 2;
            sbuf3 = to_string(peak1);
            sbuf4 = to_string(peak2);
            sbuf5 = to_string(peak_1);
            sbuf6 = to_string(peak_2);
            otbcopar_[felix].insert(stoi(sbuf3));
            otbcopar_[felix].insert(stoi(sbuf4));
            otbcopar_[felix].insert(stoi(sbuf5));
            otbcopar_[felix].insert(stoi(sbuf6));
          }
        }
        idx++;
      }
    }
    file.close();
  }

  cout << "BCO peaks " << '\n';
  for (int i = 0; i < 8; i++)
  {
    cout << "    felix " << i << " : ";
    for (int const itr : bcopar_[i])
    {
      cout << itr << " ";
    }
    cout << '\n';
  }

  cout << "Around BCO peaks " << '\n';
  for (int i = 0; i < 8; i++)
  {
    cout << "    felix " << i << " : ";
    for (int const itr : otbcopar_[i])
    {
      cout << itr << " ";
    }
    cout << '\n';
  }
}

/*void InttMixupQA::Hotchancut()
{
  //tf_hotchan_=new TFile(hotchan_file.c_str(),"READ");
  ifstream tf_hotchan_(hotchan_file.c_str());
 // cout <<tf_hotchan_<< endl;

  //map<int,int> hotmap;
  string sline;
  int lineNumber=0;

  getline(tf_hotchan_, sline);

  while (getline(tf_hotchan_,sline))
  {
    istringstream ss(sline);
    int felix,felix_ch,chip,channel;

    ss>>felix>>felix_ch>>chip>>channel;

    int hotchanID=10000000 * (felix+1) + 100000 * felix_ch + 1000 * chip + channel;

    hotmap[lineNumber++]=hotchanID;

  }

  tf_hotchan_.close();

  for(const auto&[key,value]:hotmap){
    cout<<"Line"<<key+1<<":"<<value<<endl;
  }
}*/

void InttMixupQA::DrawHists()
{
  std::cout << "output pdf: " << output_pdf_ << '\n';

  /////////////////////////Mixup plot Draw//////////////////////////////
  TCanvas *c0 = new TCanvas("c0", "c0");
  c0->Print((output_pdf_ + "[").c_str());
  TCanvas *c1[8];
  TCanvas *c2[8];
  TCanvas *c3[8];
  // TCanvas *c4[8];
  TCanvas *c5[8];
  TCanvas *c6[8];
  TCanvas *c10[8];
  TCanvas *c7 = new TCanvas("canvas7", "interval", 1600, 1200);
  TCanvas *c8 = new TCanvas("canvas8", "fraction", 1600, 1200);
  TCanvas *c9 = new TCanvas("canvas9", "bco bco_full", 1600, 1200);

  gStyle->SetOptStat(0);

  for (int felix = 0; felix < kFelix_num_; felix++)
  {
    if (felix != felix_num_)
    {
      continue;
    }
    c1[felix] = new TCanvas(Form("canvas1_%d", felix), "bcofull&bco", 1600, 1200);
    c1[felix]->cd();
    gStyle->SetOptStat(0);
    c1[felix]->Divide(2, 2);

    c1[felix]->cd(1);
    h_vsfull_bco[felix]->Draw("colz");
    c1[felix]->cd(2);
    h_vsprefull_bco[felix]->Draw("colz");
    c1[felix]->cd(3);
    h_full_bco[felix]->SetMinimum(0);
    h_full_bco[felix]->Draw();
    c1[felix]->cd(4);
    h_prefull_bco[felix]->SetMinimum(0);
    h_prefull_bco[felix]->Draw();

    c2[felix] = new TCanvas(Form("canvas2_%d", felix), "Multiplicity", 1600, 1200);
    c2[felix]->cd();
    gStyle->SetOptStat(0);
    gPad->SetRightMargin(0.15);
    TLegend *legend = new TLegend(0.7, 0.45, 0.9, 0.9);
    gPad->SetLogy();
    h_allmulti_[felix]->GetXaxis()->SetNdivisions(405);
    h_allmulti_[felix]->GetYaxis()->SetNdivisions(405);
    // h_allmulti[i]->SetXTitle("Hit Multiplcity/event");
    int const a = h_allmulti_[felix]->GetMaximum();
    h_allmulti_[felix]->GetYaxis()->SetRangeUser(1, a * 1.6);
    h_allmulti_[felix]->SetXTitle("Number of hits of previous event");
    h_allmulti_[felix]->SetYTitle("Entry");
    h_allmulti_[felix]->Draw();
    h_allmulti_[felix]->SetLineColor(1);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_allmulti_[felix], "Number of", "");
    legend->AddEntry(h_allmulti_[felix], "Mixup hits", "");
    legend->AddEntry(h_allmulti_[felix], "All event", "l");

    for (int p = 0; p < 10; p++)
    {
      h_mixupmulti[felix][p]->Draw("same");
      int const color = p + 2;
      if (color == 10)
      {
        h_mixupmulti[felix][p]->SetLineColor(46);
      }
      else
      {
        h_mixupmulti[felix][p]->SetLineColor(color);
      }
      legend->AddEntry(h_mixupmulti[felix][p], Form("%d", p + 1), "l");
    }
    legend->Draw();

    c3[felix] = new TCanvas(Form("canvas3_%d", felix), "divide multiplicity", 1600, 1200);
    c3[felix]->cd();
    gStyle->SetOptStat(0);
    gPad->SetRightMargin(0.15);
    TLegend *leg = new TLegend(0.75, 0.5, 0.9, 0.9);
    gPad->SetLogy();
    for (int p = 0; p < 10; p++)
    {
      h_divmul[felix][p]->SetXTitle("Mixup/all hit Multiplicity");
      if (p == 0)
      {
        h_divmul[felix][p]->Draw();
      }
      else
      {
        h_divmul[felix][p]->Draw("same");
      }
      int const color2 = p + 2;
      if (color2 == 10)
      {
        h_divmul[felix][p]->SetLineColor(46);
      }
      else
      {
        h_divmul[felix][p]->SetLineColor(color2);
      }
      leg->AddEntry(h_divmul[felix][p], Form("%d", p + 1), "l");
    }

    /*c4[felix]= new TCanvas( Form("canvas4_%d",felix), "copy hit",1600, 1200 );
    c4[felix]->cd();
    gStyle->SetOptStat( 0 );
    c4[felix]->Divide(2,2);
    c4[felix]->cd(1);
    h_nocopyhit[felix]->Draw();
    c4[felix]->cd(2);
    h_copyhit[felix]->Draw();
    c4[felix]->cd(3);
    h_mixcopy[felix]->Draw();
    c4[felix]->cd(4);
    h_mixvscopy[felix]->Draw();*/

    c5[felix] = new TCanvas(Form("canvas5_%d", felix), "previous event hit vs Mixup hit", 1600, 1200);
    c5[felix]->cd();
    gStyle->SetOptStat(0);
    gPad->SetRightMargin(0.15);
    h_prevsNmix[felix]->Draw("colz");

    c6[felix] = new TCanvas(Form("canvas6_%d", felix), "Mixup", 1600, 1200);
    c6[felix]->cd();
    c6[felix]->Divide(2, 2);
    c6[felix]->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    h_mixup[felix]->Draw();

    c6[felix]->cd(4);
    h_prefull_bco_all[felix]->Draw();
    c6[felix]->cd(2);
    h_vsprefull_bco_all[felix]->Draw();

    c10[felix] = new TCanvas(Form("canvas10_%d", felix), "for fraction");
    c10[felix]->cd();
    c10[felix]->Divide(2, 1);
    c10[felix]->cd(1);
    h_hitfra[felix]->Draw("colz");
    c10[felix]->cd(2);
    h_bghit[felix]->Draw("colz");

    c1[felix]->Print(output_pdf_.c_str());
    c2[felix]->Print(output_pdf_.c_str());
    c3[felix]->Print(output_pdf_.c_str());
    // c4[felix]->Print(output_pdf_.c_str());
    c5[felix]->Print(output_pdf_.c_str());
    c6[felix]->Print(output_pdf_.c_str());
    c10[felix]->Print(output_pdf_.c_str());
    tf_output_->cd();
    /*c1[felix]->Write();
    c2[felix]->Write();
    c3[felix]->Write();
    c4[felix]->Write();
    c5[felix]->Write();
    c6[felix]->Write();*/

    h_vsprefull_bco[felix]->Write();
    h_vsfull_bco[felix]->Write();
    h_prefull_bco[felix]->Write();
    h_full_bco[felix]->Write();
    /*for (int p = 0; p < 10; p++)
    {
      h_allmulti_[felix]->Add(h_mixupmulti[felix][p]);
    }
    h_allmulti_[felix]->Write();*/
    c2[felix]->Write();
    h_mixup[felix]->Write();
    h_prevsNmix[felix]->Write();
    h_hitfra[felix]->Write();
    h_bghit[felix]->Write();
    h_NmixEv->Write();
    h_AllEv->Write();
    h_divinter->Write();
    // tf_output_->Close();
  }
  // h_interval->Write();
  // h_mixinterval->Write();

  c7->cd();
  gStyle->SetOptStat(0);
  TLegend *legsph3 = new TLegend(.55, .4, .65, .6);
  gPad->SetRightMargin(0.15);
  h_interval->Draw();
  h_mixinterval->Draw("same");
  h_mixinterval->SetLineColor(2);
  legsph3->AddEntry(h_interval, "All event", "l");
  legsph3->AddEntry(h_mixinterval, "Mixup event", "l");

  /*c8->cd();
  gStyle->SetOptStat( 0 );
  c8->Divide(2, 2);
  c8->cd(1);
  g_evfraction->SetMarkerStyle(20);
  g_evfraction->SetMarkerSize(1.1);
  g_evfraction->SetMinimum(0);
  g_evfraction->SetTitle(Form("Mixup event fraction Run%d;Felix;Mixup Event fraction ",run_num_));
  g_evfraction->Draw("AP");

  c8->cd(2);
  g_hitfraction->SetMarkerStyle(20);
  g_hitfraction->SetMarkerSize(1.1);
  g_hitfraction->SetMinimum(0);
  g_hitfraction->SetTitle(Form("Mixup Hit fraction Run%d;Felix;Mixup Hit fraction ",run_num_));
  g_hitfraction->Draw("AP");

  c8->cd(3);
  g_cloevfraction->SetMarkerStyle(20);
  g_cloevfraction->SetMarkerSize(1.1);
  g_cloevfraction->SetMinimum(0);
  g_cloevfraction->SetTitle(Form("Mixclone event fraction Run%d/Mixup;Felix;Mixclone event fraction ",run_num_));
  g_cloevfraction->Draw("AP");

  c8->cd(4);
  g_copyfraction->SetMarkerStyle(20);
  g_copyfraction->SetMarkerSize(1.1);
  g_copyfraction->SetMinimum(0);
  g_copyfraction->SetTitle(Form("Mixclone hit fraction Run%d/Mixup;Felix;clone hit from prev/Mixup ",run_num_));
  g_copyfraction->Draw("AP");

  fgraph->cd();
  g_hitfraction->Write("g_hitfraction");
  g_evfraction->Write("g_evfraction");
  g_copyfraction->Write("g_copyfraction");
  g_cloevfraction->Write("g_cloevfraction");
  fgraph->Close();*/

  c9->cd();
  gStyle->SetOptStat(0);
  c9->Divide(2, 1);
  c9->cd(1);
  h_bcofull_7->Draw();
  c9->cd(2);
  h_bco->Draw();

  tf_output_->cd();
  c7->Write();
  /*c8->Write();
  c9->Write();*/

  /*g_hitfraction->Write("g_hitfraction");
  g_evfraction->Write("g_evfraction");
  g_copyfraction->Write("g_copyfraction");
  g_cloevfraction->Write("g_cloevfraction");*/
  tf_output_->Close();

  c7->Print((output_pdf_).c_str());
  c8->Print((output_pdf_).c_str());
  c9->Print((output_pdf_).c_str());

  c9->Print((output_pdf_ + "]").c_str());
}
