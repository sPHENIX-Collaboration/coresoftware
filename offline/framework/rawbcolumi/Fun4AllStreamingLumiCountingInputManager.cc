#include "Fun4AllStreamingLumiCountingInputManager.h"

#include <fun4allraw/InputManagerType.h>
#include "SingleStreamingInputv2.h"

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1Packetv2.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllInputManager.h>  // for Fun4AllInputManager
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>

#include <ffaobjects/SyncObject.h>  // for SyncObject
#include <ffaobjects/SyncObjectv1.h>

#include <frog/FROG.h>
#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <boost/format.hpp>

#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>  // for max
#include <cassert>
#include <cstdint>  // for uint64_t, uint16_t
#include <cstdlib>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <utility>   // for pair

Fun4AllStreamingLumiCountingInputManager::Fun4AllStreamingLumiCountingInputManager(const std::string &name, const std::string &dstnodename, const std::string &topnodename)
  : Fun4AllInputManager(name, dstnodename, topnodename)
  , m_SyncObject(new SyncObjectv1())
{
  Fun4AllServer *se = Fun4AllServer::instance();
  m_topNode = se->topNode(TopNodeName());

  createLuminosityHistos();
  return;
}

Fun4AllStreamingLumiCountingInputManager::~Fun4AllStreamingLumiCountingInputManager()
{
  if (IsOpen())
  {
    fileclose();
  }
  //  std::cout<<"----Write? files to output.root"<<std::endl;

  delete m_SyncObject;
  // clear leftover raw event maps and vectors with poolreaders
  // GL1
  for (auto iter : m_Gl1InputVector)
  {
    delete iter;
  }

  m_Gl1InputVector.clear();
}

int Fun4AllStreamingLumiCountingInputManager::run(const int /*nevents*/)
{
  int iret = 0;
  if (m_gl1_registered_flag)  // Gl1 first to get the reference
  {
    iret += FillGl1();
  }
  /*   std::cout<<"----Storing files to output.root"<<std::endl;
    tfile = TFile::Open(m_outputFileName.c_str(), "RECREATE", "");//UPDATE
    ttree->Write("", TObject::kOverwrite);
    h_lumibco->Write("", TObject::kOverwrite);
    h_bunchnumber->Write("", TObject::kOverwrite);
    h_bunchnumber_occur->Write("", TObject::kOverwrite);
    tfile->Close();
    delete tfile;
  */
  return iret;
}

void Fun4AllStreamingLumiCountingInputManager::SetOutputFileName(const std::string &fileName)
{
  m_outputFileName = fileName;  // Update the filename
}

int Fun4AllStreamingLumiCountingInputManager::fileclose()
{
  //  std::cout<<"----fileclose()"<<std::endl;
  //  QAHistManagerDef::saveQARootFile(m_output_file);
  return 0;
}

void Fun4AllStreamingLumiCountingInputManager::Print(const std::string &what) const
{
  if (what == "ALL" || what == "INPUTFILES")
  {
    std::cout << "-----------------------------" << std::endl;
    for (const auto &iter : m_Gl1InputVector)
    {
      std::cout << "Single Streaming Input Manager " << iter->Name() << " reads run "
                << iter->RunNumber()
                << " from file " << iter->FileName()
                << std::endl;
    }
  }
  Fun4AllInputManager::Print(what);
  return;
}

int Fun4AllStreamingLumiCountingInputManager::ResetEvent()
{
  // zhiwan
  // m_RefBCO = 0;
  return 0;
}

int Fun4AllStreamingLumiCountingInputManager::PushBackEvents(const int /*i*/)
{
  return 0;
}

int Fun4AllStreamingLumiCountingInputManager::GetSyncObject(SyncObject **mastersync)
{
  // here we copy the sync object from the current file to the
  // location pointed to by mastersync. If mastersync is a 0 pointer
  // the syncobject is cloned. If mastersync allready exists the content
  // of syncobject is copied
  if (!(*mastersync))
  {
    if (m_SyncObject)
    {
      *mastersync = dynamic_cast<SyncObject *>(m_SyncObject->CloneMe());
      assert(*mastersync);
    }
  }
  else
  {
    *(*mastersync) = *m_SyncObject;  // copy syncobject content
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

int Fun4AllStreamingLumiCountingInputManager::SyncIt(const SyncObject *mastersync)
{
  if (!mastersync)
  {
    std::cout << PHWHERE << Name() << " No MasterSync object, cannot perform synchronization" << std::endl;
    std::cout << "Most likely your first file does not contain a SyncObject and the file" << std::endl;
    std::cout << "opened by the Fun4AllDstInputManager with Name " << Name() << " has one" << std::endl;
    std::cout << "Change your macro and use the file opened by this input manager as first input" << std::endl;
    std::cout << "and you will be okay. Fun4All will not process the current configuration" << std::endl
              << std::endl;
    return Fun4AllReturnCodes::SYNC_FAIL;
  }
  int iret = m_SyncObject->Different(mastersync);
  if (iret)
  {
    std::cout << "big problem" << std::endl;
    exit(1);
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

std::string Fun4AllStreamingLumiCountingInputManager::GetString(const std::string &what) const
{
  std::cout << PHWHERE << " called with " << what << " , returning empty string" << std::endl;
  return "";
}

void Fun4AllStreamingLumiCountingInputManager::registerStreamingInput(SingleStreamingInputv2 *evtin, InputManagerType::enu_subsystem system)
{
  evtin->StreamingLumiInputManager(this);
  // if the streaming flag is set, we only want the first event from the GL1 to
  // get the starting BCO of that run which enables us to dump all the junk which
  // is taken before the run starts in the streaming systems. But we don't want the
  // GL1 in the output, so we do not create its dst node if running in streaming
  if (system == InputManagerType::GL1)
  {
    if (!m_StreamingFlag)
    {
      evtin->CreateDSTNode(m_topNode);
    }
  }
  else
  {
    evtin->CreateDSTNode(m_topNode);
  }
  evtin->ConfigureStreamingInputManager();
  if (system == InputManagerType::GL1)
  {
    m_gl1_registered_flag = true;
    m_Gl1InputVector.push_back(evtin);
  }
  else
  {
    std::cout << "invalid subsystem flag " << system << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  if (Verbosity() > 3)
  {
    std::cout << "registering " << evtin->Name()
              << " number of registered inputs: "
              << m_Gl1InputVector.size()
              << std::endl;
  }
  std::cout << m_Gl1InputVector.size() << std::endl;
}

void Fun4AllStreamingLumiCountingInputManager::AddGl1RawHit(uint64_t bclk, Gl1Packet *hit)
{
  m_Gl1RawHitMap[bclk].Gl1RawHitVector.push_back(hit);
}

void Fun4AllStreamingLumiCountingInputManager::AddGl1Window(uint64_t bco_trim, int negative_window, int positive_window)
{
  m_BCOWindows[bco_trim] = std::make_pair(bco_trim - negative_window, bco_trim + positive_window);
}

void Fun4AllStreamingLumiCountingInputManager::AddGl1BunchNumber(uint64_t bco_trim, int bunch_number)
{
  m_BCOBunchNumber[bco_trim] = bunch_number;
}

int Fun4AllStreamingLumiCountingInputManager::FillGl1()
{
  // unsigned int alldone = 0;
  for (auto iter : m_Gl1InputVector)
  {
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllStreamingLumiCountingInputManager::FillGl1 - fill pool for " << iter->Name() << std::endl;
      std::cout << "Run number " << iter->RunNumber() << std::endl;
    }
    iter->FillPool();

    if (m_RunNumber == 0)
    {
      m_RunNumber = iter->RunNumber();
      SetRunNumber(m_RunNumber);
    }
    else
    {
      if (m_RunNumber != iter->RunNumber())
      {
        std::cout << PHWHERE << " Run Number mismatch, run is "
                  << m_RunNumber << ", " << iter->Name() << " reads "
                  << iter->RunNumber() << std::endl;
        std::cout << "You are likely reading files from different runs, do not do that" << std::endl;
        Print("INPUTFILES");
        gSystem->Exit(1);
        exit(1);
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "Here BCO " << m_BCOWindows.begin()->first << " left " << m_BCOWindows.begin()->second.first << " right " << m_BCOWindows.begin()->second.second << std::endl;
  }
  /*
  for (const auto &entry : m_BCOWindows) {
      uint64_t key = entry.first;
      uint64_t valueFirst = entry.second.first;
      uint64_t valueSecond = entry.second.second;
      std::cout << "Key: " << key
                << ", Value First: " << valueFirst
                << ", Value Second: " << valueSecond
                << std::endl;
  }

  for (const auto& [bco_trim, bunch_number] : m_BCOBunchNumber) {
      std::cout << "Here BCO " << bco_trim << " Bunch Number " << bunch_number << std::endl;
  }
  */
  // std::cout << "Here BCO " <<m_BCOWindows.begin()->first <<std::endl;

  if (m_BCOWindows.size() > 1)
  {
    auto first_element = m_BCOWindows.begin();
    auto second_element = std::next(m_BCOWindows.begin());
    //	std::cout<<"Key 1: "<<first_element->first<<" Value ( "<<first_element->second.first<<" , "<<first_element->second.second<<std::endl;
    //	std::cout<<"Key 2: "<<second_element->first<<" Value ( "<<second_element->second.first<<" , "<<second_element->second.second<<std::endl;
    // 	std::cout<<second_element->first - first_element->first<<" compared with window "<< m_negative_bco_window+m_positive_bco_window <<std::endl;
    // constexpr uint64_t MAX_40BIT_VALUE = 0xFFFFFFFFFF;

    // special case for overflow: second_element - first_element > 1099511000000, then switch them
    m_diffBCO = second_element->first - first_element->first;

    if (second_element->first - first_element->first > 1099510000000)
    {
      flat_overflow = true;
      // int temp_m_diffBCO=first_element->first+1099511627775+1-second_element->first;
      bco_temp = first_element->first;
      m_BCOWindows.erase(m_BCOWindows.begin());
      bco_temp += 1099511627775 + 1;
      m_BCOWindows[bco_temp] = std::make_pair(bco_temp - m_negative_bco_window, bco_temp + m_positive_bco_window);
      first_element = m_BCOWindows.begin();
      second_element = std::next(m_BCOWindows.begin());
      m_diffBCO = second_element->first - first_element->first;
      std::cout << "overflow new diff " << m_diffBCO << " new first element " << first_element->first << " new second element " << second_element->first << std::endl;
    }
    h_diffbco->Fill(m_diffBCO);
    if (m_diffBCO < static_cast<int>(m_negative_bco_window + m_positive_bco_window))
    {
      m_BCOWindows.begin()->second.second = second_element->second.first;
      std::cout << "*** new Key 1 BCO " << m_BCOWindows.begin()->first << " left " << m_BCOWindows.begin()->second.first << " right " << m_BCOWindows.begin()->second.second << std::endl;
    }
  }

  m_bco_trim = m_BCOWindows.begin()->first;
  m_lower_bound = m_BCOWindows.begin()->second.first;
  m_upper_bound = m_BCOWindows.begin()->second.second;
  m_bunch_number = m_BCOBunchNumber[m_BCOWindows.begin()->first];
  // ttree->Fill();
  h_bunchnumber->Fill(m_BCOBunchNumber[m_BCOWindows.begin()->first]);
  h_lumibco->Fill(m_BCOWindows.begin()->second.second - m_BCOWindows.begin()->second.first);

  int lower = -1 * static_cast<int>(m_bco_trim - m_lower_bound);
  int upper = (m_upper_bound > m_bco_trim) ? static_cast<int>(m_upper_bound - m_bco_trim) : -1 * static_cast<int>(m_bco_trim - m_upper_bound);  // it is possible that upper is <0
                                                                                                                                                // std::cout<<"lower="<<lower<<", upper = "<<upper<<std::endl;//<<" or upper2 = "<<lower+(m_BCOWindows.begin()->second.second - m_BCOWindows.begin()->second.first)<<std::endl;
  for (int i = lower; i < upper; i++)
  {
    int adjusted_bunch = m_bunch_number + i;
    while (adjusted_bunch < 0)
    {
      adjusted_bunch += 120;
    }
    while (adjusted_bunch > 119)
    {
      adjusted_bunch -= 120;
    }
    if (i != 0)
    {
      h_bunchnumber_occur->Fill(adjusted_bunch);
    }  // else{std::cout<<"same gl1 removed"<<std::endl;}
  }

  if (!m_BCOBunchNumber.empty())
  {
    m_BCOBunchNumber.erase(m_BCOWindows.begin()->first);
    // m_BCOBunchNumber.erase(m_BCOBunchNumber.begin());
  }
  if (!m_BCOWindows.empty())
  {
    m_BCOWindows.erase(m_BCOWindows.begin());
  }
  if (flat_overflow)
  {
    m_BCOWindows.erase(m_BCOWindows.begin());
    bco_temp -= 1099511627775 + 1;
    m_BCOWindows[bco_temp] = std::make_pair(bco_temp - m_negative_bco_window, bco_temp + m_positive_bco_window);
    std::cout << " Change back, new bco window map  " << m_BCOBunchNumber.begin()->first << std::endl;
    flat_overflow = false;
  }

  // mow use new

  Gl1Packet *gl1packet = findNode::getClass<Gl1Packet>(m_topNode, "GL1RAWHIT");
  for (auto gl1hititer : m_Gl1RawHitMap.begin()->second.Gl1RawHitVector)
  {
    if (!m_StreamingFlag)  // if streaming flag is set, the gl1packet is a nullptr
    {
      gl1packet->FillFrom(gl1hititer);
      MySyncManager()->CurrentEvent(gl1packet->getEvtSequence());
    }
  }

  // add for mbd p_gl1
  Gl1Packet *p_gl1 = findNode::getClass<Gl1Packetv2>(m_topNode, "GL1RAWHIT");  //"GL1Packet");
  if (!p_gl1)
  {
    std::cout << "CAN not find this Gl1Packetv2" << std::endl;
  }
  else
  {
    int bunchnumber = p_gl1->getBunchNumber();
    //	uint64_t evtBCO_gl1 = p_gl1->getBCO() & 0xFFFFFFFFFFU;
    //        for (int i = 0; i <9;i++)// int(GL1PScaler_raw_vec.size()); i++)
    //        {
    if (p_gl1->lValue(0, "GL1PRAW"))  // 0-8, 0 is MBDSN
    {
      //    GL1PScaler_raw_vec[i][bunchnumber] = p_gl1->lValue(i, "GL1PRAW");
      //		std::cout<<"evtBCO: "<<evtBCO_gl1<<" bunchnumber ="<<bunchnumber<<" i = "<<i<<" ,gl1praw = " <<p_gl1->lValue(i, "GL1PRAW")<<std::endl;
      m_bunchnumber_MBDNS_raw[bunchnumber] = p_gl1->lValue(0, "GL1PRAW");
      m_bunchnumber_MBDNS_live[bunchnumber] = p_gl1->lValue(0, "GL1PLIVE");
      m_bunchnumber_MBDNS_scaled[bunchnumber] = p_gl1->lValue(0, "GL1PSCALED");
      m_bunchnumber_ZDCCoin_raw[bunchnumber] = p_gl1->lValue(5, "GL1PRAW");  // zdc coincidence
      // h_gl1p_MBDSN_bunchid->Fill(bunchnumber, p_gl1->lValue(0, "GL1PRAW"));
      // std::cout<<" bunchnumber ="<<bunchnumber<<" ,gl1praw = " <<p_gl1->lValue(0, "GL1PRAW")<<std::endl;
    }
    //      }
    if (p_gl1->lValue(0, 0))
    {
      //		m_bunchnumber_rawgl1scaler[bunchnumber] = p_gl1->lValue(0, 0);
      //	std::cout<<" bunchnumber ="<<bunchnumber<<" ,gl1rawscaler = "<< p_gl1->lValue(0, 0)<<std::endl;
      m_rawgl1scaler = p_gl1->lValue(0, 0);
    }
  }
  ttree->Fill();

  if (m_lastevent_flag)
  {
    for (const auto &[bunchnumber, mbdns_value] : m_bunchnumber_MBDNS_raw)
    {
      h_gl1p_MBDSN_bunchid_raw->Fill(bunchnumber, mbdns_value);
    }
    for (const auto &[bunchnumber, mbdns_value] : m_bunchnumber_MBDNS_live)
    {
      h_gl1p_MBDSN_bunchid_live->Fill(bunchnumber, mbdns_value);
    }
    for (const auto &[bunchnumber, mbdns_value] : m_bunchnumber_MBDNS_scaled)
    {
      h_gl1p_MBDSN_bunchid_scaled->Fill(bunchnumber, mbdns_value);
    }
    // for (const auto &[bunchnumber, mbdns_value] : m_bunchnumber_rawgl1scaler) {
    h_gl1p_rawgl1scaler->Fill(1, m_rawgl1scaler);
    //}
    for (const auto &[bunchnumber, mbdns_value] : m_bunchnumber_ZDCCoin_raw)
    {
      h_gl1p_ZDCCoin_bunchid_raw->Fill(bunchnumber, mbdns_value);
    }
  }
  // if we run streaming, we only need the first gl1 bco to skip over all the junk
  // which is taken before the daq actually starts. But once we have the first event
  // and set the refBCO to the beginning of the run, we don't want the gl1 anymore
  // so we delete its input manager(s) and unregister it
  // deleting it also deletes all its allocated memory, so we don't have to worry
  // about clearing all gl1 related maps
  if (m_StreamingFlag)
  {
    for (auto iter : m_Gl1InputVector)
    {
      delete iter;
    }
    m_gl1_registered_flag = false;
    m_Gl1InputVector.clear();
  }
  else
  {
    for (auto iter : m_Gl1InputVector)
    {
      iter->CleanupUsedPackets(m_Gl1RawHitMap.begin()->first);
    }
    m_Gl1RawHitMap.begin()->second.Gl1RawHitVector.clear();
    m_Gl1RawHitMap.erase(m_Gl1RawHitMap.begin());
  }
  // std::cout << "size  m_Gl1RawHitMap: " <<  m_Gl1RawHitMap.size()
  // 	    << std::endl;

  if (Verbosity() > 0)
  {
    if (m_alldone_flag)
    {
      std::cout << "all done is true" << std::endl;
    }
  }

  if (m_alldone_flag)
  {
    std::cout << m_event_number << " Events -- Storing files to output.root" << std::endl;
    std::string updatedFileName = m_outputFileName + "_" + std::to_string(m_event_number) + ".root";
    if (TFile::Open(updatedFileName.c_str(), "READ"))
    {
      updatedFileName = m_outputFileName + "_" + std::to_string(m_event_number + 1) + ".root";
    }
    tfile = TFile::Open(updatedFileName.c_str(), "RECREATE", "");
    ttree->Write("", TObject::kOverwrite);
    h_lumibco->Write("", TObject::kOverwrite);
    h_bunchnumber->Write("", TObject::kOverwrite);
    h_bunchnumber_occur->Write("", TObject::kOverwrite);
    h_diffbco->Write("", TObject::kOverwrite);
    h_gl1p_MBDSN_bunchid_raw->Write("", TObject::kOverwrite);
    h_gl1p_MBDSN_bunchid_live->Write("", TObject::kOverwrite);
    h_gl1p_MBDSN_bunchid_scaled->Write("", TObject::kOverwrite);
    h_gl1p_rawgl1scaler->Write("", TObject::kOverwrite);
    h_gl1p_ZDCCoin_bunchid_raw->Write("", TObject::kOverwrite);
    tfile->Close();
    delete tfile;

    ttree->Reset();
    h_lumibco->Reset();
    h_bunchnumber->Reset();
    h_bunchnumber_occur->Reset();
    h_diffbco->Reset();
    h_gl1p_MBDSN_bunchid_raw->Reset();
    h_gl1p_MBDSN_bunchid_live->Reset();
    h_gl1p_MBDSN_bunchid_scaled->Reset();
    h_gl1p_rawgl1scaler->Reset();
    h_gl1p_ZDCCoin_bunchid_raw->Reset();
  }

  return 0;
}

void Fun4AllStreamingLumiCountingInputManager::SetNegativeWindow(const unsigned int i)
{
  m_negative_bco_window = std::max(i, m_negative_bco_window);
}

void Fun4AllStreamingLumiCountingInputManager::SetPositiveWindow(const unsigned int i)
{
  m_positive_bco_window = std::max(i, m_positive_bco_window);
}

void Fun4AllStreamingLumiCountingInputManager::createLuminosityHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  // zhiwan
  {
    auto tr = new TTree("BCOWindowTree", "BCO Window Data");
    tr->Branch("bco_trim", &m_bco_trim);
    tr->Branch("lower_bound", &m_lower_bound);
    tr->Branch("upper_bound", &m_upper_bound);
    tr->Branch("bunch_number", &m_bunch_number);
    //  tr->Branch("rawgl1scaler", &m_rawgl1scaler);
    tr->SetAutoFlush(100000);
    hm->registerHisto(tr);
  }

  {
    auto h = new TH1I("h_LumiBCO", "Lumi BCO", 500, 0, 500);
    h->GetXaxis()->SetTitle(" Lumi BCO per event");
    h->SetTitle("Number of BCO matched");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1I("h_BunchNumber", "Bunch Number Lumi BCO", 121, -0.5, 120.5);
    h->GetXaxis()->SetTitle("Bunch Number per event");
    h->SetTitle("Number of crossing");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D("h_BunchNumberOccurance", "Bunch Number Lumi BCO", 120, -0.5, 119.5);
    h->GetXaxis()->SetTitle("Bunch Number per time window");
    h->SetTitle("Number of crossing");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1I("h_diffBCO", "gl1 bco 1-2", 3500, 0, 3500);
    h->GetXaxis()->SetTitle("GL1 BCO difference");
    h->SetTitle("Number of crossing");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D("h_MBDSNraw_BunchID", "Bunch Number Lumi BCO", 121, -0.5, 120.5);
    h->GetXaxis()->SetTitle("Bunch Number per event");
    h->SetTitle("MBDSN Number of crossing");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D("h_MBDSNlive_BunchID", "Bunch Number Lumi BCO", 121, -0.5, 120.5);
    h->GetXaxis()->SetTitle("Bunch Number per event");
    h->SetTitle("MBDSN Number of crossing");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D("h_MBDSNscaled_BunchID", "Bunch Number Lumi BCO", 121, -0.5, 120.5);
    h->GetXaxis()->SetTitle("Bunch Number per event");
    h->SetTitle("MBDSN Number of crossing");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D("h_rawgl1scalerBunchID", "Bunch Number Lumi BCO", 10, -0.5, 9.5);
    h->GetXaxis()->SetTitle("Bunch Number per event");
    h->SetTitle("raw GL1 scaler");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D("h_gl1p_ZDCCoin_BunchID", "Bunch Number Lumi BCO", 121, -0.5, 120.5);
    h->GetXaxis()->SetTitle("Bunch Number per event");
    h->SetTitle("raw GL1 scaler");
    hm->registerHisto(h);
  }
  // Get the global pointers
  h_lumibco = dynamic_cast<TH1 *>(hm->getHisto("h_LumiBCO"));
  h_bunchnumber = dynamic_cast<TH1 *>(hm->getHisto("h_BunchNumber"));
  h_bunchnumber_occur = dynamic_cast<TH1 *>(hm->getHisto("h_BunchNumberOccurance"));
  ttree = dynamic_cast<TTree *>(hm->getHisto("BCOWindowTree"));
  h_diffbco = dynamic_cast<TH1 *>(hm->getHisto("h_diffBCO"));
  h_gl1p_MBDSN_bunchid_raw = dynamic_cast<TH1 *>(hm->getHisto("h_MBDSNraw_BunchID"));
  h_gl1p_MBDSN_bunchid_live = dynamic_cast<TH1 *>(hm->getHisto("h_MBDSNlive_BunchID"));
  h_gl1p_MBDSN_bunchid_scaled = dynamic_cast<TH1 *>(hm->getHisto("h_MBDSNscaled_BunchID"));
  h_gl1p_rawgl1scaler = dynamic_cast<TH1 *>(hm->getHisto("h_rawgl1scalerBunchID"));
  h_gl1p_ZDCCoin_bunchid_raw = dynamic_cast<TH1 *>(hm->getHisto("h_gl1p_ZDCCoin_BunchID"));
}
