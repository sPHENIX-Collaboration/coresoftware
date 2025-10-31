// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef RAWBCOLUMI_FUN4ALLSTREAMINGLUMICOUNTINGINPUTMANAGER_H
#define RAWBCOLUMI_FUN4ALLSTREAMINGLUMICOUNTINGINPUTMANAGER_H

#include <fun4allraw/InputManagerType.h>
// #include <fun4allraw/Fun4AllStreamingInputManager.h>
#include <fun4all/Fun4AllInputManager.h>

#include <TTree.h>
#include <map>
#include <set>
#include <string>
class SingleStreamingInputv2;
class Gl1Packet;
class PHCompositeNode;
class SyncObject;
class TH1;
class TTree;
class Fun4AllStreamingLumiCountingInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllStreamingLumiCountingInputManager(const std::string &name = "DUMMY", const std::string &dstnodename = "DST", const std::string &topnodename = "TOP");
  ~Fun4AllStreamingLumiCountingInputManager() override;
  int fileopen(const std::string & /*filenam*/) override { return 0; }
  // cppcheck-suppress virtualCallInConstructor
  int fileclose() override;
  int run(const int nevents = 0) override;

  void Print(const std::string &what = "ALL") const override;
  int ResetEvent() override;
  int PushBackEvents(const int i) override;
  int GetSyncObject(SyncObject **mastersync) override;
  int SyncIt(const SyncObject *mastersync) override;
  int HasSyncObject() const override { return 1; }
  std::string GetString(const std::string &what) const override;
  void registerStreamingInput(SingleStreamingInputv2 *evtin, InputManagerType::enu_subsystem);
  int FillGl1();
  void AddGl1RawHit(uint64_t bclk, Gl1Packet *hit);
  void AddGl1Window(uint64_t bclk, int negative_window, int positive_window);
  void AddGl1BunchNumber(uint64_t bclk, int bunch_number);
  void SetNegativeWindow(const unsigned int i);
  void SetPositiveWindow(const unsigned int i);
  void Streaming(bool b = true) { m_StreamingFlag = b; }
  void SetOutputFileName(const std::string &fileName);
  void SetEndofEvent(bool flag = false, bool flag2 = false)
  {
    m_alldone_flag = flag;
    m_lastevent_flag = flag2;
  }
  void SetEventNumber(int num) { m_event_number = num; }

 private:
  struct Gl1RawHitInfo
  {
    std::vector<Gl1Packet *> Gl1RawHitVector;
    unsigned int EventFoundCounter{0};
  };

  void createLuminosityHistos();

  SyncObject *m_SyncObject{nullptr};
  PHCompositeNode *m_topNode{nullptr};

  int m_RunNumber{0};
  unsigned int m_negative_bco_window{0};
  unsigned int m_positive_bco_window{0};
  uint64_t m_rawgl1scaler{0};
  //  std::string m_output_file="output.root";
  bool m_alldone_flag = {false};
  bool m_lastevent_flag = {false};
  int m_event_number{0};
  int m_diffBCO{0};
  bool m_gl1_registered_flag{false};
  bool m_StreamingFlag{false};
  bool flat_overflow{false};
  uint64_t bco_temp = 0;

  std::vector<SingleStreamingInputv2 *> m_Gl1InputVector;
  std::map<uint64_t, Gl1RawHitInfo> m_Gl1RawHitMap;
  std::map<uint64_t, std::pair<uint64_t, uint64_t>> m_BCOWindows;
  std::map<uint64_t, int> m_BCOBunchNumber;
  std::map<int, long> m_bunchnumber_MBDNS_raw;
  std::map<int, long> m_bunchnumber_MBDNS_live;
  std::map<int, long> m_bunchnumber_MBDNS_scaled;
  std::map<int, long> m_bunchnumber_ZDCCoin_raw;
  // std::map<int, long> m_bunchnumber_rawgl1scaler;

  // QA histos
  TH1 *h_lumibco{nullptr};
  TH1 *h_bunchnumber{nullptr};
  TH1 *h_bunchnumber_occur{nullptr};
  TH1 *h_diffbco{nullptr};
  TH1 *h_gl1p_MBDSN_bunchid_raw{nullptr};
  TH1 *h_gl1p_MBDSN_bunchid_live{nullptr};
  TH1 *h_gl1p_MBDSN_bunchid_scaled{nullptr};
  TH1 *h_gl1p_rawgl1scaler{nullptr};
  TH1 *h_gl1p_ZDCCoin_bunchid_raw{nullptr};
  uint64_t m_bco_trim{};
  uint64_t m_lower_bound{};
  uint64_t m_upper_bound{};
  int m_bunch_number{};
  TTree *ttree = nullptr;
  TFile *tfile = nullptr;
  std::string m_outputFileName = "/sphenix/user/xuzhiwan/luminosity/streaming-macro/macro/output.root";  // Default value
};

#endif /* RAWBCOLUMI_FUN4ALLSTREAMINGLUMICOUNTINGINPUTMANAGER_H */
