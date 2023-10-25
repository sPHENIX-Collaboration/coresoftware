// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLSTREAMINGINPUTMANAGER_H
#define FUN4ALLRAW_FUN4ALLSTREAMINGINPUTMANAGER_H

#include <fun4all/Fun4AllInputManager.h>

#include <Event/phenixTypes.h>

#include <map>
#include <string>

class SingleStreamingInput;
class ospEvent;
class MvtxRawHit;
class InttRawHit;
class MicromegasRawHit;
class Packet;
class PHCompositeNode;
class SyncObject;
class TpcRawHit;

class Fun4AllStreamingInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllStreamingInputManager(const std::string &name = "DUMMY", const std::string &dstnodename = "DST", const std::string &topnodename = "TOP");
  ~Fun4AllStreamingInputManager() override;

  enum enu_subsystem { MVTX = 1, INTT = 2, TPC = 3, MICROMEGAS = 4};

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
  void registerStreamingInput(SingleStreamingInput *evtin, enu_subsystem);
  int FillMvtx();
  int FillIntt();
  int FillMicromegas();
  int FillTpc();
  void AddMvtxRawHit(uint64_t bclk, MvtxRawHit *hit);
  void AddInttRawHit(uint64_t bclk, InttRawHit *hit);
  void AddMicromegasRawHit(uint64_t /* bclk */, MicromegasRawHit* /* hit */);
  void AddTpcRawHit(uint64_t bclk, TpcRawHit *hit);
  void SetTpcBcoRange(const unsigned int i);

 private:
  struct MvtxRawHitInfo
  {
    std::vector<MvtxRawHit *> MvtxRawHitVector;
    unsigned int EventFoundCounter = 0;
  };

  struct InttRawHitInfo
  {
    std::vector<InttRawHit *> InttRawHitVector;
    unsigned int EventFoundCounter = 0;
  };

  struct MicromegasRawHitInfo
  {
    std::vector<MicromegasRawHit *> MicromegasRawHitVector;
    unsigned int EventFoundCounter = 0;
  };

  struct TpcRawHitInfo
  {
    std::vector<TpcRawHit *> TpcRawHitVector;
    unsigned int EventFoundCounter = 0;
  };

  int m_RunNumber = 0;
  bool m_mvtx_registered_flag = false;
  bool m_intt_registered_flag = false;
  bool m_micromegas_registered_flag = false;
  bool m_tpc_registered_flag = false;
  unsigned int m_tpc_bco_range {0};
  std::vector<SingleStreamingInput *> m_EvtInputVector;
  SyncObject *m_SyncObject = nullptr;
  PHCompositeNode *m_topNode = nullptr;
  std::map<uint64_t, MvtxRawHitInfo> m_MvtxRawHitMap;
  std::map<uint64_t, InttRawHitInfo> m_InttRawHitMap;
  std::map<uint64_t, MicromegasRawHitInfo> m_MicromegasRawHitMap;
  std::map<uint64_t, TpcRawHitInfo> m_TpcRawHitMap;
  std::map<int, std::map<int, uint64_t>> m_InttPacketFeeBcoMap;
};

#endif /* FUN4ALL_FUN4ALLSTREAMINGINPUTMANAGER_H */
