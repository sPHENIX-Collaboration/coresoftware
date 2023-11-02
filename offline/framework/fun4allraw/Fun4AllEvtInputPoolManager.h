// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLEVTINPUTPOOLMANAGER_H
#define FUN4ALLRAW_FUN4ALLEVTINPUTPOOLMANAGER_H

#include <fun4all/Fun4AllInputManager.h>

#include <Event/phenixTypes.h>

#include <map>
#include <set>
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

class Fun4AllEvtInputPoolManager : public Fun4AllInputManager
{
 public:
  Fun4AllEvtInputPoolManager(const std::string &name = "DUMMY", const std::string &dstnodename = "DST", const std::string &topnodename = "TOP");
  ~Fun4AllEvtInputPoolManager() override;

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
//  SingleEvtInput *AddEvtInputList(const std::string &listfile);
//  SingleEvtInput *AddEvtInputFile(const std::string &filename);
  void registerStreamingInput(SingleStreamingInput *evtin, enu_subsystem);
  void AddPacket(uint64_t bclk, Packet *p);
  void UpdateEventFoundCounter(const int evtno);
  int FillMvtx();
  int FillIntt();
  int FillMicromegas();
  int FillTpc();
  void AddMvtxFeeId(uint64_t bclk, uint16_t feeid);
  void AddMvtxL1TrgBco(uint64_t bclk, uint64_t lv1Bco);
  void AddMvtxRawHit(uint64_t bclk, MvtxRawHit *hit);
  void AddInttRawHit(uint64_t bclk, InttRawHit *hit);
  void AddMicromegasRawHit(uint64_t /* bclk */, MicromegasRawHit* /* hit */);
  void AddTpcRawHit(uint64_t bclk, TpcRawHit *hit);

 private:
  struct PacketInfo
  {
    std::vector<Packet *> PacketVector;
    unsigned int EventFoundCounter = 0;
  };

  struct MvtxRawHitInfo
  {
    std::set<uint16_t> MvtxFeeIds;
    std::set<uint64_t> MvtxL1TrgBco;
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
  std::vector<SingleStreamingInput *> m_EvtInputVector;
  std::vector<int> m_MvtxEvtInputList;
  SyncObject *m_SyncObject = nullptr;
  PHCompositeNode *m_topNode = nullptr;
  std::map<uint64_t, PacketInfo> m_PacketInfoMap;
  std::map<uint64_t, MvtxRawHitInfo> m_MvtxRawHitMap;
  std::map<uint64_t, InttRawHitInfo> m_InttRawHitMap;
  std::map<uint64_t, MicromegasRawHitInfo> m_MicromegasRawHitMap;
  std::map<uint64_t, TpcRawHitInfo> m_TpcRawHitMap;
};

#endif /* FUN4ALL_FUN4ALLEVTINPUTPOOLMANAGER_H */
