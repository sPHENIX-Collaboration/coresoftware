// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLSTREAMINGINPUTMANAGER_H
#define FUN4ALLRAW_FUN4ALLSTREAMINGINPUTMANAGER_H

#include "InputManagerType.h"

#include <fun4all/Fun4AllInputManager.h>

#include <map>
#include <set>
#include <string>
#include <vector>
#include <cinttypes>

class SingleStreamingInput;
class Gl1Packet;
class InttRawHit;
class MicromegasRawHit;
class MvtxRawHit;
class MvtxFeeIdInfo;
class PHCompositeNode;
class SyncObject;
class TpcRawHit;
class TH1;
class TH2;

class Fun4AllStreamingInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllStreamingInputManager(const std::string &name = "DUMMY", const std::string &dstnodename = "DST", const std::string &topnodename = "TOP");
  ~Fun4AllStreamingInputManager() override;
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
  void registerStreamingInput(SingleStreamingInput *evtin, InputManagerType::enu_subsystem);
  int FillGl1();
  int FillIntt();
  int FillMicromegas();
  int FillMvtx();
  int FillTpc();
  void AddGl1RawHit(uint64_t bclk, Gl1Packet *hit);
  void AddInttRawHit(uint64_t bclk, InttRawHit *hit);
  void AddMicromegasRawHit(uint64_t /* bclk */, MicromegasRawHit * /* hit */);
  void AddMvtxFeeIdInfo(uint64_t bclk, uint16_t feeid, uint32_t detField);
  void AddMvtxL1TrgBco(uint64_t bclk, uint64_t lv1Bco);
  void AddMvtxRawHit(uint64_t bclk, MvtxRawHit *hit);
  void AddTpcRawHit(uint64_t bclk, TpcRawHit *hit);
  void SetInttBcoRange(const unsigned int i);
  void SetInttNegativeBco(const unsigned int i);
  void SetMicromegasBcoRange(const unsigned int i);
  void SetMicromegasNegativeBco(const unsigned int i);
  void SetMvtxBcoRange(const unsigned int i);
  void SetMvtxNegativeBco(const unsigned int i);
  void SetTpcBcoRange(const unsigned int i);
  void SetTpcNegativeBco(const unsigned int i);
  int FillInttPool();
  int FillMicromegasPool();
  int FillMvtxPool();
  int FillTpcPool();
  void Streaming(bool b = true) { m_StreamingFlag = b; }

  void runMvtxTriggered(bool b = true) { m_mvtx_is_triggered = b; }

  // configuration for INTT hit carry-over issue mitigation (hit duplication)
  void EnableInttHitDuplication(bool b = true) { m_InttHitDuplication = b; }
  void SetIsRejectInttNoiseCrossings(bool b = true) {m_IsRejectInttNoiseCrossings = b;}
  void SetIsInttStreaming(bool b = true) {m_IsInttStreaming = b;}
  void SetIsDuplicateInttFPHXBCOResetHit(bool b = true) {m_IsDuplicateInttFPHXBCOResetHit = b;}
  void SetInttResetFphxBcoVec(std::vector<int> input_vec) {m_InttResetFphxBcoVec = input_vec;}
  
  void SetInttStreamingSignalCrossing(std::pair<int,int> input_pair);
  void InttHitCarryOverShiftMaxMultiple(const int i);


 private:
  struct MvtxRawHitInfo
  {
    std::set<uint64_t> MvtxL1TrgBco;
    std::vector<MvtxFeeIdInfo *> MvtxFeeIdInfoVector;
    std::vector<MvtxRawHit *> MvtxRawHitVector;
    unsigned int EventFoundCounter{0};
  };

  struct Gl1RawHitInfo
  {
    std::vector<Gl1Packet *> Gl1RawHitVector;
    unsigned int EventFoundCounter{0};
  };

  struct InttRawHitInfo
  {
    std::vector<InttRawHit *> InttRawHitVector;
    unsigned int EventFoundCounter{0};
  };

  struct MicromegasRawHitInfo
  {
    std::vector<MicromegasRawHit *> MicromegasRawHitVector;
    unsigned int EventFoundCounter{0};
  };

  struct TpcRawHitInfo
  {
    std::vector<TpcRawHit *> TpcRawHitVector;
    unsigned int EventFoundCounter{0};
  };

  void createQAHistos();

  SyncObject *m_SyncObject{nullptr};
  PHCompositeNode *m_topNode{nullptr};

  uint64_t m_RefBCO{0};

  int m_RunNumber{0};
  unsigned int m_intt_bco_range{0};
  unsigned int m_intt_negative_bco{0};
  unsigned int m_micromegas_bco_range{0};
  unsigned int m_micromegas_negative_bco{0};
  unsigned int m_mvtx_bco_range{0};
  unsigned int m_mvtx_negative_bco{0};
  unsigned int m_tpc_bco_range{0};
  unsigned int m_tpc_negative_bco{0};

  bool m_gl1_registered_flag{false};
  bool m_intt_registered_flag{false};
  bool m_micromegas_registered_flag{false};
  bool m_mvtx_registered_flag{false};
  bool m_StreamingFlag{false};
  bool m_tpc_registered_flag{false};
  bool m_mvtx_is_triggered{false};

  std::vector<SingleStreamingInput *> m_Gl1InputVector;
  std::vector<SingleStreamingInput *> m_InttInputVector;
  std::vector<SingleStreamingInput *> m_MicromegasInputVector;
  std::vector<SingleStreamingInput *> m_MvtxInputVector;
  std::vector<SingleStreamingInput *> m_TpcInputVector;
  std::map<uint64_t, Gl1RawHitInfo> m_Gl1RawHitMap;
  std::map<uint64_t, InttRawHitInfo> m_InttRawHitMap;
  std::map<uint64_t, MicromegasRawHitInfo> m_MicromegasRawHitMap;
  std::map<uint64_t, MvtxRawHitInfo> m_MvtxRawHitMap;
  std::map<uint64_t, TpcRawHitInfo> m_TpcRawHitMap;
  std::map<int, std::map<int, uint64_t>> m_InttPacketFeeBcoMap;

  // QA histos
  TH1 *h_refbco_mvtx[12]{nullptr};
  TH1 *h_taggedAllFelixes_mvtx{nullptr};
  TH1 *h_taggedAllFelixesAllFees_mvtx{nullptr};
  TH1 *h_tagBcoFelix_mvtx[12]{nullptr};

  TH1 *h_bcoGL1LL1diff[12]{nullptr};
  TH1 *h_bcoLL1Strobediff[12]{nullptr};
  TH1 *h_tagStBcoFelix_mvtx[12]{nullptr};
  TH1 *h_tagBcoFelixAllFees_mvtx[12]{nullptr};
  TH1 *h_tagStBcoFEE_mvtx[12]{nullptr};
  TH1 *h_tagL1BcoFEE_mvtx[12]{nullptr};

  TH1 *h_refbco_intt[8]{nullptr};
  TH1 *h_taggedAll_intt{nullptr};
  TH1 *h_taggedAllFee_intt{nullptr};
  TH1 *h_gl1tagged_intt[8]{nullptr};
  TH1 *h_taggedAllFees_intt[8]{nullptr};
  TH1 *h_gl1taggedfee_intt[8][14]{{nullptr}};
  TH2 *h_bcodiff_intt[8]{nullptr};

  // for INTT hit carry-over issue mitigation (hit duplication)
  bool m_InttHitDuplication{false}; // default to false; Should be set to true when running streaming data
  const unsigned int m_InttHitCarryOverShift{120}; // 120 BCOs as the default shift. Fixed value 
  int m_InttHitCarryOverShiftMaxMultiple{4}; // the max multiple of the shift. For a max multiple of M, duplicate hits from N + [1..M] * shift BCOs to N
  bool m_IsRejectInttNoiseCrossings{false}; // not allow hits in the abort-gap crossings being saved to the INTTRawHit container
  bool m_IsInttStreaming{true}; // is INTT in the streaming readout mode
  bool m_IsDuplicateInttFPHXBCOResetHit{true}; // Allow duplicating hits with FPHXBCO in the range given by std::vector<int>m_InttResetFphxBcoVec

  std::pair<int,int> m_InttStreamingSignalCrossing{6,116};
  std::map<std::string, int> m_InttRawHitCount_FEE;
  int m_Intt_print_count{0};
  std::vector<int> m_InttResetFphxBcoVec{0,1,2,3,4,5};

};

#endif /* FUN4ALL_FUN4ALLSTREAMINGINPUTMANAGER_H */
