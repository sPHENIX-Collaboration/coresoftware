// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLTRIGGEREDINPUTMANAGER_H
#define FUN4ALLRAW_FUN4ALLTRIGGEREDINPUTMANAGER_H

#include "InputManagerType.h"

#include <fun4all/Fun4AllInputManager.h>

#include <Event/phenixTypes.h>

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

class Event;
class SinglePrdfInput;
class CaloPacket;
class Gl1Packet;
class LL1Packet;
class PHCompositeNode;
class SingleTriggeredInput;
class SyncObject;
class OfflinePacket;

class Fun4AllTriggeredInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllTriggeredInputManager(const std::string &name = "DUMMY", const std::string &prdfnodename = "PRDF", const std::string &topnodename = "TOP");
  ~Fun4AllTriggeredInputManager() override;

  int fileopen(const std::string & /* filenam */) override { return 0; }
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
  void registerTriggeredInput(SingleTriggeredInput *prdfin);
  void registerGl1TriggeredInput(SingleTriggeredInput *prdfin);

 private:
 
  int m_RunNumber{0};
  int m_RefEventNo{std::numeric_limits<int>::min()};
  bool m_gl1_registered_flag{false};
  bool m_mbd_registered_flag{false};
  bool m_cemc_registered_flag{false};
  bool m_hcal_registered_flag{false};
  bool m_ll1_registered_flag{false};
  bool m_zdc_registered_flag{false};
  bool m_resync_flag{false};
  unsigned int m_InitialPoolDepth = 10;
  unsigned int m_DefaultPoolDepth = 10;
  unsigned int m_PoolDepth{m_InitialPoolDepth};
  std::set<int> m_Gl1DroppedEvent;
  SingleTriggeredInput *m_Gl1TriggeredInput{nullptr};
  std::vector<SingleTriggeredInput *> m_TriggeredInputVector;
  SyncObject *m_SyncObject{nullptr};
  PHCompositeNode *m_topNode{nullptr};
};

#endif
