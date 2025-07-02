// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLTRIGGEREDINPUTMANAGER_H
#define FUN4ALLRAW_FUN4ALLTRIGGEREDINPUTMANAGER_H

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
  void EventNumber(const int i) { m_EventNumber = i; }
  int EventNumber() const { return m_EventNumber; }
  int FillPools();

 private:
  int m_RunNumber{0};
  int m_EventNumber{0};
  bool m_OnlyGl1Flag{false};
  std::set<int> m_Gl1DroppedEvent;
  SingleTriggeredInput *m_Gl1TriggeredInput{nullptr};
  std::vector<SingleTriggeredInput *> m_TriggeredInputVector;
  SyncObject *m_SyncObject{nullptr};
  PHCompositeNode *m_topNode{nullptr};
};

#endif
