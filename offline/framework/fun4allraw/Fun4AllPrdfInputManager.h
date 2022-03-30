// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLPRDFINPUTMANAGER_H
#define FUN4ALLRAW_FUN4ALLPRDFINPUTMANAGER_H

#include <fun4all/Fun4AllInputManager.h>

#include <string>

class Event;
class Eventiterator;
class PHCompositeNode;
class SyncObject;

class Fun4AllPrdfInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllPrdfInputManager(const std::string &name = "DUMMY", const std::string &prdfnodename = "PRDF", const std::string &topnodename = "TOP");
  ~Fun4AllPrdfInputManager() override;
  int fileopen(const std::string &filenam) override;
  int fileclose() override;
  int run(const int nevents = 0) override;

  void Print(const std::string &what = "ALL") const override;
  int ResetEvent() override;
  int PushBackEvents(const int i) override;
  int GetSyncObject(SyncObject **mastersync) override;
  int SyncIt(const SyncObject *mastersync) override;
  int HasSyncObject() const  override {return 1;}
  std::string GetString(const std::string &what) const override;

 private:
  int m_Segment = -999;
  int m_EventsTotal = 0;
  int m_EventsThisFile = 0;
  PHCompositeNode *m_topNode = nullptr;
  Event *m_Event = nullptr;
  Event *m_SaveEvent = nullptr;
  Eventiterator *m_EventIterator = nullptr;
  SyncObject *m_SyncObject = nullptr;
  std::string m_PrdfNodeName;
};

#endif /* FUN4ALL_FUN4ALLPRDFINPUTMANAGER_H */
