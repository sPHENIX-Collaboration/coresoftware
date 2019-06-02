// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLPRDFINPUTMANAGER_H
#define FUN4ALL_FUN4ALLPRDFINPUTMANAGER_H

#include "Fun4AllInputManager.h"

#include <string>

class Event;
class Eventiterator;
class PHCompositeNode;
class SyncObject;

class Fun4AllPrdfInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllPrdfInputManager(const std::string &name = "DUMMY", const std::string &prdfnodename = "PRDF", const std::string &topnodename = "TOP");
  virtual ~Fun4AllPrdfInputManager();
  int fileopen(const std::string &filenam);
  int fileclose();
  int run(const int nevents = 0);

  void Print(const std::string &what = "ALL") const;
  int ResetEvent();
  int PushBackEvents(const int i);
  int GetSyncObject(SyncObject **mastersync);
  int SyncIt(const SyncObject *mastersync);

 protected:
  int segment;
  int events_total;
  int events_thisfile;
  PHCompositeNode *topNode;
  Event *evt;
  Event *save_evt;
  Eventiterator *eventiterator;
  SyncObject *syncobject;
  std::string m_PrdfNodeName;
};

#endif /* FUN4ALL_FUN4ALLPRDFINPUTMANAGER_H */
