#ifndef FUN4ALLPRDFINPUTMANAGER_H__
#define FUN4ALLPRDFINPUTMANAGER_H__

#include "Fun4AllInputManager.h"

#include <string>
#include <map>

class Event;
class Eventiterator;
class PHCompositeNode;
class SyncObject;

class Fun4AllPrdfInputManager : public Fun4AllInputManager
{
 public:
   Fun4AllPrdfInputManager(const std::string &name = "DUMMY", const std::string &topnodename = "TOP");
  virtual ~Fun4AllPrdfInputManager();
  int fileopen(const std::string &filenam);
  int fileclose();
  int run(const int nevents = 0);
  int isOpen() {return isopen;}

  void Print(const std::string &what = "ALL") const;
  int ResetEvent();
  int PushBackEvents(const int i);
  int GetSyncObject(SyncObject **mastersync);
  int SyncIt(const SyncObject *mastersync);

 protected:
  int OpenNextFile();
  int segment;
  int isopen;
  int events_total;
  int events_thisfile;
  PHCompositeNode *topNode;
  Event *evt;
  Event *save_evt;
  Eventiterator *eventiterator;
  SyncObject* syncobject;
};

#endif /* __FUN4ALLPRDFINPUTMANAGER_H__ */
