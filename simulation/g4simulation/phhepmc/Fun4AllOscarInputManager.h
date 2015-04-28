#ifndef FUN4ALLOSCARINPUTMANAGER_H__
#define FUN4ALLOSCARINPUTMANAGER_H__

#include <fun4all/Fun4AllInputManager.h>

#include <string>
#include <map>
#include <fstream>
#include <iostream>

// forward declaration of classes in namespace
namespace HepMC
{
    class IO_GenEvent;
    class GenEvent;
};

class PHCompositeNode;
class PHHepMCGenEvent;


class Fun4AllOscarInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllOscarInputManager(const std::string &name = "DUMMY", const std::string &topnodename = "TOP");
  virtual ~Fun4AllOscarInputManager();
  int fileopen(const std::string &filenam);
  int fileclose();
  int run(const int nevents = 0);
  int isOpen() {return isopen;}
  void Print(const std::string &what = "ALL") const;
  int ResetEvent();
  int PushBackEvents(const int i);
  int skip(const int i){ return PushBackEvents(i);}

  // Effectivly turn off the synchronization checking
  //
  int SyncIt(const SyncObject* /*mastersync*/) {return Fun4AllReturnCodes::SYNC_OK;}
  int GetSyncObject(SyncObject** /*mastersync*/) {return Fun4AllReturnCodes::SYNC_NOOBJECT;}
  int NoSyncPushBackEvents(const int nevt) {return PushBackEvents(nevt);}
  int ConvertFromOscar();


 protected:
  int OpenNextFile();
  int isopen;
  int events_total;
  int events_thisfile;
  std::string filename;
  std::string topNodeName;
  PHCompositeNode *topNode;
  HepMC::GenEvent *evt;
  //HepMC::GenEvent *tmpEvt;
  int skipEvents, skippedEvents;

  // some pointers for use in decompression handling
  std::ifstream *filestream; // holds compressed filestream
  std::istream *unzipstream; // feed into HepMc
  std::ifstream theOscarFile;

  PHHepMCGenEvent *phhepmcgenevt;
  bool isCompressed;
};

#endif /* __FUN4ALLOSCARINPUTMANAGER_H__ */
