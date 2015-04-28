#ifndef FUN4ALLHEPMCINPUTMANAGER_H__
#define FUN4ALLHEPMCINPUTMANAGER_H__

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

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

class Fun4AllHepMCInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllHepMCInputManager(const std::string &name = "DUMMY", const std::string &nodename = "DST", const std::string &topnodename = "TOP");
  virtual ~Fun4AllHepMCInputManager();
  int fileopen(const std::string &filenam);
  int fileclose();
  int run(const int nevents = 0);
  int isOpen() {return isopen;}
  void ReadOscar(const int i=1) {readoscar = i;}
  void MomentumUnit(const int i) {momentumunit = i;}
  void LengthUnit(const int i) {lengthunit=i;}
  void Print(const std::string &what = "ALL") const;
  int PushBackEvents(const int i);

  // Effectivly turn off the synchronization checking
  //
  int SyncIt(const SyncObject* /*mastersync*/) {return Fun4AllReturnCodes::SYNC_OK;}
  int GetSyncObject(SyncObject** /*mastersync*/) {return Fun4AllReturnCodes::SYNC_NOOBJECT;}
  int NoSyncPushBackEvents(const int nevt) {return PushBackEvents(nevt);}
  HepMC::GenEvent *ConvertFromOscar();

 protected:
  int OpenNextFile();
  int isopen;
  int events_total;
  int events_thisfile;
  int readoscar;
  int momentumunit;
  int lengthunit;
  std::string filename;
  std::string topNodeName;
  PHCompositeNode *topNode;
  HepMC::IO_GenEvent *ascii_in;
  HepMC::GenEvent *evt;
  HepMC::GenEvent *save_evt;

  // some pointers for use in decompression handling
  std::ifstream *filestream; // holds compressed filestream
  std::istream *unzipstream; // feed into HepMc
  std::ifstream theOscarFile;
};

#endif /* __FUN4ALLHEPMCINPUTMANAGER_H__ */
