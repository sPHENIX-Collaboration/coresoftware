#ifndef PHHEPMC_FUN4ALLOSCARINPUTMANAGER_H
#define PHHEPMC_FUN4ALLOSCARINPUTMANAGER_H

#include "PHHepMCGenHelper.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>  // for SYNC_NOOBJECT, SYNC_OK

#include <fstream>
#include <string>

// forward declaration of classes in namespace
namespace HepMC
{
  class GenEvent;
}  // namespace HepMC

class PHCompositeNode;
class SyncObject;

class Fun4AllOscarInputManager : public Fun4AllInputManager, public PHHepMCGenHelper
{
 public:
  Fun4AllOscarInputManager(const std::string &name = "DUMMY", const std::string &topnodename = "TOP");
  ~Fun4AllOscarInputManager() override;
  int fileopen(const std::string &filenam) override;
  int fileclose() override;
  int run(const int nevents = 0) override;
  void Print(const std::string &what = "ALL") const override;
  int ResetEvent() override;
  int PushBackEvents(const int i) override;
  int skip(const int i) override { return PushBackEvents(i); }

  // Effectivly turn off the synchronization checking
  //
  int SyncIt(const SyncObject * /*mastersync*/) override { return Fun4AllReturnCodes::SYNC_OK; }
  int GetSyncObject(SyncObject ** /*mastersync*/) override { return Fun4AllReturnCodes::SYNC_NOOBJECT; }
  int NoSyncPushBackEvents(const int nevt) override { return PushBackEvents(nevt); }
  int ConvertFromOscar();

 protected:
  int events_total;
  int events_thisfile;
  std::string filename;
  std::string topNodeName;
  PHCompositeNode *topNode;
  HepMC::GenEvent *evt;
  //HepMC::GenEvent *tmpEvt;
  int skipEvents, skippedEvents;

  // some pointers for use in decompression handling
  std::ifstream *filestream;  // holds compressed filestream
  std::istream *unzipstream;  // feed into HepMc
  std::ifstream theOscarFile;

  bool isCompressed;
};

#endif /* PHHEPMC_FUN4ALLOSCARINPUTMANAGER_H */
