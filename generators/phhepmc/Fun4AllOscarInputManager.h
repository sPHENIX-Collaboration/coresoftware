#ifndef PHHEPMC_FUN4ALLOSCARINPUTMANAGER_H__
#define PHHEPMC_FUN4ALLOSCARINPUTMANAGER_H__
#include "PHHepMCGenHelper.h"

#include <fun4all/Fun4AllInputManager.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>

// forward declaration of classes in namespace
namespace HepMC
{
class IO_GenEvent;
class GenEvent;
};  // namespace HepMC

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
  int isOpen() { return isopen; }
  void Print(const std::string &what = "ALL") const;
  int ResetEvent();
  int PushBackEvents(const int i);
  int skip(const int i) { return PushBackEvents(i); }

  // Effectivly turn off the synchronization checking
  //
  int SyncIt(const SyncObject * /*mastersync*/) { return Fun4AllReturnCodes::SYNC_OK; }
  int GetSyncObject(SyncObject ** /*mastersync*/) { return Fun4AllReturnCodes::SYNC_NOOBJECT; }
  int NoSyncPushBackEvents(const int nevt) { return PushBackEvents(nevt); }
  int ConvertFromOscar();

  //! toss a new vertex according to a Uniform or Gaus distribution
  void set_vertex_distribution_function(PHHepMCGenHelper::VTXFUNC x, PHHepMCGenHelper::VTXFUNC y, PHHepMCGenHelper::VTXFUNC z, PHHepMCGenHelper::VTXFUNC t)
  {
    hepmc_helper.set_vertex_distribution_function(x, y, z, t);
  }

  //! set the mean value of the vertex distribution, use PHENIX units of cm, ns
  void set_vertex_distribution_mean(const double x, const double y, const double z, const double t)
  {
    hepmc_helper.set_vertex_distribution_mean(x, y, z, t);
  }

  //! set the width of the vertex distribution function about the mean, use PHENIX units of cm, ns
  void set_vertex_distribution_width(const double x, const double y, const double z, const double t)
  {
    hepmc_helper.set_vertex_distribution_width(x, y, z, t);
  }
  //
  //! reuse vertex from another PHHepMCGenEvent with embedding_id = src_embedding_id Additional smearing and shift possible with set_vertex_distribution_*()
  void set_reuse_vertex(int src_embedding_id)
  {
    hepmc_helper.set_reuse_vertex(src_embedding_id);
  }

  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int get_embedding_id() const { return hepmc_helper.get_embedding_id(); }
  //
  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void set_embedding_id(int id) { hepmc_helper.set_embedding_id(id); }

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
  std::ifstream *filestream;  // holds compressed filestream
  std::istream *unzipstream;  // feed into HepMc
  std::ifstream theOscarFile;

  //! helper for insert HepMC event to DST node and add vertex smearing
  PHHepMCGenHelper hepmc_helper;
  bool isCompressed;
};

#endif /* PHHEPMC_FUN4ALLOSCARINPUTMANAGER_H */
