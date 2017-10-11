#ifndef FUN4ALLHEPMCINPUTMANAGER_H__
#define FUN4ALLHEPMCINPUTMANAGER_H__

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#ifndef __CINT__
#include <boost/iostreams/filtering_streambuf.hpp>
#endif

#include "PHHepMCGenEvent.h"

#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

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
  //! supported function distributions
  enum VTXFUNC
  {
    Uniform,
    Gaus
  };

  Fun4AllHepMCInputManager(const std::string &name = "DUMMY", const std::string &nodename = "DST", const std::string &topnodename = "TOP");
  virtual ~Fun4AllHepMCInputManager();
  virtual int fileopen(const std::string &filenam);
  virtual int fileclose();
  virtual int run(const int nevents = 0);
  int isOpen() { return isopen; }
  void ReadOscar(const int i = 1) { readoscar = i; }
  virtual void Print(const std::string &what = "ALL") const;
  virtual int PushBackEvents(const int i);

  // Effectivly turn off the synchronization checking
  //
  int SyncIt(const SyncObject * /*mastersync*/) { return Fun4AllReturnCodes::SYNC_OK; }
  int GetSyncObject(SyncObject ** /*mastersync*/) { return Fun4AllReturnCodes::SYNC_NOOBJECT; }
  int NoSyncPushBackEvents(const int nevt) { return PushBackEvents(nevt); }
  HepMC::GenEvent *ConvertFromOscar();

  //! toss a new vertex according to a Uniform or Gaus distribution
  void set_vertex_distribution_function(VTXFUNC x, VTXFUNC y, VTXFUNC z, VTXFUNC t);

  //! set the mean value of the vertex distribution, use PHENIX units of cm, ns
  void set_vertex_distribution_mean(const double x, const double y, const double z, const double t);

  //! set the width of the vertex distribution function about the mean, use PHENIX units of cm, ns
  void set_vertex_distribution_width(const double x, const double y, const double z, const double t);

  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int get_embedding_id() const { return _embedding_id; }

  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void set_embedding_id(int id) { _embedding_id = id; }

 protected:
  int OpenNextFile();

  bool shift_vertex(PHHepMCGenEvent *event) const;
  double smear(const double position, const double width, VTXFUNC dist) const;

  int isopen;
  int events_total;
  int events_thisfile;
  int readoscar;

  std::string filename;
  std::string topNodeName;
  PHCompositeNode *topNode;

  HepMC::IO_GenEvent *ascii_in;
  HepMC::GenEvent *evt;
  HepMC::GenEvent *save_evt;

  // some pointers for use in decompression handling
  std::ifstream *filestream;  // holds compressed filestream
  std::istream *unzipstream;  // feed into HepMc
  std::ifstream theOscarFile;

#ifndef __CINT__
  boost::iostreams::filtering_streambuf<boost::iostreams::input> zinbuffer;
#endif

  VTXFUNC _vertex_func_x;
  VTXFUNC _vertex_func_y;
  VTXFUNC _vertex_func_z;
  VTXFUNC _vertex_func_t;

  double _vertex_x;
  double _vertex_y;
  double _vertex_z;
  double _vertex_t;

  double _vertex_width_x;
  double _vertex_width_y;
  double _vertex_width_z;
  double _vertex_width_t;

  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int _embedding_id;

#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif
};

#endif /* __FUN4ALLHEPMCINPUTMANAGER_H__ */
