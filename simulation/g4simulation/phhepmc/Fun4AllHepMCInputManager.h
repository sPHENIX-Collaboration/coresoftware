#ifndef FUN4ALLHEPMCINPUTMANAGER_H__
#define FUN4ALLHEPMCINPUTMANAGER_H__

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <string>
#include <map>
#include <fstream>
#include <iostream>

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
  enum VTXFUNC {Uniform,Gaus};

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

  //! toss a new vertex according to a Uniform or Gaus distribution
  void set_vertex_distribution_function(VTXFUNC x, VTXFUNC y, VTXFUNC z);

  //! set the mean value of the vertex distribution
  void set_vertex_distribution_mean(const double x, const double y, const double z);

  //! set the width of the vertex distribution function about the mean
  void set_vertex_distribution_width(const double x, const double y, const double z);

 protected:
  int OpenNextFile();

  bool shift_vertex(PHHepMCGenEvent* event) const;
  double smear(const double position, const double width, VTXFUNC dist) const;

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

  VTXFUNC _vertex_func_x;
  VTXFUNC _vertex_func_y;
  VTXFUNC _vertex_func_z;
  double _vertex_x;
  double _vertex_y;
  double _vertex_z;
  double _vertex_width_x;
  double _vertex_width_y;
  double _vertex_width_z;

#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif
};

#endif /* __FUN4ALLHEPMCINPUTMANAGER_H__ */
