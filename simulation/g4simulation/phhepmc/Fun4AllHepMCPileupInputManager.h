#ifndef FUN4ALLHEPMCPILEUPINPUTMANAGER_H__
#define FUN4ALLHEPMCPILEUPINPUTMANAGER_H__

#include "Fun4AllHepMCInputManager.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <string>
#include <map>
#include <fstream>
#include <iostream>

// rootcint barfs with this header so we need to hide it
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

class Fun4AllHepMCPileupInputManager : public Fun4AllHepMCInputManager
{
 public:
  Fun4AllHepMCPileupInputManager(const std::string &name = "DUMMY",
                                 const std::string &nodename = "DST",
                                 const std::string &topnodename = "TOP");
  virtual ~Fun4AllHepMCPileupInputManager();

  int run(const int nevents = 0);

  /// past times are negative, future times are positive
  void set_time_window(double past_nsec,double future_nsec) {
    _min_integration_time = past_nsec;
    _max_integration_time = future_nsec;    
  }
  void set_collision_rate(double kHz) {_collision_rate = kHz;}
  void set_time_between_crossings(double nsec) {_time_between_crossings = nsec;}
  
 private:

  double _min_integration_time;
  double _max_integration_time;
  double _collision_rate;
  double _time_between_crossings;

  double   _ave_coll_per_crossing;
  int      _min_crossing;
  int      _max_crossing;
  
  unsigned int seed;
#ifndef __CINT__
//  gsl_rng *RandomGenerator;
#endif
};

#endif /* __FUN4ALLHEPMCINPUTMANAGER_H__ */
