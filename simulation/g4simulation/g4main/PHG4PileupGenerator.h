// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PILEUPGENERATOR_H
#define G4MAIN_PHG4PILEUPGENERATOR_H

#include "PHG4ParticleGeneratorBase.h"

#include <string>                       // for string

class PHCompositeNode;

class PHG4PileupGenerator : public PHG4ParticleGeneratorBase {

public:

  PHG4PileupGenerator(const std::string &name="PILEUPGENERATOR");
  virtual ~PHG4PileupGenerator();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  int EndRun(const int runnumber);
  int End(PHCompositeNode *topNode);

  void set_generator(PHG4ParticleGeneratorBase* generator) {_generator = generator;}

  /// past times are negative, future times are positive
  void set_time_window(double past_nsec,double future_nsec) {
    _min_integration_time = past_nsec;
    _max_integration_time = future_nsec;    
  }
  void set_collision_rate(double kHz) {_collision_rate = kHz;}
  void set_time_between_crossings(double nsec) {_time_between_crossings = nsec;}

private:

  PHG4ParticleGeneratorBase* _generator;

  double _min_integration_time;
  double _max_integration_time;
  double _collision_rate;
  double _time_between_crossings;

  double   _ave_coll_per_crossing;
  int      _min_crossing;
  int      _max_crossing;
};

#endif
