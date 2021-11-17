// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PILEUPGENERATOR_H
#define G4MAIN_PHG4PILEUPGENERATOR_H

#include "PHG4ParticleGeneratorBase.h"

#include <string>  // for string

class PHCompositeNode;

class PHG4PileupGenerator : public PHG4ParticleGeneratorBase
{
 public:
  PHG4PileupGenerator(const std::string &name = "PILEUPGENERATOR");
  ~PHG4PileupGenerator() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int Reset(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode *topNode) override;

  void set_generator(PHG4ParticleGeneratorBase *generator) { _generator = generator; }

  /// past times are negative, future times are positive
  void set_time_window(double past_nsec, double future_nsec)
  {
    _min_integration_time = past_nsec;
    _max_integration_time = future_nsec;
  }
  void set_collision_rate(double kHz) { _collision_rate = kHz; }
  void set_time_between_crossings(double nsec) { _time_between_crossings = nsec; }

 private:
  PHG4ParticleGeneratorBase *_generator = nullptr;

  double _min_integration_time = -1000.;
  double _max_integration_time = 1000.;
  double _collision_rate = 100.;  // kHz
  double _time_between_crossings = 106.;

  double _ave_coll_per_crossing = 1.;  // recalculated
  int _min_crossing = 0;               // recalculated
  int _max_crossing = 0;               // recalculated
};

#endif
