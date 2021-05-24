#ifndef PHHEPMC_FUN4ALLHEPMCPILEUPINPUTMANAGER_H
#define PHHEPMC_FUN4ALLHEPMCPILEUPINPUTMANAGER_H

#include "Fun4AllHepMCInputManager.h"

#include <gsl/gsl_rng.h>

#include <map>
#include <string>

//! Generate pile up collisions based on beam parameter
//! If set_embedding_id(i) with a negative number or 0, the pile up event will be inserted with increasing positive embedding_id. This is the default operation mode.
//! If set_embedding_id(i) with a positive number, the pile up event will be inserted with increasing positive embedding_id. This would be a strange way to use pile up.
class Fun4AllHepMCPileupInputManager : public Fun4AllHepMCInputManager
{
 public:
  Fun4AllHepMCPileupInputManager(const std::string &name = "DUMMY",
                                 const std::string &nodename = "DST",
                                 const std::string &topnodename = "TOP");
  ~Fun4AllHepMCPileupInputManager() override;

  int run(const int nevents = 0) override { return run(nevents, false); }

  int run(const int nevents, const bool skip);
  int ResetEvent() override;
  /// past times are negative, future times are positive
  void set_time_window(double past_nsec, double future_nsec)
  {
    _min_integration_time = past_nsec;
    _max_integration_time = future_nsec;
  }

  /// collision rate in Hz
  void set_collision_rate(double Hz) { _collision_rate = Hz; }
  /// time between bunch crossing in ns
  void set_time_between_crossings(double nsec) { _time_between_crossings = nsec; }

  int SkipForThisManager(const int nevents) override;
  void SignalInputManager(Fun4AllHepMCInputManager *in) { m_SignalInputManager = in; }
  int PushBackEvents(const int i) override;

 private:
  int InsertEvent(HepMC::GenEvent *evt, const double crossing_time);

  Fun4AllHepMCInputManager *m_SignalInputManager = nullptr;
  gsl_rng *RandomGenerator = nullptr;

  int m_SignalEventNumber = 0;
  /// past times are negative, future times are positive
  double _min_integration_time = -17500.0;
  double _max_integration_time = 17500.0;
  /// collision rate in Hz
  double _collision_rate = 100.0e3;
  /// time between bunch crossing in ns
  double _time_between_crossings = 106.0;

  //derived parameters
  double _ave_coll_per_crossing = 1.;
  int _min_crossing = 0;
  int _max_crossing = 0;

  bool _first_run = true;

  std::map<int, double> m_EventNumberMap;
};

#endif /* PHHEPMC_FUN4ALLHEPMCINPUTMANAGER_H */
