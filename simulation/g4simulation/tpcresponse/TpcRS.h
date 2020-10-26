// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCRS_H
#define TPCRS_H

#include <fun4all/SubsysReco.h>

//#if defined(__CLING__)
//namespace tpcrs {
//template< typename Base_t, typename Chair_t, typename Struct_t >
//  std::string ConfigStruct< Base_t, Chair_t, Struct_t >::name{};
//}
//#endif
#include <tpcrs/tpcrs.h>

#include <string>
#include <vector>

class PHCompositeNode;

class TpcRS : public SubsysReco
{
 public:
  TpcRS(const std::string &name = "TpcRS");

  virtual ~TpcRS();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  void SetupConfigurator(const std::string &filename);

 private:
  tpcrs::Configurator *cfg = nullptr;
  tpcrs::Simulator *simulator = nullptr;
  std::vector<tpcrs::SimulatedHit> simu_hits;
  std::vector<tpcrs::DistortedHit> dist_hits;
  std::vector<tpcrs::DigiHit> digi_hits;
};

#endif  // TPCRS_H
