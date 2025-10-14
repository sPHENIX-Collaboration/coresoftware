// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HEPMCJETTRIGGER_H
#define HEPMCJETTRIGGER_H

#include <fun4all/SubsysReco.h>

#include <fastjet/PseudoJet.hh>

#include <string>
#include <vector>

class PHCompositeNode;
namespace HepMC
{
  class GenEvent;
}

class HepMCJetTrigger : public SubsysReco
{
 public:
  HepMCJetTrigger(float trigger_thresh = 10., int n_incom = 1000, bool up_lim = false, const std::string& name = "HepMCJetTrigger");

  ~HepMCJetTrigger() override = default;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode* topNode) override;

  /// Clean up internals after each event.

  /// Called at the end of each run.

  /// Called at the end of all processing.

  /// Reset

 private:
  bool isGoodEvent(HepMC::GenEvent* e1);
  std::vector<fastjet::PseudoJet> findAllJets(HepMC::GenEvent* e1);
  int jetsAboveThreshold(const std::vector<fastjet::PseudoJet>& jets) const;
  float threshold{0.};
  int goal_event_number{1000};
  int n_evts{0};
  int n_good{0};
  bool set_event_limit{false};
};

#endif  // HEPMCJETTRIGGER_H
