// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HEPMCJETTRIGGER_H
#define HEPMCJETTRIGGER_H
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <string>
#include <vector>
#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHObject.h>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <HepMC/GenEvent.h>


class PHCompositeNode;

class HepMCJetTrigger : public SubsysReco
{
 public:

  HepMCJetTrigger(float trigger_thresh=10., int n_incom=1000, bool up_lim=false, const std::string &name = "HepMCJetTrigger");

  ~HepMCJetTrigger() override;

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
  int n_evts=0;
  int n_good=0;
 private:
	bool isGoodEvent(HepMC::GenEvent* e1);
	std::vector<fastjet::PseudoJet> findAllJets(HepMC::GenEvent* e1);
	int jetsAboveThreshold(std::vector<fastjet::PseudoJet> jets);
	float threshold=0.;
	int goal_event_number=1000;
	bool set_event_limit=false;
};

#endif // HEPMCJETTRIGGER_H
