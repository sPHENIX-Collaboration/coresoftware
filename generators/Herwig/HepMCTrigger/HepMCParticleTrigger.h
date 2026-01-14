// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HEPMCPARTICLETRIGGER_H
#define HEPMCPARTICLETRIGGER_H

#include <fun4all/SubsysReco.h>

#include <fastjet/PseudoJet.hh>

#include <cmath>
#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
namespace HepMC
{
  class GenEvent;
}

class HepMCParticleTrigger : public SubsysReco
{
 public:
  HepMCParticleTrigger(float trigger_thresh = 10., int n_incom = 1000, bool up_lim = false, const std::string& name = "HepMCParticleTrigger");

  ~HepMCParticleTrigger() override = default;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode* topNode) override;

  /// Clean up internals after each event.

  /// Called at the end of each run.

  /// Called at the end of all processing.

  /// Reset
  void AddParticles(const std::vector<int>&); //exclusively take input in the form of a pdg_ids (22 for photon, primary use case)
  void AddParticle(int);

  /*  void AddParents(const std::string &parents);
    void AddParents(int parent);
    void AddParents(std::vector<int> parents);
    void AddParentspID(std::vector<int> parents);
  */
  void SetPtHigh(double);
  void SetPtLow(double);
  void SetPtHighLow(double, double);

  void SetPHigh(double);
  void SetPLow(double);
  void SetPHighLow(double, double);

  void SetEtaHigh(double);
  void SetEtaLow(double);
  void SetEtaHighLow(double, double);

  void SetAbsEtaHigh(double);
  void SetAbsEtaLow(double);
  void SetAbsEtaHighLow(double, double);

  void SetPzHigh(double);
  void SetPzLow(double);
  void SetPzHighLow(double, double);

  void SetStableParticleOnly(bool b) { m_doStableParticleOnly = b; }

 private:
  bool isGoodEvent(HepMC::GenEvent* e1);
  std::vector<int> getParticles(HepMC::GenEvent* e1);
  int particleAboveThreshold(const std::map<int, int>& n_particles, int particle);
  //  std::vector<int> _theParentsi {};
  std::vector<int> _theParticles{};
  bool m_doStableParticleOnly{true};
  float threshold{0.};
  int goal_event_number{1000};
  int n_evts{0};
  /**
 * Number of events that passed the configured trigger selection.
 *
 * Tracks how many processed events met the selection criteria during the run.
 */
int n_good{0};
  bool set_event_limit{false};

  float _theEtaHigh{1.1};
  float _theEtaLow{-1.1};
  float _thePtHigh{999.9};
  float _thePtLow{-999.9};
  float _thePHigh{999.9};
  float _thePLow{-999.9};
  float _thePzHigh{999.9};
  float _thePzLow{-999.9};

  bool _doEtaHighCut{true};
  bool _doEtaLowCut{true};
  bool _doBothEtaCut{true};

  bool _doAbsEtaHighCut{false};
  bool _doAbsEtaLowCut{false};
  bool _doBothAbsEtaCut{false};

  bool _doPtHighCut{false};
  bool _doPtLowCut{false};
  bool _doBothPtCut{false};

  bool _doPHighCut{false};
  bool _doPLowCut{false};
  bool _doBothPCut{false};

  bool _doPzHighCut{false};
  bool _doPzLowCut{false};
  bool _doBothPzCut{false};
};

#endif  // HEPMCPARTICLETRIGGER_H