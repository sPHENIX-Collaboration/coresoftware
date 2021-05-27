#ifndef CALORECO_RAWTOWERCOMBINER_H
#define CALORECO_RAWTOWERCOMBINER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawTowerContainer;

//! \brief RawTowerCombiner module that joints multiple RawTower together to form a single readout in a separate node
//! Use this class to simulate ganged readout or trigger channels that combined multiple towers
//! Should be called after RawTowerBuilder but before RawTowerDigitizer
/*
## Introduction

A new tower analysis module to enable analysis of sPHENIX simulation to use multi-tower ganging of EMCal readout. For example 2x2-ganging as proposed for some CEMC de-scoping options.

Technically, a new module ```RawTowerCombiner``` is introduced to combine MxN towers into one readout channel on an eta-phi tower grid exiting on the DST tree. Therefore, it needs to be called after ```RawTowerBuilder``` and before ```RawTowerDigitizer```. During the merging, truth structure (Cells and showers) is maintained. In the default mode, both RawTower and RawTowerGeometry are edited on the DST tree, rather than introducing a new node for combined towers. I found this is least intrusive to our analysis code base.

If one need to use this module for ALD charge (e.g. photon position resolution, jet finding, etc.), please contact me directly, in order to speed up verification and feedback.

## Verification

The single particle simulation in the 2016-02-01 Geant4 production (```/sphenix/sim/sim01/production/2016-02-01```) was used for testing.

One example is distance between best CEMC cluster from the trajectory projection of 24 GeV/c electrons. Green reference plot is default towering and blue curve is after 2x2 ganging. The ganged distribution is roughly twice wider while the central value remain zeroed.
<img width="194" alt="2x2test" src="https://cloud.githubusercontent.com/assets/7947083/15116941/546e3fd2-15d3-11e6-8a19-3e720c302ddd.png">
Note: here cluster position is calculated with simple energy weighted average, without more sophisticated discretization corrections.

## Example macros

By default, this new module is NOT used in analysis and no ganging is applied.

Example macro to enabling 2x2 ganging as proposed for some de-scoping otpions: https://github.com/blackcathj/macros/blob/EMCal2x2/macros/g4simulations/G4_CEmc_Spacal.C
Or specifically these lines added:

\code

void CEMC_Towers(int verbosity = 0) {
 ...
  // TowerBuilder
 ...

  // Make ganged output for CEMC
  if (combin_CEMC_tower_2x2)
  {
    // group CEMC RawTower to CEMC2x2
    RawTowerCombiner * TowerCombiner = new RawTowerCombiner("RawTowerCombiner_CEMC");
    TowerCombiner->Detector("CEMC");
    TowerCombiner->set_combine_eta(2);
    TowerCombiner->set_combine_phi(2);
//    TowerCombiner->Verbosity(RawTowerCombiner::VERBOSITY_SOME);
    se->registerSubsystem( TowerCombiner );
  }

 // TowerDigitizer
 ...
}

\endcode

 * */

class RawTowerCombiner : public SubsysReco
{
 public:
  RawTowerCombiner(const std::string &name = "RawTowerCombiner");

  ~RawTowerCombiner() override
  {
  }

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void
  Detector(const std::string &d)
  {
    detector = d;
  }

  //! prefix to the tower node
  std::string
  get_tower_node_prefix() const
  {
    return _tower_node_prefix;
  }

  //! prefix to the tower node
  void
  set_tower_node_prefix(const std::string &simTowerNodePrefix)
  {
    _tower_node_prefix = simTowerNodePrefix;
  }

  //! number of eta towers to be merged into a new tower
  unsigned int
  get_combine_eta() const
  {
    return _n_combine_eta;
  }

  //! number of eta towers to be merged into a new tower
  void
  set_combine_eta(unsigned int combineEta)
  {
    _n_combine_eta = combineEta;
  }

  //! number of eta towers to be merged into a new tower
  unsigned int
  get_combine_phi() const
  {
    return _n_combine_phi;
  }

  //! number of eta towers to be merged into a new tower
  void
  set_combine_phi(unsigned int combinePhi)
  {
    _n_combine_phi = combinePhi;
  }

 protected:
  //! prefix to the tower node
  std::string _tower_node_prefix;

  //! number of eta towers to be merged into a new tower
  unsigned int _n_combine_eta;
  //! number of phi towers to be merged into a new tower
  unsigned int _n_combine_phi;

  //! get the new ieta from the old
  inline int get_output_bin_eta(int input_bin) const { return input_bin / _n_combine_eta; }
  //! get the new iphi from the old
  inline int get_output_bin_phi(int input_bin) const { return input_bin / _n_combine_phi; }

  void
  CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer *_towers;

  std::string detector;
};

#endif
