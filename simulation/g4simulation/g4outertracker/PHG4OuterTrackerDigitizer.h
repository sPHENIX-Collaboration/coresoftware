// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4MVTX_PHG4OUTERTRACKERDIGITIZER_H
#define G4MVTX_PHG4OUTERTRACKERDIGITIZER_H

#include <fun4all/SubsysReco.h>

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

#include <map>
#include <string>                // for string
#include <utility>               // for pair, make_pair
#include <vector>

class PHCompositeNode;

class PHG4OuterTrackerDigitizer : public SubsysReco
{
 public:
  PHG4OuterTrackerDigitizer(const std::string &name = "PHG4OuterTrackerDigitizer");
  virtual ~PHG4OuterTrackerDigitizer();

  //! module initialization
  int Init(PHCompositeNode *topNode) { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode *topNode) { return 0; };

  void DigitizeOuterTracker(PHCompositeNode *topNode);
  void CalculateOuterTrackerADCScale(PHCompositeNode *topNode);

 private:

  // settings
  double _max_adc;
  double _energy_scale;

#ifndef __CINT__
  //! random generator that conform with sPHENIX standard
  gsl_rng *RandomGenerator;
#endif
};

#endif
