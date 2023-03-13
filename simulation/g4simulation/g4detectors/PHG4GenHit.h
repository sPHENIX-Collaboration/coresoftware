// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4GENHIT_H
#define G4DETECTORS_PHG4GENHIT_H

#include <fun4all/SubsysReco.h>

#include <cmath>
#include <string>  // for string

class PHCompositeNode;

class PHG4GenHit : public SubsysReco
{
 public:
  PHG4GenHit(const std::string &name = "PHG4GenHit");
  ~PHG4GenHit() override {}

  int process_event(PHCompositeNode *topNode) override;

  void set_phi(const double d) { phi = d; }
  void set_theta(const double d) { theta = d; }
  void set_eloss(const double d) { eloss = d; }
  void set_layer(const int i) { layer = i; }
  void Detector(const std::string &n) { detector = n; }

 protected:
  double phi = NAN;
  double theta = NAN;
  double eloss = NAN;
  int layer = -9999;
  std::string detector;
};

#endif
