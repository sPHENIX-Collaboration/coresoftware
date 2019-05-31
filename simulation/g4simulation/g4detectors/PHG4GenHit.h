#ifndef G4DETECTORS_PHG4GENHIT_H
#define G4DETECTORS_PHG4GENHIT_H

#include <fun4all/SubsysReco.h>

#include <string>                // for string

class PHCompositeNode;

class PHG4GenHit: public SubsysReco
{
 public:
  PHG4GenHit(const std::string &name = "PHG4GenHit");
  virtual ~PHG4GenHit() {}

  int process_event(PHCompositeNode *topNode);

  void set_phi(const double d) {phi = d;}
  void set_theta(const double d) {theta = d;}
  void set_eloss(const double d) {eloss = d;}
  void set_layer(const int i) {layer = i;}
  void Detector(const std::string &n) {detector = n;}

 protected:
  double phi;
  double theta;
  double eloss;
  int layer;
  std::string detector;

};

#endif
