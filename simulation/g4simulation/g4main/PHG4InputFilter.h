// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4INPUTFILTER_H
#define G4MAIN_PHG4INPUTFILTER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class PHG4InputFilter : public SubsysReco
{
 public:
  PHG4InputFilter(const std::string &name = "G4INPUTFILTER");
  ~PHG4InputFilter() override {}

  int process_event(PHCompositeNode *topNode) override;

  void set_eta_range(const double min, const double max) {etamin = min; etamax = max;}
  void set_etamin(const double min) {etamin = min;}
  void set_etamax(const double max) {etamax = max;}
  void set_ptmin(const double min) {ptmin = min;}
  void set_ptmax(const double max) {ptmax = max;}

 protected:
  double get_eta(const double x, const double y, const double z);
  double etamin;
  double etamax;
  double ptmin;
  double ptmax;

};

#endif

