#ifndef PHG4INPUTFILTER_H__
#define PHG4INPUTFILTER_H__

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class PHG4InputFilter : public SubsysReco
{
 public:
  PHG4InputFilter(const std::string &name = "G4INPUTFILTER");
  virtual ~PHG4InputFilter() {}

  int process_event(PHCompositeNode *topNode);

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

