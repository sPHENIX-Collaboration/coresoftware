// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4HCALCELLRECO_H
#define G4DETECTORS_PHG4HCALCELLRECO_H

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <cmath>
#include <string>

class PHCompositeNode;

class PHG4HcalCellReco : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4HcalCellReco(const std::string &name = "HcalCellReco");

  ~PHG4HcalCellReco() override {}

  //! module initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  void SetDefaultParameters() override;

  void Detector(const std::string &d) { detector = d; }
  void checkenergy(const int i = 1) { chkenergyconservation = i; }

  void set_timing_window(const double tmi, const double tma);

  void set_fixed_energy(const double efix) { m_FixedEnergy = efix; }

 protected:
  int CheckEnergy(PHCompositeNode *topNode);
  std::string detector;
  std::string hitnodename;
  std::string cellnodename;

  int chkenergyconservation = 0;

  double tmin = NAN;
  double tmax = NAN;
  double m_DeltaT = NAN;
  double m_FixedEnergy = NAN;
};

#endif
