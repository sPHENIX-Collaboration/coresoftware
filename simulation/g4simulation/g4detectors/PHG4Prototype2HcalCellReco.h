// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4DETECTORS_PHG4PROTOTYPE2HCALCELLRECO_H
#define G4DETECTORS_PHG4PROTOTYPE2HCALCELLRECO_H

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>
#include <string>

class PHCompositeNode;

class PHG4Prototype2HcalCellReco : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4Prototype2HcalCellReco(const std::string &name = "Prototype2HcalCELLRECO");

  virtual ~PHG4Prototype2HcalCellReco() {}
  //! module initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  void SetDefaultParameters();

  void Detector(const std::string &d) { m_Detector = d; }
  void checkenergy(const int i = 1) { m_CheckEnergyConservationFlag = i; }
  void set_timing_window(const double tmi, const double tma);

 private:
  int CheckEnergy(PHCompositeNode *topNode);
  std::string m_Detector;
  std::string m_HitNodeName;
  std::string m_CellNodeName;
  int m_CheckEnergyConservationFlag;

  double m_Tmin;
  double m_Tmax;
};

#endif  // G4DETECTORS_PHG4PROTOTYPE2HCALCELLRECO_H
