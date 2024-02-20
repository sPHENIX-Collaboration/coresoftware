// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FULLPROJSPACALCELLRECO_H
#define G4DETECTORS_PHG4FULLPROJSPACALCELLRECO_H

#include "LightCollectionModel.h"

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <cmath>
#include <map>
#include <string>

class LightCollectionModel;
class PHCompositeNode;
class PHG4Cell;
class TH2;
class TH1;

class PHG4FullProjSpacalCellReco : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4FullProjSpacalCellReco(const std::string &name = "SPACALCELLRECO");

  ~PHG4FullProjSpacalCellReco() override {}

  //! module initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! reset after event processing
  int ResetEvent(PHCompositeNode *topNode) override;

  void SetDefaultParameters() override;

  void Detector(const std::string &d) { detector = d; }

  void checkenergy(const int i = 1) { chkenergyconservation = i; }

  void set_timing_window(const double tmin, const double tmax);

  LightCollectionModel &get_light_collection_model() { return light_collection_model; }

 protected:
  int CheckEnergy(PHCompositeNode *topNode);

  std::string detector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;

  double sum_energy_g4hit = 0;
  int chkenergyconservation = 0;
  std::map<unsigned int, PHG4Cell *> celllist;

  //! timing window size in ns. This is for a simple simulation of the ADC integration window starting from 0ns to this value. Default to infinity, i.e. include all hits
  double tmin = NAN;
  double tmax = NAN;
  double m_DeltaT = NAN;
  LightCollectionModel light_collection_model;
};

#endif
