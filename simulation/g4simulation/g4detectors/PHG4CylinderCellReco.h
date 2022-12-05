// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLRECO_H
#define G4DETECTORS_PHG4CYLINDERCELLRECO_H

#include <phparameter/PHParameterContainerInterface.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <utility>  // for pair

class PHCompositeNode;
class PHG4Cell;

class PHG4CylinderCellReco : public SubsysReco, public PHParameterContainerInterface
{
 public:
  explicit PHG4CylinderCellReco(const std::string &name = "CYLINDERRECO");

  ~PHG4CylinderCellReco() override {}

  //! module initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  void SetDefaultParameters() override;

  void Detector(const std::string &d);
  void cellsize(const int i, const double sr, const double sz);
  void etaphisize(const int i, const double deltaeta, const double deltaphi);
  void checkenergy(const int i = 1) { chkenergyconservation = i; }
  void OutputDetector(const std::string &d) { outdetector = d; }

  double get_timing_window_min(const int i) { return tmin_max[i].first; }
  double get_timing_window_max(const int i) { return tmin_max[i].second; }
  void set_timing_window(const int detid, const double tmin, const double tmax);

 protected:
  void set_size(const int i, const double sizeA, const double sizeB);
  int CheckEnergy(PHCompositeNode *topNode);

  std::map<int, int> binning;
  std::map<int, std::pair<double, double> > cell_size;  // cell size in phi/z
  std::map<int, std::pair<double, double> > zmin_max;   // zmin/zmax for each layer for faster lookup
  std::map<int, double> phistep;
  std::map<int, double> etastep;
  std::set<int> implemented_detid;
  std::string detector;
  std::string outdetector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;
  std::map<int, std::pair<int, int>> n_phi_z_bins;
  std::map<unsigned long long, PHG4Cell *> cellptmap;  // This map holds the hit cells
  std::map<unsigned long long, PHG4Cell *>::iterator it;
  std::map<int, std::pair<double, double>> tmin_max;
  std::map<int, double> m_DeltaTMap;

  int nbins[2];
  int chkenergyconservation = 0;

  double sum_energy_before_cuts = 0.;
  double sum_energy_g4hit = 0.;
};

#endif
