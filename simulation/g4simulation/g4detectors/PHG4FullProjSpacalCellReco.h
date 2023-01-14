// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FULLPROJSPACALCELLRECO_H
#define G4DETECTORS_PHG4FULLPROJSPACALCELLRECO_H

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <cmath>
#include <map>
#include <string>

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

  class LightCollectionModel
  {
   public:
    LightCollectionModel();

    //! delete copy ctor and assignment opertor (cppcheck)
    explicit LightCollectionModel(const LightCollectionModel &) = delete;
    LightCollectionModel &operator=(const LightCollectionModel &) = delete;

    virtual ~LightCollectionModel();

    //! input data file
    void load_data_file(const std::string &input_file, const std::string &histogram_light_guide_model, const std::string &histogram_fiber_model);

    //! load from CDB
    void load_data_from_CDB(const std::string &domain, const std::string &histogram_light_guide_model, const std::string &histogram_fiber_model);

    //! Whether use light collection model
    bool use_light_guide_model() const { return data_grid_light_guide_efficiency != nullptr; }

    //! Whether use Light Transmission Efficiency model for the fiber
    bool use_fiber_model() const { return data_grid_fiber_trans != nullptr; }

    //! get Light Collection Efficiency for the light guide as function of x,y position in fraction of tower width
    double get_light_guide_efficiency(const double x_fraction, const double y_fraction);

    //! get Light Transmission Efficiency for the fiber as function of z position (cm) in the fiber. Z=0 is at the middle of the fiber
    double get_fiber_transmission(const double z_distance);

   private:
    //! 2-D data grid for Light Collection Efficiency for the light guide as function of x,y position in fraction of tower width
    TH2 *data_grid_light_guide_efficiency = nullptr;

    //! 1-D data grid for the light transmission efficiency in the fiber as function of distance to location in the fiber. Z=0 is at the middle of the fiber
    TH1 *data_grid_fiber_trans = nullptr;

    // These two histograms are handed off to Fun4All and will be deleted there
    // this suppresses the cppcheck warning
    // cppcheck-suppress unsafeClassCanLeak
    TH2 *data_grid_light_guide_efficiency_verify = nullptr;
    // cppcheck-suppress unsafeClassCanLeak
    TH1 *data_grid_fiber_trans_verify = nullptr;
  };

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
