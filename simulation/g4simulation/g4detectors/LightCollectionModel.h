// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_LIGHTCOLLECTIONMODEL_H
#define G4DETECTORS_LIGHTCOLLECTIONMODEL_H
#include <string>

class TH2;
class TH1;

class LightCollectionModel
{
 public:
  LightCollectionModel() = default;

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
  TH2 *data_grid_light_guide_efficiency{nullptr};

  //! 1-D data grid for the light transmission efficiency in the fiber as function of distance to location in the fiber. Z=0 is at the middle of the fiber
  TH1 *data_grid_fiber_trans{nullptr};

  // These two histograms are handed off to Fun4All and will be deleted there
  // this suppresses the cppcheck warning
  // cppcheck-suppress unsafeClassCanLeak
  TH2 *data_grid_light_guide_efficiency_verify{nullptr};
  // cppcheck-suppress unsafeClassCanLeak
  TH1 *data_grid_fiber_trans_verify{nullptr};
};

#endif
