#ifndef RawTowerDigitizer_H__
#define RawTowerDigitizer_H__

#include <fun4all/SubsysReco.h>
#include <string>

#include <phool/PHTimeServer.h>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeom;

class RawTower;

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

//! simple tower digitizer which sum all cell to produce photon yield and pedstal noises
class RawTowerDigitizer : public SubsysReco
{

public:
  RawTowerDigitizer(const std::string& name = "RawTowerDigitizer");
  virtual
  ~RawTowerDigitizer();

  int
  InitRun(PHCompositeNode *topNode);
  int
  process_event(PHCompositeNode *topNode);
  int
  End(PHCompositeNode *topNode);
  void
  Detector(const std::string &d)
  {
    detector = d;
  }

  void
  set_seed(const unsigned int iseed);
  unsigned int
  get_seed() const
  {
    return seed;
  }

  enum enu_digi_algorithm
  {
    //! simple digitalization with photon statistics, ADC conversion and pedstal
    ksimple_photon_digitalization

  };

  enu_digi_algorithm
  get_digi_algorithm() const
  {
    return _digi_algorithm;
  }

  void
  set_digi_algorithm(enu_digi_algorithm digiAlgorithm)
  {
    _digi_algorithm = digiAlgorithm;
  }

  double
  get_pedstal_central_ADC() const
  {
    return _pedstal_central_ADC;
  }

  void
  set_pedstal_central_ADC(double pedstalCentralAdc)
  {
    _pedstal_central_ADC = pedstalCentralAdc;
  }

  double
  get_pedstal_width_ADC() const
  {
    return _pedstal_width_ADC;
  }

  void
  set_pedstal_width_ADC(double pedstalWidthAdc)
  {
    _pedstal_width_ADC = pedstalWidthAdc;
  }

  double
  get_photonelec_ADC() const
  {
    return _photonelec_ADC;
  }

  void
  set_photonelec_ADC(double photonelecAdc)
  {
    _photonelec_ADC = photonelecAdc;
  }

  double
  get_photonelec_yield_visible_Ge_V() const
  {
    return _photonelec_yield_visible_GeV;
  }

  void
  set_photonelec_yield_visible_Ge_V(double photonelecYieldVisibleGeV)
  {
    _photonelec_yield_visible_GeV = photonelecYieldVisibleGeV;
  }

  double
  get_zero_suppression_ADC() const
  {
    return _zero_suppression_ADC;
  }

  void
  set_zero_suppression_ADC(double zeroSuppressionAdc)
  {
    _zero_suppression_ADC = zeroSuppressionAdc;
  }

protected:
  void
  CreateNodes(PHCompositeNode *topNode);

  enu_digi_algorithm _digi_algorithm;

  //! simple digitalization with photon statistics, ADC conversion and pedstal
  //! \param  sim_tower simulation tower input
  //! \return a new RawTower object contain digitalized value of ADC output in RawTower::get_energy()
  RawTower *
  simple_photon_digitalization(RawTower * sim_tower);

  RawTowerContainer* _sim_towers;
  RawTowerContainer* _raw_towers;
  RawTowerGeom *rawtowergeom;

  std::string detector;
  std::string SimTowerNodeName;
  std::string RawTowerNodeName;
  std::string TowerGeomNodeName;

  //! photon electron yield per GeV of visible energy
  double _photonelec_yield_visible_GeV;

  //! photon electron per ADC unit
  double _photonelec_ADC;

  //! pedstal central in unit of ADC
  double _pedstal_central_ADC;

  //! pedstal width in unit of ADC
  double _pedstal_width_ADC;

  //! zero suppression in unit of ADC
  double _zero_suppression_ADC;

  PHTimeServer::timer _timer;

  unsigned int seed;
#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif
};

#endif /* RawTowerDigitizer_H__ */
