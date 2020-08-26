#ifndef G4CALO_RAWTOWERDIGITIZER_H
#define G4CALO_RAWTOWERDIGITIZER_H

#include <fun4all/SubsysReco.h>

#include <phparameter/PHParameters.h>

#include <string>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;
class RawTowerDeadMap;

class RawTower;

// rootcint barfs with this header so we need to hide it
#if !defined(__CINT__) || defined(__CLING__)
#include <gsl/gsl_rng.h>
#endif

//! simple tower digitizer which sum all cell to produce photon yield and pedstal noises
//! default input DST node is TOWER_SIM_DETECTOR
//! default output DST node is TOWER_RAW_DETECTOR
class RawTowerDigitizer : public SubsysReco
{
 public:
  RawTowerDigitizer(const std::string &name = "RawTowerDigitizer");
  virtual ~RawTowerDigitizer();

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  void Detector(const std::string &d) { m_Detector = d; _tower_params.set_name(d);}
  void TowerType(const int type) { m_TowerType = type; }
  void set_seed(const unsigned int iseed);
  unsigned int get_seed() const { return m_Seed; }
  enum enu_digi_algorithm
  {
    //! directly pass the energy of sim tower to digitized tower
    kNo_digitization = 0,
    //! wrong spelling, kept for macro compatibility
    kNo_digitalization = 0,

    //! simple digitization with photon statistics, single amplitude ADC conversion and pedestal
    kSimple_photon_digitization = 1,
    //! wrong spelling, kept for macro compatibility
    kSimple_photon_digitalization = 1,

    //! digitization with photon statistics on SiPM with an effective pixel N, ADC conversion and pedestal
    kSiPM_photon_digitization = 2
  };

  enu_digi_algorithm
  get_digi_algorithm() const
  {
    return m_DigiAlgorithm;
  }

  void
  set_digi_algorithm(enu_digi_algorithm digiAlgorithm)
  {
    m_DigiAlgorithm = digiAlgorithm;
  }

  double
  get_pedstal_central_ADC() const
  {
    return m_PedstalCentralADC;
  }

  void
  set_pedstal_central_ADC(const double pedstalCentralAdc)
  {
    m_PedstalCentralADC = pedstalCentralAdc;
  }

  double
  get_pedstal_width_ADC() const
  {
    return m_PedstalWidthADC;
  }

  void
  set_pedstal_width_ADC(const double pedstalWidthAdc)
  {
    m_PedstalWidthADC = pedstalWidthAdc;
  }

  double
  get_photonelec_ADC() const
  {
    return m_PhotonElecADC;
  }

  void
  set_photonelec_ADC(const double photonelecAdc)
  {
    m_PhotonElecADC = photonelecAdc;
  }

  double
  get_photonelec_yield_visible_GeV() const
  {
    return m_PhotonElecYieldVisibleGeV;
  }

  void
  set_photonelec_yield_visible_GeV(const double photonelecYieldVisibleGeV)
  {
    m_PhotonElecYieldVisibleGeV = photonelecYieldVisibleGeV;
  }

  double
  get_zero_suppression_ADC() const
  {
    return m_ZeroSuppressionADC;
  }

  void
  set_zero_suppression_ADC(const double zeroSuppressionAdc)
  {
    m_ZeroSuppressionADC = zeroSuppressionAdc;
  }

  void
  set_variable_zero_suppression(const bool value)
  {
    m_ZeroSuppressionFile = value;
  }

  void
  set_variable_pedestal(const bool value)
  {
    m_pedestalFile = value;
  }

  PHParameters &
  GetParameters()
  {
    return _tower_params;
  }

  std::string
  get_raw_tower_node_prefix() const
  {
    return m_RawTowerNodePrefix;
  }

  void
  set_raw_tower_node_prefix(const std::string &rawTowerNodePrefix)
  {
    m_RawTowerNodePrefix = rawTowerNodePrefix;
  }

  std::string
  get_sim_tower_node_prefix() const
  {
    return m_SimTowerNodePrefix;
  }

  void
  set_sim_tower_node_prefix(const std::string &simTowerNodePrefix)
  {
    m_SimTowerNodePrefix = simTowerNodePrefix;
  }

  // ! SiPM effective pixel per tower, only used with kSiPM_photon_digitalization
  void set_sipm_effective_pixel(const unsigned int &d) { m_SiPMEffectivePixel = d; }

  // ! SiPM effective pixel per tower, only used with kSiPM_photon_digitalization
  unsigned int get_sipm_effective_pixel() { return m_SiPMEffectivePixel; }

 private:
  void CreateNodes(PHCompositeNode *topNode);

  //! simple digitization with photon statistics, ADC conversion and pedstal
  //! \param  sim_tower simulation tower input
  //! \return a new RawTower object contain digitalized value of ADC output in RawTower::get_energy()
  RawTower *simple_photon_digitization(RawTower *sim_tower);

  //! digitization with photon statistics on SiPM with an effective pixel N, ADC conversion and pedestal
  //! this function use the effective pixel to count for the effect that the sipm is not evenly lit
  RawTower *sipm_photon_digitization(RawTower *sim_tower);

  enu_digi_algorithm m_DigiAlgorithm;

  RawTowerContainer *m_SimTowers;
  RawTowerContainer *m_RawTowers;
  RawTowerGeomContainer *m_RawTowerGeom;
  RawTowerDeadMap *m_DeadMap;

  std::string m_Detector;

  std::string m_SimTowerNodePrefix;
  std::string m_RawTowerNodePrefix;

  //! photon electron yield per GeV of visible energy
  double m_PhotonElecYieldVisibleGeV;

  //! photon electron per ADC unit
  double m_PhotonElecADC;

  //! pedstal central in unit of ADC
  double m_PedstalCentralADC;

  //! pedstal width in unit of ADC
  double m_PedstalWidthADC;

  //! pedestal from file
  bool m_pedestalFile;

  //! zero suppression in unit of ADC
  double m_ZeroSuppressionADC;

  //! zero suppression from file
  bool m_ZeroSuppressionFile;

  //! tower type to act on
  int m_TowerType;

  unsigned int m_Seed;

  // ! SiPM effective pixel per tower, only used with kSiPM_photon_digitalization
  unsigned int m_SiPMEffectivePixel;

  PHParameters _tower_params;

#if !defined(__CINT__) || defined(__CLING__)
  gsl_rng *m_RandomGenerator;
#endif
};

#endif /* G4CALO_RAWTOWERDIGITIZER_H */
