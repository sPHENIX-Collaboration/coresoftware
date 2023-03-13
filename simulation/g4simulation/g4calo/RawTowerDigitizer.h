#ifndef G4CALO_RAWTOWERDIGITIZER_H
#define G4CALO_RAWTOWERDIGITIZER_H

#include <fun4all/SubsysReco.h>

#include <phparameter/PHParameters.h>
#include <calobase/TowerInfoContainer.h>

#include <gsl/gsl_rng.h>

#include <cmath>
#include <string>

class CaloCalibSimpleCorrFile;
class PHCompositeNode;
class RawTower;
class TowerInfo;
class RawTowerContainer;
class RawTowerGeomContainer;
class RawTowerDeadMap;

//! simple tower digitizer which sum all cell to produce photon yield and pedstal noises
//! default input DST node is TOWER_SIM_DETECTOR
//! default output DST node is TOWER_RAW_DETECTOR
class RawTowerDigitizer : public SubsysReco
{
 public:
  RawTowerDigitizer(const std::string &name = "RawTowerDigitizer");
  ~RawTowerDigitizer() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void Detector(const std::string &d)
  {
    m_Detector = d;
    _tower_params.set_name(d);
  }
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

  enum ProcessTowerType
  {
    kRawTowerOnly= 0,
    kTowerInfoOnly = 1,
    kBothTowers =2
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

  // calo calib decal stuff JEF Feb 2022
  void set_DoTowerDecal(const bool doTowerDecal,
                        const char *decalFileName = "",
                        const bool doInverse = false)
  {
    m_DoDecal = doTowerDecal;
    m_DecalInverse = doInverse;
    set_DecalFileName(decalFileName);
  }

  void set_DecalFileName(const char *inCalFname)
  {
    m_DecalFileName = inCalFname;
  }

  void set_UseConditionsDB(const bool setUseCondDB)
  {
    m_UseConditionsDB = setUseCondDB;
  }

  void set_towerinfo(RawTowerDigitizer::ProcessTowerType UseTowerInfo )
  {
    m_UseTowerInfo = UseTowerInfo;
  }

 private:
  void CreateNodes(PHCompositeNode *topNode);

  //! simple digitization with photon statistics, ADC conversion and pedstal
  //! \param  sim_tower simulation tower input
  //! \return a new RawTower object contain digitalized value of ADC output in RawTower::get_energy()
  RawTower *simple_photon_digitization(RawTower *sim_tower);
  TowerInfo *simple_photon_digitization(TowerInfo *sim_tower);

  //! digitization with photon statistics on SiPM with an effective pixel N, ADC conversion and pedestal
  //! this function use the effective pixel to count for the effect that the sipm is not evenly lit
  RawTower *sipm_photon_digitization(RawTower *sim_tower);
  TowerInfo *sipm_photon_digitization(TowerInfo *sim_tower);

  enu_digi_algorithm m_DigiAlgorithm = kNo_digitization;

  RawTowerContainer *m_SimTowers = nullptr;
  RawTowerContainer *m_RawTowers = nullptr;

  TowerInfoContainer *m_SimTowerInfos = nullptr;
  TowerInfoContainer *m_RawTowerInfos = nullptr;



  RawTowerGeomContainer *m_RawTowerGeom = nullptr;
  RawTowerDeadMap *m_DeadMap = nullptr;

  std::string m_Detector = "NONE";

  std::string m_SimTowerNodePrefix = "SIM";
  std::string m_RawTowerNodePrefix = "RAW";

  //! photon electron yield per GeV of visible energy
  double m_PhotonElecYieldVisibleGeV = NAN;

  //! photon electron per ADC unit
  double m_PhotonElecADC = NAN;

  //! pedstal central in unit of ADC
  double m_PedstalCentralADC = NAN;

  //! pedstal width in unit of ADC
  double m_PedstalWidthADC = NAN;

  //! pedestal from file
  bool m_pedestalFile = false;

  //! zero suppression in unit of ADC
  double m_ZeroSuppressionADC = -1000;

  //! zero suppression from file
  bool m_ZeroSuppressionFile = false;

  //! tower type to act on
  int m_TowerType = -1;

  unsigned int m_Seed = 0;

  // ! SiPM effective pixel per tower, only used with kSiPM_photon_digitalization
  // ! sPHENIX EMCal default, 4x Hamamatsu S12572-015P MPPC [sPHENIX TDR]
  unsigned int m_SiPMEffectivePixel = 40000 * 4;

  PHParameters _tower_params;

  gsl_rng *m_RandomGenerator = nullptr;

  // calo calibs decal stuff JEF Feb 2022
  bool m_DoDecal = false;
  bool m_DecalInverse = false;
  bool m_Decal = true;
  std::string m_DecalFileName;
  bool m_UseConditionsDB = false;
  CaloCalibSimpleCorrFile *m_CalDBFile = nullptr;


  RawTowerDigitizer::ProcessTowerType m_UseTowerInfo = RawTowerDigitizer::ProcessTowerType::kBothTowers;  // 0 just produce RawTowers, 1 just produce TowerInfo objects, and 2 produce both


};

#endif /* G4CALO_RAWTOWERDIGITIZER_H */
