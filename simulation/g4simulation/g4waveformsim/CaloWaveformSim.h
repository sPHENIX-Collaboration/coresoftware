// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!  
 * \file CaloWaveformSim.h
 * \brief create waveform from hits, could also use for event overlay
 * \author Shuhang Li <sli7@bnl.gov>
 * \version $Revision:   $
 * \date    $Date: $
 */
#ifndef CALOWAVEFORMSIM_H
#define CALOWAVEFORMSIM_H

#include <calobase/TowerInfoDefs.h>
#include <caloreco/CaloTowerDefs.h>
#include <fun4all/SubsysReco.h>
#include <g4detectors/LightCollectionModel.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <vector>

class PHCompositeNode;
class TProfile;
class PHG4Hit;
class PHG4CylinderCellGeom_Spacalv1;
class PHG4CylinderGeom_Spacalv3;
class TRandom3;
class TTree;
class CDBTTree;
class TowerInfoContainer;

class CaloWaveformSim : public SubsysReco
{
public:
  CaloWaveformSim(const std::string &name = "CaloWaveformSim");
  ~CaloWaveformSim() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  // Detector configuration
  void set_detector_type(CaloTowerDefs::DetectorSystem dettype) { m_dettype = dettype; }
  void set_detector(const std::string &detector) { m_detector = detector; }

  // Calibration settings (data energy)
  void set_fieldname(const std::string &fieldname) { m_fieldname = fieldname; m_overrideFieldName = true; }
  void set_calibName(const std::string &calibName) { m_calibName = calibName; m_overrideCalibName = true; }
  void set_directURL_calib(const std::string &url) { m_giveDirectURL = true; m_directURL = url; }
  void set_overrideCalibName(bool overrideCalib) { m_overrideCalibName = overrideCalib; }
  void set_overrideFieldName(bool overrideField) { m_overrideFieldName = overrideField; }

  // Calibration settings (MC energy)
  void set_MC_fieldname(const std::string &MC_fieldname) { m_MC_fieldname = MC_fieldname; m_overrideMCFieldName = true; }
  void set_MC_calibName(const std::string &MC_calibName) { m_MC_calibName = MC_calibName; m_overrideMCCalibName = true; }
  void set_directURL_MCcalib(const std::string &url) { m_giveDirectURL_MC = true; m_directURL_MC = url; }
  void set_overrideMCFieldName(bool overrideField) { m_overrideMCFieldName = overrideField; }
  void set_overrideMCCalibName(bool overrideCalib) { m_overrideMCCalibName = overrideCalib; }

  // Time calibration (data)
  void set_fieldname_time(const std::string &fieldname_time) { m_fieldname_time = fieldname_time; m_overrideTimeFieldName = true; }
  void set_calibName_time(const std::string &calibName_time) { m_calibName_time = calibName_time; m_overrideTimeCalibName = true; }
  void set_directURL_timecalib(const std::string &url) { m_giveDirectURL_time = true; m_directURL_time = url; }
  void set_dotimecalib(bool dotimecalib) { m_dotimecalib = dotimecalib; }
  void set_overrideTimeFieldName(bool overrideField) { m_overrideTimeFieldName = overrideField; }
  void set_overrideTimeCalibName(bool overrideCalib) { m_overrideTimeCalibName = overrideCalib; }

  // Time calibration (MC)
  void set_MC_fieldname_time(const std::string &MC_fieldname_time) { m_MC_fieldname_time = MC_fieldname_time; m_overrideMCTimeFieldName = true; }
  void set_MC_calibName_time(const std::string &MC_calibName_time) { m_MC_calibName_time = MC_calibName_time; m_overrideMCTimeCalibName = true; }
  void set_directURL_MCtimecalib(const std::string &url) { m_giveDirectURL_MC_time = true; m_directURL_MC_time = url; }
  void set_overrideMCTimeFieldName(bool overrideField) { m_overrideMCTimeFieldName = overrideField; }
  void set_overrideMCTimeCalibName(bool overrideCalib) { m_overrideMCTimeCalibName = overrideCalib; }

  // Waveform template & sampling
  void set_templatefile(const std::string &templatefile) { m_templatefile = templatefile; }
  void set_nsamples(int nsamples) { m_nsamples = nsamples; }
  void set_pedestalsamples(int pedestalsamples) { m_pedestalsamples = pedestalsamples; }
  void set_sampletime(float sampletime) { m_sampletime = sampletime; }
  void set_nchannels(int nchannels) { m_nchannels = nchannels; }
  void set_sampling_fraction(float fraction) { m_sampling_fraction = fraction; }

  // Signal shaping parameters
  void set_deltaT(float deltaT) { m_deltaT = deltaT; }
  void set_timewidth(float timewidth) { m_timeshiftwidth = timewidth; }
  void set_peakpos(float peakpos) { m_peakpos = peakpos; }
  void set_highgain(bool highgain = true) { m_highgain = highgain; }
  void set_gain(int gain) { m_gain = gain; }
  void set_pedestal_scale(float scale) { m_pedestal_scale = scale; }

  // Noise configuration
  enum NoiseType { NOISE_NONE = 0, NOISE_GAUSSIAN = 1, NOISE_TREE = 2 };
  void set_noise_type(NoiseType noiseType) { m_noiseType = noiseType; }
  void set_fixpedestal(int fixpedestal) { m_fixpedestal = fixpedestal; }
  void set_gaussian_noise(int gaussian_noise) { m_gaussian_noise = gaussian_noise; }

  // Light collection model access
  LightCollectionModel &get_light_collection_model() { return light_collection_model; }

private:
  CaloTowerDefs::DetectorSystem m_dettype{CaloTowerDefs::CEMC};
  std::string m_detector{"CEMC"};

  // Data energy calibration
  std::string m_fieldname{"Femc_datadriven_qm1_correction"};
  std::string m_calibName{"cemc_pi0_twrSlope_v1"};
  bool        m_overrideCalibName{false};
  bool        m_overrideFieldName{false};
  bool        m_giveDirectURL{false};
  std::string m_directURL{""};
  // MC energy calibration
  std::string m_MC_fieldname{"Femc_datadriven_qm1_correction"};
  std::string m_MC_calibName{"cemc_pi0_twrSlope_v1"};
  bool        m_overrideMCFieldName{false};
  bool        m_overrideMCCalibName{false};
  bool        m_giveDirectURL_MC{false};
  std::string m_directURL_MC{""};

  // Data time calibration
  std::string m_fieldname_time{"time"};
  std::string m_calibName_time{"CEMC_meanTime"};
  bool        m_overrideTimeFieldName{false};
  bool        m_overrideTimeCalibName{false};
  bool        m_dotimecalib{true};
  bool        m_giveDirectURL_time{false};
  std::string m_directURL_time{""};
  // MC time calibration
  std::string m_MC_fieldname_time{"time"};
  std::string m_MC_calibName_time{"CEMC_meanTime"};
  bool        m_overrideMCTimeFieldName{false};
  bool        m_overrideMCTimeCalibName{false};
  bool        m_giveDirectURL_MC_time{false};
  std::string m_directURL_MC_time{""};

  // Waveform settings
  std::string m_templatefile{"waveformtemptempohcalcosmic.root"};
  int         m_nsamples{31};
  int         m_pedestalsamples{31};
  float       m_sampletime{50. / 3.};
  int         m_nchannels{24576};
  float       m_sampling_fraction{1.0f};

  //containers
  TowerInfoContainer *m_CaloWaveformContainer{nullptr};
  TowerInfoContainer *m_PedestalContainer{nullptr};

  // Shaping & noise
  int   m_fixpedestal{1500};
  int   m_gaussian_noise{3};
  float m_deltaT{100.};
  float m_timeshiftwidth{0.};
  bool  m_highgain{false};
  int   m_gain{1};
  float m_peakpos{6.};
  float m_pedestal_scale{1.};

  gsl_rng *m_RandomGenerator{nullptr};
  PHG4CylinderCellGeom_Spacalv1 *geo{nullptr};
  const PHG4CylinderGeom_Spacalv3  *layergeom{nullptr};
  std::vector<std::vector<float>> m_waveforms;
  int m_runNumber{0};

  unsigned int (*encode_tower)(unsigned int, unsigned int){TowerInfoDefs::encode_emcal};
  unsigned int (*decode_tower)(unsigned int){TowerInfoDefs::decode_emcal};

  CDBTTree *cdbttree{nullptr}, *cdbttree_MC{nullptr};
  CDBTTree *cdbttree_time{nullptr}, *cdbttree_MC_time{nullptr};
  TProfile *h_template{nullptr};
  LightCollectionModel light_collection_model;

  NoiseType m_noiseType{NOISE_TREE};

  void CreateNodeTree(PHCompositeNode *topNode);
  void maphitetaphi(PHG4Hit *g4hit,
                    unsigned short &etabin,
                    unsigned short &phibin,
                    float &correction);
  double template_function(double *x, double *par);
};

#endif  // CALOWAVEFORMSIM_H