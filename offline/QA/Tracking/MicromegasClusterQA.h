// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_TRACKING_MICROMEGASCLUSTERQA_H_
#define QA_TRACKING_MICROMEGASCLUSTERQA_H_

#include <fun4all/SubsysReco.h>
#include <micromegas/MicromegasCalibrationData.h>
#include <micromegas/MicromegasDefs.h>
#include <micromegas/MicromegasMapping.h>

#include <array>
#include <map>
#include <set>
#include <string>

class ActsGeometry;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;

class TH1;
class TH2;

class PHCompositeNode;

class MicromegasClusterQA : public SubsysReco
{
 public:
  MicromegasClusterQA(const std::string& name = "MicromegasClusterQA");

  ~MicromegasClusterQA() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;

  /// set default pedestal
  void set_default_pedestal(double value)
  {
    m_default_pedestal = value;
  }

  /// set whether default pedestal is used or not
  void set_use_default_pedestal(bool value)
  {
    m_use_default_pedestal = value;
  }

  /// calibration file
  void set_calibration_file(const std::string& value)
  {
    m_calibration_filename = value;
  }

  /// set min sample for signal hits
  void set_sample_min(uint16_t value) { m_sample_min = value; }

  /// set max sample for signal hits
  void set_sample_max(uint16_t value) { m_sample_max = value; }


 private:
  void create_histograms();

  std::string get_histogram_prefix() const;

  /// micromegas mapping
  MicromegasMapping m_mapping;

  /// per detector cluster multiplicity distribution
  TH2* m_h_cluster_multiplicity = nullptr;

  /// per detector cluster size distribution
  TH2* m_h_cluster_size = nullptr;

  /// per detector cluster charge distribution
  TH2* m_h_cluster_charge = nullptr;

  /// per detector reference cluster count histogram
  /*! used for standalone efficiency calculation */
  TH1* m_h_cluster_count_ref = nullptr;

  /// per detector found cluster count histogram
  /** used for standalone efficiency calculation */
  TH1* m_h_cluster_count_found = nullptr;

  /// Acts tracking geometry for surface lookup
  ActsGeometry* m_tGeometry = nullptr;

  //! hits
  TrkrHitSetContainer* m_hitsetcontainer = nullptr;

  //! clusters
  TrkrClusterContainer* m_cluster_map = nullptr;

  //! cluster to hit association
  TrkrClusterHitAssoc* m_cluster_hit_map = nullptr;

  /// first micromegas layer
  /* this is updated on the fly from geometry object */
  int m_firstlayer = 55;

  /// number of layers
  int m_nlayers = 2;

  /// keep track of detector names
  std::vector<std::string> m_detector_names;

  /// min sample for signal
  uint16_t m_sample_min = 0;

  /// max sample for signal
  uint16_t m_sample_max = 1024;

  ///@name calibration filename
  //@{

  /// if true, use default pedestal to get hit charge. Relies on calibration data otherwise
  bool m_use_default_pedestal = true;

  /// default pedestal
  double m_default_pedestal = 74.6;

  /// calibration filename
  std::string m_calibration_filename;

  /// calibration data
  MicromegasCalibrationData m_calibration_data;

  //@}
};

#endif  // MicromegasClusterQA_H
