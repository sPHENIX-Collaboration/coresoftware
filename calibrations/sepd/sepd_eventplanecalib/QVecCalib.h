#ifndef QVECCALIB_H
#define QVECCALIB_H

#include "QVecDefs.h"

#include <fun4all/SubsysReco.h>

#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

class PHCompositeNode;
class EventPlaneData;
class EpdGeom;

/**
 * @class QVecCalib
 * @brief Orchestrates the multi-pass anisotropy calibration for the sEPD Q-vectors.
 *
 * This class implements a three-pass correction procedure designed to remove
 * detector-induced biases from the sEPD event plane reconstruction:
 * * 1. **ComputeRecentering**: Calculates the first-order <Q> vector offsets (re-centering)
 * per centrality bin.
 * 2. **ApplyRecentering**: Applies the first-order offsets and computes the
 * second-order whitening/flattening matrix.
 * 3. **ApplyFlattening**: Applies the full correction (re-centering + flattening)
 * to produce final validated event planes.
 * * The class manages event-level selections based on charge-centrality correlations
 * and handles the exclusion of "bad" (hot/cold/dead) sEPD channels.
 */
class QVecCalib : public SubsysReco
{
 public:
  explicit QVecCalib(const std::string& name = "QVecCalib");

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode* topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode* topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode* topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode* topNode) override;

  enum class Pass
  {
    ComputeRecentering,
    ApplyRecentering,
    ApplyFlattening
  };

  void set_pass(int pass)
  {
    m_pass = validate_pass(pass);
  }

  void set_input_hist(std::string_view file)
  {
    m_input_hist = file;
  }

  void set_input_Q_calib(std::string_view file)
  {
    m_input_Q_calib = file;
  }

  void set_dst_tag(std::string_view tag)
  {
    m_dst_tag = tag;
  }

  void set_cdb_output_dir(std::string_view cdb_dir)
  {
    m_cdb_output_dir = cdb_dir;
  }

 private:
  static Pass validate_pass(int pass)
  {
    switch (pass)
    {
    case 0:
      return Pass::ComputeRecentering;
    case 1:
      return Pass::ApplyRecentering;
    case 2:
      return Pass::ApplyFlattening;
    default:
      throw std::invalid_argument("Invalid pass value");
    }
  }

  struct CorrectionData
  {
    QVecShared::QVec avg_Q{};

    double avg_Q_xx{0.0};
    double avg_Q_yy{0.0};
    double avg_Q_xy{0.0};

    std::array<std::array<double, 2>, 2> X_matrix{};
  };

  static constexpr size_t m_cent_bins = QVecShared::CENT_BINS;
  static constexpr auto m_harmonics = QVecShared::HARMONICS;

  static constexpr float SIGMA_HOT = 6.0F;
  static constexpr float SIGMA_COLD = -6.0F;

  double m_cent_low{-0.5};
  double m_cent_high{79.5};

  std::string m_input_hist;
  std::string m_input_Q_calib;
  std::string m_dst_tag;
  std::string m_cdb_output_dir{"."};
  Pass m_pass{Pass::ComputeRecentering};
  EventPlaneData* m_evtdata{nullptr};

  int m_event{0};
  int m_runnumber{0};

  struct EventCounters
  {
    int bad_centrality_sepd_correlation{0};
    int zero_sepd_total_charge{0};
    int total_processed{0};
  };

  EventCounters m_event_counters;

  std::array<std::array<QVecShared::QVec, 2>, m_harmonics.size()> m_q_vectors{};

  static constexpr int PROGRESS_REPORT_INTERVAL = 1000;

  // Holds all correction data
  // key: [Cent][Harmonic][Subdetector]
  // Harmonics {2,3,4} -> 3 elements
  // Subdetectors {S,N,NS} -> 3 elements
  std::array<std::array<std::array<CorrectionData, 3>, m_harmonics.size()>, m_cent_bins> m_correction_data;

  // Store harmonic orders and subdetectors for easy iteration
  static constexpr std::array<QVecShared::Subdetector, 2> m_subdetectors = {QVecShared::Subdetector::S, QVecShared::Subdetector::N};
  static constexpr std::array<QVecShared::QComponent, 2> m_components = {QVecShared::QComponent::X, QVecShared::QComponent::Y};

  // [Harmonic Index][Channel Index] -> {cos, sin}
  std::vector<std::vector<std::pair<double, double>>> m_trig_cache;

  struct AverageHists
  {
    TProfile* S_x_avg{nullptr};
    TProfile* S_y_avg{nullptr};
    TProfile* N_x_avg{nullptr};
    TProfile* N_y_avg{nullptr};

    TH2* Psi_S{nullptr};
    TH2* Psi_N{nullptr};
    TH2* Psi_NS{nullptr};
  };

  struct RecenterHists
  {
    TProfile* S_x_corr_avg{nullptr};
    TProfile* S_y_corr_avg{nullptr};
    TProfile* N_x_corr_avg{nullptr};
    TProfile* N_y_corr_avg{nullptr};

    TProfile* S_xx_avg{nullptr};
    TProfile* S_yy_avg{nullptr};
    TProfile* S_xy_avg{nullptr};
    TProfile* N_xx_avg{nullptr};
    TProfile* N_yy_avg{nullptr};
    TProfile* N_xy_avg{nullptr};

    TProfile* NS_xx_avg{nullptr};
    TProfile* NS_yy_avg{nullptr};
    TProfile* NS_xy_avg{nullptr};

    TH2* Psi_S_corr{nullptr};
    TH2* Psi_N_corr{nullptr};
    TH2* Psi_NS_corr{nullptr};
  };

  struct FlatteningHists
  {
    TProfile* S_x_corr2_avg{nullptr};
    TProfile* S_y_corr2_avg{nullptr};
    TProfile* N_x_corr2_avg{nullptr};
    TProfile* N_y_corr2_avg{nullptr};

    TProfile* S_xx_corr_avg{nullptr};
    TProfile* S_yy_corr_avg{nullptr};
    TProfile* S_xy_corr_avg{nullptr};

    TProfile* N_xx_corr_avg{nullptr};
    TProfile* N_yy_corr_avg{nullptr};
    TProfile* N_xy_corr_avg{nullptr};

    TProfile* NS_xx_corr_avg{nullptr};
    TProfile* NS_yy_corr_avg{nullptr};
    TProfile* NS_xy_corr_avg{nullptr};

    TH2* Psi_S_corr2{nullptr};
    TH2* Psi_N_corr2{nullptr};
    TH2* Psi_NS_corr2{nullptr};
  };

  // sEPD Bad Channels
  std::unordered_set<int> m_bad_channels;

  double m_sEPD_min_avg_charge_threshold{1};
  double m_sEPD_sigma_threshold{3};

  // Hists
  TH1* hCentrality{nullptr};

  TH2* h2SEPD_Charge{nullptr};
  TH2* h2SEPD_Chargev2{nullptr};

  TH2* h2SEPD_South_Charge_rbin{nullptr};
  TH2* h2SEPD_North_Charge_rbin{nullptr};

  TH2* h2SEPD_South_Charge_rbinv2{nullptr};
  TH2* h2SEPD_North_Charge_rbinv2{nullptr};

  TProfile* hSEPD_Charge_Min{nullptr};
  TProfile* hSEPD_Charge_Max{nullptr};

  TProfile* hSEPD_Bad_Channels{nullptr};

  std::map<std::string, TH2*> m_hists2D;
  std::map<std::string, TProfile*> m_profiles;

  std::vector<AverageHists> m_average_hists;
  std::vector<RecenterHists> m_recenter_hists;
  std::vector<FlatteningHists> m_flattening_hists;

  /**
   * @brief Initializes all output histograms and profiles.
   * * Dynamically generates histogram names using the shared naming helper based on
   * the current calibration pass (e.g., adding "_corr" or "_corr2" suffixes).
   */
  void init_hists();

  /**
   * @brief Safely retrieves a ROOT object from a file and returns a managed pointer.
   * * Performs a dynamic_cast to verify the requested type T and Clones the object
   * to ensure it remains valid after the source file is closed.
   * * @tparam T The ROOT class type (e.g., TProfile).
   * @param file Pointer to the source TFile.
   * @param name The name of the object within the file.
   * @return T* A managed pointer to the cloned object.
   * @throws std::runtime_error If the object is not found or type mismatch occurs.
   */
  template <typename T>
  T* load_and_clone(TFile* file, const std::string& name);

  /**
   * @brief Loads the results of previous passes from a calibration file.
   * * Populates the internal correction data structure with averages and/or
   * matrices required for the current processing pass.
   */
  int load_correction_data();

  /**
   * @brief Validates events based on sEPD total charge vs. centrality correlation.
   * * Compares the current event's total charge against the 3-sigma bounds derived
   * from the QA histograms to reject pile-up or background-dominated events.
   * @return True if the event falls within the acceptable charge window.
   */
  bool process_event_check();

  /**
   * @brief Performs the primary tower-by-tower Q-vector calculation and normalization.
   * * Loops through sEPD channels, excludes bad channels, calculates the raw Q-vector
   * for all harmonics, and normalizes the results by the total arm charge.
   * @return True if both South and North arms have non-zero total charge.
   */
  bool process_sEPD();

  /**
   * @brief Calculates and fills profiles for the initial Q-vector averages.
   * @param cent The event centrality.
   * @param q_S The South arm normalized Q-vector.
   * @param q_N The North arm normalized Q-vector.
   * @param h Reference to the cache of profiles for the first pass.
   */
  static void process_averages(double cent, const QVecShared::QVec& q_S, const QVecShared::QVec& q_N, const AverageHists& h);

  /**
   * @brief Applies re-centering offsets and fills profiles for second-moment calculation.
   * @param cent The event centrality.
   * @param h_idx Harmonic index.
   * @param q_S The South arm normalized Q-vector.
   * @param q_N The North arm normalized Q-vector.
   * @param h Reference to the cache of profiles for the second pass.
   */
  void process_recentering(double cent, size_t h_idx, const QVecShared::QVec& q_S, const QVecShared::QVec& q_N, const RecenterHists& h);

  /**
   * @brief Applies the full correction (re-centering + flattening) for validation.
   * @param cent The event centrality.
   * @param h_idx Harmonic index.
   * @param q_S The South arm normalized Q-vector.
   * @param q_N The North arm normalized Q-vector.
   * @param h Reference to the cache of profiles for the third pass.
   */
  void process_flattening(double cent, size_t h_idx, const QVecShared::QVec& q_S, const QVecShared::QVec& q_N, const FlatteningHists& h);

  /**
   * @brief Calculates the 2x2 anisotropy correction (whitening) matrix.
   * * This matrix transforms the elliptical Q-vector distribution into a circularly
   * symmetric (isotropic) distribution. It effectively corrects for detector
   * acceptance effects and gain non-uniformities by normalizing the second-order
   * moments of the Q-vector.
   * * @param xx The <Qx^2> second moment.
   * @param yy The <Qy^2> second moment.
   * @param xy The <QxQy> cross-moment.
   * @param n Harmonic order (used for error logging context).
   * @param cent_bin Centrality bin (used for error logging context).
   * @param det_label Detector label ("S" or "N").
   * @return std::array<std::array<double, 2>, 2> The 2x2 correction matrix.
   */
  std::array<std::array<double, 2>, 2> calculate_flattening_matrix(double xx, double yy, double xy, int n, int cent_bin, const std::string& det_label);

  /**
   * @brief Computes 1st-order re-centering offsets for a specific centrality bin.
   * * Extracts average Q-vector components from histograms and stores them in the
   * correction data matrix for use in subsequent processing passes.
   * * @param cent_bin The index of the centrality bin.
   * @param h_idx The index of the harmonic order in the harmonics array.
   */
  void compute_averages(size_t cent_bin, int h_idx);

  /**
   * @brief Computes re-centering parameters and solves the flattening matrices.
   * * Extracts the re-centered second moments from the profiles and populates the
   * internal CorrectionData matrix with calculated flattening coefficients.
   * @param cent_bin The centrality bin index.
   * @param h_idx The harmonic index.
   */
  void compute_recentering(size_t cent_bin, int h_idx);

  /**
   * @brief Logs the final corrected moments to verify successful flattening.
   * @param cent_bin The centrality bin index.
   * @param n The harmonic order.
   */
  void print_flattening(size_t cent_bin, int n) const;

  void prepare_hists();

  /**
   * @brief Prepares a vector of pointers to histograms used in the first pass.
   * @return A vector of AverageHists structs, indexed by harmonic.
   */
  void prepare_average_hists();

  /**
   * @brief Prepares a vector of pointers to histograms used in the second pass.
   * @return A vector of RecenterHists structs, indexed by harmonic.
   */
  void prepare_recenter_hists();

  /**
   * @brief Prepares a vector of pointers to histograms used in the third pass.
   * @return A vector of FlatteningHists structs, indexed by harmonic.
   */
  void prepare_flattening_hists();

  /**
   * @brief Top-level driver for processing Quality Assurance histograms.
   * * Loads the reference histogram file to identify bad channels and establish
   * event-level charge thresholds as a function of centrality.
   */
  int process_QA_hist();

  /**
   * @brief Identifies and catalogs "Bad" (Hot, Cold, or Dead) sEPD channels.
   * * Uses a reference charge histogram to compute Z-scores based on mean charge
   * per radial bin. Channels exceeding the sigma threshold are added to the internal exclusion set.
   * * @param file Pointer to the open TFile containing QA histograms.
   */
  int process_bad_channels(TFile* file);

  /**
   * @brief Establishes sEPD charge-cut thresholds for event selection.
   * * Uses the 2D total charge vs. centrality distribution to derive mean and
   * sigma values, generating a 1D profile of the selection window.
   * @param file Pointer to the open QA histogram file.
   */
  int process_sEPD_event_thresholds(TFile* file);

  void write_cdb();

  /**
   * @brief Writes the Event Plane calibration constants to a CDB-formatted TTree.
   * * Formats the re-centering and flattening moments into a CDBTTree payload
   * indexed by centrality bin for sPHENIX database integration.
   * * @param output_dir The filesystem directory where the .root payload will be saved.
   */
  void write_cdb_EventPlane();

  /**
   * @brief Writes the Hot/Cold tower status map to a CDB-formatted TTree.
   * * Encodes sEPD channel indices into TowerInfo keys and maps status codes (1=Dead,
   * 2=Hot, 3=Cold) to the final database payload.
   * * @param output_dir The filesystem directory where the .root payload will be saved.
   */
  void write_cdb_BadTowers();
};

#endif  // QVECCALIB_H
