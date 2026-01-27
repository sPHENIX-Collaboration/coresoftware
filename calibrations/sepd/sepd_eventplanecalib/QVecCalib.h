#ifndef QVECCALIB_H
#define QVECCALIB_H

#include "QVecDefs.h"

// ====================================================================
// ROOT Includes
// ====================================================================

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

// ====================================================================
// Standard C++ Includes
// ====================================================================
#include <map>
#include <memory>
#include <unordered_set>

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
class QVecCalib
{
 public:
  // The constructor takes the configuration
  QVecCalib(std::string input_file, std::string input_hist, std::string input_Q_calib, int pass, long long events, std::string output_dir)
    : m_input_file(std::move(input_file))
    , m_input_hist(std::move(input_hist))
    , m_input_Q_calib(std::move(input_Q_calib))
    , m_pass(static_cast<Pass>(pass))
    , m_events_to_process(events)
    , m_output_dir(std::move(output_dir))
  {
  }

  void run()
  {
    setup_chain();
    process_QA_hist();
    init_hists();
    process_events();
    save_results();
  }

  enum class Pass
  {
    ComputeRecentering,
    ApplyRecentering,
    ApplyFlattening
  };

 private:

  struct CorrectionData
  {
    QVecShared::QVec avg_Q{};
    std::array<std::array<double, 2>, 2> X_matrix{};
  };

  static constexpr size_t m_cent_bins = QVecShared::CENT_BINS;
  static constexpr auto m_harmonics = QVecShared::HARMONICS;
  double m_cent_low = -0.5;
  double m_cent_high = 79.5;

  // Holds all correction data
  // key: [Cent][Harmonic][Subdetector]
  // Harmonics {2,3,4} -> 3 elements
  // Subdetectors {S,N,NS} -> 3 elements
  std::array<std::array<std::array<CorrectionData, 3>, m_harmonics.size()>, m_cent_bins> m_correction_data;

  static constexpr size_t IDX_NS = 2;

  // Store harmonic orders and subdetectors for easy iteration
  static constexpr std::array<QVecShared::Subdetector, 2> m_subdetectors = {QVecShared::Subdetector::S, QVecShared::Subdetector::N};
  static constexpr std::array<QVecShared::QComponent, 2> m_components = {QVecShared::QComponent::X, QVecShared::QComponent::Y};

  struct EventData
  {
    int event_id{0}; // NOLINT(misc-non-private-member-variables-in-classes)
    double event_zvertex{0.0}; // NOLINT(misc-non-private-member-variables-in-classes)
    double event_centrality{0.0}; // NOLINT(misc-non-private-member-variables-in-classes)
    double sepd_totalcharge{0.0}; // NOLINT(misc-non-private-member-variables-in-classes)

    std::array<std::array<QVecShared::QVec, 2>, m_harmonics.size()> q_vectors; // NOLINT(misc-non-private-member-variables-in-classes)

    void reset()
    {
      for (auto& q_vec_harmonic : q_vectors)
      {
        for (auto& q_vec : q_vec_harmonic)
        {
          q_vec.x = 0.0;
          q_vec.y = 0.0;
        }
      }
    }

    std::vector<int>* sepd_channel{nullptr}; // NOLINT(misc-non-private-member-variables-in-classes)
    std::vector<double>* sepd_charge{nullptr}; // NOLINT(misc-non-private-member-variables-in-classes)
    std::vector<double>* sepd_phi{nullptr}; // NOLINT(misc-non-private-member-variables-in-classes)
  };

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

  // --- Member Variables ---
  EventData m_event_data;
  std::unique_ptr<TChain> m_chain;

  // Configuration stored as members
  std::string m_input_file;
  std::string m_input_hist;
  std::string m_input_Q_calib;
  Pass m_pass{0};
  long long m_events_to_process;
  std::string m_output_dir;

  // Hists
  std::map<std::string, std::unique_ptr<TH1>> m_hists1D;
  std::map<std::string, std::unique_ptr<TH2>> m_hists2D;
  std::map<std::string, std::unique_ptr<TProfile>> m_profiles;

  // sEPD Bad Channels
  std::unordered_set<int> m_bad_channels;

  double m_sEPD_min_avg_charge_threshold{1};
  double m_sEPD_sigma_threshold{3};

  // --- Private Helper Methods ---

/**
 * @brief Sets up a TChain and performs structural validation of the input ROOT file.
 * * Verifies file existence, ensures the requested TTree exists, and checks for
 * non-zero entries before returning a configured chain.
 * * @param input_filepath Path to the input .root file.
 * @param tree_name_in_file Name of the TTree inside the file.
 * @return std::unique_ptr<TChain> A configured TChain, or nullptr if validation fails.
 */
  std::unique_ptr<TChain> setupTChain(const std::string& input_filepath, const std::string& tree_name_in_file);

/**
 * @brief Orchestrates the TChain initialization and branch configuration.
 * * Sets the branch statuses and addresses for event-level data (ID, centrality, charge)
 * and sEPD tower-level data (channel, charge, phi) needed for the calibration.
 */
  void setup_chain();

/**
 * @brief Initializes all output histograms and profiles.
 * * Dynamically generates histogram names using the shared naming helper based on
 * the current calibration pass (e.g., adding "_corr" or "_corr2" suffixes).
 */
  void init_hists();

/**
 * @brief Safely retrieves a ROOT object from a file and returns a managed unique_ptr.
 * * Performs a dynamic_cast to verify the requested type T and Clones the object
 * to ensure it remains valid after the source file is closed.
 * * @tparam T The ROOT class type (e.g., TProfile).
 * @param file Pointer to the source TFile.
 * @param name The name of the object within the file.
 * @return std::unique_ptr<T> A managed pointer to the cloned object.
 * @throws std::runtime_error If the object is not found or type mismatch occurs.
 */
  template <typename T>
  std::unique_ptr<T> load_and_clone(TFile* file, const std::string& name);

/**
 * @brief Loads the results of previous passes from a calibration file.
 * * Populates the internal correction data structure with averages and/or
 * matrices required for the current processing pass.
 */
  void load_correction_data();

/**
 * @brief High-level orchestrator for the event processing phase.
 * * If the current pass requires existing calibration data (Recentering or Flattening),
 * it triggers the data loading sequence before starting the main event loop.
 */
  void process_events();

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
 * @brief Primary event loop orchestrator.
 * * Iterates through the TChain entries, performs event selection, executes
 * normalization/re-centering/flattening logic based on the current pass,
 * and fills the output histograms.
 */
  void run_event_loop();

/**
 * @brief Finalizes the analysis by writing all histograms to the output ROOT file.
 * * Creates the output directory if it does not exist and ensures all 1D, 2D
 * and TProfiles are safely persisted to disk.
 */
  void save_results() const;

/**
 * @brief Calculates and fills profiles for the initial Q-vector averages.
 * @param cent The event centrality.
 * @param q_S The South arm normalized Q-vector.
 * @param q_N The North arm normalized Q-vector.
 * @param h Reference to the cache of profiles for the first pass.
 */
  static void process_averages(double cent, QVecShared::QVec q_S, QVecShared::QVec q_N, const AverageHists& h);

/**
 * @brief Applies re-centering offsets and fills profiles for second-moment calculation.
 * @param cent The event centrality.
 * @param h_idx Harmonic index.
 * @param q_S The South arm normalized Q-vector.
 * @param q_N The North arm normalized Q-vector.
 * @param h Reference to the cache of profiles for the second pass.
 */
  void process_recentering(double cent, size_t h_idx, QVecShared::QVec q_S, QVecShared::QVec q_N, const RecenterHists& h);

/**
 * @brief Applies the full correction (re-centering + flattening) for validation.
 * @param cent The event centrality.
 * @param h_idx Harmonic index.
 * @param q_S The South arm normalized Q-vector.
 * @param q_N The North arm normalized Q-vector.
 * @param h Reference to the cache of profiles for the third pass.
 */
  void process_flattening(double cent, size_t h_idx, QVecShared::QVec q_S, QVecShared::QVec q_N, const FlatteningHists& h);

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

/**
 * @brief Prepares a vector of pointers to histograms used in the first pass.
 * @return A vector of AverageHists structs, indexed by harmonic.
 */
  std::vector<AverageHists> prepare_average_hists();

/**
 * @brief Prepares a vector of pointers to histograms used in the second pass.
 * @return A vector of RecenterHists structs, indexed by harmonic.
 */
  std::vector<RecenterHists> prepare_recenter_hists();

/**
 * @brief Prepares a vector of pointers to histograms used in the third pass.
 * @return A vector of FlatteningHists structs, indexed by harmonic.
 */
  std::vector<FlatteningHists> prepare_flattening_hists();

/**
 * @brief Top-level driver for processing Quality Assurance histograms.
 * * Loads the reference histogram file to identify bad channels and establish
 * event-level charge thresholds as a function of centrality.
 */
  void process_QA_hist();

/**
 * @brief Identifies and catalogs "Bad" (Hot, Cold, or Dead) sEPD channels.
 * * Uses a reference charge histogram to compute Z-scores based on mean charge
 * per radial bin. Channels exceeding the sigma threshold are added to the internal exclusion set.
 * * @param file Pointer to the open TFile containing QA histograms.
 */
  void process_bad_channels(TFile* file);

/**
 * @brief Establishes sEPD charge-cut thresholds for event selection.
 * * Uses the 2D total charge vs. centrality distribution to derive mean and
 * sigma values, generating a 1D profile of the selection window.
 * @param file Pointer to the open QA histogram file.
 */
  void process_sEPD_event_thresholds(TFile* file);
};

#endif // QVECCALIB_H
