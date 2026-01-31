#ifndef QVECCDB_H
#define QVECCDB_H

#include "QVecDefs.h"

// ====================================================================
// ROOT Includes
// ====================================================================
#include <TFile.h>

// ====================================================================
// Standard C++ Includes
// ====================================================================
#include <memory>
#include <utility>

/**
 * @class QVecCDB
 * @brief Generates sPHENIX Calibration Database (CDB) payloads for sEPD calibrations.
 *
 * QVecCDB is responsible for consolidating the correction parameters derived
 * during the calibration stage into standardized database formats:
 * * - **SEPD_EventPlaneCalib**: Encapsulates re-centering offsets (<Qx>, <Qy>) and
 * flattening moments (<Qx^2>, <Qy^2>, <QxQy>) indexed by centrality bin.
 * - **SEPD_HotMap**: Maps sEPD tower statuses (Dead, Hot, or Cold) to their
 * respective channel IDs for use in reconstruction.
 * * The class interfaces with the `CDBTTree` object to commit these payloads
 * for a specific run number and database tag.
 */
class QVecCDB
{
 public:
  // The constructor takes the configuration
  QVecCDB(std::string input_file, int runnumber, std::string output_dir, std::string cdb_tag)
    : m_input_file(std::move(input_file))
    , m_runnumber(runnumber)
    , m_output_dir(std::move(output_dir))
    , m_cdb_tag(std::move(cdb_tag))
  {
  }

  void run()
  {
    load_data();
    write_cdb();
  }

 private:

  static constexpr size_t m_cent_bins = QVecShared::CENT_BINS;
  static constexpr auto m_harmonics = QVecShared::HARMONICS;

  static constexpr float SIGMA_HOT = 6.0F;
  static constexpr float SIGMA_COLD = -6.0F;

  // Holds all correction data
  // key: [Harmonic][Cent][Subdetector]
  // Harmonics {2,3,4} -> 3 elements
  // Subdetectors {S,N,NS} -> 3 elements
  std::array<std::array<std::array<QVecShared::CorrectionMoments, 3>, m_cent_bins>, m_harmonics.size()> m_correction_data;

  // --- Member Variables ---
  std::string m_input_file;
  int m_runnumber;
  std::string m_output_dir;
  std::string m_cdb_tag;

  std::unique_ptr<TFile> m_tfile;

  // --- Private Helper Methods ---
/**
 * @brief Safely retrieves a ROOT object from the internal TFile and returns a managed unique_ptr.
 * * Utilizes the internal m_tfile member to locate the object. Performs a
 * dynamic_cast for type safety and Clones the object for persistent use.
 * * @tparam T The ROOT class type (e.g., TProfile, TH3).
 * @param name The name of the object within the internal file.
 * @return std::unique_ptr<T> A managed pointer to the cloned object.
 * @throws std::runtime_error If the object is not found or type mismatch occurs.
 */
  template <typename T>
  std::unique_ptr<T> load_and_clone(const std::string& name);

/**
 * @brief Provides safe access to the internal correction data storage.
 * * This accessor handles the mapping between the physics-based Subdetector enum
 * and the zero-based indexing of the underlying multi-dimensional array.
 * * @param h_idx The index of the harmonic order in the m_harmonics array.
 * @param cent_bin The index of the centrality bin.
 * @param sub The subdetector arm (South or North) using the QVecShared enum.
 * @return QVecShared::CorrectionMoments& A reference to the specific data entry.
 */
  QVecShared::CorrectionMoments& getData(size_t h_idx, size_t cent_bin, QVecShared::Subdetector sub);

/**
 * @brief High-level orchestrator for loading calibration input from a ROOT file.
 * * Opens the input file specified in the constructor and validates its integrity
 * before iteratively calling load_correction_data() for every defined harmonic.
 * * Throws a std::runtime_error if the file cannot be opened or is found to be
 * a "zombie" file.
 */
  void load_data();

/**
 * @brief Loads 1st and 2nd order correction parameters from a calibration file.
 * * Iterates through harmonics and centrality bins to populate the internal
 * correction matrix using the centralized QVecShared naming scheme.
 * * @param h_idx The index of the harmonic to load.
 */
  void load_correction_data(size_t h_idx);

/**
 * @brief Top-level orchestrator for the CDB writing phase.
 *
 * This method manages the creation of the run-specific directory structure
 * (e.g., [output_dir]/[runnumber]) within the base output path. Once the
 * filesystem is prepared, it delegates the generation and commitment of
 * specific calibration payloads to write_cdb_EventPlane() and
 * write_cdb_BadTowers().
 */
  void write_cdb();

/**
 * @brief Writes the Event Plane calibration constants to a CDB-formatted TTree.
 * * Formats the re-centering and flattening moments into a CDBTTree payload
 * indexed by centrality bin for sPHENIX database integration.
 * * @param output_dir The filesystem directory where the .root payload will be saved.
 */
  void write_cdb_EventPlane(const std::string &output_dir);

/**
 * @brief Writes the Hot/Cold tower status map to a CDB-formatted TTree.
 * * Encodes sEPD channel indices into TowerInfo keys and maps status codes (1=Dead,
 * 2=Hot, 3=Cold) to the final database payload.
 * * @param output_dir The filesystem directory where the .root payload will be saved.
 */
  void write_cdb_BadTowers(const std::string &output_dir);
};

#endif // QVECCDB_H
