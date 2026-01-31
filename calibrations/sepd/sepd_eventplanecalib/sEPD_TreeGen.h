#ifndef SEPD_TREEGEN_H
#define SEPD_TREEGEN_H

// -- sPHENIX
#include <fun4all/SubsysReco.h>

// -- c++
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

// -- ROOT
#include <TFile.h>
#include <TProfile.h>
#include <TH2.h>
#include <TTree.h>

class PHCompositeNode;

/**
 * @class sEPD_TreeGen
 * @brief SubsysReco module to produce flat TTrees and QA histograms for sEPD calibration.
 *
 * This module extracts event-level info (vertex, centrality) and sEPD tower-level info
 * (charge, phi, channel ID), applying basic event selections (Minimum Bias, Z-vertex)
 * and tower-level cuts (charge threshold, zero-suppression).
 */
class sEPD_TreeGen : public SubsysReco
{
 public:
 /**
   * @brief Constructor for sEPD_TreeGen.
   * @param name The name assigned to this SubsysReco module.
   */
  explicit sEPD_TreeGen(const std::string &name = "sEPD_TreeGen");

 /**
   * @brief Initializes the module and creates the output TTree and QA histograms.
   * @param topNode Pointer to the node tree.
   * @return Fun4All return code.
   */
  int Init(PHCompositeNode *topNode) override;

 /**
   * @brief Main event-by-event processing method.
   * @details Orchestrates event checks, centrality retrieval, sEPD tower processing,
   * and fills the TTree for valid events.
   * @param topNode Pointer to the node tree.
   * @return Fun4All return code.
   */
  int process_event(PHCompositeNode *topNode) override;

 /**
   * @brief Resets event-level data structures before the next event.
   * @param topNode Pointer to the node tree.
   * @return Fun4All return code.
   */
  int ResetEvent(PHCompositeNode *topNode) override;

 /**
   * @brief Finalizes the module, writing all histograms and the TTree to disk.
   * @param topNode Pointer to the node tree.
   * @return Fun4All return code.
   */
  int End(PHCompositeNode *topNode) override;

 /**
   * @brief Sets the filename for the QA histograms ROOT file.
   * @param file Output path for histograms.
   */
  void set_filename(const std::string &file)
  {
    m_outfile_name = file;
  }

 /**
   * @brief Sets the filename for the flat TTree ROOT file.
   * @param file Output path for the TTree.
   */
  void set_tree_filename(const std::string &file)
  {
    m_outtree_name = file;
  }

 /**
   * @brief Sets the maximum allowed Z-vertex position for event selection.
   * @param zvtx_max Maximum vertex Z in cm.
   */
  void set_zvtx_max(double zvtx_max)
  {
    m_cuts.m_zvtx_max = zvtx_max;
  }

 /**
   * @brief Sets the minimum charge threshold for individual sEPD towers.
   * @param charge_min Minimum charge to include a tower in the TTree.
   */
  void set_sepd_charge_threshold(double charge_min)
  {
    m_cuts.m_sepd_charge_min = charge_min;
  }

 /**
   * @brief Sets the maximum centrality centile allowed for processing.
   * @param cent_max Maximum centile (e.g., 80 for 0-80%).
   */
  void set_cent_max(double cent_max)
  {
    m_cuts.m_cent_max = cent_max;
  }

 private:
 /**
   * @brief Validates event-level conditions (GlobalVertex, Minimum Bias).
   * @param topNode Pointer to the node tree.
   * @return Fun4All return code.
   */
  int process_event_check(PHCompositeNode *topNode);

 /**
   * @brief Processes individual sEPD towers and calculates total charge.
   * @details Applies tower cuts, fills QA histograms, and stores tower data in vectors.
   * @param topNode Pointer to the node tree.
   * @return Fun4All return code.
   */
  int process_sEPD(PHCompositeNode *topNode);

 /**
   * @brief Retrieves and validates centrality information.
   * @param topNode Pointer to the node tree.
   * @return Fun4All return code.
   */
  int process_centrality(PHCompositeNode *topNode);

  int m_event{0};

  std::string m_outfile_name{"test.root"};
  std::string m_outtree_name{"tree.root"};

  static constexpr int PROGRESS_PRINT_INTERVAL = 20;

  // Cuts
  struct Cuts
  {
    double m_zvtx_max{10}; /*cm*/
    double m_sepd_charge_min{0.2};
    double m_cent_max{80};
  };

  Cuts m_cuts;

  struct EventData
  {
    int event_id{0};
    double event_zvertex{9999};
    double event_centrality{9999};

    double sepd_totalcharge{-9999};

    std::vector<int> sepd_channel;
    std::vector<double> sepd_charge;
    std::vector<double> sepd_phi;
  };

  EventData m_data;

  std::unique_ptr<TFile> m_output;
  std::unique_ptr<TTree> m_tree;

  std::unique_ptr<TProfile> hSEPD_Charge;
  std::unique_ptr<TH2> h2SEPD_totalcharge_centrality;
};


#endif // SEPD_TREEGEN_H
