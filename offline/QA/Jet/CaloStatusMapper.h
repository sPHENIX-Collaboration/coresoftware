/// ===========================================================================
/*! \file   CaloStatusMapper.h
 *  \author Derek Anderson
 *  \date   05.22.2024
 *
 *  A Fun4All QA module to plot no. of towers per event
 *  and vs. eta, phi as a function of status.
 */
/// ===========================================================================

#ifndef CLUSTERSTATUSMAPPER_H
#define CLUSTERSTATUSMAPPER_H

// jet qa
// FIXME change to local include when ready to merge
#include <jetqa/JetQADefs.h>

// module definitions
#include "CaloStatusMapperDefs.h"

// calo base
#include <calobase/TowerInfoContainerv2.h>

// f4a libraries
#include <fun4all/SubsysReco.h>

// c++ utilities
#include <map>
#include <string>
#include <vector>

// forward declarations
class PHCompositeNode;
class Fun4AllHistoManager;
class TH1;
class TriggerAnalyzer;



// ============================================================================
//! Map Status of Calo Towers
// ============================================================================
/*! This Fun4all modules ingests calorimeter towers, and then plots
 *  the no. of towers per event and vs. eta/phi for a given status.
 */
class CaloStatusMapper : public SubsysReco
{

  public:

    // ========================================================================
    //! Options for CaloStatusMapper module
    // ========================================================================
    struct Config
    {

      ///! turn debug messages on/off
      bool debug {true};

      ///! module name
      std::string moduleName {"CaloStatusMapper"};

      ///! histogram tag
      std::string histTag {""};

      ///! input nodes and what type of calo they are
      std::vector<CaloStatusMapperDefs::NodeDef> inNodeNames
      {
        {"TOWERINFO_CALIB_CEMC",    CaloStatusMapperDefs::Calo::EMCal},
        {"TOWERINFO_CALIB_HCALIN",  CaloStatusMapperDefs::Calo::HCal},
        {"TOWERINFO_CALIB_HCALOUT", CaloStatusMapperDefs::Calo::HCal}
      };

     ///! turn trigger selection on/off
     bool doTrgSelect {false};

     ///! trigger to select
     uint32_t trgToSelect {JetQADefs::GL1::MBDNSJet1};

    };  // end Config

    // ctor/dtor
    CaloStatusMapper(const std::string& modulename = "CaloStatusMapper", const bool debug = false);
    CaloStatusMapper(const Config& config);
    ~CaloStatusMapper() override;

    // setters
    void SetConfig(const Config& config) {m_config = config;}

    // getters
    Config GetConfig() {return m_config;}

    // f4a methods
    int Init(PHCompositeNode* topNode) override;
    int process_event(PHCompositeNode* topNode) override;
    int End(PHCompositeNode* topNode) override;

  private:

    // private methods
    void InitHistManager();
    void BuildHistograms();
    void GrabNodes(PHCompositeNode* topNode);
    std::string MakeBaseName(const std::string& base, const std::string& node, const std::string& stat = "") const;

    ///! module configuration
    Config m_config;

    ///! histogram manager
    Fun4AllHistoManager* m_manager {nullptr};

    ///! for checking which trigger fired
    TriggerAnalyzer* m_analyzer {nullptr};

    ///! status labels
    std::map<CaloStatusMapperDefs::Stat, std::string> m_mapStatLabels {CaloStatusMapperDefs::StatLabels()};

    ///! output histograms
    std::map<std::string, TH1*> m_hists;

    ///! input nodes
    std::vector<TowerInfoContainer*> m_inNodes;

    ///! no. of events processed
    uint64_t m_nEvent {0};

};  // end CaloStatusMapper

#endif

// end ========================================================================
