// ----------------------------------------------------------------------------
/*! \file    TriggerClusterMaker.h'
 *  \authors Derek Anderson
 *  \date    05.15.2024
 *
 *  A Fun4All QA module to construct trigger clusters,
 *  jet patches stored in RawCluster objects, for
 *  downstream analysis
 */
// ----------------------------------------------------------------------------

#ifndef TRIGGER_TRIGGERCLUSTERMAKER_H
#define TRIGGER_TRIGGERCLUSTERMAKER_H

// calo base
#include <calobase/RawClusterContainer.h>

// f4a libraries
#include <fun4all/SubsysReco.h>

// root libraries
#include <TFile.h>
#include <TTree.h>

#include <array>
#include <string>
#include <utility>
#include <vector>

// forward declarations
class LL1Out;
class PHCompositeNode;
class RawCluster;
class TriggerPrimitiveContainer;



// ----------------------------------------------------------------------------
//! Options for TriggerClusterMaker module
// ----------------------------------------------------------------------------
struct TriggerClusterMakerConfig {

  // general options
  bool debug = true;

  // output options
  bool        saveToNode  = false;
  bool        saveToFile  = true;
  std::string outNodeName = "TriggerCluster";
  std::string outFileName = "test.root";

  // input nodes
  std::array<std::string, 2> inLL1Nodes = {
    "LL1OUT_RAW_JET",
    "LL1OUT_JET"
  };
  std::array<std::string, 8> inPrimNodes = {
    "TRIGGERPRIMITIVES_RAW_EMCAL",
    "TRIGGERPRIMITIVES_RAW_EMCAL_LL1",
    "TRIGGERPRIMITIVES_RAW_JET",
    "TRIGGERPRIMITIVES_EMCAL",
    "TRIGGERPRIMITIVES_EMCAL_LL1",
    "TRIGGERPRIMITIVES_HCALIN",
    "TRIGGERPRIMITIVES_HCALOUT",
    "TRIGGERPRIMITIVES_HCAL_LL1"
  };

};



// ----------------------------------------------------------------------------
//! Makes Trigger Cluster
// ----------------------------------------------------------------------------
/*! This Fun4all modules ingests calorimeter triggers and
 *  trigger primitives to turn them into RawCluster objects,
 *  i.e. "Trigger Clusters", for downstream analysis. Output
 *  clusters can be placed on the node tree, or saved to
 *  a TTree in a specified output file.
 */
class TriggerClusterMaker : public SubsysReco {

  public:

    // ctor
    TriggerClusterMaker(const std::string& name = "TriggerClusterMaker");
    ~TriggerClusterMaker() override;

    // setters
    void SetConfig(const TriggerClusterMakerConfig& config) {m_config = config;}

    // getters
    TriggerClusterMakerConfig GetConfig() {return m_config;}

    // f4a methods
    int Init(PHCompositeNode* topNode)          override;
    int process_event(PHCompositeNode* topNode) override;
    int End(PHCompositeNode* topNode)           override;

  private:

    // LL1 accessors
    enum LL1 {
      Raw,
      Jet
    };

    // trigger primitive accessors
    enum Primitive {
      RawEM,
      RawEMLL1,
      RawTrg,
      EMCal,
      EMCalLL1,
      IHCal,
      OHCal,
      HCalLL1
    };

    // private methods
    void InitOutFile();
    void InitOutTree();
    void InitOutNode(PHCompositeNode* topNode);
    void GrabNodes(PHCompositeNode* topNode);

    // f4a members
    std::vector<LL1Out*>                    m_inLL1Nodes;
    std::vector<TriggerPrimitiveContainer*> m_inPrimNodes;

    // output members
    std::unique_ptr<TFile>               m_outFile       = NULL;
    std::unique_ptr<TTree>               m_outTree       = NULL;
/* TODO implement
    std::unique_ptr<RawClusterContainer> m_outClustStore = NULL;
*/

    // module configuration
    TriggerClusterMakerConfig m_config;

};

#endif

// end ------------------------------------------------------------------------
