// ----------------------------------------------------------------------------
/*! \file    TriggerClusterMaker.cc'
 *  \authors Derek Anderson
 *  \date    05.15.2024
 *
 *  A Fun4All QA module to construct trigger clusters,
 *  jet patches stored in RawCluster objects, for
 *  downstream analysis
 */
// ----------------------------------------------------------------------------

// module definition
#include "TriggerClusterMaker.h"

// calo base
#include <calobase/RawCluster.h>

// trigger includes
#include <calotrigger/LL1Out.h>
#include <calotrigger/LL1Outv1.h>
#include <calotrigger/TriggerPrimitiveContainer.h>
#include <calotrigger/TriggerPrimitiveContainerv1.h>

// f4a includes
#include <fun4all/Fun4AllReturnCodes.h>

// phool includes
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>

// c++ includes
#include <algorithm>
#include <cassert>
#include <iostream>



// ctor/dtor ==================================================================

// ----------------------------------------------------------------------------
//! Module constructor
// ----------------------------------------------------------------------------
TriggerClusterMaker::TriggerClusterMaker(const std::string &name) : SubsysReco(name) {

  // print debug message
  if (m_config.debug && (Verbosity() > 0)) {
    std::cout << "TriggerClusterMaker::TriggerClusterMaker(const std::string &name) Calling ctor" << std::endl;
  }

  // initialize input nodes to null
  std::fill(m_inLL1Nodes.begin(), m_inLL1Nodes.end(), nullptr);
  std::fill(m_inPrimNodes.begin(), m_inPrimNodes.end(), nullptr);

}  // end ctor



// ----------------------------------------------------------------------------
//! Module destructor
// ----------------------------------------------------------------------------
TriggerClusterMaker::~TriggerClusterMaker() {

  // print debug message
  if (m_config.debug && (Verbosity() > 0)) {
    std::cout << "TriggerClusterMaker::~TriggerClusterMaker() Calling dtor" << std::endl;
  }

  /* nothing to do */

}  // end dtor



// fun4all methods ============================================================

// ----------------------------------------------------------------------------
//! Initialize module
// ----------------------------------------------------------------------------
int TriggerClusterMaker::Init(PHCompositeNode* topNode) {

  if (m_config.debug) {
    std::cout << "TriggerClusterMaker::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  }

  // initialize outputs
  if (m_config.saveToNode) {
    InitOutNode(topNode);
  }
  if (m_config.saveToFile) {
    InitOutFile();
    InitOutTree();
  }
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'Init(PHCompositeNode*)'



// ----------------------------------------------------------------------------
//! Grab inputs and build trigger clusters
// ----------------------------------------------------------------------------
int TriggerClusterMaker::process_event(PHCompositeNode* topNode) {

  if (m_config.debug) {
    std::cout << "TriggerClusterMaker::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  }

  // grab input nodes
  GrabNodes(topNode);

  /* TODO fill in */

  return Fun4AllReturnCodes::EVENT_OK;


}  // end 'process_event(PHCompositeNode*)'



// ----------------------------------------------------------------------------
//! Run final calculations
// ----------------------------------------------------------------------------
int TriggerClusterMaker::End(PHCompositeNode * /*topNode*/) {

  if (m_config.debug) {
    std::cout << "TriggerClusterMaker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }

  /* TODO fill in */

  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'End(PHCompositeNode*)'



// private methods ============================================================

// ----------------------------------------------------------------------------
//! Initialize output file
// ----------------------------------------------------------------------------
void TriggerClusterMaker::InitOutFile() {

  // print debug message
  if (m_config.debug && (Verbosity() > 0)) {
    std::cout << "TriggerClusterMaker::InitOutFile() Creating output file" << std::endl;
  }

  m_outFile = std::make_unique<TFile>(m_config.outFileName.data(), "recreate");
  if (!m_outFile) {
    std::cerr << PHWHERE << ": PANIC! Couldn't output file!" << std::endl;
    assert(m_outFile);
  }
  return;

}  // end 'InitOutFile()'



// ----------------------------------------------------------------------------
//! Initialize output TTree
// ----------------------------------------------------------------------------
void TriggerClusterMaker::InitOutTree() {

  // print debug message
  if (m_config.debug && (Verbosity() > 0)) {
    std::cout << "TriggerClusterMaker::InitOutTree() Creating output tree" << std::endl;
  }

  /* TODO fill in */
  return;

}  // end 'InitOutTree()'



// ----------------------------------------------------------------------------
//! Create output node on node tree
// ----------------------------------------------------------------------------
void TriggerClusterMaker::InitOutNode(PHCompositeNode*  /*topNode*/) {

  // print debug message
  if (m_config.debug && (Verbosity() > 0)) {
    std::cout << "TriggerClusterMaker::InitOutNode(PHCompositeNode*) Creating output node" << std::endl;
  }

  /* TODO fill in */
  return;

}  // end 'InitOutNode(PHCompositeNode*)'



// ----------------------------------------------------------------------------
//! Grab input nodes
// ----------------------------------------------------------------------------
void TriggerClusterMaker::GrabNodes(PHCompositeNode* topNode) {

  // print debug message
  if (m_config.debug && (Verbosity() > 0)) {
    std::cout << "TriggerClusterMaker::GrabNodes(PHCompositeNode*) Grabbing input nodes" << std::endl;
  }

  // get LL1 nodes
  for (const std::string& inLL1Node : m_config.inLL1Nodes) {
    m_inLL1Nodes.push_back(
      findNode::getClass<LL1Out>(topNode, inLL1Node)
    );
    if (!m_inLL1Nodes.back()) {
      std::cerr << PHWHERE << ": PANIC! Couldn't grab LL1Out node '" << inLL1Node << "'!" << std::endl;
      assert(m_inLL1Nodes.back());
    }
  }

  // get trigger primitive nodes
  for (const std::string& inPrimNode : m_config.inPrimNodes) {
    m_inPrimNodes.push_back(
      findNode::getClass<TriggerPrimitiveContainer>(topNode, inPrimNode)
    );
    if (!m_inPrimNodes.back()) {
      std::cerr << PHWHERE << ": PANIC! Couldn't grab TriggerPrimitive node '" << inPrimNode << "'!" << std::endl;
      assert(m_inPrimNodes.back());
    }
  }
  return;

}  // end 'GrabNodes(PHCompositeNode*)'

// end ------------------------------------------------------------------------
