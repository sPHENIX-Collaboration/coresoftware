#include "ClusterSmearingTool.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree

#include <ffamodules/CDBInterface.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <gsl/gsl_randist.h>

#include <TSystem.h>

#include <algorithm>  // for max
#include <cmath>
#include <exception>  // for exception
#include <iostream>   // for operator<<, basic_ostream
#include <numbers>
#include <stdexcept>  // for runtime_error

//____________________________________________________________________________..
ClusterSmearingTool::ClusterSmearingTool(const std::string &name)
  : SubsysReco(name)
{
  if (Verbosity() > 0)
  {
    std::cout << "ClusterSmearingTool::ClusterSmearingTool(const std::string &name) Calling ctor" << std::endl;
  }
}

//____________________________________________________________________________..
int ClusterSmearingTool::InitRun(PHCompositeNode *topNode)
{
  if (!m_overrideSeed)
  {
    m_seed = PHRandomSeed();
  }
  m_rng.reset(gsl_rng_alloc(gsl_rng_mt19937));
  gsl_rng_set(m_rng.get(), m_seed);

  try
  {
    LoadCalib();
    CreateNodeTree(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::InitRun: seed " << m_seed
              << ", input node " << m_inputNodeName
              << ", output node " << (m_writeNewNode ? m_outputNodeName : m_inputNodeName + " (in place)")
              << std::endl;
  }
  if (Verbosity() > 3)
  {
    topNode->print();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void ClusterSmearingTool::LoadCalib()
{
  // reset per run, a previous run may have had no calibration
  m_doSmearing = true;

  std::string calibdir;
  if (m_giveDirectURL)
  {
    calibdir = m_directURL;
    std::cout << "ClusterSmearingTool::InitRun: Using setted url " << calibdir << std::endl;
  }
  else
  {
    calibdir = CDBInterface::instance()->getUrl(m_calibName);
  }

  if (calibdir.empty())
  {
    if (m_doAbortMissingCalib)
    {
      std::cout << "ClusterSmearingTool::InitRun: No smearing parameters found for " << m_calibName << " and abort mode is set. Exiting." << std::endl;
      gSystem->Exit(1);
    }
    std::cout << "ClusterSmearingTool::InitRun: No smearing parameters found for " << m_calibName << ", clusters will NOT be smeared" << std::endl;
    m_doSmearing = false;
    return;
  }

  CDBTTree cdbttree(calibdir);
  m_parA = cdbttree.GetSingleFloatValue("smearE_A");
  m_parB = cdbttree.GetSingleFloatValue("smearE_B");
  m_parC = cdbttree.GetSingleFloatValue("smearE_C");
  m_parD = cdbttree.GetSingleFloatValue("smearPos_D");
  m_parF = cdbttree.GetSingleFloatValue("smearPos_F");

  if (!std::isfinite(m_parA) || !std::isfinite(m_parB) || !std::isfinite(m_parC) ||
      !std::isfinite(m_parD) || !std::isfinite(m_parF))
  {
    std::cout << "ClusterSmearingTool::InitRun: missing smearing parameter in " << calibdir << std::endl;
    throw std::runtime_error("Failed to read smearing parameters in ClusterSmearingTool::LoadCalib");
  }

  std::cout << "ClusterSmearingTool::InitRun: sigma_E/E = " << m_parA << "/E (+) " << m_parB << "/sqrt(E) (+) " << m_parC
            << ", sigma_pos = " << m_parD << "/sqrt(E) (+) " << m_parF << std::endl;
}

//____________________________________________________________________________..
int ClusterSmearingTool::process_event(PHCompositeNode *topNode)
{
  RawClusterContainer *_input_clusters = findNode::getClass<RawClusterContainer>(topNode, m_inputNodeName);
  if (!_input_clusters)
  {
    std::cout << Name() << "::process_event: " << m_inputNodeName << " node missing, abort event" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  RawClusterContainer *_smear_clusters = _input_clusters;
  if (m_writeNewNode)
  {
    _smear_clusters = findNode::getClass<RawClusterContainer>(topNode, m_outputNodeName);
    if (!_smear_clusters)
    {
      std::cout << Name() << "::process_event: " << m_outputNodeName << " node missing, abort run" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    _smear_clusters->Reset();
    // keep the cluster keys identical to the input node
    for (const auto &[key, cluster] : _input_clusters->getClustersMap())
    {
      RawCluster *clone = dynamic_cast<RawCluster *>(cluster->CloneMe());
      if (!clone)
      {
        std::cout << Name() << "::process_event: could not clone cluster " << key << std::endl;
        continue;
      }
      _smear_clusters->getClustersMap()[key] = clone;
    }
  }

  if (!m_doSmearing)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  float vtx_z = 0;
  if (m_useVertex)
  {
    GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (vertexmap && !vertexmap->empty())
    {
      GlobalVertex *vtx = vertexmap->begin()->second;
      if (vtx)
      {
        vtx_z = vtx->get_z();
      }
    }
  }

  for (const auto &[key, cluster] : _smear_clusters->getClustersMap())
  {
    SmearCluster(cluster, vtx_z);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void ClusterSmearingTool::SmearCluster(RawCluster *cluster, float vtx_z)
{
  const float E = cluster->get_energy();
  if (!(E > 0))
  {
    return;
  }
  const float sqrtE = std::sqrt(E);

  if (m_doEnergySmearing)
  {
    const float termA = m_parA / E;
    const float termB = m_parB / sqrtE;
    const float sigmaE = std::sqrt(termA * termA + termB * termB + m_parC * m_parC);
    const float Efac = std::max(0.0F, 1.0F + static_cast<float>(gsl_ran_gaussian(m_rng.get(), sigmaE)));
    cluster->set_energy(E * Efac);
    cluster->set_ecore(cluster->get_ecore() * Efac);
  }

  if (m_doPositionSmearing)
  {
    const float termD = m_parD / sqrtE;
    const float sigmaPos = std::sqrt(termD * termD + m_parF * m_parF);
    const float deta = static_cast<float>(gsl_ran_gaussian(m_rng.get(), sigmaPos));
    const float dphi = static_cast<float>(gsl_ran_gaussian(m_rng.get(), sigmaPos));

    // eta is smeared relative to the vertex, the cluster radius is kept fixed
    const float r = cluster->get_r();
    if (r > 0)
    {
      const float eta = std::asinh((cluster->get_z() - vtx_z) / r);
      cluster->set_z(vtx_z + r * std::sinh(eta + deta));
    }

    float phi = cluster->get_phi() + dphi;
    while (phi > std::numbers::pi_v<float>)
    {
      phi -= 2 * std::numbers::pi_v<float>;
    }
    while (phi <= -std::numbers::pi_v<float>)
    {
      phi += 2 * std::numbers::pi_v<float>;
    }
    cluster->set_phi(phi);
  }
}

void ClusterSmearingTool::CreateNodeTree(PHCompositeNode *topNode)
{
  if (!m_writeNewNode)
  {
    return;
  }

  if (m_outputNodeName == m_inputNodeName)
  {
    std::cout << Name() << "::CreateNodeTree: output node name " << m_outputNodeName
              << " is the same as the input node name, smearing in place" << std::endl;
    m_writeNewNode = false;
    return;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in ClusterSmearingTool::CreateNodeTree");
  }

  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", m_detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_detector);
    dstNode->addNode(DetNode);
  }

  RawClusterContainer *_smear_clusters = findNode::getClass<RawClusterContainer>(dstNode, m_outputNodeName);
  if (!_smear_clusters)
  {
    _smear_clusters = new RawClusterContainer();
    PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_smear_clusters, m_outputNodeName, "PHObject");
    DetNode->addNode(clusterNode);
  }
}
