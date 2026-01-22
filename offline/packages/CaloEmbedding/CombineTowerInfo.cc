#include "CombineTowerInfo.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <iostream>
#include <stdexcept>

//____________________________________________________________________________
CombineTowerInfo::CombineTowerInfo(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________
int CombineTowerInfo::InitRun(PHCompositeNode* topNode)
{
  if (m_inputNodeA.empty() || m_inputNodeB.empty() || m_outputNode.empty())
  {
    throw std::runtime_error("CombineTowerInfo: input/output node names not set");
  }

  CreateNodes(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
void CombineTowerInfo::CreateNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode =
      dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    throw std::runtime_error("CombineTowerInfo: DST node not found");
  }

  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", m_detector));

  m_towersA = findNode::getClass<TowerInfoContainer>(topNode, m_inputNodeA);
  m_towersB = findNode::getClass<TowerInfoContainer>(topNode, m_inputNodeB);

  if (!m_towersB)
  {
    std::cout << "CombineTowerInfo: " <<m_inputNodeB << " not found" << std::endl;
    throw std::runtime_error("CombineTowerInfo: input TowerInfoContainer missing");
  }
  if (!m_towersA)
  {
    std::cout << "CombineTowerInfo: " <<m_inputNodeA << " not found" << std::endl;
    throw std::runtime_error("CombineTowerInfo: input TowerInfoContainer missing");
  }

  m_towersOut = findNode::getClass<TowerInfoContainer>(dstNode, m_outputNode);
  if (!m_towersOut)
  {
    m_towersOut =
        dynamic_cast<TowerInfoContainer*>(m_towersA->CloneMe());

    auto* node = new PHIODataNode<PHObject>(
        m_towersOut, m_outputNode, "PHObject");

    DetNode->addNode(node);
  }

  if (m_towersA->size() != m_towersB->size())
  {
    throw std::runtime_error("CombineTowerInfo: input containers have different sizes");
  }
}

//____________________________________________________________________________
int CombineTowerInfo::process_event(PHCompositeNode* /*topNode*/)
{
  const unsigned int ntowers = m_towersA->size();

  for (unsigned int ich = 0; ich < ntowers; ++ich)
  {
    TowerInfo* towerA = m_towersA->get_tower_at_channel(ich);
    TowerInfo* towerB = m_towersB->get_tower_at_channel(ich);
    TowerInfo* towerO = m_towersOut->get_tower_at_channel(ich);

    towerO->copy_tower(towerA);

    const float e_sum = towerA->get_energy() + towerB->get_energy();
    towerO->set_energy(e_sum);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

