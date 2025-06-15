#include "LaserEventRejecter.h"

#include "LaserEventInfo.h"

#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

LaserEventRejecter::LaserEventRejecter(const std::string &name)
  : SubsysReco(name)
{
}

int LaserEventRejecter::process_event(PHCompositeNode *topNode)
{
  m_laserEventInfo = findNode::getClass<LaserEventInfo>(topNode, "LaserEventInfo");
  if (!m_laserEventInfo)
  {
    std::cout << PHWHERE << "ERROR: Can't find node LaserEventInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if(m_laserEventInfo->isGl1LaserPileupEvent() || m_laserEventInfo->isLaserEvent())
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}
