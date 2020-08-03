


#include "PHActsVertexFitter.h"


#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"

#include <iostream>

PHActsVertexFitter::PHActsVertexFitter()
{

}


int PHActsVertexFitter::Init(PHCompositeNode *topNode)
{
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTRUN;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFitter::End(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFitter::process_event(PHCompositeNode *topNode)
{
  

  return Fun4AllReturnCodes::EVENT_OK;
}


int PHActsVertexFitter::getNodes(PHCompositeNode *topNode)
{

  m_svtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  
  if(!m_svtxTrackMap)
    {
      std::cout << PHWHERE << "Unable to find SvtxTrackMap. Exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
      
    }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "ActsTrackingGeometry not on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;

}
