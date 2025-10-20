/*
 * \file TpcGlobalPositionWrapper.cc
 * \brief provides the tpc 3d global position with all distortion and crossing corrections applied
 * \author Joe Osborn <josborn1@bnl.gov>, Hugo Pereira Da Costa <hugo.pereira-da-costa@lanl.gov>
 */

#include "TpcGlobalPositionWrapper.h"

#include "TpcClusterZCrossingCorrection.h"
#include "TpcDistortionCorrectionContainer.h"

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>

//____________________________________________________________________________________________________________________
void TpcGlobalPositionWrapper::loadNodes( PHCompositeNode* topNode )
{
  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  // tpc distortion corrections
  m_dcc_module_edge = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerModuleEdge");
  if (m_dcc_module_edge && m_verbosity > 0)
  {
    std::cout << "TpcGlobalPositionWrapper::loadNodes - found module edge TPC distortion correction container" << std::endl;
  }

  m_dcc_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
  if (m_dcc_static && m_verbosity > 0)
  {
    std::cout << "TpcGlobalPositionWrapper::loadNodes - found static TPC distortion correction container" << std::endl;
  }
  m_dcc_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerAverage");
  if (m_dcc_average && m_verbosity > 0)
  {
    std::cout << "TpcGlobalPositionWrapper::loadNodes - found average TPC distortion correction container" << std::endl;
  }
  m_dcc_fluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerFluctuation");
  if (m_dcc_fluctuation && m_verbosity > 0)
  {
    std::cout << "TpcGlobalPositionWrapper::loadNodes - found fluctuation TPC distortion correction container" << std::endl;
  }
}

//____________________________________________________________________________________________________________________
Acts::Vector3 TpcGlobalPositionWrapper::applyDistortionCorrections(Acts::Vector3 global) const
{
  // apply distortion corrections
  if (m_enable_module_edge_corr && m_dcc_module_edge)
  {
    global = m_distortionCorrection.get_corrected_position(global, m_dcc_module_edge);
  }

  if (m_enable_static_corr && m_dcc_static)
  {
    global = m_distortionCorrection.get_corrected_position(global, m_dcc_static);
  }

  if (m_enable_average_corr && m_dcc_average)
  {
    global = m_distortionCorrection.get_corrected_position(global, m_dcc_average);
  }

  if (m_enable_fluctuation_corr && m_dcc_fluctuation)
  {
    global = m_distortionCorrection.get_corrected_position(global, m_dcc_fluctuation);
  }

  return global;
}

//____________________________________________________________________________________________________________________
Acts::Vector3 TpcGlobalPositionWrapper::getGlobalPositionDistortionCorrected(const TrkrDefs::cluskey& key, TrkrCluster* cluster, short int crossing ) const
{

  if( !m_tGeometry )
  {
    std::cout << "TpcGlobalPositionWrapper::getGlobalPositionDistortionCorrected - m_tGeometry not set" << std::endl;
    return {0,0,0};
  }

  // get global position from acts
  Acts::Vector3 global = m_tGeometry->getGlobalPosition(key, cluster);

  // make sure cluster is from TPC
  if( TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::tpcId )
  {

    // verify crossing validity
    if(crossing == SHRT_MAX)
    {
      if(!m_suppressCrossing)
      {
        std::cout << "TpcGlobalPositionWrapper::getGlobalPositionDistortionCorrected - invalid crossing." << std::endl;
      }
      return global;
    }

    // apply crossing correction
    global.z() = TpcClusterZCrossingCorrection::correctZ(global.z(), TpcDefs::getSide(key), crossing);
    // std::cout << "Global: " << global.x() << "  " << global.y() << "  " << global.z() << std::endl;

    // apply distortion corrections
    global = applyDistortionCorrections(global);
    //std::cout << "Global after dist corr: " << global.x() << "  " << global.y() << "  " << global.z() << std::endl;
  }

  return global;
}
