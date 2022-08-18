/**
 * @file trackbase/TrkrHitTruthAssocv1.cc
 * @author R. Boucher
 * @date August 2022 
 * @brief Implementation of alignmentTransformationContainer
 */

#include "alignmentTransformationContainer.h"

#include "TrkrDefs.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <ostream>  
 
// for operator<<, endl, basic_ostream, bas...

void alignmentTransformationContainer::Reset()
{ transformMap.clear(); }

void alignmentTransformationContainer::identify(std::ostream &os)
{
  os << "-----alignmentTransformationContainer-----" << std::endl;

  for( const auto& entry : transformMap )
  {
    int layer = TrkrDefs::getLayer(entry.first);
    os << "   hitset key: "  << entry.first << " layer " << layer
       << " transform: " << entry.second
       << std::endl;
  }

  os << "------------------------------" << std::endl;

  return;
}

void alignmentTransformationContainer::addTransform(const TrkrDefs::hitsetkey hitsetkey, Eigen::Matrix4d transform)
{
   transformMap.insert(std::make_pair(hitsetkey, transform));
}


void alignmentTransformationContainer::removeTransform(const TrkrDefs::hitsetkey hitsetkey)
{
  const auto hitsetrange = transformMap.equal_range(hitsetkey);
  for( auto mapiter = hitsetrange.first; mapiter != hitsetrange.second; ++mapiter)
  {
    if(mapiter->first == hitsetkey)
    {
      transformMap.erase(mapiter);
      return;
    }
  }
}

Eigen::Matrix4d alignmentTransformationContainer::getTransform(const TrkrDefs::hitsetkey hitsetkey)
{
  Eigen::Matrix4d transform;
  const auto hitsetrange = transformMap.equal_range(hitsetkey);
  for( auto mapiter = hitsetrange.first; mapiter != hitsetrange.second; ++mapiter)
  {
    if(mapiter->first == hitsetkey)
    {
      return mapiter->second;
    }
  }
  std::cout << "Unable to find hitsetkey: "<< hitsetkey << " in alignmentTransformationContainer" << std::endl;
  exit(1);
}
