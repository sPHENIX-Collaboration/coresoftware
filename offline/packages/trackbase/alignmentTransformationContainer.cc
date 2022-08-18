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
    int layer = TrkrDefs::getLayer(entry.first.first);
    os << "   hitset key: "  << entry.first.first<< " subsurface key: " << entry.first.second << " layer " << layer
       << " transform: " << entry.second
       << std::endl;
  }

  os << "------------------------------" << std::endl;

  return;
}

void alignmentTransformationContainer::addTransform(const TrkrDefs::hitsetkey hitsetkey, unsigned int subsurfkey, Eigen::Matrix4d transform)
{
  std::pair<TrkrDefs::hitsetkey,unsigned int> key = std::make_pair(hitsetkey,subsurfkey);
  transformMap.insert(std::make_pair(key, transform));
}


void alignmentTransformationContainer::removeTransform(const TrkrDefs::hitsetkey hitsetkey, unsigned int subsurfkey)
{
  std::pair<TrkrDefs::hitsetkey, unsigned int> key = std::make_pair(hitsetkey,subsurfkey);

  const auto range = transformMap.equal_range(key);
  for( auto mapiter = range.first; mapiter != range.second; ++mapiter)
  {
    if(mapiter->first == key)
    {
      transformMap.erase(mapiter);
      return;
    }
  }
}

Eigen::Matrix4d alignmentTransformationContainer::getTransform(const TrkrDefs::hitsetkey hitsetkey, unsigned int subsurfkey)
{
  Eigen::Matrix4d transform;
  std::pair<TrkrDefs::hitsetkey, unsigned int> key = std::make_pair(hitsetkey,subsurfkey);

  const auto range = transformMap.equal_range(key);
  for( auto mapiter = range.first; mapiter != range.second; ++mapiter)
  {
    if(mapiter->first == key)
    {
      return mapiter->second;
    }
  }
  std::cout << "Unable to find hitsetkey: "<< key.first<< " subsurface key:"<< key.second << " in alignmentTransformationContainer" << std::endl;
  exit(1);
}
