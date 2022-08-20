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
 


void alignmentTransformationContainer::Reset()
{ transformMap.clear(); }

void alignmentTransformationContainer::identify(std::ostream &os)
{
  os << "-----alignmentTransformationContainer-----" << std::endl;
  for( const auto& entry : transformMap )
  {
    os << " Acts Id: "  << entry.first 
       << " Transform: " << entry.second
       << std::endl;
  }
  os << "------------------------------" << std::endl;

  return;
}

void alignmentTransformationContainer::addTransform(const Acts::GeometryIdentifier id, Eigen::Matrix4d transform)
{
  transformMap.insert(std::make_pair(id, transform));
}

void alignmentTransformationContainer::removeTransform(const Acts::GeometryIdentifier id)
{
  const auto range = transformMap.equal_range(id);
  for( auto mapiter = range.first; mapiter != range.second; ++mapiter)
  {
    if(mapiter->first == id)
    {
      transformMap.erase(mapiter);
      return;
    }
  }
}

Eigen::Matrix4d alignmentTransformationContainer::getTransform(const Acts::GeometryIdentifier id)
{
  Eigen::Matrix4d transform;

  const auto range = transformMap.equal_range(id);
  for( auto mapiter = range.first; mapiter != range.second; ++mapiter)
  {
    if(mapiter->first == id)
    {
      return mapiter->second;
    }
  }
  std::cout << "Unable to find Acts Id: "<< id<<  " in alignmentTransformationContainer" << std::endl;
  exit(1);
}


std::map<const Acts::GeometryIdentifier, Eigen::Matrix4d> alignmentTransformationContainer::getContainer()
{
  return transformMap;
}


void alignmentTransformationContainer::setContainer(){}
