/**
 * @file trackbase/TrkrHitTruthAssocv1.cc
 * @author R. Boucher
 * @date August 2022 
 * @brief Implementation of alignmentTransformationContainer
 */

#include "alignmentTransformationContainer.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <algorithm>
#include <ostream>  
 
bool alignmentTransformationContainer::use_alignment = false;

void alignmentTransformationContainer::Reset()
{ transformMap.clear(); }

void alignmentTransformationContainer::identify(std::ostream &os)
{

  os << "-----alignmentTransformationContainer-----" << std::endl;
  for( const auto& entry : transformMap )
  {
    os << " Acts Id: "  << entry.first 
       << " Transform: " << entry.second.matrix()
       << std::endl;
  }
  os << "------------------------------" << std::endl;
 
  return;
}

void alignmentTransformationContainer::addTransform(const Acts::GeometryIdentifier id, Acts::Transform3 transform)
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

Acts::Transform3& alignmentTransformationContainer::getTransform(const Acts::GeometryIdentifier id)
{
  Acts::Transform3 transform; // intialize this or return something instead of exit

  const auto range = transformMap.equal_range(id);
  for( auto mapiter = range.first; mapiter != range.second; ++mapiter)
  {
    if(mapiter->first == id)
    {
      return mapiter->second;
    }
  }
  std::cout << "Unable to find Acts Id: "<< id<<  " in alignmentTransformationContainer" << std::endl;
  exit(1); //later find way to detect null transform that you must return 
}


const std::map<Acts::GeometryIdentifier, Acts::Transform3> alignmentTransformationContainer::getMap()
{
  return transformMap;
}


void alignmentTransformationContainer::set(){}
