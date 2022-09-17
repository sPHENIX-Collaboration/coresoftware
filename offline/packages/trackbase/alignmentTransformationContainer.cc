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
  for( auto& entry : transformMap )
  {
    os << " Layer: "  << entry.first << std::endl; 

    auto layerMap = entry.second;
    for(auto lyr_entry : layerMap)
      {
	os << " Acts Id: "  << lyr_entry.first 
	   << " Transform: " << lyr_entry.second.matrix()
	   << std::endl;
      }
  }
  os << "------------------------------" << std::endl;
  
  return;
}

void alignmentTransformationContainer::addTransform(const Acts::GeometryIdentifier id, Acts::Transform3 transform)
{
  unsigned int sphlayer = getsphlayer(id);
  // find or create the relevant map
  auto it = transformMap.find(sphlayer);
  if(it == transformMap.end())
    {
      std::map<Acts::GeometryIdentifier, Acts::Transform3> newMap; 
      newMap.insert(std::make_pair(id, transform));
      transformMap.insert(std::make_pair(sphlayer, newMap));
    }
  else
    {
      it->second.insert(std::make_pair(id, transform));
    }
}

void alignmentTransformationContainer::removeTransform(const Acts::GeometryIdentifier id)
{
  unsigned int sphlayer = getsphlayer(id);
  auto it = transformMap.find(sphlayer);
  if(it != transformMap.end())
    {
      auto lyr_mapit = it->second.find(id);

	if(lyr_mapit != it->second.end())
	  {
	    it->second.erase(lyr_mapit);
	    return;
	  }
    }
}

Acts::Transform3& alignmentTransformationContainer::getTransform(const Acts::GeometryIdentifier id)
{
  unsigned int sphlayer = getsphlayer(id);
  auto it = transformMap.find(sphlayer);
  if(it != transformMap.end())
    {
      auto lyr_mapit = it->second.find(id);

	if(lyr_mapit != it->second.end())
	  {
	    return lyr_mapit->second;
	  }
    }

  std::cout << "Unable to find Acts Id: "<< id<<  " in alignmentTransformationContainer" << std::endl;
  exit(1); //later find way to detect null transform that you must return 
}

const std::map<unsigned int, std::map<Acts::GeometryIdentifier, Acts::Transform3>>& alignmentTransformationContainer::getMap()
{
  return transformMap;
}

unsigned int alignmentTransformationContainer::getsphlayer(Acts::GeometryIdentifier id)
{
  unsigned int layer = id.layer(); 
  unsigned int volume = id.volume(); 
  unsigned int sphlayer = base_layer_map.find(volume)->second + layer / 2 -1;
  
  return sphlayer;
}

void alignmentTransformationContainer::set(){}
