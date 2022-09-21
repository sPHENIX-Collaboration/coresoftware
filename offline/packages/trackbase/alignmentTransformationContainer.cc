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
  //  for( auto& entry : transformMap )
  for(unsigned int i=0; i< transformMap.size(); ++i)
    {
      auto& layerMap = transformMap[i];
      if(layerMap.empty()) continue;

      os << " Layer: "  << i << std::endl; 

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


  if(sphlayer < transformMap.size())
    {
      auto& layerMap = transformMap[sphlayer]; 
      if(layerMap.empty()) 
	{
	  std::map<Acts::GeometryIdentifier, Acts::Transform3> newmap;
	  transformMap[sphlayer] = newmap;
	}

      transformMap[sphlayer].insert(std::make_pair(id, transform));
    }
    else if(sphlayer == transformMap.size())
      {
	std::map<Acts::GeometryIdentifier, Acts::Transform3> newmap;
	newmap.insert(std::make_pair(id, transform));
	transformMap.push_back(newmap);
      }
    else
      {
	std::map<Acts::GeometryIdentifier, Acts::Transform3> newmap;
	newmap.insert(std::make_pair(id, transform));
	transformMap.resize(sphlayer+1, newmap);
      }
}

void alignmentTransformationContainer::removeTransform(const Acts::GeometryIdentifier id)
{
  unsigned int sphlayer = getsphlayer(id);
  auto& layerMap = transformMap[sphlayer];
  
  auto lyr_mapit = layerMap.find(id);
  
  if(lyr_mapit != layerMap.end())
    {
      layerMap.erase(lyr_mapit);
      return;
    }
    
}

Acts::Transform3& alignmentTransformationContainer::getTransform(const Acts::GeometryIdentifier id)
{
  unsigned int sphlayer = getsphlayer(id);

  auto& layerMap = transformMap.at(sphlayer);
  if(!layerMap.empty())
    {
      auto lyr_mapit = layerMap.find(id);
      if(lyr_mapit != layerMap.end())
	{
	  return lyr_mapit->second;
	}
    }
  
  std::cout << "Unable to find Acts Id: "<< id<<  " in alignmentTransformationContainer" << std::endl;
  exit(1); 
}

//
//const std::map<unsigned int, std::map<Acts::GeometryIdentifier, Acts::Transform3>>& alignmentTransformationContainer::getMap()
const std::vector<std::map<Acts::GeometryIdentifier, Acts::Transform3>>& alignmentTransformationContainer::getMap()
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
