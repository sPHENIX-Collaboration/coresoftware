#ifndef TRACKBASE_ALIGNMENTTRANSFORMATIONCONTAINER_H
#define TRACKBASE_ALIGNMENTTRANSFORMATIONCONTAINER_H

/**
 * @file trackbase/alignmentTranformationContainer
 * @author R. Boucher
 * @date August 2022
 * @brief Storage class for alignment transformation to node tree
 */

#include "TrkrDefs.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <iostream>              // for cout, ostream
#include <map>
#include <utility>               // for pair

#include <phool/PHObject.h>


/**
 * @brief Container class for Alignment transformations
 *
 * Association object holding transformations associated with given tracker hitset
 */
class alignmentTransformationContainer
{
  
  public:

  alignmentTransformationContainer() = default;

  virtual ~alignmentTransformationContainer(){}
 
  void Reset(); 

  void identify(std::ostream &os = std::cout);  

  void addTransform(const TrkrDefs::hitsetkey, unsigned int, Eigen::Matrix4d); 

  void removeTransform(const TrkrDefs::hitsetkey, unsigned int); 
  
  Eigen::Matrix4d getTransform(const TrkrDefs::hitsetkey hitsetkey, unsigned int subsurfkey);

  private:
  
  std::map<const std::pair<TrkrDefs::hitsetkey,unsigned int>, Eigen::Matrix4d> transformMap;
  
  ClassDef(alignmentTransformationContainer,1);

};

#endif //TRACKBASE_ALIGNMENTTRANSFORMATIONCONTAINER_H
