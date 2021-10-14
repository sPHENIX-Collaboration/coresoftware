#ifndef DECAYFINDERCONTAINERV1_H
#define DECAYFINDERCONTAINERV1_H

#include "DecayFinderContainerBase.h"

//#include <phool/PHObject.h>
//
//#include <cstddef>   // for size_t
//#include <iostream>  // for cout, ostream
//#include <map>
//#include <utility>  // for pair

/**
 * @brief DecayFinder container object
 *
 * Container for DecayFinder objects, based off KFParticle_Container 
 */

class PHObject;

class DecayFinderContainer_v1 : public DecayFinderContainerBase 
{
 public:
  DecayFinderContainer_v1();
  DecayFinderContainer_v1(const DecayFinderContainer_v1& decayfindermap);
  DecayFinderContainer_v1& operator=(const DecayFinderContainer_v1& decayfindermap);
  ~DecayFinderContainer_v1() override;

 private:
  ClassDefOverride(DecayFinderContainer_v1, 1)
};

#endif  //DECAYFINDERCONTAINERV1_H
