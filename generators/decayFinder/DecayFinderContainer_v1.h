#ifndef DECAYFINDER_DECAYFINDERCONTAINERV1_H
#define DECAYFINDER_DECAYFINDERCONTAINERV1_H

#include "DecayFinderContainerBase.h"

/**
 * @brief DecayFinder container object
 *
 * Container for DecayFinder objects, based off KFParticle_Container 
 */

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

#endif  //DECAYFINDER_DECAYFINDERCONTAINERV1_H
