#ifndef DECAYFINDER_DECAYFINDERCONTAINERV1_H
#define DECAYFINDER_DECAYFINDERCONTAINERV1_H

#include "DecayFinderContainerBase.h"

#include <cstddef>        // for size_t
#include <iostream>        // for cout, ostream

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

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new DecayFinderContainerBase(*this); }

  bool empty() const override { return m_decaymap.empty(); }
  size_t size() const override { return m_decaymap.size(); }
  size_t count(unsigned int key) const override { return m_decaymap.count(key); }
  void clear() override { Reset(); }

  const Decay get(unsigned int key) const override;
  Decay get(unsigned int key) override;

  ConstIter begin() const override { return m_decaymap.begin(); }
  ConstIter find(unsigned int key) const override { return m_decaymap.find(key); }
  ConstIter end() const override { return m_decaymap.end(); }

  Iter begin() override { return m_decaymap.begin(); }
  Iter find(unsigned int key) override { return m_decaymap.find(key); }
  Iter end() override { return m_decaymap.end(); }

  Decay insert(const Decay decay) override;

  //Use the PDG MC ID to return a subset of the DecayFinder container, if those particle exist in the container
  Map returnDecaysByPDGid(int PDGid) override;
  
  size_t erase(unsigned int key) override;
 
 private:
  Map m_decaymap;

  ClassDefOverride(DecayFinderContainer_v1, 1);
};

#endif  //DECAYFINDER_DECAYFINDERCONTAINERV1_H
