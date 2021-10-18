#ifndef DECAYFINDER_DECAYFINDERCONTAINERBASE_H
#define DECAYFINDER_DECAYFINDERCONTAINERBASE_H

#include <phool/PHObject.h>

#include <cstddef>   // for size_t
#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for pair
#include <vector>

/**
 * @brief DecayFinder container object
 *
 * Container for DecayFinder objects, based off KFParticle_Container 
 */

class DecayFinderContainerBase : public PHObject
{
 public:
  typedef std::vector<std::pair<int, int>> Decay;
  typedef std::map<unsigned int, Decay> Map;
  typedef std::map<unsigned int, Decay>::const_iterator ConstIter;
  typedef std::map<unsigned int, Decay>::iterator Iter;

  DecayFinderContainerBase();
  DecayFinderContainerBase(const DecayFinderContainerBase& decayfindermap);
  DecayFinderContainerBase& operator=(const DecayFinderContainerBase& decayfindermap);
  virtual ~DecayFinderContainerBase();

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new DecayFinderContainerBase(*this); }

  bool empty() const { return m_decaymap.empty(); }
  size_t size() const { return m_decaymap.size(); }
  size_t count(unsigned int key) const { return m_decaymap.count(key); }
  void clear() { Reset(); }

  const Decay get(unsigned int key) const;
  Decay get(unsigned int key);

  ConstIter begin() const { return m_decaymap.begin(); }
  ConstIter find(unsigned int key) const { return m_decaymap.find(key); }
  ConstIter end() const { return m_decaymap.end(); }

  Iter begin() { return m_decaymap.begin(); }
  Iter find(unsigned int key) { return m_decaymap.find(key); }
  Iter end() { return m_decaymap.end(); }

  Decay insert(const Decay decay);

  ///Use the PDG MC ID to return a subset of the DecayFinder container, if those particle exist in the container
  Map returnDecaysByPDGid(int PDGid);

  size_t erase(unsigned int key);

  Map m_decaymap;

 private:
  ClassDefOverride(DecayFinderContainerBase, 1)
};

#endif  //DECAYFINDER_DECAYFINDERCONTAINERBASE_H
