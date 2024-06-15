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
  // Each decay is stored as an initial pair of embedding ID and barcode.
  // This pair matches with another int, the PDGID
  using Decay = std::vector<std::pair<std::pair<int, int>, int>>;
  using Map = std::map<unsigned int, Decay>;
  using ConstIter = Map::const_iterator;
  using Iter = Map::iterator;

  ~DecayFinderContainerBase() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "DecayFinderContainer base class" << std::endl;
  }
  void Reset() override;
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int /*unused*/) const { return 0; }
  virtual void clear();

  virtual const Decay get(unsigned int) const;
  virtual Decay get(unsigned int);

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int);
  virtual Iter end();

  virtual Decay insert(const Decay&);

  virtual Map returnDecaysByPDGid(int);

  virtual size_t erase(unsigned int key);

 protected:
  DecayFinderContainerBase() = default;

 private:
  ClassDefOverride(DecayFinderContainerBase, 1);
};

#endif  // DECAYFINDER_DECAYFINDERCONTAINERBASE_H
