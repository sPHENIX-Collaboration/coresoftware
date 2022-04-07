#ifndef TRACKBASEHISTORIC_PHG4PARTICLESVTXMAP_H
#define TRACKBASEHISTORIC_PHG4PARTICLESVTXMAP_H

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <set>

class PHG4ParticleSvtxMap : public PHObject
{
 public:
  /// Truth->reco map with structure <g4part id, std::map< weight, std::set<reco track id>>>
  typedef std::map<float, std::set<unsigned int>> WeightedRecoTrackMap;
  typedef std::map<int, WeightedRecoTrackMap> Map;
  typedef std::map<int, WeightedRecoTrackMap>::const_iterator ConstIter;
  typedef std::map<int, WeightedRecoTrackMap>::iterator Iter;

  ~PHG4ParticleSvtxMap() override {}

  void identify(std::ostream& os = std::cout) const override
  {
    os << "PHG4ParticleSvtxMap base class " << std::endl;
  }
  
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }
  
  virtual bool empty() const { return true; }
  virtual std::size_t size() const { return 0; }
  virtual std::size_t count() const { return 0; }
  virtual void clear() {}

  virtual const WeightedRecoTrackMap get(unsigned int) const;
  virtual WeightedRecoTrackMap get(unsigned int);
  virtual ConstIter insert(const WeightedRecoTrackMap);
  virtual std::size_t erase(unsigned int) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int);
  virtual Iter end();

 protected:
  PHG4ParticleSvtxMap() {}

 private:
  ClassDefOverride(PHG4ParticleSvtxMap, 1);
  
};

#endif
