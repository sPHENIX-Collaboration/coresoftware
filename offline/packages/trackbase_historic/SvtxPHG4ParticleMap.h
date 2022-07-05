#ifndef TRACKBASEHISTORIC_SVTXPHG4PARTICLEMAP_H
#define TRACKBASEHISTORIC_SVTXPHG4PARTICLEMAP_H

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <set>

class SvtxPHG4ParticleMap : public PHObject
{
 public:
  /// Reco->Truth map with structure <reco track key, std::map< weight, std::set<g4part id>>>
  typedef std::map<float, std::set<int>> WeightedTruthTrackMap;
  typedef std::map<unsigned int, WeightedTruthTrackMap> Map;
  typedef std::map<unsigned int, WeightedTruthTrackMap>::const_iterator ConstIter;
  typedef std::map<unsigned int, WeightedTruthTrackMap>::iterator Iter;

  ~SvtxPHG4ParticleMap() override {}

  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxPHG4ParticleMap base class " << std::endl;
  }
  
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }
  void Reset() override {}
  virtual bool empty() const { return true; }
  virtual std::size_t size() const { return 0; }
  virtual std::size_t count(const unsigned int) const { return 0; }
  virtual void clear() {}

  virtual bool processed() const { return false; }
  virtual void setProcessed(const bool) {}

  virtual const WeightedTruthTrackMap & get(const unsigned int) const;
  virtual WeightedTruthTrackMap & get(const unsigned int);
  virtual WeightedTruthTrackMap insert(const unsigned int, const WeightedTruthTrackMap);
  virtual std::size_t erase(const unsigned int) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(const unsigned int) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(const unsigned int);
  virtual Iter end();

 protected:
  SvtxPHG4ParticleMap() {}

 private:
  ClassDefOverride(SvtxPHG4ParticleMap, 1);
  
};

#endif
