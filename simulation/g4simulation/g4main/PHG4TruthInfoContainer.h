#ifndef __PHG4TRUTHINFOCONTAINER_H__
#define __PHG4TRUTHINFOCONTAINER_H__


#include <phool/PHObject.h>
#include <map>
#include <set>

class PHG4Particle;
class PHG4VtxPoint;

class PHG4TruthInfoContainer: public PHObject
{
  public:
  typedef std::map<int,PHG4Particle *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  typedef std::map<int,PHG4VtxPoint *> VtxMap;
  typedef VtxMap::iterator VtxIterator;
  typedef VtxMap::const_iterator ConstVtxIterator;
  typedef std::pair<VtxIterator, VtxIterator> VtxRange;
  typedef std::pair<ConstVtxIterator, ConstVtxIterator> ConstVtxRange;

  PHG4TruthInfoContainer();

  virtual ~PHG4TruthInfoContainer();

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  //! Add a hit that the user has created (use with caution)
  ConstIterator AddHit(const int detid, PHG4Particle *newhit);

  //! Add a hit and return an iterator to the user
  //  Iterator AddHit(const int detid);

  PHG4Particle* GetHit(const int detid);
  PHG4Particle* GetPrimaryHit(const int detid);

  //! Add a vertex and return an iterator to the user
  VtxIterator AddVertex(const int detid);

  //! Add a vertex and return an iterator to the user
  ConstVtxIterator AddVertex(const int detid, PHG4VtxPoint *);

  //! Add a primary vertex and return index to the user
  int AddPrimaryVertex(PHG4VtxPoint *);

  ConstIterator AddPrimaryParticle(PHG4Particle *newparticle);

  PHG4VtxPoint* GetVtx(const int detid);

  PHG4VtxPoint* GetPrimaryVtx(const int vtxid);

  //! Get a range of iterators covering the entire container
  Range GetHitRange();
  ConstRange GetHitRange() const;
  
  //! Get a range of iterators covering the entire vertex container
  VtxRange GetVtxRange();
  ConstVtxRange GetVtxRange() const;

  //! hit size
  unsigned int size( void ) const
  { return hitmap.size(); }

  //! Get the number of vertices stored
  unsigned int GetNumVertices() const { return vtxmap.size(); }

  //! Get the map itself
  const Map& GetMap() const { return hitmap; }
  const Map& GetPrimaryMap() const { return primary_particle_map; }
  const VtxMap& GetVtxMap() const { return vtxmap; }

  int maxtrkindex() const;
  int mintrkindex() const;

  int maxvtxindex() const;
  int minvtxindex() const;
 
  void delete_hit(Iterator hiter);
  void delete_vtx(VtxIterator viter);

  int GetPrimaryVertexIndex() {return (primary_vtxmap.begin())->first;}
  int GetNumPrimaryVertexParticles() {return primary_particle_map.size();}
  
  int GetLastParticleIndex() {return (primary_particle_map.rbegin())->first;}

  std::pair< std::set<int>::const_iterator, std::set<int>::const_iterator > GetEmbeddedTrkIds() const
    {return std::make_pair(embedded_trkid.begin(), embedded_trkid.end());}
  void AddEmbededTrkId(const int i) {embedded_trkid.insert(i);}

  int isEmbeded(const int trackid) const;

 protected:

  //! generate a key
  static int get_key(const int detid);

  Map hitmap;
  VtxMap vtxmap;

  Map primary_particle_map;
  VtxMap primary_vtxmap;
  // track ids of embedded particles
  std::set<int> embedded_trkid;

  ClassDef(PHG4TruthInfoContainer,1)
};

#endif
