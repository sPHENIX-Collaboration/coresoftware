#include "PHG4TruthInfoContainer.h"
#include "PHG4Particlev1.h"
#include "PHG4VtxPointv1.h"

#include <phool/phool.h>

#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>

#include <algorithm>
#include <stdexcept>
#include <cstdlib>

using namespace std;

int 
PHG4TruthInfoContainer::get_key(const int detid)
{
  return detid;
}

PHG4TruthInfoContainer::PHG4TruthInfoContainer()
{
  //std::cout << __PRETTY_FUNCTION__ << std::endl;
}

PHG4TruthInfoContainer::~PHG4TruthInfoContainer()
{
  //std::cout << __PRETTY_FUNCTION__ << std::endl;
}

// Deleter object to that takes a pair and deletes the second member
template<class T>
struct SecondDeleter : std::unary_function<T,void>
{
  void operator()(T v) { delete v.second; }
};

void
PHG4TruthInfoContainer::Reset()
{
  std::for_each(hitmap.begin(),hitmap.end(),SecondDeleter<Map::value_type>());
  hitmap.clear();

  std::for_each(vtxmap.begin(),vtxmap.end(),SecondDeleter<VtxMap::value_type>());
  vtxmap.clear();

  std::for_each(primary_particle_map.begin(),primary_particle_map.end(),SecondDeleter<Map::value_type>());
  primary_particle_map.clear();

  std::for_each(primary_vtxmap.begin(),primary_vtxmap.end(),SecondDeleter<VtxMap::value_type>());
  primary_vtxmap.clear();

  embedded_trkid.clear();

  return;
}

void
PHG4TruthInfoContainer::identify(ostream& os) const
{
   ConstIterator iter;
   cout << "---particlemap--------------------------" << endl;
   for (iter = hitmap.begin(); iter != hitmap.end(); ++iter)
     {
       cout << "hit key " <<  iter->first << endl;
       (iter->second)->identify();
     }
   ConstVtxIterator vter;
   cout << "---vtxmap-------------------------------" << endl;
   for (vter = vtxmap.begin(); vter != vtxmap.end(); ++vter)
     {
       cout << "vtx id: " << vter ->first << endl;
       (vter ->second)->identify();
     }
   cout << "---primaryparticlemap-------------------" << endl;
   for (iter = primary_particle_map.begin(); iter != primary_particle_map.end(); ++iter)
     {
       cout << "primary_particle_ key " <<  iter->first << endl;
       (iter->second)->identify();
     }
   cout << "---primary vertex emap-------------------" << endl;
   for (vter = primary_vtxmap.begin(); vter != primary_vtxmap.end(); ++vter)
     {
       cout << "vtx id: " << vter ->first << endl;
       (vter ->second)->identify();
     }
   cout << "---list of embeded tracks-------------------" << endl;
   for (std::set<int>::const_iterator eter = embedded_trkid.begin(); eter != embedded_trkid.end(); ++eter)
     {
       cout << "embeded track ID " << *eter << endl;
     }
   
  return;
}

PHG4TruthInfoContainer::ConstIterator
PHG4TruthInfoContainer::AddHit(const int trackid, PHG4Particle *newhit)
{
  int key = get_key(trackid);
  ConstIterator it;
  bool added = false;
  boost::tie(it,added) = hitmap.insert(std::make_pair(key,newhit));
  if ( added ) return it;
  cerr << "PHG4TruthInfoContainer::AddHit - Attempt to add hit with existing trackid "
      << trackid <<": "<<newhit ->get_name()<<" id "
      << newhit->get_track_id()
      <<", p = ["<<newhit->get_px()<<", "<<newhit->get_py()<<", "<<newhit->get_pz()<<"], "
      <<" parent ID "<<newhit->get_parent_id()
      << std::endl;
  return hitmap.end();
}

PHG4TruthInfoContainer::ConstIterator
PHG4TruthInfoContainer::AddPrimaryParticle(PHG4Particle *newhit)
{
  int key = primary_particle_map.size()+1;
  if ( primary_particle_map.find(key) !=  primary_particle_map.end())
    {
      cout << PHWHERE << "particle with id " << key << " already inserted, fix indexing" << endl;
      exit(1);
    }

  ConstIterator it;
  bool added = false;
  boost::tie(it,added) = primary_particle_map.insert(std::make_pair(key,newhit));
  if ( added ) return it;
  cerr << "PHG4TruthInfoContainer::AddHit - Attempt to add hit with existing trackid " << key << std::endl;
  return hitmap.end();
}

PHG4Particle*
PHG4TruthInfoContainer::GetHit(const int trackid)
{
  int key = get_key(trackid);
  Iterator it = hitmap.find(key);
  if ( it != hitmap.end() ) return it->second;
  return 0;
}

PHG4Particle*
PHG4TruthInfoContainer::GetPrimaryHit(const int trackid)
{
  int key = get_key(trackid);
  Iterator it = primary_particle_map.find(key);
  if ( it != primary_particle_map.end() ) return it->second;
  return 0;
}

PHG4TruthInfoContainer::Range
PHG4TruthInfoContainer::GetHitRange()
{
  return Range(hitmap.begin(),hitmap.end());
}

PHG4TruthInfoContainer::ConstRange
PHG4TruthInfoContainer::GetHitRange() const
{
  return ConstRange(hitmap.begin(),hitmap.end());
}

PHG4VtxPoint*
PHG4TruthInfoContainer::GetVtx(const int vtxid)
{
  int key = get_key(vtxid);
  VtxIterator it = vtxmap.find(key);
  if ( it != vtxmap.end() )
    {
      return it->second;
    }
  return NULL;
}

PHG4VtxPoint*
PHG4TruthInfoContainer::GetPrimaryVtx(const int vtxid)
{
  int key = get_key(vtxid);
  VtxIterator it = primary_vtxmap.find(key);
  if ( it != primary_vtxmap.end() )
    {
      return it->second;
    }
  return NULL;
}

PHG4TruthInfoContainer::VtxRange
PHG4TruthInfoContainer::GetVtxRange()
{
  return VtxRange(vtxmap.begin(),vtxmap.end());
}

PHG4TruthInfoContainer::ConstVtxRange
PHG4TruthInfoContainer::GetVtxRange() const
{
  return ConstVtxRange(vtxmap.begin(),vtxmap.end());
}

PHG4TruthInfoContainer::VtxIterator
PHG4TruthInfoContainer::AddVertex(const int id)
{
  int key = get_key(id);
  VtxIterator it = vtxmap.find(key);
  if ( it != vtxmap.end() ) 
    {
      std::cerr << "PHG4TruthInfoContainer::AddVertex - Attempt to add vertex with existing id " << id << std::endl;
      exit(2);
    }
  bool tmp = false;
  boost::tie(it,tmp) = vtxmap.insert( std::make_pair(key,new PHG4VtxPointv1(NAN,NAN,NAN,NAN,key)) );
  return it;
}

PHG4TruthInfoContainer::ConstVtxIterator
PHG4TruthInfoContainer::AddVertex(const int id, PHG4VtxPoint *newvtx)
{
  int key = get_key(id);
  ConstVtxIterator it;
  bool added = false;
  if (vtxmap.find(id) != vtxmap.end())
    {
      //PHG4VtxPoint *pt = 
      cout << "trying to add existing vtx " << id 
	   << " vtx pos: " << endl;
 (vtxmap.find(id)->second)->identify();
 identify();
    }
  boost::tie(it,added) = vtxmap.insert(std::make_pair(key,newvtx));
  if ( added ) {
      newvtx->set_id(key);
      return it;
  }
  cerr << "PHG4TruthInfoContainer::AddVertex - Attempt to add vertex with existing id " << id << std::endl;
  return vtxmap.end();
}

int
PHG4TruthInfoContainer::AddPrimaryVertex(PHG4VtxPoint *newvtx)
{
// initial guess of key, if we are the first index (map.size() == 0)
//  no need to check anything just add the vertex
  int key = primary_vtxmap.size()+1;  
  if (key == 1)
    {
      primary_vtxmap.insert(make_pair(key,newvtx));
      newvtx->set_id(key);
    }
  else
    {
      // we already have a vertex, check if the newly added vertex is identical to existing vertex
      ConstVtxIterator vtxiter;
      for (vtxiter = primary_vtxmap.begin(); vtxiter != primary_vtxmap.end(); ++vtxiter)
	{
	  PHG4VtxPoint *savedvtx = vtxiter->second;
	  if (*savedvtx == *newvtx)
	    {
	      key = vtxiter->first;
	      delete newvtx; // since we do not store it we have to delete it here
	      return key;
	    }
	}
      // if not found among existing vertices, add it
       primary_vtxmap.insert(make_pair(key,newvtx));
       newvtx->set_id(key);
    }
  return key;
}

int
PHG4TruthInfoContainer::maxprimarytrkindex() const
{
  int key = 0;
  if (!primary_particle_map.empty())
    {
      key = primary_particle_map.rbegin()->first;
    }
  if (key < 0)
    {
      key = 0;
    }
  return key;
}

int
PHG4TruthInfoContainer::minprimarytrkindex() const
{
  int key = 0;
  if (!primary_particle_map.empty())
    {
       key = primary_particle_map.begin()->first;
    }
  if (key > 0)
    {
      key = 0;
    }
  return key;
}

int
PHG4TruthInfoContainer::maxtrkindex() const
{
  int key = 0;
  if (!hitmap.empty())
    {
      key = hitmap.rbegin()->first;
    }
  if (key < 0)
    {
      key = 0;
    }
  return key;
}

int
PHG4TruthInfoContainer::mintrkindex() const
{
  int key = 0;
  if (!hitmap.empty())
    {
       key = hitmap.begin()->first;
    }
  if (key > 0)
    {
      key = 0;
    }
  return key;
}

int
PHG4TruthInfoContainer::maxvtxindex() const
{
  int key = 0;
  if (!vtxmap.empty())
    {
      key = vtxmap.rbegin()->first;
    }
  if (key < 0)
    {
      key = 0;
    }
  return key;
}

int
PHG4TruthInfoContainer::minvtxindex() const
{
  int key = 0;
  if (!vtxmap.empty())
    {
       key = vtxmap.begin()->first;
    }
  if (key > 0)
    {
      key = 0;
    }
  return key;
}

void
PHG4TruthInfoContainer::delete_hit(Iterator hiter)
{
  delete hiter->second;
  hitmap.erase(hiter);
  return;
}

void
PHG4TruthInfoContainer::delete_vtx(VtxIterator viter)
{
  delete viter->second;
  vtxmap.erase(viter);
  return;
}

int
PHG4TruthInfoContainer::isEmbeded(const int trackid) const
{
  if (embedded_trkid.find(trackid) != embedded_trkid.end())
    {
      return true;
    }
  return false;
}
