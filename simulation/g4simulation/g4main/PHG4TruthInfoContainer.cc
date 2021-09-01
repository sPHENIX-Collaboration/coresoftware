#include "PHG4TruthInfoContainer.h"

#include "PHG4Particle.h"
#include "PHG4Shower.h"
#include "PHG4VtxPoint.h"

#include <boost/tuple/tuple.hpp>

#include <limits>
#include <string>

using namespace std;

PHG4TruthInfoContainer::PHG4TruthInfoContainer()
  : particlemap()
  , vtxmap()
  , particle_embed_flags()
  , vertex_embed_flags()
{
}

PHG4TruthInfoContainer::~PHG4TruthInfoContainer() { Reset(); }

void PHG4TruthInfoContainer::Reset()
{
  for (Iterator iter = particlemap.begin(); iter != particlemap.end(); ++iter)
  {
    delete iter->second;
  }
  particlemap.clear();

  for (VtxIterator iter = vtxmap.begin(); iter != vtxmap.end(); ++iter)
  {
    delete iter->second;
  }
  vtxmap.clear();

  for (ShowerIterator iter = showermap.begin(); iter != showermap.end(); ++iter)
  {
    delete iter->second;
  }
  showermap.clear();

  particle_embed_flags.clear();
  vertex_embed_flags.clear();

  return;
}

void PHG4TruthInfoContainer::identify(ostream& os) const
{
  os << "---particlemap--------------------------" << endl;
  for (ConstIterator iter = particlemap.begin(); iter != particlemap.end(); ++iter)
  {
    os << "particle id " << iter->first << endl;
    (iter->second)->identify();
  }

  os << "---vtxmap-------------------------------" << endl;
  for (ConstVtxIterator vter = vtxmap.begin(); vter != vtxmap.end(); ++vter)
  {
    os << "vtx id: " << vter->first << endl;
    (vter->second)->identify();
  }

  os << "---showermap-------------------------------" << endl;
  for (ConstShowerIterator ster = showermap.begin(); ster != showermap.end(); ++ster)
  {
    os << "shower id: " << ster->first << endl;
    (ster->second)->identify();
  }

  os << "---list of embeded track flags-------------------" << endl;
  for (std::map<int, int>::const_iterator eter = particle_embed_flags.begin();
       eter != particle_embed_flags.end();
       ++eter)
  {
    os << "embeded track id: " << eter->first
         << " flag: " << eter->second << endl;
  }

  os << "---list of embeded vtx flags-------------------" << endl;
  for (std::map<int, int>::const_iterator eter = vertex_embed_flags.begin();
       eter != vertex_embed_flags.end();
       ++eter)
  {
    os << "embeded vertex id: " << eter->first
         << " flag: " << eter->second << endl;
  }

  os << "---primary vertex-------------------" << endl;
  os << "Vertex " << GetPrimaryVertexIndex() << " is identified as the primary vertex" << endl;

  return;
}

PHG4TruthInfoContainer::ConstIterator
PHG4TruthInfoContainer::AddParticle(const int trackid, PHG4Particle* newparticle)
{
  int key = trackid;
  ConstIterator it;
  bool added = false;
  boost::tie(it, added) = particlemap.insert(std::make_pair(key, newparticle));
  if (added) return it;

  cerr << "PHG4TruthInfoContainer::AddParticle"
       << " - Attempt to add particle with existing trackid "
       << trackid << ": " << newparticle->get_name() << " id "
       << newparticle->get_track_id()
       << ", p = [" << newparticle->get_px()
       << ", " << newparticle->get_py()
       << ", " << newparticle->get_pz() << "], "
       << " parent ID " << newparticle->get_parent_id()
       << std::endl;
  return particlemap.end();
}

PHG4Particle* PHG4TruthInfoContainer::GetParticle(const int trackid)
{
  int key = trackid;
  Iterator it = particlemap.find(key);
  if (it != particlemap.end()) return it->second;
  return nullptr;
}

PHG4Particle* PHG4TruthInfoContainer::GetPrimaryParticle(const int trackid)
{
  if (trackid <= 0) return nullptr;
  Iterator it = particlemap.find(trackid);
  if (it != particlemap.end()) return it->second;
  return nullptr;
}

PHG4VtxPoint* PHG4TruthInfoContainer::GetVtx(const int vtxid)
{
  int key = vtxid;
  VtxIterator it = vtxmap.find(key);
  if (it != vtxmap.end()) return it->second;
  return nullptr;
}

PHG4VtxPoint* PHG4TruthInfoContainer::GetPrimaryVtx(const int vtxid)
{
  if (vtxid <= 0) return nullptr;
  VtxIterator it = vtxmap.find(vtxid);
  if (it != vtxmap.end()) return it->second;
  return nullptr;
}

PHG4Shower* PHG4TruthInfoContainer::GetShower(const int showerid)
{
  int key = showerid;
  ShowerIterator it = showermap.find(key);
  if (it != showermap.end()) return it->second;
  return nullptr;
}

PHG4Shower* PHG4TruthInfoContainer::GetPrimaryShower(const int showerid)
{
  if (showerid <= 0) return nullptr;
  ShowerIterator it = showermap.find(showerid);
  if (it != showermap.end()) return it->second;
  return nullptr;
}

PHG4TruthInfoContainer::ConstVtxIterator
PHG4TruthInfoContainer::AddVertex(const int id, PHG4VtxPoint* newvtx)
{
  int key = id;
  ConstVtxIterator it;
  bool added = false;

  if (vtxmap.find(id) != vtxmap.end())
  {
    cout << "trying to add existing vtx " << id
         << " vtx pos: " << endl;
    (vtxmap.find(id)->second)->identify();
    identify();
  }

  boost::tie(it, added) = vtxmap.insert(std::make_pair(key, newvtx));
  if (added)
  {
    newvtx->set_id(key);
    return it;
  }

  cerr << "PHG4TruthInfoContainer::AddVertex"
       << " - Attempt to add vertex with existing id " << id << std::endl;
  return vtxmap.end();
}

PHG4TruthInfoContainer::ConstShowerIterator
PHG4TruthInfoContainer::AddShower(const int id, PHG4Shower* newshower)
{
  int key = id;
  ConstShowerIterator it;
  bool added = false;

  if (showermap.find(id) != showermap.end())
  {
    cout << "trying to add existing shower " << id
         << " shower pos: " << endl;
    (showermap.find(id)->second)->identify();
    identify();
  }

  boost::tie(it, added) = showermap.insert(std::make_pair(key, newshower));
  if (added)
  {
    newshower->set_id(key);
    return it;
  }

  cerr << "PHG4TruthInfoContainer::AddShower"
       << " - Attempt to add shower with existing id " << id << std::endl;
  return showermap.end();
}

int PHG4TruthInfoContainer::maxtrkindex() const
{
  int key = 0;
  if (!particlemap.empty()) key = particlemap.rbegin()->first;
  if (key < 0) key = 0;
  return key;
}

int PHG4TruthInfoContainer::mintrkindex() const
{
  int key = 0;
  if (!particlemap.empty()) key = particlemap.begin()->first;
  if (key > 0) key = 0;
  return key;
}

int PHG4TruthInfoContainer::maxvtxindex() const
{
  int key = 0;
  if (!vtxmap.empty()) key = vtxmap.rbegin()->first;
  if (key < 0) key = 0;
  return key;
}

int PHG4TruthInfoContainer::minvtxindex() const
{
  int key = 0;
  if (!vtxmap.empty()) key = vtxmap.begin()->first;
  if (key > 0) key = 0;
  return key;
}

int PHG4TruthInfoContainer::maxshowerindex() const
{
  int key = 0;
  if (!showermap.empty()) key = showermap.rbegin()->first;
  if (key < 0) key = 0;
  return key;
}

int PHG4TruthInfoContainer::minshowerindex() const
{
  int key = 0;
  if (!showermap.empty()) key = showermap.begin()->first;
  if (key > 0) key = 0;
  return key;
}

void PHG4TruthInfoContainer::delete_particle(Iterator piter)
{
  delete piter->second;
  particlemap.erase(piter);
  return;
}

void PHG4TruthInfoContainer::delete_particle(int trackid)
{
  Iterator it = particlemap.find(trackid);
  if (it != particlemap.end())
    delete_particle(it);
}

void PHG4TruthInfoContainer::delete_vtx(VtxIterator viter)
{
  delete viter->second;
  vtxmap.erase(viter);
  return;
}

void PHG4TruthInfoContainer::delete_vtx(int vtxid)
{
  VtxIterator it = vtxmap.find(vtxid);
  if (it != vtxmap.end())
    delete_vtx(it);
}

void PHG4TruthInfoContainer::delete_shower(ShowerIterator siter)
{
  delete siter->second;
  showermap.erase(siter);
  return;
}

int PHG4TruthInfoContainer::isEmbeded(const int trackid) const
{
  std::map<int, int>::const_iterator iter = particle_embed_flags.find(trackid);
  if (iter == particle_embed_flags.end())
  {
    return 0;
  }

  return iter->second;
}

int PHG4TruthInfoContainer::isEmbededVtx(const int vtxid) const
{
  std::map<int, int>::const_iterator iter = vertex_embed_flags.find(vtxid);
  if (iter == vertex_embed_flags.end())
  {
    return 0;
  }

  return iter->second;
}

bool PHG4TruthInfoContainer::is_primary_vtx(const PHG4VtxPoint* v) const
{
  return (v->get_id() > 0);
}

bool PHG4TruthInfoContainer::is_primary(const PHG4Particle* p) const
{
  return (p->get_track_id() > 0);
}

int PHG4TruthInfoContainer::GetPrimaryVertexIndex() const
{
  ConstVtxRange vrange = GetPrimaryVtxRange();

  int highest_embedding_ID = numeric_limits<int>::min();
  int vtx_id_for_highest_embedding_ID = 0;
  for (auto iter = vrange.first; iter != vrange.second; ++iter)
  {
    //        cout <<"PHG4TruthInfoContainer::GetPrimaryVertexIndex - vertex ID "<<iter->first<<" embedding ID "<< isEmbededVtx(iter->first) <<": ";
    // iter->second->identify();
    const int embedding_ID = isEmbededVtx(iter->first);

    if (embedding_ID >= highest_embedding_ID)
    {
      highest_embedding_ID = embedding_ID;
      vtx_id_for_highest_embedding_ID = iter->first;
    }
  }

  if (highest_embedding_ID == numeric_limits<int>::min())
  {
    cout << "PHG4TruthInfoContainer::GetPrimaryVertexIndex - "
         << "WARNING: no valid primary vertex. Return an invalid ID of 0"
         << endl;
    return 0;
  }

  return vtx_id_for_highest_embedding_ID;
}

bool operator==(const PHG4TruthInfoContainer::Map::value_type& lhs, const PHG4TruthInfoContainer::Map::value_type& rhs)
{
  return *lhs.second == *rhs.second;
}

bool operator==(const PHG4TruthInfoContainer::VtxMap::value_type& lhs, const PHG4TruthInfoContainer::VtxMap::value_type& rhs)
{
  return *lhs.second == *rhs.second;
}

bool operator==(const PHG4TruthInfoContainer::ShowerMap::value_type& lhs, const PHG4TruthInfoContainer::ShowerMap::value_type& rhs)
{
  return *lhs.second == *rhs.second;
}
