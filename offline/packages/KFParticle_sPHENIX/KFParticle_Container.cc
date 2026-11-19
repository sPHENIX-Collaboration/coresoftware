/*****************/
/* Cameron Dean  */
/*   LANL 2020   */
/* cdean@bnl.gov */
/*****************/

// Ideas taken from SvtxTrackMap_v1 & TrkrClusterContainer

#include "KFParticle_Container.h"

#include <KFParticle.h>

#include <cstdlib>
#include <iterator>  // for reverse_iterator
#include <map>       // for _Rb_tree_const_iterator, _Rb_tree_iterator
#include <ostream>   // for operator<<, endl, ostream, basic_ostream, bas...
#include <utility>   // for pair, make_pair

KFParticle_Container::KFParticle_Container() = default;

// NOLINTNEXTLINE(bugprone-copy-constructor-init)
KFParticle_Container::KFParticle_Container(const KFParticle_Container& kfparticlemap)
{
  for (const auto& [id, particle] : kfparticlemap)
  {
    m_kfpmap.emplace(id, dynamic_cast<KFParticle*>(particle->Clone()));
  }
}

KFParticle_Container& KFParticle_Container::operator=(const KFParticle_Container& kfparticlemap)
{
  Reset();
  for (const auto& [id, particle] : kfparticlemap)
  {
    m_kfpmap.emplace(id, dynamic_cast<KFParticle*>(particle->Clone()));
  }
  return *this;
}

KFParticle_Container::~KFParticle_Container()
{
  Reset();
}

void KFParticle_Container::Reset()
{
  for (auto& iter : m_kfpmap)
  {
    KFParticle* particle = iter.second;
    delete particle;
  }
  m_kfpmap.clear();
}

void KFParticle_Container::identify(std::ostream& os) const
{
  os << "KFParticle_Container: size = " << m_kfpmap.size() << std::endl;
  return;
}

const KFParticle* KFParticle_Container::get(unsigned int id) const
{
  const auto iter = m_kfpmap.find(id);
  return iter == m_kfpmap.end() ? nullptr : iter->second;
}

KFParticle* KFParticle_Container::get(unsigned int id)
{
  const auto iter = m_kfpmap.find(id);
  return iter == m_kfpmap.end() ? nullptr : iter->second;
}

KFParticle* KFParticle_Container::insert(const KFParticle* particle)
{
  const unsigned int index = m_kfpmap.empty() ? 0 : m_kfpmap.rbegin()->first + 1;
  // m_kfpmap[index]->SetId(index);
  return m_kfpmap.emplace(index, dynamic_cast<KFParticle*>(particle->Clone())).first->second;
}

KFParticle_Container::ConstIter
KFParticle_Container::addParticle(KFParticle* particle)
{
  return addParticleSpecifyKey(particle->Id(), particle);
}

KFParticle_Container::ConstIter
KFParticle_Container::addParticleSpecifyKey(unsigned int id, KFParticle* particle)
{
  auto ret = m_kfpmap.insert(std::make_pair(id, particle));
  if (!ret.second)
  {
    std::cout << "KFParticle_Container::AddParticleSpecifyKey: duplicate id: " << id << " exiting now" << std::endl;
    exit(1);
  }
  return ret.first;
}

KFParticle_Container::Map
KFParticle_Container::returnParticlesByPDGid(int PDGid)
{
  Map requiredParticles;

  for (const auto& iter : m_kfpmap)
  {
    if (iter.second->GetPDG() == PDGid)
    {
      requiredParticles.insert(iter);
    }
  }

  return requiredParticles;
}

size_t KFParticle_Container::erase(unsigned int key)
{
  delete m_kfpmap[key];
  return m_kfpmap.erase(key);
}
