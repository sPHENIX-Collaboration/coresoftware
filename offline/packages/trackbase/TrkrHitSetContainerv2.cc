/**
 * @file trackbase/TrkrHitSetContainerv2.cc
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 * @brief Implementation for TrkrHitSetContainerv2
 */
#include "TrkrHitSetContainerv2.h"

#include "TrkrDefs.h"
#include "TrkrHitSetv1.h"

#include <TClass.h>

#include <cassert>
#include <cstdlib>

TrkrHitSetContainerv2::
    TrkrHitSetContainerv2(const std::string& hitsetclass, const size_t estimated_size)
  : m_hitArray(hitsetclass.c_str(), estimated_size)
{
}

void TrkrHitSetContainerv2::Reset()
{
  // force rebuild of indexing map
  m_hitmap.clear();

  //! fast clear without calling destructor and without marking hitsets removed.
  //! This is inspired by but even faster than TClonesArray::Clear()
  //! We just reset hitset values and reuse the hitset objects
  Int_t n = m_hitArray.GetEntriesFast();
  for (Int_t i = 0; i < n; i++) {
     TObject *obj = m_hitArray.UncheckedAt(i);
     if (obj) {
        obj->Clear();
        obj->ResetBit( kHasUUID );
        obj->ResetBit( kIsReferenced );
        obj->SetUniqueID( 0 );
     }
  }

  // alternative is to also marking hitset removed, which is not used here but optional for a v3 container
  // m_hitArray.Clear("C");
}

void TrkrHitSetContainerv2::identify(std::ostream& os) const
{
  syncMapArray();

  os << "TrkrHitSetContainerv2 with class "
     << m_hitArray.GetClass()->GetName()
     << ": Number of hits: " << size() << " index map size = " << m_hitmap.size() << std::endl;
  ConstIterator iter;
  for (const auto& pair : m_hitmap)
  {
    int layer = TrkrDefs::getLayer(pair.first);
    os << "hitsetkey " << pair.first << " layer " << layer << std::endl;
    pair.second->identify();
  }
  return;
}

TrkrHitSetContainerv2::ConstIterator
TrkrHitSetContainerv2::addHitSet(TrkrHitSet* newhit)
{
  std::cout << __PRETTY_FUNCTION__
            << " : deprecated. Use findOrAddHitSet()." << std::endl;
  return addHitSetSpecifyKey(newhit->getHitSetKey(), newhit);
}

TrkrHitSetContainerv2::ConstIterator
TrkrHitSetContainerv2::addHitSetSpecifyKey(const TrkrDefs::hitsetkey key, TrkrHitSet* newhit)
{
  std::cout << __PRETTY_FUNCTION__
            << " : deprecated. Use findOrAddHitSet()." << std::endl;

  exit(1);

  return TrkrHitSetContainer::addHitSetSpecifyKey(key, newhit);
}

void TrkrHitSetContainerv2::removeHitSet(TrkrDefs::hitsetkey )
{
  std::cout << __PRETTY_FUNCTION__
            << " : deprecated. This function still works but slows down operation." << std::endl;

  exit(1);
}

void TrkrHitSetContainerv2::removeHitSet(TrkrHitSet* hitset)
{
  removeHitSet(hitset->getHitSetKey());
}

TrkrHitSetContainerv2::ConstRange
TrkrHitSetContainerv2::getHitSets(const TrkrDefs::TrkrId trackerid) const
{
  syncMapArray();
  const TrkrDefs::hitsetkey keylo = TrkrDefs::getHitSetKeyLo(trackerid);
  const TrkrDefs::hitsetkey keyhi = TrkrDefs::getHitSetKeyHi(trackerid);
  return std::make_pair(m_hitmap.lower_bound(keylo), m_hitmap.upper_bound(keyhi));
}

TrkrHitSetContainerv2::ConstRange
TrkrHitSetContainerv2::getHitSets(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const
{
  syncMapArray();
  TrkrDefs::hitsetkey keylo = TrkrDefs::getHitSetKeyLo(trackerid, layer);
  TrkrDefs::hitsetkey keyhi = TrkrDefs::getHitSetKeyHi(trackerid, layer);
  return std::make_pair(m_hitmap.lower_bound(keylo), m_hitmap.upper_bound(keyhi));
}

TrkrHitSetContainerv2::ConstRange
TrkrHitSetContainerv2::getHitSets() const
{
  syncMapArray();
  return std::make_pair(m_hitmap.cbegin(), m_hitmap.cend());
}

TrkrHitSetContainerv2::Iterator
TrkrHitSetContainerv2::findOrAddHitSet(TrkrDefs::hitsetkey key)
{
  syncMapArray();
  auto it = m_hitmap.lower_bound(key);
  if (it == m_hitmap.end() || (key < it->first))
  {
    TrkrHitSet* hitset = (TrkrHitSet*) m_hitArray.ConstructedAt(m_hitArray.GetLast() + 1);
    assert(hitset);
    hitset -> setHitSetKey(key);
    it = m_hitmap.insert(it, std::make_pair(key, hitset));
  }
  return it;
}

TrkrHitSet*
TrkrHitSetContainerv2::findHitSet(TrkrDefs::hitsetkey key)
{
  syncMapArray();
  auto it = m_hitmap.find(key);
  if (it != m_hitmap.end())
  {
    return it->second;
  }
  else
  {
    return nullptr;
  }
}

void TrkrHitSetContainerv2::syncMapArray(void) const
{
  if (m_hitmap.size() == (size_t) size()) return;

  if (m_hitmap.size() > 0)
  {
    std::cout
        << __PRETTY_FUNCTION__ << " Error: m_hitmap and m_hitArray get out of sync, which should not happen unless DST readback. "
        << " size() = " << size()
        << " m_hitmap.size( ) = " << m_hitmap.size()
        << " m_hitArray.GetSize() = " << m_hitArray.GetSize()
        << " m_hitArray.GetLast() = " << m_hitArray.GetLast()
        << " m_hitArray.GetEntries() = " << m_hitArray.GetEntries()
        << std::endl;

    assert(m_hitmap.size() == 0);
  }

  for (unsigned int i = 0; i < size(); ++i)
  {
    TrkrHitSet* hitset = dynamic_cast<TrkrHitSet*>(m_hitArray[i]);

    if (hitset == nullptr)
    {
      std::cout << __PRETTY_FUNCTION__ << " : fatal error, invalid hitset in m_hitArray at position " << i << ". "
                << " size() = " << size()
                << " m_hitmap.size( ) = " << m_hitmap.size()
                << " m_hitArray.GetSize() = " << m_hitArray.GetSize()
                << " m_hitArray.GetLast() = " << m_hitArray.GetLast()
                << " m_hitArray.GetEntries() = " << m_hitArray.GetEntries()
                << std::endl;
      assert(hitset);
    }
    else
      m_hitmap[hitset->getHitSetKey()] = hitset;
  }
}
