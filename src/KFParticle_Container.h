#ifndef KFParticle_Container_H
#define KFParticle_Container_H

#include "KFParticle.h"

#include <phool/PHObject.h>

#include <map>

#include <cstddef>        // for size_t
#include <iostream>       // for cout, ostream
#include <utility>        // for pair

using namespace std;

class KFParticle;

/**
 * @brief KFParticle container object
 *
 * Container for KFParticle objects
 */

class KFParticle_Container : public PHObject
{
 public:

  typedef map<unsigned int, KFParticle*> Map;
  typedef map<unsigned int, KFParticle*>::const_iterator ConstIter;
  typedef map<unsigned int, KFParticle*>::iterator Iter;

  KFParticle_Container();
  KFParticle_Container(const KFParticle_Container& kfparticlemap);
  KFParticle_Container& operator=(const KFParticle_Container& kfparticlemap);
  virtual ~KFParticle_Container();

  void identify(ostream& os = cout) const;
  void Reset();
  int isValid() const { return 1; }
  PHObject* CloneMe() const { return new KFParticle_Container(*this); }

  bool empty() const { return m_kfpmap.empty(); }
  size_t size() const { return m_kfpmap.size(); }
  size_t count(unsigned int key) const { return m_kfpmap.count(key); }
  void clear() { Reset(); }

  const KFParticle* get(unsigned int key) const;
  KFParticle* get(unsigned int key);

  ConstIter begin() const { return m_kfpmap.begin(); }
  ConstIter find(unsigned int key) const { return m_kfpmap.find(key); }
  ConstIter end() const { return m_kfpmap.end(); }

  Iter begin() { return m_kfpmap.begin(); }
  Iter find(unsigned int key) { return m_kfpmap.find(key); }
  Iter end() { return m_kfpmap.end(); }

  KFParticle* insert(const KFParticle *particle);
  ConstIter addParticle(KFParticle *particle);
  ConstIter addParticleSpecifyKey(unsigned int key, KFParticle *particle);

  Map returnParticlesByPDGid( int PDGid );

  size_t erase(unsigned int key)
  {
    delete m_kfpmap[key];
    return m_kfpmap.erase(key);
  }

 protected:

  Map m_kfpmap;
  //ClassDef(KFParticle_Container, 1)

};

#endif //KFParticle_Container_H
