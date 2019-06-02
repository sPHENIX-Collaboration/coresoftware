#ifndef PHHEPMC_PHGENEVENTLISTV1_H
#define PHHEPMC_PHGENEVENTLISTV1_H

#include "PHGenEventList.h"
#include "PHGenEventv1.h"

#include <cstddef>          // for size_t
#include <iostream>
#include <map>
#include <vector>

class PHGenEventListv1 : public PHGenEventList
{
 public:
  typedef PHGenEventv1 PHGenEventVersion;

  PHGenEventListv1();
  virtual ~PHGenEventListv1();

  size_t size() const { return _genevents.size(); }
  const PHGenEvent* at(size_t i) const;
  PHGenEvent* at(size_t i);

  bool has(unsigned int id) const;
  size_t find(unsigned int id) const;
  const PHGenEvent* fetch(unsigned int id) const;

  void insert(const PHGenEvent* genevent);
  unsigned int generate_id() const;
  void remove(size_t i);
  void clear();

  void identify(std::ostream& out = std::cout) const;
  void print(std::ostream& out = std::cout) const;
  void Reset();

 private:
  bool stale() const { return _stale; }
  void refresh() const;

  std::vector<PHGenEventVersion> _genevents;

#ifndef __CINT____                                                     // hide from dictionary generation
  mutable bool _stale;                                                 //! exclude from ROOT I/O
  mutable std::map<unsigned int, const PHGenEvent*> _id_genevent_map;  //! exclude from ROOT I/O
#endif                                                                 // __CINT__

  ClassDef(PHGenEventListv1, 1)
};

#endif  // PHHEPMC_PHGENEVENTLISTV1_H
