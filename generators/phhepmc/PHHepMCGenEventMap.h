#ifndef PHHEPMC_PHHEPMCGENEVENTMAP_H
#define PHHEPMC_PHHEPMCGENEVENTMAP_H

#include "PHHepMCGenEvent.h"

#include <phool/PHObject.h>

#include <cstddef>           // for size_t
#include <iostream>           // for cout, ostream
#include <map>

//! \brief PHHepMCGenEventMap is collection of HEPMC events input into this simulation
//! map of embedding ID -> PHHepMCGenEvent
//! positive ID is the embedded event of interest, e.g. jetty event from pythia
//! negative IDs are backgrounds, .e.g out of time pile up collisions
//! Usually, ID = 0 means the primary Au+Au collision background
class PHHepMCGenEventMap : public PHObject
{
 public:
  //! map of embedding ID -> PHHepMCGenEvent
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  typedef std::map<int, PHHepMCGenEvent*> HepMCGenEventMap;
  typedef std::map<int, PHHepMCGenEvent*>::const_iterator ConstIter;
  typedef std::map<int, PHHepMCGenEvent*>::iterator Iter;
  typedef std::map<int, PHHepMCGenEvent*>::const_reverse_iterator ConstReverseIter;
  typedef std::map<int, PHHepMCGenEvent*>::reverse_iterator ReverseIter;

  PHHepMCGenEventMap() = default;
  PHHepMCGenEventMap(const PHHepMCGenEventMap& eventmap);
  PHHepMCGenEventMap& operator=(const PHHepMCGenEventMap& eventmap);

  ~PHHepMCGenEventMap() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return 1; }
  PHHepMCGenEventMap* CloneMe() const override { return new PHHepMCGenEventMap(*this); }
  //! container service
  bool empty() const { return _map.empty(); }
  size_t size() const { return _map.size(); }
  size_t count(int idkey) const { return _map.count(idkey); }
  void clear() { Reset(); }
  //! fetch event
  const PHHepMCGenEvent* get(int idkey) const;
  PHHepMCGenEvent* get(int idkey);

  //! insert a event of interest, e.g. jetty event from pythia
  PHHepMCGenEvent* insert(const PHHepMCGenEvent* event) { return insert_active_event(event); }
  //! insert a event of interest, e.g. jetty event from pythia
  PHHepMCGenEvent* insert_active_event(const PHHepMCGenEvent* event = nullptr);
  //! insert a event of background, e.g. Au+Au collision background. First event has embedding ID = 0, which is usually the primary Au+Au collision in the case of HI embedding
  PHHepMCGenEvent* insert_background_event(const PHHepMCGenEvent* event = nullptr);
  //! insert a event with specific embedding ID
  PHHepMCGenEvent* insert_event(const int embedding_id, const PHHepMCGenEvent* event = nullptr);

  size_t erase(int idkey)
  {
    delete _map[idkey];
    return _map.erase(idkey);
  }

  //! find
  ConstIter find(unsigned int idkey) const { return _map.find(idkey); }
  Iter find(int idkey) { return _map.find(idkey); }
  //! iterator from lowest ID to highest, i.e. background to signal
  ConstIter begin() const { return _map.begin(); }
  ConstIter end() const { return _map.end(); }
  Iter begin() { return _map.begin(); }
  Iter end() { return _map.end(); }
  //! iterator from lowest ID to highest, i.e. signal to background
  ConstReverseIter rbegin() const { return _map.rbegin(); }
  ConstReverseIter rend() const { return _map.rend(); }
  ReverseIter rbegin() { return _map.rbegin(); }
  ReverseIter rend() { return _map.rend(); }
  //
  //! for c++11 range-based for loop
  const HepMCGenEventMap& get_map() const { return _map; }
  HepMCGenEventMap& get_map() { return _map; }

 private:
  HepMCGenEventMap _map;

  ClassDefOverride(PHHepMCGenEventMap, 4);
};

#endif
