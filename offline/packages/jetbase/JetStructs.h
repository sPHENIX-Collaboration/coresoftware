#ifndef JETBASE_JETSTRUCTS_H
#define JETBASE_JETSTRUCTS_H

// A few convenience structures for either iterating over data, or referencing between classes

#include <TClonesArray.h>
#include "Jet.h"

class Jetv2;

// ---------------------------------------------------------------------------------------
// Convenience class for iterating over jets in TClonesArray in a JetContainer
// ---------------------------------------------------------------------------------------
struct IterJetv2TCA {
  TClonesArray* tca   { nullptr };
  Jetv2*& current_jet ;
  int index { 0 };
  int size;

  // build Iterator -- capture reference to current_jet pointer from JetContainer
  IterJetv2TCA  (TClonesArray* _tca, Jetv2*& _in_jet) 
      : tca{_tca}, current_jet{_in_jet}, size { tca->GetEntriesFast() }
  {
    current_jet = (Jetv2*) tca->UncheckedAt(0);
  }

  void operator++() {
      current_jet = (Jetv2*) tca->UncheckedAt(++index);
  };

  Jetv2* operator*() { return current_jet; };

  bool operator!=(const IterJetv2TCA& rhs) { 
    if (index == rhs.size) {
      current_jet = (Jetv2*) tca->UncheckedAt(0);
      return false;
    } else {
      return true;
    }
  };
};

struct JetV2SortingCriteria {
  Jet::SORT       criteria   { Jet::SORT::PT };
  Jet::PROPERTY   property   { Jet::PROPERTY::no_property }; // when sorted by property
  unsigned int    prop_index { 0 }; // for use when sorting by criteria
  Jet::SORT_ORDER order      { Jet::SORT_ORDER::DESCENDING }; // 
};

#endif
