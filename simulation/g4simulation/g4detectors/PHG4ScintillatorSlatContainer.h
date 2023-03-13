// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SCINTILLATORSLATCONTAINER_H
#define G4DETECTORS_PHG4SCINTILLATORSLATCONTAINER_H

#include "PHG4ScintillatorSlatDefs.h"  // for keytype

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <map>
#include <set>
#include <utility>  // for pair

class PHG4ScintillatorSlat;

class PHG4ScintillatorSlatContainer : public PHObject
{
 public:
  typedef std::map<PHG4ScintillatorSlatDefs::keytype, PHG4ScintillatorSlat *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::set<int>::const_iterator LayerIter;
  typedef std::pair<LayerIter, LayerIter> LayerRange;

  PHG4ScintillatorSlatContainer() {}

  ~PHG4ScintillatorSlatContainer() override {}

  // from PHObject
  void identify(std::ostream &os = std::cout) const override;
  void Reset() override;

  ConstIterator AddScintillatorSlat(const PHG4ScintillatorSlatDefs::keytype key, PHG4ScintillatorSlat *newscintillatorSlat);

  //! preferred removal method, key is currently the slat id
  void RemoveScintillatorSlat(PHG4ScintillatorSlatDefs::keytype key)
  {
    slatmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveScintillatorSlat(PHG4ScintillatorSlat *slat)
  {
    Iterator its = slatmap.begin();
    Iterator last = slatmap.end();
    for (; its != last;)
    {
      if (its->second == slat)
      {
        slatmap.erase(its++);
      }
      else
      {
        ++its;
      }
    }
  }

  //! return all scintillatorSlats matching a given detid
  ConstRange getScintillatorSlats(const short icolumn) const;

  //! return all hist
  ConstRange getScintillatorSlats(void) const;

  PHG4ScintillatorSlat *findScintillatorSlat(PHG4ScintillatorSlatDefs::keytype key);

  unsigned int size(void) const
  {
    return slatmap.size();
  }

  double getTotalEdep() const;

 protected:
  Map slatmap;

  ClassDefOverride(PHG4ScintillatorSlatContainer, 1)
};

#endif
