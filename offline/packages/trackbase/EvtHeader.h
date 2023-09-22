/**
 * @file trackbase/EvtHeader.h
 * @author Joseph Bertaux
 * @date 21 Sept 2023
 * @brief Base class for Event Header object
**/
#ifndef TRACKBASE_EVT_HEADER_H
#define TRACKBASE_EVT_HEADER_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <cstdint>

#include <iostream>
#include <map>

/**
 * @brief Base class for event header object
 *
 * This is the empty virtual base class for an event header object.
 * Each subsystem should implement an inherited version
 * which contains the actual storage information.
 */
class EvtHeader : public PHObject
{
 public:
  using FEE_t = int;
  using BCO_t = uint64_t; //Alternatives are Long64_t or long long
  using Map = std::map<std::pair<TrkrDefs::TrkrId, FEE_t>, BCO_t>;
  using Interator = Map::iterator;
  using ConstIterator = Map::const_iterator;

  ~EvtHeader() override {}
  void identify(std::ostream& = std::cout) const override;
  void Reset() override;

  virtual void AddBCO(TrkrDefs::TrkrId const&, FEE_t const&, BCO_t const&);
  virtual ConstIterator GetBCO(TrkrDefs::TrkrId const&, FEE_t const&) const;

 protected:
  ClassDefOverride(EvtHeader, 1);
};

#endif//TRACKBASE_EVT_HEADER_H
