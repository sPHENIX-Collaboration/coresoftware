#ifndef TRUTHNEUTRALMESONCONTAINER_H
#define TRUTHNEUTRALMESONCONTAINER_H

#include "TruthNeutralMeson.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <utility>

class TruthNeutralMesonContainer : public PHObject
{
 public:
  using Map = std::map<unsigned int, TruthNeutralMeson*>;
  using Iterator = Map::iterator;
  using ConstIterator = Map::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  TruthNeutralMesonContainer() = default;
  ~TruthNeutralMesonContainer() override = default;

  void Reset() override;
  int isValid() const override { return 1; }
  void identify(std::ostream& os = std::cout) const override;

  ConstIterator AddMeson(TruthNeutralMeson* meson);

  TruthNeutralMeson* getMeson(unsigned int id);
  const TruthNeutralMeson* getMeson(unsigned int id) const;

  ConstRange getMesons() const;
  Range getMesons();

  const Map& getMesonMap() const { return _mesons; }
  Map& getMesonMap() { return _mesons; }

  unsigned int size() const { return _mesons.size(); }

 private:
  Map _mesons{};

  ClassDefOverride(TruthNeutralMesonContainer, 1)
};

#endif
