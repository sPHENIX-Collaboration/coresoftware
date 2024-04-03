#include "JetAlgo.h"

typedef std::map<Jet::PROPERTY, unsigned int> PropMap;

PropMap DummyPropMap;

PropMap& JetAlgo::property_indices() { return DummyPropMap; }
