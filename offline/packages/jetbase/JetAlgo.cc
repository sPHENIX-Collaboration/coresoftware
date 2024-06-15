#include "JetAlgo.h"

using PropMap = std::map<Jet::PROPERTY, unsigned int>;

PropMap DummyPropMap;

PropMap& JetAlgo::property_indices() { return DummyPropMap; }
