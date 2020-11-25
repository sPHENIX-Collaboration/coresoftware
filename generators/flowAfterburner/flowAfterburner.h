#ifndef FLOWAFTERBURNER_FLOWAFTERBURNER_H
#define FLOWAFTERBURNER_FLOWAFTERBURNER_H

#include <string>

namespace CLHEP
{
  class HepRandomEngine;
}
namespace HepMC
{
  class GenEvent;
}

enum flowAfterburnerAlgorithm
{
  minbias_algorithm,
  minbias_v2_algorithm,
  custom_algorithm
};

int flowAfterburner(HepMC::GenEvent *inEvent,
                    CLHEP::HepRandomEngine *engine,
                    std::string algorithmName,
                    float mineta, float maxeta,
                    float minpt, float maxpt);

#endif
