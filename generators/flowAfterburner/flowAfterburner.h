#include <HepMC/GenEvent.h>
#include <CLHEP/Random/RandFlat.h>

#include <string>

enum flowAfterburnerAlgorithm 
  {
    minbias_algorithm,
    minbias_v2_algorithm,
    custom_algorithm
  };

  
int
flowAfterburner(HepMC::GenEvent *inEvent, 
		CLHEP::HepRandomEngine *engine, 
		std::string algorithmName,
		float mineta, float maxeta,
		float minpt, float maxpt);

