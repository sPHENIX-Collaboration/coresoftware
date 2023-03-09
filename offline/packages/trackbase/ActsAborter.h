
#ifndef TRACKBASE_ACTSABORTER_H
#define TRACKBASE_ACTSTRACKFITTINGALGORITHM_H


#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>

struct ActsAborter {
  
  unsigned int m_abortlayer = std::numeric_limits<unsigned int>::max();

  template <typename propagator_state_t, typename stepper_t>
    bool operator()(propagator_state_t& state, const stepper_t& /*stepper*/) const {

    if(state.navigation.targetReached) 
      {
	return true;
      }
    
    auto layerno = state.navigation.currentSurface->geometryId().layer();
    if(layerno == m_abortlayer)
      {
	state.navigation.targetReached = true;
	return true;
      }

    return false;

  }
};


#endif
