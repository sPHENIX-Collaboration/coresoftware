
#ifndef TRACKBASE_ACTSABORTER_H
#define TRACKBASE_ACTSABORTER_H

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>

struct ActsAborter
{
  unsigned int abortlayer = std::numeric_limits<unsigned int>::max();
  unsigned int abortvolume = std::numeric_limits<unsigned int>::max();

  template <typename propagator_state_t, typename stepper_t,
    typename navigator_t>
    bool operator()(propagator_state_t& state, const stepper_t& /*stepper*/,
		    const navigator_t& navigator, const Acts::Logger&/*logger*/) const
  {
    if (navigator.targetReached(state.navigation))
    {
      return true;
    }

    //if (!state.navigation.currentSurface)
    if(!navigator.currentSurface(state.navigation))
    {
      return false;
    }

    auto volumeno = state.navigation.currentSurface->geometryId().volume();
    auto layerno = state.navigation.currentSurface->geometryId().layer();
    auto sensitive = state.navigation.currentSurface->geometryId().sensitive();

    /// Check that we are in the proper layer and that we've also reached
    /// a sensitive surface
    if (layerno == abortlayer and volumeno == abortvolume and sensitive != 0)
    {
      navigator.targetReached(state.navigation, true);
      return true;
    }

    return false;
  }
};

#endif
