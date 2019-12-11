/*!
 *  \file		Tools.cc
 *  \brief		tools
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

//GenFit
#include <GenFit/RKTrackRep.h>
#include <GenFit/AbsTrackRep.h>           // for AbsTrackRep
#include <GenFit/Exception.h>             // for Exception
#include <GenFit/MeasuredStateOnPlane.h>  // for MeasuredStateOnPlane


#include <TVector3.h>                     // for TVector3

//STL
#include <cassert>                       // for assert
#include <iostream>                       // for operator<<, basic_ostream
#include <limits>

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl

//#define _DEBUG_

namespace PHGenFit
{
double extrapolateToCylinder(
    genfit::MeasuredStateOnPlane* state,
    double radius, TVector3 line_point, TVector3 line_direction,
    const int pdg_code, const int direction,
    const int verbosity)
{
  assert(direction == 1 or direction == -1);

  assert(state);

  double pathlenth = std::numeric_limits<double>::quiet_NaN();

  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg_code);
  assert(rep);

  state->setRep(rep);

  try
  {
    pathlenth = rep->extrapolateToCylinder(*state, radius, line_point, line_direction);
  }
  catch (genfit::Exception& e)
  {
    if (verbosity > 1)
    {
      LogWarning("Can't extrapolate track!");
      std::cerr << e.what();
    }
    return pathlenth;
  }

  return pathlenth;
}

}  // namespace PHGenFit
