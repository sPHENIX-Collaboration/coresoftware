/*!
 *  \file		Track.h
 *  \brief		tools
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef PHGENFIT_TOOLS_H
#define PHGENFIT_TOOLS_H

//STL
#include <memory>
#include <vector>

namespace genfit
{
class AbsTrackRep;
class StateOnPlane;
class Track;
class MeasuredStateOnPlane;

}  // namespace genfit

namespace PHGenFit
{
double extrapolateToCylinder(
    const genfit::MeasuredStateOnPlane* state,
    double radius, TVector3 line_point, TVector3 line_direction,
    const int pdg_code = 211, const int direction = 1,
    const int verbosity = 0);

}  // namespace PHGenFit

#endif  //__PHGenFit_Tools__
