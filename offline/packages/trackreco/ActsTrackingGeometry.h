#ifndef TRACKRECO_ACTSTRACKINGGEOMETRY_H
#define TRACKRECO_ACTSTRACKINGGEOMETRY_H

#include <memory>
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <ActsExamples/Plugins/BField/BFieldOptions.hpp>

/**
 * A struct to carry around Acts geometry on node tree, so as to not put 
 * all of the MakeActsGeometry tree
 */
class ActsTrackingGeometry{
 public:
  ActsTrackingGeometry(){}
  ActsTrackingGeometry(std::shared_ptr<const Acts::TrackingGeometry> tGeo,
		       ActsExamples::Options::BFieldVariant mag,
		       Acts::CalibrationContext calib,
		       Acts::GeometryContext geoCtxt,
		       Acts::MagneticFieldContext magFieldCtxt)
  : tGeometry(tGeo)
  , magField(mag)
  , calibContext(calib)
  , geoContext(geoCtxt)
  , magFieldContext(magFieldCtxt)
  {}

  /// Tracking geometry and magnetic field, for fitter function
  std::shared_ptr<const Acts::TrackingGeometry> tGeometry;
  ActsExamples::Options::BFieldVariant magField;

  /// Acts context, for Kalman options
  Acts::CalibrationContext calibContext;
  Acts::GeometryContext geoContext;
  Acts::MagneticFieldContext magFieldContext;
};


#endif
