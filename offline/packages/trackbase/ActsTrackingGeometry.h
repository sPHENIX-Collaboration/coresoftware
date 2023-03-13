#ifndef TRACKBASE_ACTSTRACKINGGEOMETRY_H
#define TRACKBASE_ACTSTRACKINGGEOMETRY_H

#include <memory>
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

/**
 * A struct to carry around Acts geometry on node tree, so as to not put 
 * all of the MakeActsGeometry tree
 */
struct ActsTrackingGeometry{
  ActsTrackingGeometry(){}
  ActsTrackingGeometry(std::shared_ptr<const Acts::TrackingGeometry> tGeo,
		       //ActsTrackingGeometry(std::shared_ptr<Acts::TrackingGeometry> tGeo,
		       std::shared_ptr<Acts::MagneticFieldProvider> mag,
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
  //std::shared_ptr<Acts::TrackingGeometry> tGeometry;

  std::shared_ptr<Acts::MagneticFieldProvider> magField;

  /// Acts context, for Kalman options
  Acts::CalibrationContext calibContext;
  Acts::GeometryContext geoContext;
  Acts::MagneticFieldContext magFieldContext;

  Acts::GeometryContext& getGeoContext() 
  {
    return geoContext;
  }

  double tpcSurfStepPhi = 0;
  double tpcSurfStepZ = 0;
};


#endif
