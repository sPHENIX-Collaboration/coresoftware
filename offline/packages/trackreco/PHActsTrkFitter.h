/*!
 *  \file		PHActsTrkFitter.h
 *  \brief		Refit SvtxTracks with Acts.
 *  \details	Refit SvtxTracks with Acts
 *  \author		Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TRACKRECO_ACTSTRKFITTER_H
#define TRACKRECO_ACTSTRKFITTER_H

#include "PHTrackFitting.h"
#include "PHActsSourceLinks.h"

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <memory>
#include <string>

namespace FW
{
  namespace Data
  {
    class TrkrClusterSourceLink;
  }
}  // namespace FW

class ActsTrack;
class MakeActsGeometry;

using SourceLink = FW::Data::TrkrClusterSourceLink;

class PHActsTrkFitter : public PHTrackFitting
{
 public:
  /// Default constructor
  PHActsTrkFitter(const std::string& name = "PHActsTrkFitter");

  /// Destructor
  ~PHActsTrkFitter();

  /// End, write and close files
  int End(PHCompositeNode*);

  /// Get and create nodes
  int Setup(PHCompositeNode* topNode);

  /// Process each event by calling the fitter
  int Process();

 private:
  /// Event counter
  int m_event;

  /// Get all the nodes
  int getNodes(PHCompositeNode*);

  /// Create new nodes
  int createNodes(PHCompositeNode*);

  /// Vector of acts tracks created by PHActsTracks
  std::vector<ActsTrack>* m_actsProtoTracks;

  /// Options that Acts::Fitter needs to run from MakeActsGeometry
  ActsTrackingGeometry *m_tGeometry;
};

#endif
