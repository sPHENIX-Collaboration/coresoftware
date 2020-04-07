#ifndef TRACKRECO_PHACTSTRKPROP_H
#define TRACKRECO_PHACTSTRKPROP_H

#include "PHTrackPropagating.h"
#include "PHActsSourceLinks.h"

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/Propagator/detail/DebugOutputActor.hpp>
#include <Acts/Propagator/detail/StandardAborters.hpp>
#include <Acts/Propagator/detail/SteppingLogger.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>

#include <memory>
#include <string>

struct ActsTrack;
class SvtxTrack;

using RecordedMaterial = Acts::MaterialInteractor::result_type;
using PropagationOutput
    = std::pair<std::vector<Acts::detail::Step>, RecordedMaterial>;

using PerigeeSurface = std::shared_ptr<const Acts::PerigeeSurface>;

class PHActsTrkProp : public PHTrackPropagating
{
 public:
  /// Default constructor
  PHActsTrkProp(const std::string& name = "PHActsTrkProp");

  /// Destructor
  ~PHActsTrkProp();

  /// End, write and close files
  int End();

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

  /// Run the propagation algorithm in acts. Returns result of propagation
  PropagationOutput propagate(FW::TrackParameters parameters);

  /// The acts geometry constructed in MakeActsGeometry
  ActsGeometry *m_actsGeometry;

  /// Minimum track pT to propagate, for acts propagator, units of GeV
  double m_minTrackPt;

  /// Propagation step size, units of mm
  double m_maxStepSize; 

  /// List of acts track fits, created by PHActsTrkFitter
  std::vector<ActsTrack>* m_actsTracks;

};

#endif
