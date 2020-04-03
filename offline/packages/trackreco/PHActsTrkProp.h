
#ifndef TRACKRECO_ACTSTRKPROP_H
#define TRACKRECO_ACTSTRKPROP_H

#include "PHTrackPropagating.h"
#include "PHActsSourceLinks.h"

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <memory>
#include <string>

namespace FW
{
  namespace Data
  {
    class TrkrClusterSourceLink;
  }
}

class SvtxTrack;

using PerigeeSurface = std::shared_ptr<const Acts::PerigeeSurface>;
using SourceLink = FW::Data::TrkrClusterSourceLink;

class PHActsTrkProp : public PHTrackPropagating
{
 public:
  /// Default constructor
  PHActsTrkProp(const std::string& name = "PHActsTrkProp");

  /// Destructor
  ~PHActsTrkProp();

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

  Acts::BoundSymMatrix getActsCovMatrix(const SvtxTrack *track);

  FitCfgOptions *m_fitCfgOptions;

};

#endif
