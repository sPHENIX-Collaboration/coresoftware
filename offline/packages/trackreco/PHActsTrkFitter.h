/*!
 *  \file		PHActsTrkFitter.h
 *  \brief		Refit SvtxTracks with Acts.
 *  \details	Refit SvtxTracks with Acts
 *  \author		Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TRACKRECO_ACTSTRKFITTER_H
#define TRACKRECO_ACTSTRKFITTER_H

#include "PHTrackFitting.h"

#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/BinnedArray.hpp>                       // for Binne...
#include <Acts/Utilities/Logger.hpp>                            // for getDe...


#include <TMatrixDfwd.h>                      // for TMatrixD

#include <map>
#include <memory>                // for shared_ptr
#include <string>
#include <vector>

namespace FW {
  namespace Data {
    class TrkrClusterSourceLink;
  }
}
struct ActsTrack;
struct FitCfgOptions;

using SourceLink = FW::Data::TrkrClusterSourceLink;

class PHActsTrkFitter : public PHTrackFitting
{
 public:

  //! Default constructor
  PHActsTrkFitter(const std::string& name = "PHActsTrkFitter");

  //! dtor
  ~PHActsTrkFitter();


  //!End, write and close files
  int End(PHCompositeNode*);

  int Setup(PHCompositeNode* topNode);

  int Process();

  


 private:
  /// Event counter
  int m_event;

  /// Get all the nodes
  int getNodes(PHCompositeNode*);

  /// Create new nodes
  int createNodes(PHCompositeNode*);

  /// Vector of acts tracks created by PHActsTracks
  std::vector<ActsTrack> *m_actsProtoTracks;

  /// Options that Acts::Fitter needs to run from MakeActsGeometry
  FitCfgOptions *m_fitCfgOptions;

};

#endif
