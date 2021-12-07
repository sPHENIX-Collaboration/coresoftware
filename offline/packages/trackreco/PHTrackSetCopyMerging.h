/*!
 *  \file		PHTrackSetMerging.h
 *  \brief		Base class for track container merging
 *  \author		Christof Roland <cer@mit.edu>
 */

#ifndef TRACKRECO_PHTRACKSETCOPYMERGING_H
#define TRACKRECO_PHTRACKSETCOPYMERGING_H

// PHENIX includes
#include "PHTrackSetMerging.h"

// forward declarations
class PHCompositeNode;

//class SvtxClusterMap;
class SvtxTrackMap;

/// \class PHTrackSetMerging

class PHTrackSetCopyMerging : public PHTrackSetMerging
{
 public:
  PHTrackSetCopyMerging(const std::string &name = "PHTrackSetCopyMerging");
  ~PHTrackSetCopyMerging() override {}

  int Process(PHCompositeNode *topNode) override;
 protected:

  //  virtual int Setup(PHCompositeNode *topNode) override;

  //int End() override;

  // private:

};

#endif
