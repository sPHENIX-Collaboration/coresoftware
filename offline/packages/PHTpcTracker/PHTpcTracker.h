/*!
 *  \file       PHTpcTracker.h
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#ifndef PHTPCTRACKER_H_
#define PHTPCTRACKER_H_

// PHENIX includes
#include <fun4all/SubsysReco.h>

#include <cmath>    // for M_PI
#include <cstddef>  // for size_t
#include <string>   // for string

class PHCompositeNode;
class PHField;
class TGeoManager;
class PHTpcSeedFinder;
class PHTpcTrackFollower;
class PHTpcVertexFinder;
class PHTpcEventExporter;
class PHTpcLookup;

namespace PHGenFit2
{
  class Fitter;
}

/// \class PHTpcTracker
///
/// \brief
///
class PHTpcTracker : public SubsysReco
{
 public:
  PHTpcTracker(const std::string& name = "PHTpcTracker");
  ~PHTpcTracker();

  int process_event(PHCompositeNode* topNode);

  void set_seed_finder_options(double maxdistance1 = 3.0, double tripletangle1 = M_PI / 8, size_t minhits1 = 10,
                               double maxdistance2 = 6.0, double tripletangle2 = M_PI / 8, size_t minhits2 = 5, size_t nthreads = 1);
  void set_seed_finder_optimization_remove_loopers(bool opt = false, double minr = 10.0, double maxr = 70.0);
  void set_track_follower_optimization_helix(bool opt = false);
  void set_track_follower_optimization_precise_fit(bool opt = true);

  PHTpcSeedFinder* get_seed_finder() { return mSeedFinder; }
  PHTpcTrackFollower* get_track_follower() { return mTrackFollower; }
  PHTpcVertexFinder* get_vertex_finder() { return mVertexFinder; }
  PHTpcEventExporter* get_event_exporter() { return mEventExporter; }

  void enable_vertexing(bool opt = false) { mEnableVertexing = opt; }
  void enable_json_export(bool opt = false) { mEnableJsonExport = opt; }

 protected:
  PHField* getMagField(PHCompositeNode* topNode, double& B);

 private:
  PHTpcSeedFinder* mSeedFinder;
  PHTpcTrackFollower* mTrackFollower;
  PHTpcVertexFinder* mVertexFinder;
  PHTpcEventExporter* mEventExporter;
  PHTpcLookup* mLookup;

  PHGenFit2::Fitter* mFitter;
  TGeoManager* mTGeoManager;
  PHField* mField;
  double mB;

  bool mEnableVertexing;
  bool mEnableJsonExport;
};

#endif /* __PHTPCTRACKER_H */
