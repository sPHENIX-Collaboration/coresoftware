/*!
 *  \file       PHTpcSeedFinder.h
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#ifndef PHTPCSEEDFINDER_H_
#define PHTPCSEEDFINDER_H_


#include <cmath>    // for M_PI
#include <cstddef>  // for size_t
#include <vector>    // for vector

class TrkrClusterContainer;

namespace kdfinder { template <class T> class TrackCandidate; }

/// \class PHTpcSeedFinder
///
/// \brief
///
class PHTpcSeedFinder
{
 public:
  PHTpcSeedFinder();
  virtual ~PHTpcSeedFinder() {}

  void set_options(double max_distance1 = 3.0, double triplet_angle1 = M_PI / 8, size_t minhits1 = 10,
                   double max_distance2 = 6.0, double triplet_angle2 = M_PI / 8, size_t minhits2 = 5,
                   size_t nthreads = 1)
  {
    mMaxDistance1 = max_distance1;
    mTripletAngle1 = triplet_angle1;
    mMinHits1 = minhits1;
    mMaxDistance2 = max_distance2;
    mTripletAngle2 = triplet_angle2;
    mMinHits2 = minhits2;
    mNThreads = nthreads;
  }

  void set_optimization_remove_loopers(bool opt = false, double minr = 10.0, double maxr = 70.0)
  {
    mRemoveLoopers = opt;
    mMinLooperRadius = minr;
    mMaxLooperRadius = maxr;
  }

  std::vector<kdfinder::TrackCandidate<double>*> findSeeds(TrkrClusterContainer* cluster_map, double B /* magfield */);

 protected:
 private:
  double mMaxDistance1;   // 3.0 /* max distance in cm*/,
  double mTripletAngle1;  // M_PI / 8 /* triplet angle */,
  size_t mMinHits1;       // 10 /* min hits to keep track */, // first iteration
  double mMaxDistance2;   //  6.0,
  double mTripletAngle2;  // M_PI / 8,
  size_t mMinHits2;       // 5, // second iteration params
  size_t mNThreads;       // 1 /* nthreads */,

  bool mRemoveLoopers;
  double mMinLooperRadius;
  double mMaxLooperRadius;
};

#endif /* PHTPCSEEDFINDER_H_ */
