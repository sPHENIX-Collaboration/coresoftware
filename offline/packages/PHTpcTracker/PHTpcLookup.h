/*!
 *  \file       PHTpcLookup.h
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#ifndef PHTPCLOOKUP_H_
#define PHTPCLOOKUP_H_

#include "externals/kdfinder.hpp"
#include "externals/nanoflann.hpp"  // for KDTreeSingleIndexAdaptor, L2_Simp...

#include <cstddef>  // for size_t
#include <vector>

class TrkrClusterContainer;
class TrkrHitSetContainer;

/// \class PHTpcLookup
///
/// \brief
///
class PHTpcLookup
{
 public:
  PHTpcLookup();
  ~PHTpcLookup();

  void init(TrkrClusterContainer* cluster_map, TrkrHitSetContainer* hitsets);
  void clear();

  std::vector<std::vector<double>*> find(double x, double y, double z, double radius, size_t& nMatches);

 protected:

  TrkrClusterContainer* mClusterMap;
  TrkrHitSetContainer* mHitsets;
  std::vector<std::vector<double> > mKDhits;
  kdfinder::KDPointCloud<double> mCloud;
  nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, kdfinder::KDPointCloud<double> >,
                                      kdfinder::KDPointCloud<double>, 3>* mKDindex;

 private:
};

#endif /* PHTPCLOOKUP_H_ */
