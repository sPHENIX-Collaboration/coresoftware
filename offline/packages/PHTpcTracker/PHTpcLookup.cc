/*!
 *  \file       PHTpcLookup.cc
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */
#include "PHTpcLookup.h"
#include "PHTpcTrackerUtil.h"

#include <phool/PHLog.h>

#include <log4cpp/CategoryStream.hh>  // for CategoryStream

#include <cmath>                     // for sqrt
#include <memory>                     // for allocator_traits<>::value_type
#include <utility>                    // for pair

class TrkrClusterContainer;

PHTpcLookup::PHTpcLookup()
  : mClusterMap(nullptr)
  , mKDindex(nullptr)
{
}

PHTpcLookup::~PHTpcLookup()
{
  delete mKDindex;
}

void PHTpcLookup::init(TrkrClusterContainer* cluster_map)
{
  mClusterMap = cluster_map;

  // convert cluster_map to kdhits
  mKDhits = PHTpcTrackerUtil::convert_clusters_to_hits(mClusterMap);

  // import kdhits into KDPointCloud
  const size_t TOTALHITS = mKDhits.size();

  mCloud.pts.resize(TOTALHITS);
  for (size_t i = 0, ilen = TOTALHITS; i < ilen; i++)
  {
    mCloud.pts[i] = mKDhits[i];
  }

  // build mKDindex

  delete mKDindex;

  mKDindex = new nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, kdfinder::KDPointCloud<double> >,
                                                     kdfinder::KDPointCloud<double>, 3>(3, mCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));

  mKDindex->buildIndex();
}

void PHTpcLookup::clear()
{
  mCloud.pts.clear();
  mKDhits.clear();
}

std::vector<std::vector<double>*> PHTpcLookup::find(double x, double y, double z, double radius, size_t& nMatches)
{
  std::vector<std::vector<double>*> matches;
  if (!mKDindex)
  {
    LOG_ERROR("tracking.PHTpcLookup.find") << "using tpc lookup before init";
    return matches;
  }
  radius *= radius;  // KD-tree is using search radius^2
  nanoflann::SearchParams params;
  double query_pt[3] = {x, y, z};
  std::vector<std::pair<size_t, double> > ret_matches;
  nMatches = mKDindex->radiusSearch(&query_pt[0], radius, ret_matches, params);

  for (size_t i = 0; i < nMatches; i++)
  {
    std::vector<double>& hit = mKDhits[ret_matches[i].first];
    matches.push_back(&hit);
    double distance = std::sqrt(ret_matches[i].second);
    LOG_DEBUG("tracking.PHTpcLookup.find") << "hit x: " << hit[0] << ", y: " << hit[1] << ", z: " << hit[2] << ", dist: " << distance;
  }
  return matches;
}
