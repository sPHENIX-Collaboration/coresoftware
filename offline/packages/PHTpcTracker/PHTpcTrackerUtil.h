/*!
 *  \file       PHTpcTrackerUtil.h
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#ifndef PHTPCTRACKERUTIL_H_
#define PHTPCTRACKERUTIL_H_

#include <vector>

class TrkrClusterContainer;
class TrkrHitSetContainer;

namespace PHTpcTrackerUtil
{
  std::vector<std::vector<double> > convert_clusters_to_hits(TrkrClusterContainer* cluster_map, TrkrHitSetContainer *hitsets);

}  // namespace PHTpcTrackerUtil

#endif  // PHTPCTRACKERUTIL_H_
