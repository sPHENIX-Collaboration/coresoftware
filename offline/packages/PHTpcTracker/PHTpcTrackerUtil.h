/*!
 *  \file       PHTpcTrackerUtil.h
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#ifndef PHTPCTRACKERUTIL_H_
#define PHTPCTRACKERUTIL_H_

#include <trackbase/TrkrClusterContainer.h>
#include <vector>

namespace PHTpcTrackerUtil
{
  std::vector<std::vector<double> > convert_clusters_to_hits(TrkrClusterContainer* cluster_map);

}  // namespace PHTpcTrackerUtil

#endif  // PHTPCTRACKERUTIL_H_
