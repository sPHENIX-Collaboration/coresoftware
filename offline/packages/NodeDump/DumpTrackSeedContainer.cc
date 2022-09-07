#include "DumpTrackSeedContainer.h"

#include <phool/PHIODataNode.h>

#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<TrackSeedContainer>;

DumpTrackSeedContainer::DumpTrackSeedContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTrackSeedContainer::process_Node(PHNode *myNode)
{
  TrackSeedContainer *trackseedcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    trackseedcontainer = thisNode->getData();
  }
  if (trackseedcontainer)
  {
    TrackSeedContainer::ConstIter hiter;
    *fout << "size: " << trackseedcontainer->size() << std::endl;
    for (hiter = trackseedcontainer->begin(); hiter != trackseedcontainer->end(); ++hiter)
    {
      if (!*hiter)
      {
        continue;
      }
      *fout << "get_pz(): " << (*hiter)->get_pz() << std::endl;
      *fout << "get_x(): " << (*hiter)->get_x() << std::endl;
      *fout << "get_y(): " << (*hiter)->get_y() << std::endl;
      *fout << "get_z(): " << (*hiter)->get_z() << std::endl;
      *fout << "get_qOverR(): " << (*hiter)->get_qOverR() << std::endl;
      *fout << "get_X0(): " << (*hiter)->get_X0() << std::endl;
      *fout << "get_Y0(): " << (*hiter)->get_Y0() << std::endl;
      *fout << "get_slope(): " << (*hiter)->get_slope() << std::endl;
      *fout << "get_Z0(): " << (*hiter)->get_Z0() << std::endl;
      *fout << "get_eta(): " << (*hiter)->get_eta() << std::endl;
      *fout << "get_theta(): " << (*hiter)->get_theta() << std::endl;
      *fout << "get_pt(): " << (*hiter)->get_pt() << std::endl;
      *fout << "get_p(): " << (*hiter)->get_p() << std::endl;
      *fout << "get_crossing(): " << (*hiter)->get_crossing() << std::endl;
      *fout << "size_cluster_keys(): " << (*hiter)->size_cluster_keys() << std::endl;
      for (TrackSeed::ConstClusterKeyIter citer = (*hiter)->begin_cluster_keys(); citer != (*hiter)->end_cluster_keys(); ++citer)
      {
        *fout << "keyid: " << *citer << std::endl;
      }
    }
  }
  return 0;
}
