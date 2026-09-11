#include "DumpTpcCrossingDecisionContainer.h"

#include <phool/PHIODataNode.h>

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<TpcCrossingDecisionContainer>;

DumpTpcCrossingDecisionContainer::DumpTpcCrossingDecisionContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
}

int DumpTpcCrossingDecisionContainer::process_Node(PHNode *myNode)
{
  TpcCrossingDecisionContainer *container{nullptr};
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  if (thisNode)
  {
    container = thisNode->getData();
  }
  if (container)
  {
    unsigned int is = container->size();
    *fout << "size: " << is << std::endl;
    for (unsigned int i = 0; i < is; i++)
    {
      TpcCrossingDecision *dec = container->get_decision(i);
      *fout << "get_assembled_track_id(" << i << "): " << dec->get_assembled_track_id() << std::endl;
      *fout << "get_selected_crossing(" << i << "): " << dec->get_selected_crossing() << std::endl;
      *fout << "get_silicon_vertex_id(" << i << "): " << dec->get_silicon_vertex_id() << std::endl;
      *fout << "get_tpc_z0(" << i << "): " << dec->get_tpc_z0() << std::endl;
      *fout << "get_silicon_vertex_z(" << i << "): " << dec->get_silicon_vertex_z() << std::endl;
      *fout << "get_delta_z(" << i << "): " << dec->get_delta_z() << std::endl;
      *fout << "get_best_abs_delta_z(" << i << "): " << dec->get_best_abs_delta_z() << std::endl;
      *fout << "get_second_best_abs_delta_z(" << i << "): " << dec->get_second_best_abs_delta_z() << std::endl;
      *fout << "get_selected_tier(" << i << "): " << dec->get_selected_tier() << std::endl;
      *fout << "get_selected_score(" << i << "): " << dec->get_selected_score() << std::endl;
      *fout << "get_number_of_available_crossings(" << i << "): " << dec->get_number_of_available_crossings() << std::endl;
      *fout << "get_number_of_allowed_crossings(" << i << "): " << dec->get_number_of_allowed_crossings() << std::endl;
      *fout << "get_number_of_tested_crossings(" << i << "): " << dec->get_number_of_tested_crossings() << std::endl;
      *fout << "get_number_of_tpc_valid_crossings(" << i << "): " << dec->get_number_of_tpc_valid_crossings() << std::endl;
      *fout << "get_number_of_vertex_compatible_crossings(" << i << "): " << dec->get_number_of_vertex_compatible_crossings() << std::endl;
      unsigned int n = dec->get_number_of_candidates();
      *fout << "get_number_of_candidates(" << i << "): " << n << std::endl;
      for (unsigned int j = 0; j < n; j++)
      {
        const TpcCrossingCandidate *cand = dec->get_candidate(j);
        *fout << "crossing( " << j << "): " << cand->crossing << std::endl;
      }
      *fout << "get_status(" << i << "): " << dec->get_status() << std::endl;
    }
  }
  return 0;
}
