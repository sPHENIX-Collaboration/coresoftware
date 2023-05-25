#include "ClusterHitsVerbose.h"

namespace
{
  std::vector<std::pair<int,float>> dummy_vec;
}

std::vector<std::pair<int,float>>& ClusterHitsVerbose::get_phi_hits() {
  return dummy_vec;
}

std::vector<std::pair<int,float>>& ClusterHitsVerbose::get_z_hits() {
  return dummy_vec;
}


