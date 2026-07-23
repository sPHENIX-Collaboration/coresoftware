#pragma once

#include "Tpc_PolyClusterContainer.h"

#include <iostream>
#include <vector>

class Tpc_PolyCluster;

class Tpc_PolyClusterContainerv1 : public Tpc_PolyClusterContainer
{
 public:
  Tpc_PolyClusterContainerv1();
  ~Tpc_PolyClusterContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_clusters.size()); }
  void add_cluster(Tpc_PolyCluster* cluster) override { m_clusters.push_back(cluster); }
  const Tpc_PolyCluster* get_cluster(unsigned int i) const override
  {
    if (i >= m_clusters.size()) return nullptr;
    return m_clusters[i];
  }
  Tpc_PolyCluster* get_cluster(unsigned int i) override
  {
    if (i >= m_clusters.size()) return nullptr;
    return m_clusters[i];
  }

 private:
  std::vector<Tpc_PolyCluster*> m_clusters;

  ClassDefOverride(Tpc_PolyClusterContainerv1, 1)
};
