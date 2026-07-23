#include "Tpc_PolyClusterContainerv1.h"
#include "Tpc_PolyCluster.h"

ClassImp(Tpc_PolyClusterContainerv1)

    Tpc_PolyClusterContainerv1::Tpc_PolyClusterContainerv1()
{
  Reset();
}

Tpc_PolyClusterContainerv1::~Tpc_PolyClusterContainerv1()
{
  Reset();
}

void Tpc_PolyClusterContainerv1::identify(std::ostream& os) const
{
  os << "Tpc_PolyClusterContainerv1 with " << m_clusters.size() << " clusters" << std::endl;
}

void Tpc_PolyClusterContainerv1::Reset()
{
  for (auto& m_cluster : m_clusters)
  {
    delete m_cluster;
  }
  m_clusters.clear();
}

int Tpc_PolyClusterContainerv1::isValid() const
{
  return m_clusters.empty() ? 0 : 1;
}

PHObject* Tpc_PolyClusterContainerv1::CloneMe() const
{
  Tpc_PolyClusterContainerv1* copy = new Tpc_PolyClusterContainerv1();
  for (auto* m_cluster : m_clusters)
  {
    if (m_cluster)
    {
      copy->m_clusters.push_back(static_cast<Tpc_PolyCluster*>(m_cluster->CloneMe()));
    }
  }
  return copy;
}
