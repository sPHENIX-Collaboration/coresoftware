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
  for (unsigned int i = 0; i < m_clusters.size(); ++i) delete m_clusters[i];
  m_clusters.clear();
}

int Tpc_PolyClusterContainerv1::isValid() const
{
  return m_clusters.empty() ? 0 : 1;
}

PHObject* Tpc_PolyClusterContainerv1::CloneMe() const
{
  Tpc_PolyClusterContainerv1* copy = new Tpc_PolyClusterContainerv1();
  for (unsigned int i = 0; i < m_clusters.size(); ++i)
  {
    if (m_clusters[i]) copy->m_clusters.push_back(static_cast<Tpc_PolyCluster*>(m_clusters[i]->CloneMe()));
  }
  return copy;
}
