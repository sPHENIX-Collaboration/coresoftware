#ifndef CLUSTERHITSVERBOSE_V1_H
#define CLUSTERHITSVERBOSE_V1_H

#include "ClusterHitsVerbose.h"

#include <phool/PHObject.h>
#include <vector>
#include <utility>
 
/* class VtxPoint; */
class ClusterHitsVerbose_v1 : public ClusterHitsVerbose
{
  public:
    ClusterHitsVerbose_v1() = default;
    ~ClusterHitsVerbose_v1() override = default;

    std::vector<std::pair<int,float>>& get_phi_hits() override;// { return m_phi_hits; };
    std::vector<std::pair<int,float>>& get_z_hits()   override;// { return m_z_hits; };

    void add_phi_hit (int i, float v) { m_phi_hits.push_back({i,v}); };
    void add_z_hit   (int i, float v) { m_z_hits.push_back({i,v});   };

  private:
    std::vector<std::pair<int,float>> m_phi_hits {};
    std::vector<std::pair<int,float>> m_z_hits {};


  public:
    // PHObject virtual overloads
    void identify(std::ostream& os = std::cout) const override
    {
      os << "ClusterHitsVerbose base class" << std::endl;
    };

  protected:
  ClassDefOverride(ClusterHitsVerbose_v1, 1)
};

#endif // CLUSTERHITSVERBOSE_V1_H
