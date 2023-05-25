#ifndef CLUSTERHITSVERBOSE_H
#define CLUSTERHITSVERBOSE_H

#include <phool/PHObject.h>
#include <vector>

class ClusterHitsVerbose : public PHObject
{
  public:
   ~ClusterHitsVerbose() override {};

  virtual std::vector<std::pair<int,float>>& get_phi_hits() = 0;
  virtual std::vector<std::pair<int,float>>& get_z_hits()  = 0;

  virtual void identify(std::ostream& os = std::cout) const override
    {
      os << "ClusterHitsVerbose" << std::endl;
    };

  protected:
  /* std::vector<std::pair<int,float>>* dummy_vec {nullptr}; */
  ClusterHitsVerbose() = default;
  ClassDefOverride(ClusterHitsVerbose, 1)
};

#endif
