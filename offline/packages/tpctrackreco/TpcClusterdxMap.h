#ifndef FUN4ALLRAW_TPCCLUSTERDXMAP_H
#define FUN4ALLRAW_TPCCLUSTERDXMAP_H

#include <phool/PHObject.h>
#include <limits>

#include <map>

class TpcClusterdxMap : public PHObject
{
 public:
  TpcClusterdxMap() = default;
  virtual ~TpcClusterdxMap() = default;

  virtual double Add_Cluster_dx(uint64_t,double) { return std::numeric_limits<double>::quiet_NaN(); }
  virtual double get_Cluster_dx(uint64_t) { return std::numeric_limits<double>::quiet_NaN(); }
  virtual std::map<uint64_t,double> get_map() {
    std::map<uint64_t,double> m;
    return m; }

 private:
  std::map<uint64_t,double> Cluster_dx;
  ClassDefOverride(TpcClusterdxMap, 0)
};

#endif
