#ifndef FUN4ALLRAW_TPCCLUSTERDXMAPV1_H
#define FUN4ALLRAW_TPCCLUSTERDXMAPV1_H

#include "TpcClusterdxMap.h"

#include <map>

class TpcClusterdxMapv1 : public TpcClusterdxMap
{
 public:
  TpcClusterdxMapv1();
  ~TpcClusterdxMapv1() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  double Add_Cluster_dx(uint64_t index, double cluster_dx) override;
  double get_Cluster_dx(uint64_t index) override;
  std::map<uint64_t,double> get_map() override;

 private:
  std::map<uint64_t,double> Cluster_dx;
  ClassDefOverride(TpcClusterdxMapv1, 1)
};

#endif
