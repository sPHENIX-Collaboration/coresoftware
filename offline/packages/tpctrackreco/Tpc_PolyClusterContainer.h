#pragma once

#include <phool/PHObject.h>

#include <iostream>

class Tpc_PolyCluster;

class Tpc_PolyClusterContainer : public PHObject
{
 public:
  Tpc_PolyClusterContainer() = default;
  ~Tpc_PolyClusterContainer() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Tpc_PolyClusterContainer base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  virtual unsigned int size() const { return 0; }
  virtual void add_cluster(Tpc_PolyCluster*) {}
  virtual const Tpc_PolyCluster* get_cluster(unsigned int) const { return nullptr; }
  virtual Tpc_PolyCluster* get_cluster(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(Tpc_PolyClusterContainer, 0)
};
