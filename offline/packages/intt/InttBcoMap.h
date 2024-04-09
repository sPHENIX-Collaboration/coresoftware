#ifndef INTTBCOMAP_H
#define INTTBCOMAP_H

#include "InttLoadable.h"
#include "InttMapping.h"

#include <array>
#include <string>

class CDBTTree;

class InttBcoMap : public InttLoadable
{
 public:
  static int const WIDTH = 1;

  InttBcoMap() = default;
  virtual ~InttBcoMap() override;

  using InttLoadable::LoadFromCDB;
  using InttLoadable::LoadFromFile;

  int LoadFromFile() override { return LoadFromFile("InttBcoMap.root"); }
  int LoadFromCDB() override { return LoadFromCDB("InttBcoMap"); }

  virtual bool IsBad(int const&,
                     int const&,
                     uint64_t const&,
                     int const&);

  virtual bool IsBad(InttNameSpace::RawData_s const&, uint64_t const&, const int&);
  virtual bool IsBad(InttNameSpace::Offline_s const&, uint64_t const&, const int&);

 protected:
  int LoadFromCDBTTree(CDBTTree&) override;

 private:
  typedef std::array<std::array<int, 14>, 8> bco_array_t;
  bco_array_t* m_bco{nullptr};  // [Felix server][Felix channel]
};

#endif  // INTTBCOMAP_H
