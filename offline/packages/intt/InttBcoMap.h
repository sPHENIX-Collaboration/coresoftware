#ifndef INTTBCOMAP_H
#define INTTBCOMAP_H

#include "InttLoadable.h"
#include "InttMap.h"

#include <array>
#include <string>

class CDBTTree;

class InttBcoMap : public InttLoadable
{
 public:
  static int const WIDTH = 1;

  InttBcoMap();
  virtual ~InttBcoMap() override = default;

  virtual bool IsBad(int const&,
                     int const&,
                     uint64_t const&,
                     int const&);

  virtual bool IsBad(InttMap::RawData_s const&, uint64_t const&, const int&);

  std::string DefaultFileName() const override { return "INTT_BCOMAP.root"; }
  std::string DefaultCDBName() const override { return "INTT_BCOMAP"; }

 protected:
  int LoadFromCDBTTree(CDBTTree&) override;

 private:
  typedef std::array<std::array<int, 14>, 8> bco_array_t;
  bco_array_t m_bco;  // [Felix server][Felix channel]
};

#endif  // INTTBCOMAP_H
