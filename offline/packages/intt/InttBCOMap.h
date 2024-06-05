#ifndef INTTBCOMAP_H
#define INTTBCOMAP_H

#include "InttLoadable.h"
#include "InttMap.h"

#include <array>
#include <string>

class CDBTTree;

class InttBCOMap : public InttLoadable
{
 public:
  InttBCOMap();
  virtual ~InttBCOMap() = default;

  virtual bool IsBad(int const &felix_server,
                     int const &fexlix_channel,
                     uint64_t const &bco_full,
                     const int &bco);

  virtual void Verbosity(const int& verbosity) { m_verbosity = verbosity; }

 protected:
  int LoadFromCdbTTree(CDBTTree &cdbttree) override;

 private:
  typedef std::array<std::array<int, 14>, 8> BCOArray;
  BCOArray m_bco{};  //[Felix server][Felix channel]

  int      m_verbosity{0};
};

#endif
