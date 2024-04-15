#ifndef INTTBCOMAP_H
#define INTTBCOMAP_H

#include "InttMapping.h"

#include <array>
#include <string>

class CDBTTree;

class InttBCOMap
{
 public:
  InttBCOMap();
  virtual ~InttBCOMap() {}

  virtual int LoadFromCDB(std::string const &calibname);
  virtual int LoadFromFile(std::string const &filename);

  virtual bool IsBad(int const &felix_server,
                     int const &fexlix_channel,
                     uint64_t const &bco_full,
                     const int &bco);

  virtual bool IsBad(InttNameSpace::RawData_s const &rawdata, uint64_t const &bco_full, const int &bco);
  virtual bool IsBad(InttNameSpace::Offline_s const &offline, uint64_t const &bco_full, const int &bco);

  virtual void Verbosity(const int& verbosity) { m_verbosity = verbosity; }

 protected:
  int LoadFromCDBTTree(CDBTTree &cdbttree);

 private:
  typedef std::array<std::array<int, 14>, 8> BCOArray;
  BCOArray m_bco{};  //[Felix server][Felix channel]

  int      m_verbosity{0};
};

#endif
