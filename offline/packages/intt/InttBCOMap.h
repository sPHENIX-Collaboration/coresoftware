#ifndef INTTBCOMAP_H
#define INTTBCOMAP_H

#include <array>
#include <string>
#include <TROOT.h>
#include "InttMapping.h"

class CDBTTree;

class InttBCOMap
{
public:
  InttBCOMap() = default;
  virtual ~InttBCOMap() {}

  virtual int LoadFromCDB(std::string const &calibname); 
  virtual int LoadFromFile(std::string const &filename); 

  virtual bool IsBad(int const &felix_server,
                     int const &fexlix_channel,
                     Long64_t const &bco_full,
                     const int &bco); 

  virtual bool IsBad(InttNameSpace::RawData_s const &rawdata, Long64_t const &bco_full, const int &bco);
  virtual bool IsBad(InttNameSpace::Offline_s const &offline, Long64_t const &bco_full, const int &bco);


protected:
  int LoadFromCDBTTree(CDBTTree &cdbttree); 
private:
  typedef std::array<std::array<int, 14>, 8> BCOArray;
  BCOArray m_bco; //[Felix server][Felix channel]
};
#endif