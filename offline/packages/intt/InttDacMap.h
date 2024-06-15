#ifndef INTT_INTTDACMAP_H
#define INTT_INTTDACMAP_H

#include "InttMapping.h"

#include <array>
#include <string>

class CDBTTree;

class InttDacMap
{
 public:
  InttDacMap();
  virtual ~InttDacMap() {}

  virtual int LoadFromCDB(std::string const& calibname);
  virtual int LoadFromFile(std::string const& filename);
  virtual int WriteToFile(std::string const& filename);

  // Access by OnlineChannel
  virtual unsigned short GetDAC(const uint& felix_server,
                                const uint& felix_channel,
                                const uint& chip,
                                const uint& channel,
                                const uint& adc);
  virtual unsigned short GetDAC(InttNameSpace::RawData_s const& rawdata, const uint& adc);
  virtual unsigned short GetDAC(InttNameSpace::Offline_s const& offline, const uint& adc);

  virtual void Verbosity(const int& verbosity) { m_verbosity = verbosity; };

  virtual void SetDefault(const uint& Adc0 = 15,
                          const uint& Adc1 = 30,
                          const uint& Adc2 = 60,
                          const uint& Adc3 = 90,
                          const uint& Adc4 = 120,
                          const uint& Adc5 = 150,
                          const uint& Adc6 = 180,
                          const uint& Adc7 = 210);

 protected:
  int LoadFromCDBTTree(CDBTTree& cdbttree);
  void FillToCDBTTree(CDBTTree& cdbttree);

 private:
  typedef std::array<std::array<std::array<std::array<int, 8>, 26>, 14>, 8> DacArray;

  DacArray m_dac{};  // [FELIX_SERVER:8][FELIX_CHANNEL:14][CHIP:26][DAC:8]

  int m_verbosity{0};
};

#endif  // INTT_INTTDACMAP_H
