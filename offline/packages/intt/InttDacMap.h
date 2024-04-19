#ifndef INTT_INTTDACMAP_H
#define INTT_INTTDACMAP_H

#include "InttLoadable.h"
#include "InttMap.h"

#include <array>
#include <string>

class CDBTTree;

class InttDacMap : public InttLoadable
{
 public:
  InttDacMap();
  virtual ~InttDacMap() override = default;


  virtual int WriteToFile(std::string const& filename);

  // Access by OnlineChannel
  virtual unsigned short GetDAC(const uint& felix_server,
                                const uint& felix_channel,
                                const uint& chip,
                                const uint& channel,
                                const uint& adc);
  virtual unsigned short GetDAC(InttMap::RawData_s const& rawdata, const uint& adc);

  virtual void SetDefault(const uint& Adc0 = 15,
                          const uint& Adc1 = 30,
                          const uint& Adc2 = 60,
                          const uint& Adc3 = 90,
                          const uint& Adc4 = 120,
                          const uint& Adc5 = 150,
                          const uint& Adc6 = 180,
                          const uint& Adc7 = 210);

  std::string DefaultFileName() const override { return "INTT_DACMAP.root"; }
  std::string DefaultCDBName() const override { return "INTT_DACMAP"; }

 protected:
  int LoadFromCDBTTree(CDBTTree& cdbttree) override;

  void FillToCDBTTree(CDBTTree& cdbttree);

  using InttLoadable::m_verbosity;

 private:
  typedef std::array<std::array<std::array<std::array<int, 8>, 26>, 14>, 8> DacArray;
  DacArray m_dac{};  // [FELIX_SERVER:8][FELIX_CHANNEL:14][CHIP:26][DAC:8]
};

#endif  // INTT_INTTDACMAP_H
