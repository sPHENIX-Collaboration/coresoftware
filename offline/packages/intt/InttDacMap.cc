#include "InttDacMap.h"

#include <cdbobjects/CDBTTree.h>

#include <filesystem>
#include <iostream>

InttDacMap::InttDacMap()
{
  InttDacMap::SetDefault();
}

int InttDacMap::WriteToFile(std::string const& filename)
{
  // std::cout<<"WriteMapToFile"<<std::endl;
  CDBTTree cdbttree = CDBTTree(filename);

  FillToCDBTTree(cdbttree);

  cdbttree.Commit();
  cdbttree.CommitSingle();
  cdbttree.WriteCDBTTree();

  // std::cout<<"WriteMapToFile completed"<<std::endl;
  return 0;
}

int InttDacMap::LoadFromCdbTTree(CDBTTree& cdbttree)
{
  uint64_t N = cdbttree.GetSingleIntValue("size");
  for (uint64_t n = 0; n < N; ++n)
  {
    int felix_server = cdbttree.GetIntValue(n, "felix_server");
    int felix_channel = cdbttree.GetIntValue(n, "felix_channel");
    int chip = cdbttree.GetIntValue(n, "chip");
    int adc = cdbttree.GetIntValue(n, "adc");
    int dac = cdbttree.GetIntValue(n, "dac");
    m_dac[felix_server][felix_channel][chip][adc] = dac;

    if (m_verbosity > 0)
    {
      std::cout << "felix_server" << felix_server << " ";
      std::cout << "felix_channel" << felix_channel << " ";
      std::cout << "chip" << chip << " ";
      std::cout << "adc" << adc << " ";
      std::cout << "dac" << dac << std::endl;
    }
  }

  return 0;
}

unsigned short InttDacMap::GetDAC(const uint& felix_server,
                                  const uint& felix_channel,
                                  const uint& chip,
                                  const uint& /*channel*/,
                                  const uint& adc)
{
  if (felix_server < m_dac.size() &&
      felix_channel < m_dac[felix_server].size() &&
      chip < m_dac[felix_server][chip].size() &&
      adc < m_dac[felix_server][chip][adc].size())
  {
    return m_dac[felix_server][felix_channel][chip][adc];
  }
  else
  {
    std::cout << "Range Invalid :"
              << " " << felix_server
              << " " << felix_channel
              << " " << chip
              << " " << adc
              << ", return -1" << std::endl;

    return -1;
  }
}

void InttDacMap::SetDefault(const uint& Adc0,
                            const uint& Adc1,
                            const uint& Adc2,
                            const uint& Adc3,
                            const uint& Adc4,
                            const uint& Adc5,
                            const uint& Adc6,
                            const uint& Adc7)
{
  // std::cout<<m_dac.size()<<std::endl;
  for (auto& felix_channel : m_dac)
  {
    // std::cout<<"\t"<<felix_channel.size()<<std::endl;

    for (auto& chip : felix_channel)
    {
      // std::cout<<"\t\t"<<chip.size()<<std::endl;

      for (auto& adc : chip)
      {
        // std::cout<<"\t\t\t"<<adc.size()<<std::endl;

        adc[0] = Adc0;
        adc[1] = Adc1;
        adc[2] = Adc2;
        adc[3] = Adc3;
        adc[4] = Adc4;
        adc[5] = Adc5;
        adc[6] = Adc6;
        adc[7] = Adc7;
      }
    }
  }
}

void InttDacMap::FillToCDBTTree(CDBTTree& cdbttree)
{
  // std::cout<<"InttDacMap::FillToCDBTTree"<<std::endl;

  // std::cout<<m_dac.size()<<std::endl;
  int n = 0;
  for (uint felix_server = 0; felix_server < m_dac.size(); felix_server++)
  {
    for (uint felix_channel = 0; felix_channel < m_dac[felix_server].size(); felix_channel++)
    {
      for (uint chip = 0; chip < m_dac[felix_server][felix_channel].size(); chip++)
      {
        for (uint adc = 0; adc < m_dac[felix_server][felix_channel][chip].size(); adc++)
        {
          int dac = m_dac[felix_server][felix_channel][chip][adc];
          cdbttree.SetIntValue(n, "felix_server", felix_server);
          cdbttree.SetIntValue(n, "felix_channel", felix_channel);
          cdbttree.SetIntValue(n, "chip", chip);
          cdbttree.SetIntValue(n, "adc", adc);
          cdbttree.SetIntValue(n, "dac", dac);
          std::cout << "felix_server" << felix_server << " ";
          std::cout << "felix_channel" << felix_channel << " ";
          std::cout << "chip" << chip << " ";
          std::cout << "adc" << adc << " ";
          std::cout << "dac" << dac << std::endl;
          n++;
        }
      }
    }
  }
  cdbttree.SetSingleIntValue("size", n);
  // std::cout<<"size "<<           n<< std::endl;

  // std::cout<<"InttDacMap::FillToCDBTTree done"<<std::endl;
}
