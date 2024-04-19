#include "InttDacMap.h"

#include <ffamodules/CDBInterface.h>

#include <cdbobjects/CDBTTree.h>

#include <filesystem>
#include <iostream>
#include <map>

InttDacMap::InttDacMap()
{
  InttDacMap::SetDefault();
}

int InttDacMap::WriteToFile(std::string const& filename)
{
  CDBTTree cdbttree = CDBTTree(filename);

  FillToCDBTTree(cdbttree);

  cdbttree.Commit();
  cdbttree.CommitSingle();
  cdbttree.WriteCDBTTree();

  return 0;
}

int InttDacMap::LoadFromCDBTTree(CDBTTree& cdbttree)
{
  int felix_server;
  int felix_channel;
  int chip;
  int adc;
  int dac;
  std::map<std::string, int*> fields = {
    {"felix_server",  &felix_server},
    {"felix_channel", &felix_channel},
    {"chip",          &chip},
    {"adc",           &adc},
    {"dac",           &dac},
  };
  for (int n = 0, N = cdbttree.GetSingleIntValue("size"); n < N; ++n)
  {
    for(auto& p : fields)
    {
      if((*p.second = cdbttree.GetIntValue(n, p.first)) != std::numeric_limits<int>::min())
      {
        continue;
      }

      std::cerr << __PRETTY_FUNCTION__ << "\n"
                << "\tCDBTTree::GetIntValue returned std::numeric_limits<int>::min()\n"
                << "\tname, channel combination does not exist in \"cdbttree\"\n" << std::endl;
      return 1;
    }

	try
	{
      m_dac.at(felix_server).at(felix_channel).at(chip).at(adc) = dac;
	}
	catch (std::out_of_range const& e)
	{
      std::cerr << __PRETTY_FUNCTION__ << "\n"
                << "\t" << e.what() << std::endl;
	  return 1;
	}

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
  try
  {
    return m_dac.at(felix_server).at(felix_channel).at(chip).at(adc);
  }
  catch (std::out_of_range const& e)
  {
    std::cerr << __PRETTY_FUNCTION__ << "\n"
              << "\t" << e.what() << std::endl;
    return 0xFFFFU;
  }
}

unsigned short InttDacMap::GetDAC(InttMap::RawData_s const& rawdata, const uint& adc)
{
  return GetDAC(rawdata.pid - 3001,
                rawdata.fee,
                rawdata.chp,
                rawdata.chn,
                adc);
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
