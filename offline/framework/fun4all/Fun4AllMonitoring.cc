#include "Fun4AllMonitoring.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/algorithm/string.hpp>
#pragma GCC diagnostic pop

#include <boost/tokenizer.hpp>

#include <unistd.h>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <vector>

Fun4AllMonitoring *Fun4AllMonitoring::mInstance = nullptr;

Fun4AllMonitoring::Fun4AllMonitoring()
  : Fun4AllBase("Fun4AllMonitoring")
{
}

void Fun4AllMonitoring::Snapshot(const std::string & /*what*/)
{
  if (mOutFileName.empty())
  {
    return;
  }
  Get_Memory();
  if (mEvent == 0)  // called for the first time, write header
  {
    std::ofstream outfile(mOutFileName, std::ios_base::trunc);
    outfile << "Event     HeapPss     mmap   OtherPss" << std::endl;
    outfile.close();
  }
  std::ofstream outfile(mOutFileName, std::ios_base::app);
  outfile << mEvent << " " << mHeapPss << " " << mMMapPSS << " " << mOtherPss << std::endl;
  outfile.close();
  mHeapPss = 0;
  mOtherPss = 0;
  mMMapPSS = 0;
  mEvent++;
  return;
}

void Fun4AllMonitoring::Get_Memory()
{
  std::ifstream smap_stat("/proc/self/smaps");
  std::string instring;
  std::string libraryname;
  int i = 0;
  while (smap_stat)
  {
    getline(smap_stat, instring);
    if (instring.empty())
    {
      continue;
    }
    boost::trim(instring);  // remove leading + trailing spaces
    std::regex reg(R"(\s+)");
    instring = std::regex_replace(instring, reg, " ");  // replace multiple spaces with one space
    std::string firststring = instring.substr(0, instring.find(' '));
    if (firststring.find('-') != std::string::npos || (firststring.find("0000000") != std::string::npos && i == 0))
    {
      std::vector<std::string> tokens;
      boost::split(tokens, instring, boost::is_any_of(" "));

      if (tokens.size() == 6)
      {
        libraryname = tokens.back();
      }
      else
      {
        libraryname = "mmap";
      }
      i++;
    }
    else
    {
      boost::char_separator<char> sep(" ");
      using tokenizer = boost::tokenizer<boost::char_separator<char>>;
      tokenizer tok(instring, sep);
      tokenizer::iterator tok_iter = tok.begin();
      if ((*tok_iter).find("Pss") != std::string::npos)
      {
        ++tok_iter;
        std::string number = *tok_iter;
        boost::trim(number);
        if (libraryname.find("[heap]") != std::string::npos)
        {
          mHeapPss += std::stol(number);
        }
        else if (libraryname == "mmap")
        {
          mMMapPSS += std::stol(number);
        }
        else
        {
          mOtherPss += std::stol(number);
        }
      }
    }
  }
}

void Fun4AllMonitoring::PrintsMaps() const
{
  static int icnt = 0;
  std::stringstream smaps;
  smaps << "/proc/" << getpid() << "/smaps" << std::ends;
  std::ifstream smap_stat(smaps.str());
  std::string fname = "smaps.list." + std::to_string(icnt);
  std::ofstream outfile(fname);
  std::string fname1 = "smapsfilt.list." + std::to_string(icnt);
  std::ofstream outfilefilt(fname1);
  std::string key_str;
  std::string value_str;
  uint64_t vmem = 0;
  uint64_t rss = 0;
  while (smap_stat)
  {
    smap_stat >> key_str >> value_str;
    outfile << smap_stat.rdbuf();
    if (key_str == "Size:")
    {
      outfilefilt << key_str << " " << value_str << std::endl;
      vmem += std::stol(value_str);
    }
    else if (key_str == "Rss:")
    {
      rss += std::stol(value_str);
      outfilefilt << key_str << " " << value_str << std::endl;
    }
  }
  std::cout << "vmem: " << vmem << ", rss: " << rss << std::endl;
  outfile.close();
  outfilefilt.close();
  //  std::string cmd = "cat /proc/" + std::to_string(getpid()) + "/smaps";
  //  gSystem->Exec(cmd.c_str());
  icnt++;
  return;
}

void Fun4AllMonitoring::OutFileName(const std::string &fname)
{
  mOutFileName = fname;
}
