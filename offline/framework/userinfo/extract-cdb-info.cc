// extract-cdb-info: the CDB payloads in the TTree UserInfo written by TreeUserInfo::Add (also
// concatenated by hadd-meta --merge-meta) as YAML, one block per cdb_timestamp:
//   51741:
//     CEMC_meanTime: /cvmfs/.../84b11f82..._CEMC_meanTime_ana509_2024p022_v001_51741.root
// A type with several payloads for one timestamp gets a list. With -r RUN it instead writes the
// '<type> <url>' lines of that run for CDBInterface::ReadCalibrationsFromFile().
// A FILE ending in .yaml/.yml is read as this YAML, so YAML files can be combined or turned into -r output.
#include "TreeUserInfo.h"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{
  void usage()
  {
    std::cerr << "usage: extract-cdb-info [-t TREE] [-r RUN] [-o OUT] FILE...\n"
              << "  -t TREE  tree to read (default: the first tree with a cdb_timestamp in each file)\n"
              << "  FILE     ROOT file, or .yaml/.yml written by extract-cdb-info\n"
              << "  -r RUN   write '<type> <url>' lines of RUN (for CDBInterface::ReadCalibrationsFromFile)\n"
              << "           instead of the YAML\n"
              << "  -o OUT   write to OUT instead of stdout" << std::endl;
  }

  // plain YAML scalar if safe, else double quoted
  std::string yaml(const std::string &s)
  {
    bool plain = !s.empty() && s[0] != '-';
    for (char c : s)
    {
      if (!isalnum(static_cast<unsigned char>(c)) && std::string("/._+-").find(c) == std::string::npos)
      {
        plain = false;
      }
    }
    if (plain)
    {
      return s;
    }
    std::string q = "\"";
    for (char c : s)
    {
      if (c == '"' || c == '\\')
      {
        q += '\\';
      }
      q += c;
    }
    return q + "\"";
  }
}  // namespace

int main(int argc, char **argv)
{
  std::string treename;
  std::string out;
  bool calib = false;
  uint64_t run = 0;
  std::vector<std::string> files;
  for (int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    if ((arg == "-t" || arg == "-r" || arg == "-o") && i + 1 < argc)
    {
      std::string val = argv[++i];
      if (arg == "-t")
      {
        treename = val;
      }
      else if (arg == "-o")
      {
        out = val;
      }
      else
      {
        char *end = nullptr;
        run = std::strtoull(val.c_str(), &end, 10);
        if (val.empty() || *end)
        {
          std::cerr << "extract-cdb-info: bad run number " << val << std::endl;
          return 1;
        }
        calib = true;
      }
    }
    else if (arg == "-h" || arg == "--help")
    {
      usage();
      return 0;
    }
    else if (!arg.empty() && arg[0] == '-')
    {
      usage();
      return 1;
    }
    else
    {
      files.push_back(arg);
    }
  }
  if (files.empty())
  {
    usage();
    return 1;
  }

  TreeUserInfo::Payloads payloads;
  for (const auto &file : files)
  {
    if (TreeUserInfo::Read(file, payloads, treename) < 0)
    {
      return 1;
    }
  }

  std::ofstream fout;
  if (!out.empty())
  {
    fout.open(out);
    if (!fout)
    {
      std::cerr << "extract-cdb-info: cannot write " << out << std::endl;
      return 1;
    }
  }
  std::ostream &os = out.empty() ? std::cout : fout;

  if (calib)
  {
    return TreeUserInfo::WriteCalibrations(os, payloads, run) < 0 ? 1 : 0;
  }

  for (const auto &[timestamp, types] : payloads)
  {
    os << timestamp << ":" << (types.empty() ? " {}" : "") << "\n";
    for (const auto &[type, urls] : types)
    {
      os << "  " << yaml(type) << ":";
      if (urls.size() == 1)
      {
        os << " " << yaml(urls.front()) << "\n";
        continue;
      }
      os << "\n";
      for (const auto &url : urls)
      {
        os << "    - " << yaml(url) << "\n";
      }
    }
  }
  return 0;
}
