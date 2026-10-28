// Reading the CDB payloads back out of a provenance tree (see TreeUserInfo.h).
// Only ROOT is used here, so extract-cdb-info is built from this file directly and runs
// without the sPHENIX libraries.
#include "TreeUserInfo.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

// Helpers used only in this file.
namespace
{
  // Add url to the list unless it is already there (keeps the order of first appearance).
  // Empty urls are skipped: CDBInterface records "" when a type has no payload in the tag
  // (e.g. HCALINnome_default), and an empty file name is no use to ReadCalibrationsFromFile.
  // Values with whitespace are CDB error messages, not files, and are skipped too.
  void addUrl(std::vector<std::string> &urls, const std::string &url)
  {
    if (url.empty() || url.find_first_of(" \t\n") != std::string::npos)
    {
      return;
    }
    if (std::find(urls.begin(), urls.end(), url) == urls.end())
    {
      urls.push_back(url);
    }
  }

  bool isYaml(const std::string &file)
  {
    auto endsWith = [&file](const std::string &ext)
    { return file.size() >= ext.size() && file.compare(file.size() - ext.size(), ext.size(), ext) == 0; };
    return endsWith(".yaml") || endsWith(".yml");
  }

  // s without leading/trailing spaces
  std::string trim(const std::string &s)
  {
    size_t b = s.find_first_not_of(" \t\r");
    return b == std::string::npos ? "" : s.substr(b, s.find_last_not_of(" \t\r") - b + 1);
  }

  // A YAML scalar as extract-cdb-info writes it: plain, or in double quotes with \" and \\ escaped.
  std::string scalar(const std::string &s)
  {
    std::string t = trim(s);
    if (t.size() < 2 || t.front() != '"' || t.back() != '"')
    {
      return t;
    }
    std::string u;
    for (size_t i = 1; i + 1 < t.size(); ++i)
    {
      if (t[i] == '\\' && i + 2 < t.size())
      {
        ++i;  // take the escaped character as is
      }
      u += t[i];
    }
    return u;
  }

  // Reads the YAML written by extract-cdb-info (only that layout, not general YAML):
  //   <timestamp>:            a block, at the start of the line
  //     <type>: <url>         one payload
  //     <type>:               several payloads
  //       - <url>
  // Returns the number of blocks, -1 on error.
  int readYaml(const std::string &file, TreeUserInfo::Payloads &payloads)
  {
    std::ifstream in(file);
    if (!in)
    {
      std::cerr << "TreeUserInfo::Read: cannot open " << file << std::endl;
      return -1;
    }
    int nblocks = 0;
    std::map<std::string, std::vector<std::string>> *types = nullptr;  // the current block
    std::vector<std::string> *urls = nullptr;                           // the current type's list
    std::string line;
    for (int lineno = 1; std::getline(in, line); ++lineno)
    {
      std::string body = trim(line);
      if (body.empty() || body[0] == '#' || body == "---")
      {
        continue;  // blank, comment, document start
      }
      size_t colon = body.find(':');
      if (line[0] != ' ' && colon != std::string::npos)
      {
        // not indented: "<timestamp>:" starts a new block
        std::string ts = trim(body.substr(0, colon));
        if (!ts.empty() && ts.find_first_not_of("0123456789") == std::string::npos)
        {
          types = &payloads[std::stoull(ts)];
          urls = nullptr;
          ++nblocks;
          continue;
        }
      }
      else if (body.compare(0, 2, "- ") == 0 && urls)
      {
        // "- <url>": one more payload of the current type
        addUrl(*urls, scalar(body.substr(2)));
        continue;
      }
      else if (types && colon != std::string::npos)
      {
        // "<type>: <url>", or "<type>:" followed by a list
        urls = &(*types)[scalar(body.substr(0, colon))];
        addUrl(*urls, scalar(body.substr(colon + 1)));
        continue;
      }
      std::cerr << "TreeUserInfo::Read: " << file << ":" << lineno << ": not extract-cdb-info YAML: " << line
                << std::endl;
      return -1;
    }
    return nblocks;
  }
}  // namespace

int TreeUserInfo::Read(TTree *tree, Payloads &payloads)
{
  // Only the three branches we need. Object branches (std::vector) are read through a pointer
  // that ROOT fills; starting with nullptr lets ROOT allocate the vectors (deleted below).
  ULong64_t cdb_timestamp = 0;
  std::vector<std::string> *cdb_type = nullptr;
  std::vector<std::string> *cdb_url = nullptr;
  // SetBranchAddress returns < 0 if the branch is missing or has another type
  if (tree->SetBranchAddress("cdb_timestamp", &cdb_timestamp) < 0 ||
      tree->SetBranchAddress("cdb_type", &cdb_type) < 0 ||
      tree->SetBranchAddress("cdb_url", &cdb_url) < 0)
  {
    std::cerr << "TreeUserInfo::Read: " << tree->GetName() << " is not a provenance tree" << std::endl;
    tree->ResetBranchAddresses();
    return -1;
  }

  // One entry per job. The payloads are grouped by the job's cdb_timestamp (the run number
  // for data), i.e. "what was used for a job of this run", which is what a re-run needs.
  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i)
  {
    tree->GetEntry(i);
    auto &types = payloads[cdb_timestamp];
    // cdb_type[j] and cdb_url[j] describe the same payload
    for (size_t j = 0; j < cdb_type->size() && j < cdb_url->size(); ++j)
    {
      addUrl(types[(*cdb_type)[j]], (*cdb_url)[j]);
    }
  }

  // Detach the tree from our local variables before they go away, and free the vectors ROOT made.
  tree->ResetBranchAddresses();
  delete cdb_type;
  delete cdb_url;

  // drop types that only had empty urls (addUrl skipped them all)
  for (auto &[timestamp, types] : payloads)
  {
    for (auto t = types.begin(); t != types.end();)
    {
      t = t->second.empty() ? types.erase(t) : std::next(t);
    }
  }
  return static_cast<int>(nentries);
}

int TreeUserInfo::Read(const std::string &file, Payloads &payloads, const std::string &treename)
{
  if (isYaml(file))
  {
    return readYaml(file, payloads);
  }
  // unique_ptr closes the file when we return
  std::unique_ptr<TFile> f(TFile::Open(file.c_str()));
  if (!f || f->IsZombie())
  {
    std::cerr << "TreeUserInfo::Read: cannot open " << file << std::endl;
    return -1;
  }
  auto *tree = f->Get<TTree>(treename.c_str());
  if (!tree)
  {
    std::cerr << "TreeUserInfo::Read: no tree " << treename << " in " << file << std::endl;
    return -1;
  }
  return Read(tree, payloads);
}

int TreeUserInfo::WriteCalibrations(std::ostream &os, const std::string &file, uint64_t run, const std::string &treename)
{
  Payloads payloads;
  if (Read(file, payloads, treename) < 0)
  {
    return -1;
  }
  return WriteCalibrations(os, payloads, run);
}

int TreeUserInfo::WriteCalibrations(std::ostream &os, const Payloads &payloads, uint64_t run)
{
  auto it = payloads.find(run);
  if (it == payloads.end())
  {
    std::cerr << "TreeUserInfo::WriteCalibrations: no cdb_timestamp " << run << std::endl;
    return -1;
  }
  int n = 0;
  for (const auto &[type, urls] : it->second)
  {
    // Several payloads of one type: the RUN/CdbUrl node first holds the payloads carried in from
    // the input DSTs, and CDBInterface::End() appends the ones this job fetched after them. So
    // the last one is what this job itself used (e.g. CEMC_calib_ADC_to_ETower _v2 from the
    // job rather than _v1 from the DST production).
    if (urls.size() > 1)
    {
      std::cerr << "TreeUserInfo::WriteCalibrations: " << urls.size() << " payloads of " << type << " for " << run
                << ", using the last: " << urls.back() << std::endl;
    }
    os << type << " " << urls.back() << "\n";
    ++n;
  }
  return n;
}

std::string TreeUserInfo::Calibrations(const std::string &file, uint64_t run, const std::string &treename)
{
  std::ostringstream os;
  return WriteCalibrations(os, file, run, treename) < 0 ? "" : os.str();
}
