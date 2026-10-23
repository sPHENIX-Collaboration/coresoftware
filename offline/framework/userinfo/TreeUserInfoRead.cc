#include "TreeUserInfo.h"

#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TNamed.h>
#include <TTree.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <utility>

namespace
{
  const std::string cdbPrefix = "cdb.";

  // one UserInfo block: cdb_timestamp and the cdb.<type> values
  struct Block
  {
    bool hasTimestamp = false;
    uint64_t timestamp = 0;
    std::vector<std::pair<std::string, std::string>> cdb;
  };

  void addUrl(std::vector<std::string> &urls, const std::string &url)
  {
    if (!url.empty() && std::find(urls.begin(), urls.end(), url) == urls.end())
    {
      urls.push_back(url);
    }
  }

  void store(const Block &block, TreeUserInfo::Payloads &payloads)
  {
    if (!block.hasTimestamp)
    {
      if (!block.cdb.empty())
      {
        std::cerr << "TreeUserInfo::Read: UserInfo block without cdb_timestamp, skipped" << std::endl;
      }
      return;
    }
    auto &types = payloads[block.timestamp];
    for (const auto &[type, value] : block.cdb)
    {
      if (value.find_first_of(" \t\n") != std::string::npos)
      {
        continue;
      }
      auto &urls = types[type];
      std::stringstream ss(value);
      std::string url;
      while (std::getline(ss, url, ','))
      {
        addUrl(urls, url);
      }
      if (urls.empty())
      {
        types.erase(type);
      }
    }
  }

  bool isYaml(const std::string &file)
  {
    auto endsWith = [&file](const std::string &ext)
    { return file.size() >= ext.size() && file.compare(file.size() - ext.size(), ext.size(), ext) == 0; };
    return endsWith(".yaml") || endsWith(".yml");
  }

  std::string trim(const std::string &s)
  {
    size_t b = s.find_first_not_of(" \t\r");
    return b == std::string::npos ? "" : s.substr(b, s.find_last_not_of(" \t\r") - b + 1);
  }

  // a scalar as extract-cdb-info writes it: plain or double quoted
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
        ++i;
      }
      u += t[i];
    }
    return u;
  }

  // the YAML written by extract-cdb-info (not general YAML):
  //   <timestamp>:            block
  //     <type>: <url>         one payload
  //     <type>:               several payloads
  //       - <url>
  int readYaml(const std::string &file, TreeUserInfo::Payloads &payloads)
  {
    std::ifstream in(file);
    if (!in)
    {
      std::cerr << "TreeUserInfo::Read: cannot open " << file << std::endl;
      return -1;
    }
    int nblocks = 0;
    std::map<std::string, std::vector<std::string>> *types = nullptr;
    std::vector<std::string> *urls = nullptr;
    std::string line;
    for (int lineno = 1; std::getline(in, line); ++lineno)
    {
      std::string body = trim(line);
      if (body.empty() || body[0] == '#' || body == "---")
      {
        continue;
      }
      size_t colon = body.find(':');
      if (line[0] != ' ' && colon != std::string::npos)
      {
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
        addUrl(*urls, scalar(body.substr(2)));
        continue;
      }
      else if (types && colon != std::string::npos)
      {
        std::string value = scalar(body.substr(colon + 1));
        urls = &(*types)[scalar(body.substr(0, colon))];
        addUrl(*urls, value);
        continue;
      }
      std::cerr << "TreeUserInfo::Read: " << file << ":" << lineno << ": not extract-cdb-info YAML: " << line
                << std::endl;
      return -1;
    }
    return nblocks;
  }

  bool hasTimestamp(TTree *tree)
  {
    return tree->GetUserInfo()->FindObject("cdb_timestamp") != nullptr;
  }

  // first tree (highest cycle) with a cdb_timestamp in dir and its subdirectories
  TTree *findTree(TDirectory *dir)
  {
    for (TObject *obj : *dir->GetListOfKeys())
    {
      auto *key = static_cast<TKey *>(obj);
      TClass *cl = TClass::GetClass(key->GetClassName());
      if (!cl || dir->GetKey(key->GetName())->GetCycle() != key->GetCycle())
      {
        continue;
      }
      if (cl->InheritsFrom(TTree::Class()))
      {
        auto *tree = key->ReadObject<TTree>();
        if (tree && hasTimestamp(tree))
        {
          return tree;
        }
      }
      else if (cl->InheritsFrom(TDirectory::Class()))
      {
        if (TTree *tree = findTree(key->ReadObject<TDirectory>()))
        {
          return tree;
        }
      }
    }
    return nullptr;
  }
}  // namespace

int TreeUserInfo::Read(TTree *tree, Payloads &payloads)
{
  int nblocks = 0;
  Block block;
  std::set<std::string> keys;
  for (TObject *obj : *tree->GetUserInfo())
  {
    auto *n = dynamic_cast<TNamed *>(obj);
    if (!n)
    {
      continue;
    }
    std::string key = n->GetName();
    std::string value = n->GetTitle();
    if (!keys.insert(key).second)
    {
      store(block, payloads);
      ++nblocks;
      block = Block();
      keys = {key};
    }
    if (key == "cdb_timestamp")
    {
      block.hasTimestamp = true;
      block.timestamp = std::stoull(value);
    }
    else if (key.compare(0, cdbPrefix.size(), cdbPrefix) == 0)
    {
      block.cdb.emplace_back(key.substr(cdbPrefix.size()), value);
    }
  }
  if (!keys.empty())
  {
    store(block, payloads);
    ++nblocks;
  }
  return nblocks;
}

int TreeUserInfo::Read(const std::string &file, Payloads &payloads, const std::string &treename)
{
  if (isYaml(file))
  {
    return readYaml(file, payloads);
  }
  std::unique_ptr<TFile> f(TFile::Open(file.c_str()));
  if (!f || f->IsZombie())
  {
    std::cerr << "TreeUserInfo::Read: cannot open " << file << std::endl;
    return -1;
  }
  TTree *tree = treename.empty() ? findTree(f.get()) : f->Get<TTree>(treename.c_str());
  if (!tree)
  {
    std::cerr << "TreeUserInfo::Read: no " << (treename.empty() ? "tree with a cdb_timestamp" : "tree " + treename)
              << " in " << file << std::endl;
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
    if (urls.size() > 1)
    {
      std::cerr << "TreeUserInfo::WriteCalibrations: " << urls.size() << " payloads of " << type << " for " << run
                << ", using " << urls.front() << std::endl;
    }
    os << type << " " << urls.front() << "\n";
    ++n;
  }
  return n;
}

std::string TreeUserInfo::Calibrations(const std::string &file, uint64_t run, const std::string &treename)
{
  std::ostringstream os;
  return WriteCalibrations(os, file, run, treename) < 0 ? "" : os.str();
}
