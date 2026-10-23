#include "TreeUserInfo.h"

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/CdbUrlSave.h>

#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <TList.h>
#include <TNamed.h>
#include <TTree.h>

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

#ifndef TREEUSERINFO_PREFIX
#define TREEUSERINFO_PREFIX ""
#endif

namespace
{
  std::string env(const char *name)
  {
    const char *v = std::getenv(name);
    return v ? v : "";
  }

  // symlinks resolved; the path as is if it does not exist
  std::string resolved(const std::string &path)
  {
    std::error_code ec;
    std::filesystem::path p = std::filesystem::canonical(path, ec);
    return ec ? path : p.string();
  }

  // contents of file, trailing whitespace removed; the reason if it cannot be read
  std::string fileContents(const std::string &file)
  {
    std::ifstream in(file);
    if (!in)
    {
      return "cannot read " + file;
    }
    std::stringstream ss;
    ss << in.rdbuf();
    std::string s = ss.str();
    s.erase(s.find_last_not_of(" \t\r\n") + 1);
    return s;
  }

  // domain -> payload(s) fetched through CDBInterface, from RUN/CdbUrl (includes entries
  // carried in from input DSTs); a domain with several distinct payloads gets them comma separated
  bool usedPayloads(std::map<std::string, std::string> &used)
  {
    PHCompositeNode *topNode = Fun4AllServer::instance()->topNode();
    // CDBInterface fills RUN/CdbUrl in its End(), which may run after ours; repeated End is de-duplicated
    CDBInterface::instance()->End(topNode);
    PHNodeIterator iter(topNode);
    auto *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
    CdbUrlSave *urls = runNode ? findNode::getClass<CdbUrlSave>(runNode, "CdbUrl") : nullptr;
    if (!urls)
    {
      return false;
    }
    for (const auto &entry : *urls)
    {
      std::string &value = used[std::get<0>(entry)];
      const std::string &url = std::get<1>(entry);
      if (value.empty())
      {
        value = url;
      }
      else if (("," + value + ",").find("," + url + ",") == std::string::npos)
      {
        value += "," + url;
      }
    }
    return true;
  }
}  // namespace

void TreeUserInfo::Set(TTree *tree, const std::string &key, const std::string &value)
{
  TList *ui = tree->GetUserInfo();
  if (TObject *old = ui->FindObject(key.c_str()))
  {
    ui->Remove(old);
    delete old;
  }
  ui->Add(new TNamed(key.c_str(), value.c_str()));
}

int TreeUserInfo::Add(TTree *tree)
{
  if (!tree)
  {
    std::cout << "TreeUserInfo::Add: no tree" << std::endl;
    return -1;
  }
  recoConsts *rc = recoConsts::instance();
  if (!rc->FlagExist("CDB_GLOBALTAG") || !rc->FlagExist("TIMESTAMP"))
  {
    std::cout << "TreeUserInfo::Add: recoConsts CDB_GLOBALTAG and TIMESTAMP must be set first" << std::endl;
    return -1;
  }
  std::map<std::string, std::string> used;
  if (!usedPayloads(used))
  {
    std::cout << "TreeUserInfo::Add: no RUN/CdbUrl node with the payloads used" << std::endl;
    return -1;
  }
  std::string globaltag = rc->get_StringFlag("CDB_GLOBALTAG");
  uint64_t timestamp = rc->get_uint64Flag("TIMESTAMP");

  std::string offline = resolved(env("OFFLINE_MAIN"));
  std::string install = env("MYINSTALL");
  if (install.empty())
  {
    install = TREEUSERINFO_PREFIX;
  }
  std::vector<std::pair<std::string, std::string>> info = {
      {"release", std::filesystem::path(offline).filename().string()},
      {"offline_main", offline},
      {"git_commit_id", fileContents(offline + "/GIT_COMMIT_ID")},
      {"local_install", resolved(install)},
      {"cdb_globaltag", globaltag},
      {"cdb_timestamp", std::to_string(timestamp)}};

  for (const auto &kv : used)
  {
    info.emplace_back("cdb." + kv.first, kv.second);
  }
  for (const auto &kv : info)
  {
    Set(tree, kv.first, kv.second);
  }
  return static_cast<int>(info.size());
}

void TreeUserInfo::Print(TTree *tree)
{
  if (!tree)
  {
    return;
  }
  for (TObject *obj : *tree->GetUserInfo())
  {
    if (auto *n = dynamic_cast<TNamed *>(obj))
    {
      std::cout << n->GetName() << "=" << n->GetTitle() << std::endl;
    }
  }
}
