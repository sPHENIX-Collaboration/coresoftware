// cdb-frozen-tag-check: check the CDB payloads recorded in a provenance tree (see MetaDataTree.h)
// against what a (frozen) global tag holds, i.e. "would a job run with TAG get the same payloads?"
//
//   cdb-frozen-tag-check [-t TREE] [-i] [-v] TAG FILE...
//     -t TREE  name of the provenance tree (default: provenance)
//     -i       also check the payloads only carried in from the input DSTs (cdb_from_input 1),
//              which the production steps that made the inputs used, not the jobs themselves
//     -v       also list the entries (file:entry) of each differing payload
//   exit code: 0 all match, 1 something differs, 2 error (no file/tree/branch, unknown tag, CDB error)
//
// Each recorded payload (cdb_type, cdb_url, cdb_url_timestamp) is looked up in TAG for the
// timestamp it was fetched for, the way CDBInterface::getUrl() would do it:
//   type X            TAG must give X the same url
//   type X_default    the job fell back to X_default because the tag it ran with had no X; so TAG
//                     must not have X either (else a job would now use X), and must give
//                     X_default the same url
// Recorded payloads with an empty url (types that had nothing in the CDB) are skipped.
// Trees older than version 4 have no cdb_from_input: all their payloads are checked, since it
// is not known which ones the jobs used.
//
// The CDB is asked once per distinct timestamp (SphenixClient::getUrlDict gives all types at
// once), so this stays fast for merged files with many entries.
//
// Output, for each differing payload: type, timestamp, the recorded url, what TAG has, in how
// many entries, e.g.
//   MBD_PILEUP @53756 (27 entries, from input DSTs)
//       recorded: /cvmfs/.../mbd_pileup_40000_99999.root
//       in tag:   /cvmfs/.../mbd_pileup-53756.root

#include <sphenixnpc/SphenixClient.h>

#include <TChain.h>
#include <TError.h>
#include <TFile.h>
#include <TROOT.h>

#include <nlohmann/json.hpp>

#include <algorithm>  // std::min
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unistd.h>  // getopt
#include <vector>

namespace
{
  void usage()
  {
    std::cerr << "usage: cdb-frozen-tag-check [-t TREE] [-i] [-v] TAG FILE...\n"
              << "  -t TREE  name of the provenance tree (default: provenance)\n"
              << "  -i       also check the payloads only carried in from the input DSTs\n"
              << "  -v       also list the entries (file:entry) of each differing payload\n"
              << "exit code: 0 all match, 1 something differs, 2 error\n";
  }

  // "file:entry" of an entry of the chain (entry number within that file)
  std::string entryLabel(TChain &chain, Long64_t entry)
  {
    int treeNumber = chain.GetTreeNumber();
    Long64_t first = chain.GetTreeOffset()[treeNumber];
    return std::string(chain.GetFile()->GetName()) + ":" + std::to_string(entry - first);
  }

  const std::string defaultSuffix = "_default";

  bool isDefault(const std::string &type)
  {
    return type.size() > defaultSuffix.size() &&
           type.compare(type.size() - defaultSuffix.size(), defaultSuffix.size(), defaultSuffix) == 0;
  }

  // One recorded payload, as seen in all the entries that had it.
  struct Recorded
  {
    std::vector<std::string> labels{};  // entries with this payload
    bool used{false};                   // used by the job itself (cdb_from_input 0) in at least one entry
  };
}  // namespace

int main(int argc, char **argv)
{
  // ---- 1. Arguments -----------------------------------------------------------------------------
  std::string treename = "provenance";
  bool withInput = false;
  bool verbose = false;
  int opt;
  while ((opt = getopt(argc, argv, "t:ivh")) != -1)
  {
    switch (opt)
    {
    case 't':
      treename = optarg;
      break;
    case 'i':
      withInput = true;
      break;
    case 'v':
      verbose = true;
      break;
    default:
      usage();
      return opt == 'h' ? 0 : 2;
    }
  }
  if (argc - optind < 2)
  {
    usage();
    return 2;
  }
  std::string tag = argv[optind++];

  gROOT->SetBatch(true);
  gErrorIgnoreLevel = kError;

  // ---- 2. Read the recorded payloads ------------------------------------------------------------
  TChain chain(treename.c_str());
  for (int i = optind; i < argc; ++i)
  {
    if (chain.Add(argv[i], 0) == 0)  // 0: open the file now, returns 0 if file or tree is missing
    {
      std::cerr << "cdb-frozen-tag-check: no tree '" << treename << "' in " << argv[i] << std::endl;
      return 2;
    }
  }
  Long64_t nentries = chain.GetEntries();
  if (nentries == 0)
  {
    std::cerr << "cdb-frozen-tag-check: tree '" << treename << "' has no entries" << std::endl;
    return 2;
  }
  std::vector<std::string> *cdb_type = nullptr;
  std::vector<std::string> *cdb_url = nullptr;
  std::vector<ULong64_t> *cdb_url_timestamp = nullptr;
  if (chain.SetBranchAddress("cdb_type", &cdb_type) < 0 ||
      chain.SetBranchAddress("cdb_url", &cdb_url) < 0 ||
      chain.SetBranchAddress("cdb_url_timestamp", &cdb_url_timestamp) < 0)
  {
    std::cerr << "cdb-frozen-tag-check: '" << treename << "' is missing cdb_type/cdb_url/cdb_url_timestamp" << std::endl;
    return 2;
  }
  // cdb_from_input (version 4 on): 1 for payloads only carried in from the input DSTs
  std::vector<Int_t> *cdb_from_input = nullptr;
  bool haveFromInput = chain.GetBranch("cdb_from_input") && chain.SetBranchAddress("cdb_from_input", &cdb_from_input) >= 0;
  if (!haveFromInput)
  {
    std::cerr << "cdb-frozen-tag-check: no cdb_from_input branch (tree older than version 4): checking all"
              << " payloads, also those only carried in from the input DSTs" << std::endl;
  }
  chain.SetBranchStatus("*", false);
  for (const char *b : {"cdb_type", "cdb_url", "cdb_url_timestamp", "cdb_from_input"})
  {
    chain.SetBranchStatus(b, true);
  }

  // (type, timestamp, url) -> entries; std::map keeps them sorted by type for the output
  std::map<std::tuple<std::string, uint64_t, std::string>, Recorded> recorded;
  for (Long64_t i = 0; i < nentries; ++i)
  {
    if (chain.GetEntry(i) <= 0)
    {
      std::cerr << "cdb-frozen-tag-check: cannot read entry " << i << std::endl;
      return 2;
    }
    std::string label = entryLabel(chain, i);
    std::set<std::tuple<std::string, uint64_t, std::string>> seenPayload;  // count each once per entry
    size_t n = std::min({cdb_type->size(), cdb_url->size(), cdb_url_timestamp->size()});
    for (size_t j = 0; j < n; ++j)
    {
      const std::string &type = (*cdb_type)[j];
      const std::string &url = (*cdb_url)[j];
      // without the branch (old tree) every payload counts as used: we cannot tell
      bool used = !haveFromInput || j >= cdb_from_input->size() || (*cdb_from_input)[j] == 0;
      if (url.empty() || (!withInput && !used))
      {
        continue;
      }
      auto key = std::make_tuple(type, static_cast<uint64_t>((*cdb_url_timestamp)[j]), url);
      if (!seenPayload.insert(key).second)
      {
        continue;
      }
      Recorded &r = recorded[key];
      r.labels.push_back(label);
      r.used = r.used || used;
    }
  }
  if (recorded.empty())
  {
    std::cerr << "cdb-frozen-tag-check: no CDB payloads recorded in the tree" << std::endl;
    return 2;
  }

  // ---- 3. Ask the CDB what TAG holds ------------------------------------------------------------
  SphenixClient client(tag);
  if (!client.existGlobalTag(tag))
  {
    std::cerr << "cdb-frozen-tag-check: global tag " << tag << " does not exist" << std::endl;
    return 2;
  }
  // timestamp -> (type -> url) in TAG, one getUrlDict call per timestamp
  std::map<uint64_t, std::map<std::string, std::string>> inTag;
  for (const auto &[key, r] : recorded)
  {
    uint64_t timestamp = std::get<1>(key);
    if (inTag.contains(timestamp))
    {
      continue;
    }
    nlohmann::json resp = client.getUrlDict(static_cast<long long>(timestamp));
    if (resp["code"] != 0)
    {
      std::cerr << "cdb-frozen-tag-check: CDB error for " << tag << " @" << timestamp << ": " << resp["msg"] << std::endl;
      return 2;
    }
    auto &urls = inTag[timestamp];
    for (const auto &[type, url] : resp["msg"].items())
    {
      urls[type] = url.get<std::string>();
    }
  }

  // ---- 4. Compare -------------------------------------------------------------------------------
  size_t nbad = 0;
  for (const auto &[key, r] : recorded)
  {
    const auto &[type, timestamp, url] = key;
    const auto &urls = inTag[timestamp];
    auto lookup = [&urls](const std::string &t)
    {
      auto it = urls.find(t);
      return it == urls.end() ? std::string() : it->second;
    };
    // what a job with TAG would get for this type, and under which type
    std::string tagType = type;
    std::string tagUrl;
    if (isDefault(type))
    {
      // recorded X_default: TAG having X itself takes precedence (CDBInterface tries X first)
      std::string base = type.substr(0, type.size() - defaultSuffix.size());
      tagUrl = lookup(base);
      if (tagUrl.empty())
      {
        tagUrl = lookup(type);
      }
      else
      {
        tagType = base;
      }
    }
    else
    {
      tagUrl = lookup(type);
    }
    if (tagType == type && tagUrl == url)
    {
      continue;
    }
    ++nbad;
    std::cout << type << " @" << timestamp << " (" << r.labels.size() << (r.labels.size() == 1 ? " entry" : " entries");
    if (!r.used)
    {
      std::cout << ", from input DSTs";
    }
    std::cout << ")\n"
              << "    recorded: " << url << "\n"
              << "    in tag:   ";
    if (tagUrl.empty())
    {
      std::cout << "(no payload)";
    }
    else
    {
      std::cout << (tagType != type ? tagType + ": " : "") << tagUrl;
    }
    std::cout << "\n";
    if (verbose)
    {
      for (const auto &label : r.labels)
      {
        std::cout << "             " << label << "\n";
      }
    }
  }

  // ---- 5. Summary -------------------------------------------------------------------------------
  std::cout << "checked " << recorded.size() << " payloads (" << nentries << " entries, " << inTag.size()
            << " timestamp" << (inTag.size() == 1 ? "" : "s") << ") against " << tag << "\n";
  if (nbad == 0)
  {
    std::cout << "OK: all payloads match " << tag << "\n";
    return 0;
  }
  std::cout << "DIFFERENT: " << nbad << " of " << recorded.size() << " payloads do not match " << tag << "\n";
  return 1;
}
