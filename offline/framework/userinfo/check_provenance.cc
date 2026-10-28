// check_provenance: check that all jobs in a provenance tree (see TreeUserInfo.h) ran with the
// same software and the same CDB payloads, like hadd-meta --check-meta does for UserInfo.
//
//   check_provenance [-t TREE] [-v] FILE...
//     -t TREE  name of the provenance tree (default: provenance)
//     -v       also list the entries (file:entry) with the less common values (all but the
//              most common git_commit_id / url of a payload)
//   exit code: 0 all the same, 1 something differs, 2 error (no file, no tree, missing branch)
//
// The tree has one entry per job (hadd appends entries), and several files are read as one
// chain, so the check runs over every job in all the files given.
//
// What is checked:
//   - git_commit_id must be the same in every entry.
//   - The CDB payloads, for each payload type: every entry must have the same set of URLs.
//     A set, because one entry can hold several URLs for a type: the RUN/CdbUrl node also
//     carries the payloads of the production that made the input DSTs, which may differ from
//     the ones this job fetched (e.g. CEMC_calib_ADC_to_ETower _v1 from the DST production,
//     _v2 from the job).
//   Payloads are only compared between entries with the same cdb_timestamp (= run number for
//   data): different runs are supposed to get different run-dependent calibrations.
// release and cdb_globaltag are printed in the summary (with counts if they differ), but they do
// not fail the check: the git commit pins the code, the payload URLs pin the calibrations.
//
// On a failure each differing payload is printed with how many different URLs it had and how
// many entries had each URL ("(missing)": entries without that payload type), e.g.
//   timestamp 53756, payload MBD_PILEUP: 3 different urls in 27 entries
//          27 entries  /cvmfs/.../mbd_pileup_40000_99999.root   <- from the input DSTs, in all
//          25 entries  /cvmfs/.../mbd_pileup-53756.root
//           2 entries  /cvmfs/.../mbd_pileup-53756_v2.root

#include <TChain.h>
#include <TError.h>
#include <TFile.h>
#include <TROOT.h>

#include <cstdint>
#include <functional>  // std::greater
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unistd.h>  // getopt
#include <vector>

namespace
{
  void usage()
  {
    std::cerr << "usage: check_provenance [-t TREE] [-v] FILE...\n"
              << "  -t TREE  name of the provenance tree (default: provenance)\n"
              << "  -v       also list the entries (file:entry) with the less common values\n"
              << "exit code: 0 all the same, 1 something differs, 2 error\n";
  }

  // "file:entry" of an entry of the chain, e.g. out.root:3 (entry number within that file),
  // so differing jobs can be found again.
  std::string entryLabel(TChain &chain, Long64_t entry)
  {
    int treeNumber = chain.GetTreeNumber();
    Long64_t first = chain.GetTreeOffset()[treeNumber];  // first chain entry of the current file
    return std::string(chain.GetFile()->GetName()) + ":" + std::to_string(entry - first);
  }

  // value -> labels ("file:entry") of the entries that had it; used for release, cdb_globaltag
  // and git_commit_id, and for the URLs of each payload type
  using ValueCount = std::map<std::string, std::vector<std::string>>;

  // Print "N entries  value" lines, most frequent first. With verbose, also print the entries of
  // the odd ones out: every value except those in all `total` entries and the most frequent of
  // the rest (usually the "right" one; listing its entries would just be noise).
  void printCounts(const ValueCount &entries, size_t total, bool verbose)
  {
    std::multimap<size_t, const std::string *, std::greater<>> byCount;  // sorted largest first
    for (const auto &[value, labels] : entries)
    {
      byCount.emplace(labels.size(), &value);
    }
    bool skippedMajority = false;
    for (const auto &[n, value] : byCount)
    {
      std::cout << "    " << std::setw(5) << n << (n == 1 ? " entry    " : " entries  ")
                << (value->empty() ? "(empty)" : *value) << "\n";
      if (!verbose || n == total)
      {
        continue;
      }
      if (!skippedMajority)
      {
        skippedMajority = true;
        continue;
      }
      for (const auto &label : entries.at(*value))
      {
        std::cout << "             " << label << "\n";
      }
    }
  }
}  // namespace

int main(int argc, char **argv)
{
  // ---- 1. Arguments -----------------------------------------------------------------------------
  std::string treename = "provenance";
  bool verbose = false;
  int opt;
  while ((opt = getopt(argc, argv, "t:vh")) != -1)
  {
    switch (opt)
    {
    case 't':
      treename = optarg;
      break;
    case 'v':
      verbose = true;
      break;
    default:
      usage();
      return opt == 'h' ? 0 : 2;
    }
  }
  if (optind >= argc)
  {
    usage();
    return 2;
  }

  gROOT->SetBatch(true);
  gErrorIgnoreLevel = kError;  // no "no dictionary" style warnings, real errors still show

  // ---- 2. Open the files as one chain -----------------------------------------------------------
  // TChain::Add normally does not open the file; with 0 as second argument it opens each file
  // now and returns 0 if the file or the tree is missing, so that is caught here.
  TChain chain(treename.c_str());
  for (int i = optind; i < argc; ++i)
  {
    if (chain.Add(argv[i], 0) == 0)
    {
      std::cerr << "check_provenance: no tree '" << treename << "' in " << argv[i] << std::endl;
      return 2;
    }
  }
  Long64_t nentries = chain.GetEntries();
  if (nentries == 0)
  {
    std::cerr << "check_provenance: tree '" << treename << "' has no entries" << std::endl;
    return 2;
  }

  // ---- 3. Connect the branches we need ----------------------------------------------------------
  // Object branches (std::string, std::vector) are read through a pointer that ROOT fills;
  // starting with nullptr lets ROOT allocate the objects. Simple types are read directly.
  std::string *release = nullptr;
  std::string *git_commit_id = nullptr;
  std::string *cdb_globaltag = nullptr;
  ULong64_t cdb_timestamp = 0;
  std::vector<std::string> *cdb_type = nullptr;
  std::vector<std::string> *cdb_url = nullptr;
  // SetBranchAddress returns < 0 if the branch is missing or has another type (e.g. not a
  // provenance tree, or a different provenance_version)
  if (chain.SetBranchAddress("release", &release) < 0 ||
      chain.SetBranchAddress("git_commit_id", &git_commit_id) < 0 ||
      chain.SetBranchAddress("cdb_globaltag", &cdb_globaltag) < 0 ||
      chain.SetBranchAddress("cdb_timestamp", &cdb_timestamp) < 0 ||
      chain.SetBranchAddress("cdb_type", &cdb_type) < 0 ||
      chain.SetBranchAddress("cdb_url", &cdb_url) < 0)
  {
    std::cerr << "check_provenance: '" << treename << "' is missing provenance branches" << std::endl;
    return 2;
  }
  // only read these branches, not the input file lists etc.
  chain.SetBranchStatus("*", false);
  for (const char *b : {"release", "git_commit_id", "cdb_globaltag", "cdb_timestamp", "cdb_type", "cdb_url"})
  {
    chain.SetBranchStatus(b, true);
  }

  // ---- 4. Collect -------------------------------------------------------------------------------
  ValueCount releases;
  ValueCount commits;
  ValueCount globaltags;
  // Per timestamp: the labels of its entries, and for each payload type the URLs with the
  // entries that had them. An entry that does not have a type at all is counted under the
  // URL "(missing)" afterwards.
  struct TimestampData
  {
    std::vector<std::string> labels{};
    std::map<std::string, ValueCount> urls{};  // type -> url -> labels
  };
  
  std::map<ULong64_t, TimestampData> timestamps;

  for (Long64_t i = 0; i < nentries; ++i)
  {
    if (chain.GetEntry(i) <= 0)
    {
      std::cerr << "check_provenance: cannot read entry " << i << std::endl;
      return 2;
    }
    std::string label = entryLabel(chain, i);
    releases[*release].push_back(label);
    commits[*git_commit_id].push_back(label);
    globaltags[*cdb_globaltag].push_back(label);

    TimestampData &ts = timestamps[cdb_timestamp];
    ts.labels.push_back(label);
    // this entry's URLs per type, as a set: the same (type, url) twice in one entry counts once
    std::map<std::string, std::set<std::string>> mine;
    for (size_t j = 0; j < cdb_type->size() && j < cdb_url->size(); ++j)
    {
      mine[(*cdb_type)[j]].insert((*cdb_url)[j]);
    }
    for (const auto &[type, urls] : mine)
    {
      for (const auto &url : urls)
      {
        ts.urls[type][url].push_back(label);
      }
    }
  }

  // ---- 5. Compare -------------------------------------------------------------------------------
  bool ok = true;

  // git_commit_id: one value for all entries
  if (commits.size() > 1)
  {
    ok = false;
    std::cout << "git_commit_id: " << commits.size() << " different values in " << nentries << " entries\n";
    printCounts(commits, nentries, verbose);
  }

  // payloads: within each timestamp, every URL of a type must be in every entry
  size_t ntypes = 0;
  size_t nbad = 0;
  for (auto &[timestamp, ts] : timestamps)
  {
    size_t n = ts.labels.size();
    for (auto &[type, urls] : ts.urls)
    {
      ++ntypes;
      // Entries without this type at all: the ones in none of its URL lists.
      std::set<std::string> have;
      for (const auto &[url, labels] : urls)
      {
        have.insert(labels.begin(), labels.end());
      }
      bool same = (have.size() == n);
      for (const auto &[url, labels] : urls)
      {
        same = same && labels.size() == n;
      }
      if (same)
      {
        continue;
      }
      ok = false;
      ++nbad;
      if (have.size() < n)
      {
        auto &missing = urls["(missing)"];
        for (const auto &label : ts.labels)
        {
          if (!have.count(label))
          {
            missing.push_back(label);
          }
        }
      }
      std::cout << "timestamp " << timestamp << ", payload " << type << ": " << urls.size()
                << " different urls in " << n << " entries\n";
      printCounts(urls, n, verbose);
    }
  }

  // ---- 6. Summary -------------------------------------------------------------------------------
  // one line per value, so a mix (e.g. two releases) is visible even though it does not fail
  auto summary = [](const char *name, const ValueCount &vc)
  {
    std::cout << name << ":";
    for (const auto &[value, labels] : vc)
    {
      std::cout << " " << (value.empty() ? "(empty)" : value);
      if (vc.size() > 1)
      {
        std::cout << " (" << labels.size() << ")";
      }
    }
    std::cout << "\n";
  };
  std::cout << "checked " << nentries << " entries, " << timestamps.size() << " timestamp"
            << (timestamps.size() == 1 ? "" : "s") << ", " << ntypes << " payload types\n";
  summary("  release", releases);
  summary("  git_commit_id", commits);
  summary("  cdb_globaltag", globaltags);
  if (ok)
  {
    std::cout << "OK: same git_commit_id and CDB payloads in all entries\n";
    return 0;
  }
  std::cout << "DIFFERENT:";
  if (commits.size() > 1)
  {
    std::cout << " git_commit_id";
  }
  if (nbad)
  {
    std::cout << " " << nbad << " payload" << (nbad == 1 ? "" : "s");
  }
  std::cout << "\n";
  return 1;
}
