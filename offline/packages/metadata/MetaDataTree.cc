#include "MetaDataTree.h"

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/CdbUrlSave.h>

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>
#include <fun4all/Fun4AllUtils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <TDirectory.h>
#include <TTree.h>

#include <cstdint>
#include <filesystem>
#include <iostream>
#include <regex>
#include <tuple>
#include <utility>
#include <vector>

// Helpers used only in this file (an unnamed namespace makes them invisible to other files).
namespace
{
  // The list of CDB payloads this job used.
  //
  // Every module gets its calibrations through CDBInterface::instance()->getUrl(type), and
  // CDBInterface remembers each (type, url, timestamp) it handed out, also when the payloads come
  // from a file given to CDBInterface::ReadCalibrationsFromFile(). In its End() it copies
  // that list into the RUN/CdbUrl node (a CdbUrlSave object), which is how it ends up in DSTs.
  // If the job reads DSTs, RUN/CdbUrl is also read in from them, so it contains the payloads
  // of the earlier production steps too.
  //
  // Problem: CDBInterface registers itself with Fun4AllServer the first time a module asks it
  // for a payload (usually in InitRun), i.e. after the modules of the macro. Fun4AllServer calls
  // End() in registration order, so when *our* module's End() runs, CDBInterface::End() has
  // usually not run yet and RUN/CdbUrl is still missing this job's payloads.
  // Fix: call CDBInterface::End() ourselves first. Calling it again later (when Fun4AllServer
  // does) is harmless: it skips entries that are already in the node.
  // (Once CDBInterface gets a public getter for its list, this call can go away.)
  CdbUrlSave *usedPayloads()
  {
    PHCompositeNode *topNode = Fun4AllServer::instance()->topNode();
    CDBInterface::instance()->End(topNode);
    // find the RUN node under TOP, then the CdbUrl object in it
    PHNodeIterator iter(topNode);
    auto *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
    return runNode ? findNode::getClass<CdbUrlSave>(runNode, "CdbUrl") : nullptr;
  }

  // Run and segment of an input file, from the -<run>-<segment>.<ext> ending of sPHENIX file
  // names, e.g. DST_CALO_run2pp_ana509_2024p022_v001-00051741-00005.root -> run 51741, segment 5.
  // Fun4AllUtils::GetRunSegment is Fun4All's own parser for that, but it assumes the ending is
  // there: it splits the name at '-' and '.' and takes tokens off the end without checking how
  // many there are, so a name like "input.root" (only two tokens) runs it past the end of its
  // token list and can crash the job. So we check the ending first, and for names without it
  // return GetRunSegment's own "unknown" values (run 0, segment -9999) ourselves.
  std::pair<int, int> runSegment(const std::string &file)
  {
    // only the file name: a directory may contain '-' and '.' too
    std::string name = std::filesystem::path(file).filename().string();
    // anything, then -digits-digits.ext at the very end
    static const std::regex convention(R"(.*-[0-9]+-[0-9]+\.[^.]+)");
    if (!std::regex_match(name, convention))
    {
      return {0, -9999};
    }
    return Fun4AllUtils::GetRunSegment(file);
  }
}  // namespace

// Convenience overload: make the tree in whatever directory is current, which is what you get
// after m_outfile->cd(). gDirectory is ROOT's "current directory" (a file or a folder in one).
TTree *MetaDataTree::Make(const std::string &name)
{
  return Make(gDirectory, name);
}

TTree *MetaDataTree::Make(TDirectory *dir, const std::string &name)
{
  // ---- 1. Check that the job is far enough along -------------------------------------------
  // recoConsts is Fun4All's global key/value store; the macro (or the input) sets the CDB
  // global tag and the timestamp (= run number for data) there before InitRun.
  recoConsts *rc = recoConsts::instance();
  if (!rc->FlagExist("CDB_GLOBALTAG") || !rc->FlagExist("TIMESTAMP"))
  {
    std::cout << "MetaDataTree::Make: recoConsts CDB_GLOBALTAG and TIMESTAMP must be set first" << std::endl;
    return nullptr;
  }
  CdbUrlSave *urls = usedPayloads();
  if (!urls)
  {
    std::cout << "MetaDataTree::Make: no RUN/CdbUrl node with the payloads used" << std::endl;
    return nullptr;
  }

  // ---- 2. Collect the values -----------------------------------------------------------------
  // These are plain local variables; the branches below point at them while we Fill() once.
  // ROOT types (Int_t, ULong64_t) are used so the branch types are the same on every platform.

  // Layout of this tree. Bump MetaDataTree::Version when branches are added or changed:
  // hadd refuses to merge trees with different branches, so readers can check this first.
  Int_t provenance_version = Version;

  // The CDB global tag and timestamp the job ran with.
  std::string cdb_globaltag = rc->get_StringFlag("CDB_GLOBALTAG");
  ULong64_t cdb_timestamp = rc->get_uint64Flag("TIMESTAMP");

  // The payloads, as three parallel vectors: element i of each describes the same payload.
  // CdbUrlSave holds std::tuple<type, url, timestamp>; std::get<N> picks element N of a tuple.
  // The timestamp is the one the payload was fetched for, which can differ from cdb_timestamp
  // for payloads carried in from input DSTs.
  std::vector<std::string> cdb_type;
  std::vector<std::string> cdb_url;
  std::vector<ULong64_t> cdb_url_timestamp;
  for (const auto &entry : *urls)
  {
    cdb_type.push_back(std::get<0>(entry));
    cdb_url.push_back(std::get<1>(entry));
    cdb_url_timestamp.push_back(std::get<2>(entry));
  }

  // The input files, again as parallel vectors (element i = one opened file).
  //
  // Fun4AllServer -> sync manager -> input managers (one per se->registerInputManager() in the
  // macro, e.g. "DSTCALO", "DSTJET"). Each input manager keeps the list of every file it opened
  // during the job, in order, so several files read one after the other are all there.
  // getSyncManager() without a name returns the default sync manager, which is where the
  // input managers of a normal macro live; it returns nullptr if there is none.
  std::vector<std::string> input_file;
  std::vector<Int_t> input_run;
  std::vector<Int_t> input_segment;
  if (Fun4AllSyncManager *sync = Fun4AllServer::instance()->getSyncManager())
  {
    for (Fun4AllInputManager *in : sync->GetInputManagers())
    {
      for (const auto &file : in->GetFileOpenedList())
      {
        // run 0 / segment -9999 for names not ending in -<run>-<segment>.<ext> (see runSegment)
        std::pair<int, int> runseg = runSegment(file);
        input_file.push_back(file);  // as given to the input manager: the lfn for catalog DSTs, else the path
        input_run.push_back(runseg.first);
        input_segment.push_back(runseg.second);
      }
    }
  }

  // ---- 3. Build the tree ----------------------------------------------------------------------
  // A new TTree attaches itself to the current directory (gDirectory). TContext makes dir the
  // current directory for this scope only and restores the previous one when it goes out of
  // scope, so the caller's gDirectory is unchanged afterwards.
  // With dir == nullptr there is no current directory: the tree lives in memory only, and
  // SetDirectory(nullptr) makes that explicit (the caller then owns and deletes it).
  TDirectory::TContext ctx(dir);
  auto *tree = new TTree(name.c_str(), "job provenance");
  if (!dir)
  {
    tree->SetDirectory(nullptr);
  }

  // One branch per value. Branch(name, &variable) works out the type from the variable:
  // Int_t and ULong64_t become simple leaves, std::string and std::vector<...> become object
  // branches (their dictionaries come with ROOT, so plain ROOT can read them back).
  tree->Branch("provenance_version", &provenance_version);
  tree->Branch("cdb_globaltag", &cdb_globaltag);
  tree->Branch("cdb_timestamp", &cdb_timestamp);
  tree->Branch("cdb_type", &cdb_type);
  tree->Branch("cdb_url", &cdb_url);
  tree->Branch("cdb_url_timestamp", &cdb_url_timestamp);
  tree->Branch("input_file", &input_file);
  tree->Branch("input_run", &input_run);
  tree->Branch("input_segment", &input_segment);

  // Fill() copies the current values of all the variables into the tree: one entry for this
  // job. (hadd later appends the entries of the merged files: one entry per input job.)
  tree->Fill();

  // The branches still hold the addresses of the local variables above, which are destroyed
  // when this function returns. Reset them so the tree never touches that memory again
  // (e.g. if someone calls GetEntry() on it later); ROOT then uses its own buffers.
  tree->ResetBranchAddresses();
  return tree;
}
