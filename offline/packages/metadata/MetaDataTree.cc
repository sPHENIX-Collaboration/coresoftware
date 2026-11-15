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
  // The RUN/CdbUrl node (a CdbUrlSave object), nullptr if there is none.
  //
  // When the job reads DSTs, RUN/CdbUrl is read in from them: the payloads used by the earlier
  // production steps that made those DSTs. In its End(), CDBInterface appends this job's own
  // payloads to it (that is how they get into DSTs this job writes). Whether that already
  // happened when Make() runs depends on the order of End() calls, so the node is only used
  // for the payloads carried in from the inputs; this job's own come from
  // CDBInterface::getUrls().
  CdbUrlSave *runNodePayloads()
  {
    PHNodeIterator iter(Fun4AllServer::instance()->topNode());
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
  // ---- 2. Collect the values -----------------------------------------------------------------
  // These are plain local variables; the branches below point at them while we Fill() once.
  // ROOT types (Int_t, ULong64_t) are used so the branch types are the same on every platform.

  // Layout of this tree. Bump MetaDataTree::Version when branches are added or changed:
  // hadd refuses to merge trees with different branches, so readers can check this first.
  Int_t provenance_version = Version;

  // The CDB global tag and timestamp the job ran with.
  std::string cdb_globaltag = rc->get_StringFlag("CDB_GLOBALTAG");
  ULong64_t cdb_timestamp = rc->get_uint64Flag("TIMESTAMP");

  // The payloads, as four parallel vectors: element i of each describes the same payload.
  // (type, url, timestamp) tuples; std::get<N> picks element N of a tuple. The timestamp is the
  // one the payload was fetched for, which can differ from cdb_timestamp for payloads carried
  // in from input DSTs.
  //   cdb_from_input 0: this job got it from CDBInterface (from the CDB, or from a file given to
  //                     ReadCalibrationsFromFile), i.e. the job used it. These come first.
  //                  1: only in the RUN/CdbUrl node read in from the input DSTs: used by the
  //                     production steps that made them, not fetched by this job.
  // A payload that is both (the job fetched what its inputs already had) is listed once, as 0.
  std::vector<std::string> cdb_type;
  std::vector<std::string> cdb_url;
  std::vector<ULong64_t> cdb_url_timestamp;
  std::vector<Int_t> cdb_from_input;
  auto addPayload = [&](const std::tuple<std::string, std::string, uint64_t> &payload, bool fromInput)
  {
    cdb_type.push_back(std::get<0>(payload));
    cdb_url.push_back(std::get<1>(payload));
    cdb_url_timestamp.push_back(std::get<2>(payload));
    cdb_from_input.push_back(fromInput ? 1 : 0);
  };
  const auto &ownPayloads = CDBInterface::instance()->getUrls();
  for (const auto &payload : ownPayloads)
  {
    addPayload(payload, false);
  }
  if (CdbUrlSave *node = runNodePayloads())
  {
    for (const auto &payload : *node)
    {
      // the job's own are already listed (they are in the node too once CDBInterface::End() ran)
      if (!ownPayloads.contains(payload))
      {
        addPayload(payload, true);
      }
    }
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
  tree->Branch("cdb_from_input", &cdb_from_input);
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
