#ifndef TREEUSERINFO_TREEUSERINFO_H
#define TREEUSERINFO_TREEUSERINFO_H

#include <cstdint>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

class TTree;

// Provenance of a TTree, stored in tree->GetUserInfo() as TNamed(key, value):
//   release        basename of $OFFLINE_MAIN (e.g. ana.579)
//   offline_main   $OFFLINE_MAIN, symlinks resolved
//   git_commit_id  contents of $OFFLINE_MAIN/GIT_COMMIT_ID (or "cannot read <file>")
//   local_install  $MYINSTALL (else the prefix this library was installed to)
//   cdb_globaltag  recoConsts CDB_GLOBALTAG
//   cdb_timestamp  recoConsts TIMESTAMP
//   cdb.<type>     payload file(s) of <type>
// The cdb.<type> entries are the payloads this job fetched through CDBInterface
// (RUN/CdbUrl node, which includes those carried in from input DSTs); several distinct
// payloads of one type are comma separated. Payloads read without CDBInterface are missing.
//
// Call it after the CDB tag and timestamp are set (InitRun or later) and before the tree is
// written, e.g. in End():  TreeUserInfo::Add(m_tree);
// Read back:  tree->GetUserInfo()->FindObject("cdb.CEMC_meanTime")->GetTitle()
//
// Reading (TreeUserInfoRead.cc, ROOT only) also handles UserInfo concatenated by
// hadd-meta --merge-meta: a block starts wherever a key repeats. Blocks are combined by
// cdb_timestamp; values with whitespace (CDB error messages) are not payload files and skipped.
namespace TreeUserInfo
{
  // timestamp -> payload type -> payload file(s), distinct, in the order found
  using Payloads = std::map<uint64_t, std::map<std::string, std::vector<std::string>>>;

  // adds all entries above; returns the number added, -1 on error (nothing added)
  int Add(TTree *tree);

  // sets one entry, replacing an existing one with the same key
  void Set(TTree *tree, const std::string &key, const std::string &value);

  // prints key=value for every TNamed in the UserInfo
  void Print(TTree *tree);

  // adds the cdb.<type> entries of tree to payloads; returns the number of UserInfo blocks read
  int Read(TTree *tree, Payloads &payloads);

  // same for treename in file (empty: the first tree with a cdb_timestamp); -1 on error.
  // A .yaml/.yml file is read as the YAML extract-cdb-info writes (treename ignored); returns its number of blocks
  int Read(const std::string &file, Payloads &payloads, const std::string &treename = "");

  // writes '<type> <url>' lines of run for CDBInterface::ReadCalibrationsFromFile(); a type with
  // several payloads gets the first. Returns the number of lines, -1 if run is not in file
  int WriteCalibrations(std::ostream &os, const std::string &file, uint64_t run, const std::string &treename = "");
  int WriteCalibrations(std::ostream &os, const Payloads &payloads, uint64_t run);

  // the same as a string, empty on error
  std::string Calibrations(const std::string &file, uint64_t run, const std::string &treename = "");
}  // namespace TreeUserInfo

#endif
