#ifndef PROVENANCETREE_TREEUSERINFO_H
#define PROVENANCETREE_TREEUSERINFO_H

#include <cstdint>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

class TDirectory;
class TTree;

// Provenance of a job as a TTree with one entry; hadd appends the entries of its inputs, so a
// merged file has one entry per input job. Plain ROOT types only, readable without sPHENIX libraries.
//   provenance_version  Int_t                     layout of this tree (1)
//   release             std::string               basename of $OFFLINE_MAIN (e.g. ana.579)
//   offline_main        std::string               $OFFLINE_MAIN, symlinks resolved
//   git_commit_id       std::string               contents of $OFFLINE_MAIN/GIT_COMMIT_ID (or "cannot read <file>")
//   local_install       std::string               $MYINSTALL (else the prefix this library was installed to)
//   cdb_globaltag       std::string               recoConsts CDB_GLOBALTAG
//   cdb_timestamp       ULong64_t                 recoConsts TIMESTAMP
//   cdb_type            std::vector<std::string>  payloads fetched through CDBInterface (RUN/CdbUrl node,
//   cdb_url             std::vector<std::string>    including those carried in from input DSTs): type,
//   cdb_url_timestamp   std::vector<ULong64_t>      payload file and the timestamp it was fetched for
//   input_manager       std::vector<std::string>  every file opened by the input managers of the default
//   input_file          std::vector<std::string>    sync manager, in the order opened: manager name, file
//   input_run           std::vector<Int_t>          name as given (lfn or path), and run and segment from
//   input_segment       std::vector<Int_t>          Fun4AllUtils::GetRunSegment (segment -9999 for names
//                                                   not following the -<run>-<segment>.<ext> convention)
//
// Call it after the CDB tag and timestamp are set (InitRun or later), e.g. in End():
//   m_outfile->cd();
//   TTree *prov = TreeUserInfo::Make();  // created in the current directory, i.e. m_outfile
//   prov->Write();
// Read back:  prov->Scan("release:cdb_timestamp");  TChain c("provenance") over files works too.
// Check that all jobs used the same git commit and CDB payloads:  check_provenance FILE...
// The payloads per timestamp as YAML, or as a file for CDBInterface::ReadCalibrationsFromFile():
//   extract-cdb-info FILE...   /   extract-cdb-info -r RUN FILE...   /   TreeUserInfo::Calibrations(file, run)
// The reading functions below (TreeUserInfoRead.cc) only need ROOT.
class TreeUserInfo
{
 public:
  static constexpr int Version = 1;

  // a new tree with this job's entry, attached to dir (nullptr: kept in memory); owned by dir, else
  // by the caller. nullptr on error (recoConsts not set, no RUN/CdbUrl node)
  static TTree *Make(TDirectory *dir, const std::string &name = "provenance");

  // the same in the current directory (gDirectory)
  static TTree *Make(const std::string &name = "provenance");

  // ---- Reading (TreeUserInfoRead.cc) ----------------------------------------------------------
  // cdb_timestamp -> payload type -> payload file(s), distinct, in the order found.
  // Several files for one type happen when the input DSTs were made with another payload than
  // the one this job fetched; the job's own comes last (see WriteCalibrations).
  using Payloads = std::map<uint64_t, std::map<std::string, std::vector<std::string>>>;

  // adds the payloads of all entries of a provenance tree; returns the number of entries, -1 on error
  static int Read(TTree *tree, Payloads &payloads);

  // the same for tree treename in file; -1 on error. A .yaml/.yml file is read as the YAML
  // extract-cdb-info writes (treename ignored); returns its number of timestamp blocks
  static int Read(const std::string &file, Payloads &payloads, const std::string &treename = "provenance");

  // writes '<type> <url>' lines of run for CDBInterface::ReadCalibrationsFromFile(); a type with
  // several payloads gets the last one (the one the job itself used). Returns the number of
  // lines, -1 if run is not in file
  static int WriteCalibrations(std::ostream &os, const std::string &file, uint64_t run,
                               const std::string &treename = "provenance");
  static int WriteCalibrations(std::ostream &os, const Payloads &payloads, uint64_t run);

  // the same as a string, empty on error
  static std::string Calibrations(const std::string &file, uint64_t run, const std::string &treename = "provenance");
};

#endif
