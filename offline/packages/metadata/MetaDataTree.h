#ifndef METADATATREE_METADATATREE_H
#define METADATATREE_METADATATREE_H

#include <cstdint>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

class TDirectory;
class TTree;

// Provenance of a job as a TTree with one entry; hadd appends the entries of its inputs, so a
// merged file has one entry per input job. Plain ROOT types only, readable without sPHENIX libraries.
//   provenance_version  Int_t                     layout of this tree (4)
//   cdb_globaltag       std::string               recoConsts CDB_GLOBALTAG
//   cdb_timestamp       ULong64_t                 recoConsts TIMESTAMP
//   cdb_type            std::vector<std::string>  CDB payloads: type, payload file and the timestamp it
//   cdb_url             std::vector<std::string>    was fetched for, and where it comes from:
//   cdb_url_timestamp   std::vector<ULong64_t>      cdb_from_input 0: this job used it (CDBInterface::getUrls(),
//   cdb_from_input      std::vector<Int_t>            incl. payloads read with ReadCalibrationsFromFile); listed first
//                                                   1: only in the RUN/CdbUrl node of the input DSTs (used by
//                                                     the production steps that made them, not by this job)
//   input_manager       std::vector<std::string>  every file opened by the input managers of the default
//   input_file          std::vector<std::string>    sync manager, in the order opened: manager name, file
//   input_run           std::vector<Int_t>          name as given (lfn or path), and run and segment from
//   input_segment       std::vector<Int_t>          Fun4AllUtils::GetRunSegment (run 0, segment -9999 for
//                                                   names not following the -<run>-<segment>.<ext> convention)
//
// Call it after the CDB tag and timestamp are set (InitRun or later), e.g. in End():
//   m_outfile->cd();
//   TTree *prov = MetaDataTree::Make();  // created in the current directory, i.e. m_outfile
//   prov->Write();
// Read back:  prov->Scan("cdb_globaltag:cdb_timestamp");  TChain c("provenance") over files works too.
// Check that all jobs used the same CDB payloads:  check_provenance FILE...
// The payloads per timestamp as YAML, or as a file for CDBInterface::ReadCalibrationsFromFile():
//   extract-cdb-info FILE...   /   extract-cdb-info -r RUN FILE...   /   MetaDataTree::Calibrations(file, run)
// The reading functions below (MetaDataTreeRead.cc) only need ROOT.
class MetaDataTree
{
 public:
  static constexpr int Version = 4;  // 3: release, git commit, local install and cdb_calibration_file removed
                                     // 4: cdb_from_input added

  // a new tree with this job's entry, attached to dir (nullptr: kept in memory); owned by dir, else
  // by the caller. nullptr on error (recoConsts not set)
  static TTree *Make(TDirectory *dir, const std::string &name = "provenance");

  // the same in the current directory (gDirectory)
  static TTree *Make(const std::string &name = "provenance");

  // ---- Reading (MetaDataTreeRead.cc) ----------------------------------------------------------
  // cdb_timestamp -> payload type -> payload file(s), distinct, in the order found.
  // Several files for one type happen when jobs of one run used different payloads, or when the
  // payloads carried in from the input DSTs are included (withInput).
  using Payloads = std::map<uint64_t, std::map<std::string, std::vector<std::string>>>;

  // adds the payloads the jobs used (cdb_from_input 0) of all entries of a provenance tree, with
  // withInput also those only carried in from the input DSTs; returns the number of entries, -1
  // on error. Trees older than version 4 have no cdb_from_input: all their payloads are added.
  static int Read(TTree *tree, Payloads &payloads, bool withInput = false);

  // the same for tree treename in file; -1 on error. A .yaml/.yml file is read as the YAML
  // extract-cdb-info writes (treename and withInput ignored); returns its number of timestamp blocks
  static int Read(const std::string &file, Payloads &payloads, const std::string &treename = "provenance",
                  bool withInput = false);

  // writes '<type> <url>' lines of run for CDBInterface::ReadCalibrationsFromFile(), from the
  // payloads the jobs used; a type with several payloads gets the first found (with a warning).
  // Returns the number of lines, -1 if run is not in file
  static int WriteCalibrations(std::ostream &os, const std::string &file, uint64_t run,
                               const std::string &treename = "provenance");
  static int WriteCalibrations(std::ostream &os, const Payloads &payloads, uint64_t run);

  // the same as a string, empty on error
  static std::string Calibrations(const std::string &file, uint64_t run, const std::string &treename = "provenance");
};

#endif
