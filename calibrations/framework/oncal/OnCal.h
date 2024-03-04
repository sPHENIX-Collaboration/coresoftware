#ifndef ONCAL_ONCAL_H
#define ONCAL_ONCAL_H

#include <fun4all/SubsysReco.h>
#include <iostream>
#include <string>
#include <vector>

class Event;
class PHCompositeNode;

class OnCal : public SubsysReco
{
 public:
  virtual ~OnCal() {}

  //  These might be overwritten by everyone...
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override = 0;  // Here you analyze and commit (if committing flag is set)

  // Thsse control committing to the database...
  virtual void CommitToPdbCal(const int value) = 0;  // Set the flag for whether EndOfAnalysis will commit or not
  virtual int VerificationOK() const = 0;            // Tell us whether the new calib is close enough to the old one
  virtual int CommitedToPdbCalOK() const = 0;        // Tell us whether committing was successful by re-reading the data

  // commit without verification, needed for bootstrap calib
  // which is too different from previous calibs (e.g. begin of new Run)
  virtual void CommitNoVerify(const int) { return; }

  //  These default behaviors from SubsysReco base class
  virtual void identify(std::ostream &out = std::cout) const { out << Name() << std::endl; }
  virtual int BeginRun(const int) { return 0; }
  int EndRun(const int) override { return 0; }
  int Reset(PHCompositeNode * /*topNode*/) override { return 0; }
  int ResetEvent(PHCompositeNode * /*topNode*/) override { return 0; }
  virtual void DumpCalib() const { return; }

  unsigned int AllDone() const { return alldone; }
  void AllDone(const int i) { alldone = i; }
  void AddComment(const std::string &adcom);
  std::string Comment() const { return m_Comment; }
  int GetPdbCalTables(std::vector<std::string> &vec) const
  {
    vec = pdbcaltables;
    return 0;
  }
  virtual int CopyTables(const int FromRun, const int ToRun, const int commit) const;
  virtual int CreateCalibration(const int /*runnumber*/, const std::string & /*what*/, std::string & /*comment*/, const int /*commit*/) { return -1; }
  virtual std::vector<std::string> GetLocalFileList() const { return localfilelist; }

 protected:
  OnCal(const std::string &Name);  // so noone can call it from outside
  std::string m_Comment;
  std::vector<std::string> pdbcaltables;
  std::vector<std::string> pdbcalclasses;
  std::vector<std::pair<int, int> > bankids;
  std::vector<std::string> localfilelist;
  unsigned int alldone;
};

#endif /* __ONCAL_H__ */
