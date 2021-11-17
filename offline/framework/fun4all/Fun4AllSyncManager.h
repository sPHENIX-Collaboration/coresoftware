// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// This Master manages the reading of the Input.
// In principle all of this could be done in the Fun4AllServer
// but this is simpler to develop and test while Fun4All is allready in use.

#ifndef FUN4ALL_FUN4ALLSYNCMANAGER_H
#define FUN4ALL_FUN4ALLSYNCMANAGER_H

#include "Fun4AllBase.h"

#include <string>  // for string
#include <vector>

class Fun4AllInputManager;
class SyncObject;

class Fun4AllSyncManager : public Fun4AllBase
{
 public:
  Fun4AllSyncManager(const std::string &name = "SYNCMANAGERNONAME");
  ~Fun4AllSyncManager() override;
  int registerInputManager(Fun4AllInputManager *InManager);
  Fun4AllInputManager *getInputManager(const std::string &name);

  //! run n events (0 means up to end of file
  int run(const int nevnts = 0);

  /*! 
    \brief skip n events (0 means up to the end of file). 
    Skip means read, don't process.
  */
  int skip(const int nevnts = 0);

  int fileopen(const std::string &managername, const std::string &filename);
  int fileclose(const std::string &managername = "NONE");
  int CurrentRun() { return m_CurrentRun; }
  void CurrentRun(const int ival) { m_CurrentRun = ival; }
  void CurrentEvent(const int evt);
  void Print(const std::string &what = "ALL") const override;
  void SegmentNumber(const int iseg) { m_PrdfSegment = iseg; }
  int SegmentNumber() const { return m_PrdfSegment; }
  int BranchSelect(const std::string &managername, const std::string &branch, int iflag);
  int BranchSelect(const std::string &branch, const int iflag);
  int setBranches(const std::string &managername);
  int setBranches();
  void TotalEvents(const int i) { m_EventsTotal = i; }
  int TotalEvents() const { return m_EventsTotal; }
  void PrdfEvents(const int i) { m_PrdfEvents = i; }
  int PrdfEvents() const { return m_PrdfEvents; }
  void GetInputFullFileList(std::vector<std::string> &fnames) const;
  void Repeat(const int i = -1) { m_Repeat = i; }
  void PushBackInputMgrsEvents(const int i);
  int ResetEvent();
  const std::vector<Fun4AllInputManager *> GetInputManagers() const { return m_InManager; }

 private:
  void PrintSyncProblem() const;
  int CheckSync(unsigned i);
  int m_PrdfSegment = 0;
  int m_PrdfEvents = 0;
  int m_EventsTotal = 0;
  int m_CurrentRun = 0;
  int m_CurrentEvent = 0;
  int m_Repeat = 0;
  SyncObject *m_MasterSync = nullptr;
  std::vector<Fun4AllInputManager *> m_InManager;
  std::vector<int> m_iretInManager;
};

#endif
