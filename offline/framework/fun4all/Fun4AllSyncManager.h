// This Master manages the reading of the Input.
// In principle all of this could be done in the Fun4AllServer 
// but this is simpler to develop and test while Fun4All is allready in use.

#ifndef FUN4ALL_FUN4ALLSYNCMANAGER_H
#define FUN4ALL_FUN4ALLSYNCMANAGER_H

#include "Fun4AllBase.h"

#include <set>
#include <vector>

class Fun4AllInputManager;
class SyncObject;

class Fun4AllSyncManager: public Fun4AllBase
{
 public:
  Fun4AllSyncManager(const std::string &name="SYNCMANAGERNONAME");
  virtual ~Fun4AllSyncManager();
  int registerInputManager(Fun4AllInputManager *InManager);
  Fun4AllInputManager *getInputManager(const std::string &name);

  //! run n events (0 means up to end of file
  int run(const int nevnts = 0);
  
  /*! 
    \brief skip n events (0 means up to the end of file). 
    Skip means read, don't process.
  */
  int skip(const int nevnts = 0);

  int fileopen(const std::string &managername = "NONE", const std::string &filename = "NONE");
  int fileclose(const std::string &managername = "NONE");
  int CurrentRun() {return currentrun;}
  void CurrentRun(const int ival) {currentrun = ival;}
  void CurrentEvent(const int evt);
  void Print(const std::string &what = "ALL") const;
  void SegmentNumber(const int iseg) {prdf_segment = iseg;}
  int SegmentNumber() const {return prdf_segment;}
  int BranchSelect(const std::string &managername, const std::string &branch, int iflag);
  int BranchSelect(const std::string &branch, const int iflag);
  int setBranches(const std::string &managername);
  int setBranches();
  void TotalEvents(const int i) {events_total = i;}
  int  TotalEvents() const {return events_total;}
  void PrdfEvents(const int i) {prdf_events = i;}
  int PrdfEvents() const {return prdf_events;}
  void GetInputFullFileList(std::vector<std::string> &fnames) const;
  void Repeat(const int i=-1) {repeat = i;}
  void PushBackInputMgrsEvents(const int i);
  int ResetEvent();
  const std::vector<Fun4AllInputManager *> GetInputManagers() const {return InManager;}

  private:
  int CheckSync(unsigned i);
  int prdf_segment;
  int prdf_events;
  int events_total;
  int currentrun;
  int currentevent;
  int repeat;
  SyncObject *MasterSync;
  std::vector<Fun4AllInputManager *> InManager;
  std::vector<int> iretInManager;
};

#endif
