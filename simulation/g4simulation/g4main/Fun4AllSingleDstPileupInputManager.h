// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_FUN4ALLSINGLEDSTPILEUPINPUTMANAGER_H
#define G4MAIN_FUN4ALLSINGLEDSTPILEUPINPUTMANAGER_H

/*!
 * \file Fun4AllSingleDstPileupInputManager.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h> 

#include <phool/PHCompositeNode.h>  // for PHCompositeNode
#include <phool/PHNodeIOManager.h>  // for PHNodeIOManager

#include <gsl/gsl_rng.h>

#include <map>
#include <memory>
#include <string>

class SyncObject;

/*!
 * dedicated input manager that merges single events into "merged" events, containing a trigger event
 * and a number of time-shifted pile-up events corresponding to a given pile-up rate
*/
class Fun4AllSingleDstPileupInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllSingleDstPileupInputManager(const std::string &name = "DUMMY", const std::string &nodename = "DST", const std::string &topnodename = "TOP");
  int fileopen(const std::string &filenam) override;
  int fileclose() override;
  int run(const int nevents = 0) override;
  int BranchSelect(const std::string &branch, const int iflag) override;
  int setBranches() override;
  void Print(const std::string &what = "ALL") const override;
  int PushBackEvents(const int i) override;

  // Effectivly turn off the synchronization checking (copy from Fun4AllNoSyncDstInputManager)
  int SyncIt(const SyncObject* /*mastersync*/) override { return Fun4AllReturnCodes::SYNC_OK; }
  int GetSyncObject(SyncObject** /*mastersync*/) override { return Fun4AllReturnCodes::SYNC_NOOBJECT; }
  int NoSyncPushBackEvents(const int nevt) override { return PushBackEvents(nevt); }

  /// collision rate in Hz
  void setCollisionRate(double Hz)
  { m_collision_rate = Hz; }

  /// time between bunch crossing in ns
  void setTimeBetweenCrossings(double nsec)
  { m_time_between_crossings = nsec; }

  //! set time window for pileup events (ns)
  void setPileupTimeWindow(double tmin, double tmax)
  {
    m_tmin = tmin;
    m_tmax = tmax;
  }

 private:

  //!@name event counters
  //@{
  bool m_ReadRunTTree = true;
  int m_ievent_total = 0;
  int m_ievent_thisfile = 0;
  //@}

  std::string m_fullfilename;
  std::string m_RunNode = "RUN";
  std::map<const std::string, int> m_branchread;

  //! dst node from TopNode
  PHCompositeNode *m_dstNode = nullptr;

  //! run node from TopNode
  PHCompositeNode *m_runNode = nullptr;

  //! internal dst node to copy background events
  std::unique_ptr<PHCompositeNode> m_dstNodeInternal;

  //!@name internal runnodes to perform the summation over multiple runs
  //@{
  std::unique_ptr<PHCompositeNode> m_runNodeCopy;
  std::unique_ptr<PHCompositeNode> m_runNodeSum;
  //@}

  //! input manager for active (trigger) events
  /*! corresponding nodes are copied directly to the topNode */
  std::unique_ptr<PHNodeIOManager> m_IManager;

  //! input manager for background (pileup) events
  /*! corresponding nodes are copied to the internal dst node, then merged to the top node */
  std::unique_ptr<PHNodeIOManager> m_IManager_background;

  //! time between crossings. This is a RHIC constant (ns)
  double m_time_between_crossings = 106;

  //! collision rate (Hz)
  double m_collision_rate = 5e4;

  //! min integration time for pileup in the TPC (ns)
  double m_tmin = -13500;

  //! max integration time for pileup in the TPC (ns)
  double m_tmax = 13500;

  //! random generator
  class Deleter
  {
    public:
    void operator() (gsl_rng* rng) const { gsl_rng_free(rng); }
  };

  std::unique_ptr<gsl_rng, Deleter> m_rng;

};

#endif /* __Fun4AllSingleDstPileupInputManager_H__ */
