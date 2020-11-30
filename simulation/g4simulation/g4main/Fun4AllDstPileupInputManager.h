// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLDSTPILEUPINPUTMANAGER_H
#define FUN4ALL_FUN4ALLDSTPILEUPINPUTMANAGER_H

/*!
 * \file Fun4AllDstPileupInputManager.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/Fun4AllInputManager.h>

#include <phool/PHCompositeNode.h>  // for PHCompositeNode
#include <phool/PHNodeIOManager.h>  // for PHNodeIOManager

#include <cstdint>  // for int64_t
#include <map>
#include <memory>
#include <string>
#include <vector>  // for vector

class EventHeader;
class PHG4HitContainer;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;
class SyncObject;

class Fun4AllDstPileupInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllDstPileupInputManager(const std::string &name = "DUMMY", const std::string &nodename = "DST", const std::string &topnodename = "TOP");
  int fileopen(const std::string &filenam);
  int fileclose();
  int run(const int nevents = 0);
  int GetSyncObject(SyncObject **mastersync);
  int SyncIt(const SyncObject *mastersync);
  int BranchSelect(const std::string &branch, const int iflag);
  int setBranches();
  virtual int setSyncBranches(PHNodeIOManager *IManager);
  void Print(const std::string &what = "ALL") const;
  int PushBackEvents(const int i);

  //! generate bunch crossing list
  /**
   * @param[in] nevents number of bunch crossing generated. This should match the number of minimum bias Hijing events to be processed
   * @param[in] collision_rate the collision rate, in Hz
   * previously set bunch crossing list is erased in the process
   */
  void generateBunchCrossingList( int nevents, float collision_rate = 5e4 );

  //! store event bunch crossing ids
  /*! bunch crossings are used to decide which pile-up events should be merged to a given "trigger" event */
  void setBunchCrossingList(const std::vector<int64_t> &value) { m_bunchCrossings = value; }

  //! get list of bunch crossing ids
  const std::vector<int64_t> &getBunchCrossingList() const { return m_bunchCrossings; }

  //! store event offset
  /*! offste is added to the current event number to look for the corresponding event timestamp */
  void setEventOffset(int value) { m_event_offset = value; }

  //! event offset
  int getEventOffset() const { return m_event_offset; }

  //! set time window for pileup events (ns)
  void setPileupTimeWindow(double tmin, double tmax)
  {
    m_tmin = tmin;
    m_tmax = tmax;
  }

 protected:
  int ReadNextEventSyncObject();
  void ReadRunTTree(const int i) { m_ReadRunTTree = i; }

 private:
  //! load nodes
  void load_nodes(PHCompositeNode *);

  //! copy background event
  void copy_background_event(PHCompositeNode *, double delta_t);

  //!@name event counters
  //@{
  int m_ReadRunTTree = 1;
  int m_ievent_total = 0;
  int m_ievent_thisfile = 0;
  int m_events_skipped_during_sync = 0;
  int m_events_accepted = 0;
  //@}

  std::string m_fullfilename;
  std::string m_RunNode = "RUN";
  std::map<const std::string, int> m_branchread;
  std::string m_syncbranchname;

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

  //! synchronization object
  SyncObject *m_syncobject = nullptr;

  //! collisions bunch crossing ids
  std::vector<int64_t> m_bunchCrossings;

  //! time between crossings. This is a RHIC constant (ns)
  static constexpr double m_time_between_crossings = 106;

  //! keep track of last accepted event bunch crossing
  /*! it is needed in order not to accept two consecutive events belonging to the same bunch crossing */
  int64_t m_last_bunchCrossing = 0;

  //! event offset
  /*! it is added to the current event number to look for the corresponding timestamp */
  int m_event_offset = 0;

  //! min integration time for pileup in the TPC (ns)
  double m_tmin = -13500;

  //! max integration time for pileup in the TPC (ns)
  double m_tmax = 13500;

  //! event header
  EventHeader *m_eventheader = nullptr;

  //! hepmc
  PHHepMCGenEventMap *m_geneventmap = nullptr;

  //! maps g4hit containers to node names
  std::map<std::string, PHG4HitContainer *> m_g4hitscontainers;

  //! truth information
  PHG4TruthInfoContainer *m_g4truthinfo = nullptr;
};

#endif /* __Fun4AllDstPileupInputManager_H__ */
