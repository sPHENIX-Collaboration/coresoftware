// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_FUN4ALLDSTPILEUPINPUTMANAGER_H
#define G4MAIN_FUN4ALLDSTPILEUPINPUTMANAGER_H

/*!
 * \file Fun4AllDstPileupInputManager.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/Fun4AllInputManager.h>

#include <phool/PHCompositeNode.h>  // for PHCompositeNode
#include <phool/PHNodeIOManager.h>  // for PHNodeIOManager

#include <gsl/gsl_rng.h>

#include <cstdint>  // for int64_t
#include <map>
#include <memory>
#include <string>
#include <vector>  // for vector

class EventHeader;
class PHG4HitContainer;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;

class Fun4AllDstPileupInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllDstPileupInputManager(const std::string &name = "DUMMY", const std::string &nodename = "DST", const std::string &topnodename = "TOP");
  int fileopen(const std::string &filenam) override;
  int fileclose() override;
  int run(const int nevents = 0) override;
  int BranchSelect(const std::string &branch, const int iflag);
  int setBranches();
  void Print(const std::string &what = "ALL") const override;
  int PushBackEvents(const int i) override;

  //! set time window for pileup events (ns)
  void setPileupTimeWindow(double tmin, double tmax)
  {
    m_tmin = tmin;
    m_tmax = tmax;
  }

  //! obsolete. Does nothing, kept for backward API compatibility.
  void generateBunchCrossingList( int nevents, float collision_rate )
  {}

 private:
  //! load nodes
  void load_nodes(PHCompositeNode *);

  //! copy background event
  void copy_background_event(PHCompositeNode *, double delta_t);

  //!@name event counters
  //@{
  bool m_ReadRunTTree = true;
  int m_ievent_total = 0;
  int m_ievent_thisfile = 0;
  int m_events_accepted = 0;
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
  static constexpr double m_time_between_crossings = 106;

  //! collision rate (Hz)
  double m_collision_rate = 5e4;

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

  //! random generator
  class Deleter
  {
    public:
    void operator() (gsl_rng* rng) const { gsl_rng_free(rng); }
  };

  std::unique_ptr<gsl_rng, Deleter> m_rng;

};

#endif /* __Fun4AllDstPileupInputManager_H__ */
