// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_Fun4AllDstPileupMerger_H
#define G4MAIN_Fun4AllDstPileupMerger_H

/*!
 * \file Fun4AllDstPileupMerger.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <map>
#include <string>

class PHCompositeNode;
class PHG4HitContainer;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;

/*!
 * utility class that can merge the relevant nodes together, once time shifted
 * in order to generate full pileup events from raw events
 * it is used internally by Fun4AllDstPileupInputManager and Fun4AllSingleDstPileupInputManager
 */
class Fun4AllDstPileupMerger final
{

  public:

  //! constructor
  Fun4AllDstPileupMerger() = default;

  //! destructor
  ~Fun4AllDstPileupMerger() = default;

  //! load destination nodes from composite
  void load_nodes(PHCompositeNode*);

  //! time-shift and copy content of source nodes to destination
  void copy_background_event(PHCompositeNode *, double delta_t) const;

  void copyDetectorActiveCrossings(const std::map<std::string,std::pair<double ,double>> &dmap) {m_DetectorTiming = dmap;}

  private:

  //! hepmc
  PHHepMCGenEventMap *m_geneventmap = nullptr;

  //! maps g4hit containers to node names
  std::map<std::string, PHG4HitContainer *> m_g4hitscontainers;

  //! truth information
  PHG4TruthInfoContainer *m_g4truthinfo = nullptr;

  std::map<std::string,std::pair<double ,double>> m_DetectorTiming;
};

#endif
