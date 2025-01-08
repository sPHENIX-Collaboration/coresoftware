// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOEMBEDDING_HEPMCCOLLISIONVERTEX_H
#define CALOEMBEDDING_HEPMCCOLLISIONVERTEX_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class HepMCCollisionVertex : public SubsysReco
{
 public:
  HepMCCollisionVertex(const std::string &name = "HepMCCollisionVertex");

  ~HepMCCollisionVertex() override = default;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

 private:
};

#endif  // CALOEMBEDDING_HEPMCCOLLISIONVERTEX_H
