// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_HEPMCNODEREADER_H
#define G4MAIN_HEPMCNODEREADER_H

#include <fun4all/SubsysReco.h>

// rootcint barfs with this header so we need to hide it
#include <gsl/gsl_rng.h>

#include <string>

class PHCompositeNode;

//! HepMCNodeReader take input from all subevents from PHHepMCGenEventMap and send them to simulation in Geant4
//! For HepMC subevent which is already simulated, they will not be simulated again in Geant4.
class HepMCNodeReader : public SubsysReco
{
 public:
  HepMCNodeReader(const std::string &name = "HEPMCREADER");
  ~HepMCNodeReader() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void pythia(const bool pythia)
  { is_pythia = pythia; }

  //! this function is depreciated.
  //! Embedding IDs are controlled for individually HEPMC subevents in Fun4AllHepMCInputManagers and event generators.
  void Embed(const int i = 1);

  //! this function is depreciated.
  //! HepMCNodeReader::VertexPosition() move all HEPMC subevents to a new vertex location.
  //! And the vertex shifts are better controlled for individually HEPMC subevents in Fun4AllHepMCInputManagers and event generators.
  void VertexPosition(const double v_x, const double v_y, const double v_z);

  //! HepMCNodeReader::SmearVertex - WARNING - this function is depreciated.
  //! HepMCNodeReader::SmearVertex() smear each HEPMC subevents to a new vertex location.
  //! And the vertex smears are better controlled for individually HEPMC subevents in Fun4AllHepMCInputManagers and event generators.
  //! Positive value is Gauss smear, and negative values are flat smear
  void SmearVertex(const double s_x, const double s_y, const double s_z);

  //! Arbitary time shift for all sub-events.
  //! And the vertex shifts are better controlled for individually HEPMC subevents in Fun4AllHepMCInputManagers and event generators.
  void SetT0(const double t0) { vertex_t0 = t0; }
  //
  //! Override seed
  void SetSeed(const unsigned int i)
  {
    seed = i;
    use_seed = 1;
  }

 private:
  double smeargauss(const double width);
  double smearflat(const double width);

  gsl_rng *RandomGenerator;
  bool is_pythia;
  int use_seed;
  unsigned int seed;
  double vertex_pos_x;
  double vertex_pos_y;
  double vertex_pos_z;
  double vertex_t0;
  double width_vx;
  double width_vy;
  double width_vz;
};

#endif
