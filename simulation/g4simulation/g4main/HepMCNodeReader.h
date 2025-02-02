// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_HEPMCNODEREADER_H
#define G4MAIN_HEPMCNODEREADER_H

#include <TF1.h>

#include <fun4all/SubsysReco.h>

// rootcint barfs with this header so we need to hide it
#include <gsl/gsl_rng.h>

#include <string>
#include <vector>

class PHCompositeNode;

//! HepMCNodeReader take input from all subevents from PHHepMCGenEventMap and send them to simulation in Geant4
//! For HepMC subevent which is already simulated, they will not be simulated again in Geant4.
class HepMCNodeReader : public SubsysReco
{
 public:
  HepMCNodeReader(const std::string &name = "HepMCNodeReader");
  ~HepMCNodeReader() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void pythia(const bool pythia) { is_pythia = pythia; }

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

  // Method to add strangeness content of the event: f is the fraction of additional strangeness particles, in percent
  void AddStrangeness(const float f) { addfraction = f; }

 private:
  double smeargauss(const double width);
  double smearflat(const double width);

  // Exponentially modified Gaussian distribution function, modelling pT
  static double EMGFunction(double *x, double *par);
  // Double Gaussian function, modelling eta
  static double DBGFunction(double *x, double *par);

  gsl_rng *RandomGenerator{nullptr};
  bool is_pythia{false};
  int use_seed{0};
  unsigned int seed{0};
  double vertex_pos_x{0.0};
  double vertex_pos_y{0.0};
  double vertex_pos_z{0.0};
  double vertex_t0{0.0};
  double width_vx{0.0};
  double width_vy{0.0};
  double width_vz{0.0};

  // Method to change the strangeness content of the event
  std::vector<int> list_strangePID = {310, 3122, -3122};                               // K_s0 (PID=310), Lambda (PID=3122)
  std::vector<double> list_strangePIDprob = {1 - 0.333, 0.333 / 2., 0.333 / 2.};       // K_s0 (PID=310) has 2/3 probability, Lambda (PID=3122) has 1/3 probability
  std::vector<std::pair<int, std::pair<double, double>>> list_strangePID_probrange{};  // list of strange particles and their probability ranges
  float addfraction{0.0};                                                              // additional strangeness particles, in percent; default 0
  int Nstrange_add{0};                                                                 // number of strange particles to be added; default 0
  TF1 *fpt{nullptr};
  TF1 *feta{nullptr};
};

#endif
