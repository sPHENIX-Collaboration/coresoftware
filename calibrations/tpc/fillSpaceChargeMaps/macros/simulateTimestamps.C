#include <TFile.h>
#include <TTree.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <Rtypes.h>

#include <format>
#include <fstream>
#include <iostream>

R__LOAD_LIBRARY(libgslcblas.so)
R__LOAD_LIBRARY(libgsl.so)

//_______________________________________________________________
void simulateTimestamps()
{
  // random generator
  auto *rng = gsl_rng_alloc(gsl_rng_mt19937);
  // gsl_rng_set( rng, time(0) );//time(0) is the seed
  gsl_rng_set(rng, 1);

  // time interval (s) between two bunch crossing
  /* value copied from generators/phhepmc/Fun4AllHepMCPileupInputManager.cc */
  static constexpr double deltat_crossing = 106e-9;

  // triggered rate (Hz)
  // static constexpr double trigger_rate = 1e5; // Hijing 100kHz
  static constexpr double trigger_rate = 5e4;  // Hijing 50kHz
  // static constexpr double trigger_rate = 3e6; // pp collisions

  // mean number of collision per crossing
  static constexpr double mu = trigger_rate * deltat_crossing;

  // print configuration
  std::cout << "SimulateCollisions - deltat_crossing: " << deltat_crossing << std::endl;
  std::cout << "SimulateCollisions - trigger_rate: " << trigger_rate << std::endl;
  std::cout << "SimulateCollisions - mu: " << mu << std::endl;

  // number of requested triggers (should correspond approximately to 1/10 second of data taking)
  // static constexpr uint ntrigtot = 1e4;
  // static constexpr uint ntrigtot = 6e6;
  static constexpr uint ntrigtot = 1e6;
  // write all timestamps to a file
  std::ofstream out("./data/timestamps_50kHz_1M.txt");
  out << "// bunchcrossin id; time (ns)" << std::endl;
  out << "// assuming 106ns between bunches" << std::endl;

  // running collision time
  int64_t bunchcrossing = 0;
  float time = 0;
  float timeD = 0;
  TTree *timestamps = new TTree("timestamps", "beamcrossings timestamps");
  timestamps->Branch("bc", &bunchcrossing, "bc/I");
  timestamps->Branch("t", &time, "t/F");
  timestamps->Branch("dt", &timeD, "dt/F");

  // keep track of the last collision time
  double previous_trigger_time = 0;

  // generate triggers
  for (int itrig = 0; itrig < ntrigtot;)
  {
    ++bunchcrossing;
    time += deltat_crossing;

    auto ntrig = gsl_ran_poisson(rng, mu);
    for (uint i = 0; i < ntrig; ++i)
    {
      out << " " << static_cast<int>(bunchcrossing) << " " << static_cast<int>(time * 1e9) << std::endl;
      timeD = (time - previous_trigger_time) * 1e9;
      timestamps->Fill();
      previous_trigger_time = time;
      ++itrig;
      if (itrig % 100000 == 0)
      {
        std::cout << "SimulateCollisions - itrig: " << itrig << std::endl;
      }
    }
  }

  // printout last trigger time
  std::cout << "SimulateCollisions - last trigger time: " << time << std::endl;

  out.close();
  // out.close();
  gsl_rng_free(rng);

  TFile *outf = new TFile("./data/timestamps_50kHz_1M_t.root", "recreate");
  timestamps->Write();

  outf->Write();
  outf->Close();
}
