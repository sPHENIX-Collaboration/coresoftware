#include <TTree.h>
#include <TFile.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h> 

#include <stdio.h>
#include <stdlib.h>

R__LOAD_LIBRARY(libgslcblas.so)
R__LOAD_LIBRARY(libgsl.so)

//_______________________________________________________________
void simulateTimestamps()
{

  // random generator
  auto rng = gsl_rng_alloc(gsl_rng_mt19937);
  //gsl_rng_set( rng, time(0) );//time(0) is the seed
  gsl_rng_set( rng, 1 );

  // time interval (s) between two bunch crossing
  /* value copied from generators/phhepmc/Fun4AllHepMCPileupInputManager.cc */
  static constexpr double deltat_crossing = 106e-9; 
  
  // triggered rate (Hz)
  // static constexpr double trigger_rate = 1e5; // Hijing 100kHz
  static constexpr double trigger_rate = 5e4; // Hijing 50kHz
  //static constexpr double trigger_rate = 3e6; // pp collisions
  
  // mean number of collision per crossing
  static constexpr double mu = trigger_rate*deltat_crossing;
  
  // print configuration
  printf("SimulateCollisions - deltat_crossing: %f \n", deltat_crossing );
  printf("SimulateCollisions - trigger_rate: %f \n" , trigger_rate);
  printf("SimulateCollisions - mu: %f \n" , mu);
  
  // number of requested triggers (should correspond approximately to 1/10 second of data taking)
  // static constexpr uint ntrigtot = 1e4;
  //static constexpr uint ntrigtot = 6e6;
  static constexpr uint ntrigtot = 1e6;
  // write all timestamps to a file
  FILE *fptr;
  fptr = fopen("./data/timestamps_50kHz_1M.txt","w");
  fprintf(fptr,"// bunchcrossin id; time (ns)\n");
  fprintf(fptr,"// assuming 106ns between bunches\n" );

    
  // running collision time
  int64_t bunchcrossing = 0;
  float time = 0;
  float timeD = 0;
  TTree *timestamps = new TTree("timestamps","beamcrossings timestamps");
  timestamps -> Branch("bc",  &bunchcrossing,   "bc/I");
  timestamps -> Branch("t",  &time,   "t/F");
  timestamps -> Branch("dt",  &timeD,   "dt/F");

  // keep track of the last collision time
  double previous_trigger_time = 0;
  
  // generate triggers
  for( int itrig = 0; itrig < ntrigtot; )
  {    
    ++ bunchcrossing;
    time += deltat_crossing;
    
    auto ntrig = gsl_ran_poisson( rng, mu );    
    for( uint i = 0; i < ntrig; ++i )
    { 
      fprintf(fptr,"%d %d\n",int(bunchcrossing), int(time*1e9) );
      timeD = (time-previous_trigger_time)*1e9;
      timestamps ->Fill();
      previous_trigger_time = time;
      ++itrig;
      if( itrig%100000==0 ) printf("SimulateCollisions - itrig: %d \n", itrig );
    }        
  }
  
  // printout last trigger time
  printf("SimulateCollisions - last trigger time: %f \n", time );
  
  fclose(fptr);
  //out.close();
  gsl_rng_free(rng);

  TFile *outf = new TFile("./data/timestamps_50kHz_1M_t.root","recreate");
  timestamps -> Write();
 
  outf -> Write();
  outf -> Close();  

}
