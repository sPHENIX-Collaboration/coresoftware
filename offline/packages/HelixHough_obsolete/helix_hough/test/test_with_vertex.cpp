#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "HelixHough.h"
#include "HelixRange.h"
#include "HelixResolution.h"
#include "VertexFitFunc.h"
#include "NewtonMinimizerGradHessian.h"
#include <sys/time.h>
#include <Eigen/LU>
#include <Eigen/Core>
#include <math.h>
#include "FourHitSeedFinder/FourHitSeedFinder.h"

using namespace std;
using namespace FitNewton;
using namespace Eigen;

int main(int argc, char** argv)
{
  bool print_truth = false;
  bool print_reco = false;
  
  TFile infile(argv[1]);
  TTree* etree=0;
  TTree* ttree=0;
  infile.GetObject("events", etree);
  etree->SetBranchAddress("tracklist", &ttree);
  unsigned int track=0;
  unsigned int nhits=0;
  int layer[12];
  double x_hits[12];
  double y_hits[12];
  double z_hits[12];
  double trk_kappa=0.;
  double trk_d=0.;
  double trk_phi=0.;
  double trk_dzdl=0.;
  double trk_z0=0.;
  unsigned char charge=0;
  
  
  int nlayers = 4;
  vector<float> radii;
  radii.assign(nlayers,0.);
  radii[0]=2.5;
  radii[1]=5.0;
  radii[2]=10.0;
  radii[3]=14.0;
  vector<float> smear_xy_layer;smear_xy_layer.assign(nlayers,0);
  vector<float> smear_z_layer;smear_z_layer.assign(nlayers,0);
  float sqrt_12 = sqrt(12.);
  smear_xy_layer[0] = (50.0e-4/sqrt_12);
  smear_z_layer[0] = (425.0e-4/sqrt_12);
  smear_xy_layer[1] = (50.0e-4/sqrt_12);
  smear_z_layer[1] = (425.0e-4/sqrt_12);
  smear_xy_layer[2] = (80.0e-4/sqrt_12);
  smear_z_layer[2] = (1000.0e-4/sqrt_12);
  smear_xy_layer[3] = (80.0e-4/sqrt_12);
  smear_z_layer[3] = (1000.0e-4/sqrt_12);
  
  //phi,d,kappa,dzdl,z0
  HelixResolution min_res(0.01, 0.5, 0.002, 0.01, 2.);
  HelixResolution max_res(0.001, 0.5, 0.001, 0.01, 2.);
  HelixRange top_range(0., 1.*M_PI,   -0.0025, 0.0025,   0., 0.03,   -0.9, 0.9,   -0.005, 0.005);
  FourHitSeedFinder tracker(radii, 4, 1, 2, 4, 1, min_res, max_res, top_range);
  tracker.setLayerResolution(smear_xy_layer, smear_z_layer);
  tracker.setVertexResolution(0.005, 0.01);
  tracker.setUsingVertex(true);
  tracker.setChi2Cut(3.0);
  tracker.setPrintTimings(true);
  unsigned int max_hits = 5;
  
  
  
  
  
  for(unsigned int ev=0;ev<etree->GetEntries();ev++)
  {
    etree->GetEntry(ev);
    ttree->SetBranchAddress("track", &track);
    ttree->SetBranchAddress("nhits", &nhits);
    ttree->SetBranchAddress("x_hits", &x_hits);
    ttree->SetBranchAddress("y_hits", &y_hits);
    ttree->SetBranchAddress("z_hits", &z_hits);
    ttree->SetBranchAddress("layer", &layer);
    ttree->SetBranchAddress("trk_kappa", &trk_kappa);
    ttree->SetBranchAddress("trk_d", &trk_d);
    ttree->SetBranchAddress("trk_phi", &trk_phi);
    ttree->SetBranchAddress("trk_dzdl", &trk_dzdl);
    ttree->SetBranchAddress("trk_z0", &trk_z0);
    ttree->SetBranchAddress("charge", &charge);
    
    cout<<"event "<<ev<<":"<<endl<<endl;
    
    vector<SimpleHit3D> hits;     // 3d space point + hit index
    vector<SimpleTrack3D> tracks; // list of space points _and_ helix parameters
    unsigned int index=0;
    cout<<"total MC tracks = "<<ttree->GetEntries()<<endl;
    cout<<endl;
    
    for(unsigned int trk=0;trk<ttree->GetEntries();trk++)
    {
      ttree->GetEntry(trk);
      if(print_truth==true)
      {
        cout<<"truth track : "<<endl;
        cout<<"phi = "<<trk_phi<<endl;
        cout<<"d = "<<trk_d<<endl;
        cout<<"kappa = "<<trk_kappa<<endl;
        cout<<"dzdl = "<<trk_dzdl<<endl;
        cout<<"z0 = "<<trk_z0<<endl;
        cout<<endl;
      }
      for(unsigned int hit=0;hit<nhits;hit++)
      {
        float phi = atan2(y_hits[hit], x_hits[hit]);
        float xy_error = smear_xy_layer[layer[hit]]*sqrt_12*0.5;
        float x_error = fabs(xy_error*sin(phi));
        float y_error = fabs(xy_error*cos(phi));
        float z_error = smear_z_layer[layer[hit]]*sqrt_12*0.5;
        
        hits.push_back(SimpleHit3D(x_hits[hit],x_error, y_hits[hit],y_error, z_hits[hit],z_error, index, layer[hit]));
        index++;
      }
    }
    
    double z_avg = 0.;
    double z_weight = 0.;
    for(unsigned int ht=0;ht<hits.size();++ht)
    {
      double temp = 1./(hits[ht].dz);
      z_avg += (hits[ht].z)*temp;
      z_weight += temp;
    }
    z_avg/=z_weight;
    cout<<"z average = "<<z_avg<<endl;
    
    HelixResolution min_res_init(0.01, 10.5, 0.01, 0.02, 0.01);
    HelixResolution max_res_init(0.001, 0.5, 0.01, 0.02, 0.01);
    HelixRange top_range_init(0., 1.*M_PI,   -0.1, 0.1,   0., 0.005,   -0.9, 0.9,   -0.1 + z_avg, 0.1 + z_avg);
    FourHitSeedFinder tracker_init(radii, 2, 1, 2, 4, 2, min_res_init, max_res_init, top_range_init);
    tracker_init.setLayerResolution(smear_xy_layer, smear_z_layer);
    tracker_init.setVertexResolution(0.001, 0.005);
    tracker_init.setUsingVertex(false);
    tracker_init.setChi2Cut(3.0);
    tracker_init.setPrintTimings(true);
    unsigned int max_hits_init = 8;
    unsigned int maxtracks = 160;
    
    
    maxtracks = 15*(hits.size())/1000;
    if(maxtracks < 40){maxtracks=40;}
    
    timeval t1,t2;
    double time1=0.;
    double time2=0.;
    gettimeofday(&t1, NULL);
    tracker_init.findHelices(hits, 4, max_hits_init, tracks, maxtracks);
    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
    cout<<"initial tracking time = "<<(time2 - time1)<<endl;
    
    unsigned int ngood = 0;
    unsigned int n_trk = tracks.size();
    for(unsigned int trk=0;trk<n_trk;++trk)
    {
      unsigned int high = 0;
      unsigned int low = (1<<31);
      unsigned int n_hit = tracks[trk].hits.size();
      for(unsigned int ht=0;ht<n_hit;++ht)
      {
        tracks[trk].hits[ht].index > high ? high = tracks[trk].hits[ht].index : high = high;
        tracks[trk].hits[ht].index < low ? low = tracks[trk].hits[ht].index : low = low;
      }
      if( (high - low) <= 3 ){ngood++;}
    }
    
    cout<<"tracks found before vertex finding = "<<tracks.size()<<endl;
    cout<<"good tracks found before vertex finding = "<<ngood<<endl;
    
    VertexFitFunc vfit;
    vfit.setTracks(&tracks);
    
    // very general minimization routine (good for other things too)
    NewtonMinimizerGradHessian minimizer;
    minimizer.setFunction(&vfit);
    
    VectorXd start_point = VectorXd::Zero(3); // input initial guess
    VectorXd min_point = VectorXd::Zero(3);   // output vector for minimize method below
    
    //first use a large sigma for expo-dca^2 calculation (same in x,y,z)
    vfit.setFixedPar(0, 3.0);
    // -> guess
    // <- best fit
    // -> relative tolerance
    // -> max number of iterations
    // -> absolute tolerance
    // return true if successful
    minimizer.minimize(start_point, min_point, 1.0e-9, 16, 1.0e-15);
    //then continue with a smaller sigma
    vfit.setFixedPar(0, 0.3);
    start_point = min_point;
    minimizer.minimize(start_point, min_point, 1.0e-9, 16, 1.0e-15);
    vfit.setFixedPar(0, 0.03);
    start_point = min_point;
    minimizer.minimize(start_point, min_point, 1.0e-9, 16, 1.0e-15);
    
    
    // store output vertex spatial point
    vector<double> vertex;
    vertex.assign(3,0.);
    vertex[0] = min_point(0);
    vertex[1] = min_point(1);
    vertex[2] = min_point(2);
    
//     vertex[0] = 0.;
//     vertex[1] = 0.;
//     vertex[2] = 0.;
    
    cout << "Found Event vertex = " << vertex[0] << " " << vertex[1] << " " << vertex[2] << endl << endl;
    
    tracks.clear();
    
    
    for(unsigned int ht=0;ht<hits.size();++ht)
    {
      hits[ht].x -= vertex[0];
      hits[ht].y -= vertex[1];
      hits[ht].z -= vertex[2];
    }
    
    gettimeofday(&t1, NULL);
    //     tracker.findHelices(hits, 4, max_hits, tracks);
    tracker.findHelices(hits, 4, max_hits, tracks);
    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec)/1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec)/1000000.);
    cout<<"primary tracking time = "<<(time2 - time1)<<endl;
    
    for(unsigned int trk=0;trk<tracks.size();++trk)
    {
      for(unsigned int ht=0;ht<tracks[trk].hits.size();++ht)
      {
        tracks[trk].hits[ht].x += vertex[0];
        tracks[trk].hits[ht].y += vertex[1];
        tracks[trk].hits[ht].z += vertex[2];
      }
    }
    
    
    
    ngood = 0;
    n_trk = tracks.size();
    for(unsigned int trk=0;trk<n_trk;++trk)
    {
      unsigned int high = 0;
      unsigned int low = (1<<31);
      unsigned int n_hit = tracks[trk].hits.size();
      for(unsigned int ht=0;ht<n_hit;++ht)
      {
        tracks[trk].hits[ht].index > high ? high = tracks[trk].hits[ht].index : high = high;
        tracks[trk].hits[ht].index < low ? low = tracks[trk].hits[ht].index : low = low;
      }
      if( (high - low) <= 3 ){ngood++;}
    }
    
    if(print_reco==true)
    {
      for(unsigned int trk=0;trk<n_trk;++trk)
      {
        cout<<"track "<<trk<<" : "<<endl;
        unsigned int n_hit = tracks[trk].hits.size();
        for(unsigned int ht=0;ht<n_hit;++ht)
        {
          cout<<"hit "<<ht<<" : ";
          cout<<tracks[trk].hits[ht].x<<" ";
          cout<<tracks[trk].hits[ht].y<<" ";
          cout<<tracks[trk].hits[ht].z<<" ";
          cout<<tracks[trk].hits[ht].index<<endl;
        }
        cout<<"phi = "<<tracks[trk].phi<<endl;
        cout<<"d = "<<tracks[trk].d<<endl;
        cout<<"kappa = "<<tracks[trk].kappa<<endl;
        cout<<"dzdl = "<<tracks[trk].dzdl<<endl;
        cout<<"z0 = "<<tracks[trk].z0<<endl;
        cout<<endl;
      }
    }
    
    
    cout<<n_trk<<" tracks found"<<endl;
    cout<<ngood<<" good tracks found"<<endl;
    cout<<endl<<endl;
  }
  
  
  return 0;

}
