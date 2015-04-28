#include "FourHitSeedFinder.h"
#include <algorithm>


using namespace std;


void FourHitSeedFinder::findTracks_3_4(vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks, const HelixRange& range)
{
  vector<vector<SimpleHit3D> > layer_sorted;
  vector<SimpleHit3D> one_vec;
  layer_sorted.assign(4, one_vec);
  for(unsigned int i=0;i<hits.size();++i)
  {
    layer_sorted[hits[i].layer].push_back(hits[i]);
  }
  
  float k = 0.5*(range.max_k + range.min_k);
  float dzdl = 0.5*(range.max_dzdl + range.min_dzdl);
  float p_inv = 3.33333333333333314e+02*k*Bfield_inv*sqrt(1. - dzdl*dzdl);
  float scatter_weight = 2.;
  
  unsigned int l0_size = layer_sorted[0].size();
  unsigned int l1_size = layer_sorted[1].size();
  unsigned int l2_size = layer_sorted[2].size();
  unsigned int l3_size = layer_sorted[3].size();
  
  vector<vector<unsigned int> > three_layers;
  vector<unsigned int> one_hit;one_hit.assign(3,0);
  
  vector<double> chi2_hits;
  vector<unsigned int> tempcomb;
  tempcomb.assign(3,0);
  vector<float> kappa_vals;
  vector<float> dzdl_vals;
  for(unsigned int i0=0;i0<l0_size;++i0)
  {
    for(unsigned int i1=0;i1<l1_size;++i1)
    {
      for(unsigned int i2=0;i2<l2_size;++i2)
      {
        // have we seen this combo before?
        tempcomb[0] = layer_sorted[0][i0].index;
        tempcomb[1] = layer_sorted[1][i1].index;
        tempcomb[2] = layer_sorted[2][i2].index;
        sort(tempcomb.begin(), tempcomb.end());
        set<vector<unsigned int> >::iterator it = combos_3_fail.find(tempcomb);
        // if we know this combo sucks, ignore it
        if(it != combos_3_fail.end()){continue;}
//         it = combos_3_pass.find(tempcomb);
//         // if we know this combo is good, attach it to the list
//         if(it != combos_3_pass.end())
//         {
//           one_hit[0] = i0;
//           one_hit[1] = i1;
//           one_hit[2] = i2;
//           three_layers.push_back(one_hit);
//           continue;
//         }
        
        // we haven't seen this combo before ... let's fit it!
        
        SimpleTrack3D temp_track;
        temp_track.hits.push_back(layer_sorted[0][i0]);
        temp_track.hits.push_back(layer_sorted[1][i1]);
        temp_track.hits.push_back(layer_sorted[2][i2]);
        
        // calculate how much this thing should scatter
        float rad_diff = detector_radii[2] - detector_radii[1];
        float scatter = scatter_weight*p_inv*detector_scatter[1]*7.07106781186547462e-01;
        float dr = scatter*rad_diff;
        layer_xy_resolution[2] += dr;
        layer_z_resolution[2] += dr;
        double chi2 = fitTrack(temp_track, chi2_hits);
        layer_xy_resolution[2] -= dr;
        layer_z_resolution[2] -= dr;
        if(fabs(chi2) > chi2_cut){combos_3_fail.insert(tempcomb);continue;}
        one_hit[0] = i0;
        one_hit[1] = i1;
        one_hit[2] = i2;
        three_layers.push_back(one_hit);
        combos_3_pass.insert(tempcomb);
        kappa_vals.push_back(temp_track.kappa);
        dzdl_vals.push_back(temp_track.dzdl);
      }
    }
  }
  
  tempcomb.resize(4,0);
  // now three_layers contains a list of three-hit combinations from the inner 3 layers which pass the chi^2 cut
  // next, add in the hits from the 4th layer
  for(unsigned int j=0;j<three_layers.size();++j)
  {
    for(unsigned int i3=0;i3<l3_size;++i3)
    {
      tempcomb[0] = layer_sorted[0][three_layers[j][0]].index;
      tempcomb[1] = layer_sorted[1][three_layers[j][1]].index;
      tempcomb[2] = layer_sorted[2][three_layers[j][2]].index;
      tempcomb[3] = layer_sorted[3][i3].index;
      sort(tempcomb.begin(), tempcomb.end());
      set<vector<unsigned int> >::iterator it = combos.find(tempcomb);
      if(it != combos.end()){continue;}
      combos.insert(tempcomb);
      
      SimpleTrack3D temp_track;
      temp_track.hits.push_back(layer_sorted[1][three_layers[j][1]]);
      temp_track.hits.push_back(layer_sorted[2][three_layers[j][2]]);
      temp_track.hits.push_back(layer_sorted[3][i3]);
      
      // calculate how much this thing should scatter
      float rad_diff = detector_radii[3] - detector_radii[2];
      float scatter = scatter_weight*p_inv*detector_scatter[2]*7.07106781186547462e-01;
      float dr = scatter*rad_diff;
      layer_xy_resolution[2] += dr;
      layer_z_resolution[2] += dr;
      double chi2 = fitTrack(temp_track, chi2_hits);
      layer_xy_resolution[2] -= dr;
      layer_z_resolution[2] -= dr;
      if(fabs(chi2) > chi2_cut){continue;}
      
//       float sm_k = temp_track.kappa;
//       float big_k = kappa_vals[j];
//       if(sm_k > big_k)
//       {
//         big_k = temp_track.kappa;
//         sm_k = kappa_vals[j];
//       }
//       float rel_diff = (big_k - sm_k)/big_k;
//       if(rel_diff > (0.006/k)*0.5){continue;}
//       
//       float sm_z = fabs(temp_track.dzdl);
//       float big_z = fabs(dzdl);
//       if(sm_z > big_z)
//       {
//         big_z = fabs(temp_track.dzdl);
//         sm_z = fabs(dzdl);
//       }
//       rel_diff = (big_z - sm_z)/big_z;
//       if(rel_diff > (0.006/k)*0.5){continue;}
      
      SimpleTrack3D good_track;
      good_track.hits.push_back(layer_sorted[0][three_layers[j][0]]);
      good_track.hits.push_back(layer_sorted[1][three_layers[j][1]]);
      good_track.hits.push_back(layer_sorted[2][three_layers[j][2]]);
      good_track.hits.push_back(layer_sorted[3][i3]);
      
      tracks.push_back(good_track);
    }
  }
}




