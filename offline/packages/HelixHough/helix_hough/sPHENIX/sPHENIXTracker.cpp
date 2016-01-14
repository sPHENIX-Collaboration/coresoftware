#include "sPHENIXTracker.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <float.h>
#include <sys/time.h>


using namespace std;
using namespace Eigen;
using namespace SeamStress;


class hitTriplet
{
public:
  hitTriplet(unsigned int h1, unsigned int h2, unsigned int h3, unsigned int t, float c) : hit1(h1), hit2(h2), hit3(h3), track(t), chi2(c) {}
  ~hitTriplet(){}
  
  bool operator<(const hitTriplet& other) const
  {
    return ( hit1 < other.hit1 ) || ( ( hit2 < other.hit2 ) && ( hit1 == other.hit1 ) ) || ( ( hit3 < other.hit3 ) && ( hit1 == other.hit1 ) && ( hit2 == other.hit2 ) );
  }
  
  bool operator==(const hitTriplet& other) const
  {
    return ( (hit1 == other.hit1) && (hit2 == other.hit2) && (hit3 == other.hit3) );
  }
  
  unsigned int hit1, hit2, hit3, track;
  float chi2;
};


void sPHENIXTracker::tripletRejection(vector<SimpleTrack3D>& input, vector<SimpleTrack3D>& output, vector<bool>& usetrack, vector<float>& next_best_chi2)
{
  vector<hitTriplet> trips;
  for(unsigned int i=0;i<input.size();++i)
  {
    for(unsigned int h1=0;h1<input[i].hits.size();++h1)
    {
      for(unsigned int h2=(h1+1);h2<input[i].hits.size();++h2)
      {
        for(unsigned int h3=(h2+1);h3<input[i].hits.size();++h3)
        {
          if(cut_on_dca == false)
          {
            trips.push_back(hitTriplet(input[i].hits[h1].index,input[i].hits[h2].index,input[i].hits[h3].index,i,track_states[i].chi2));
          }
        }
      }
    }
  }
  if(trips.size() == 0){return;}
  sort(trips.begin(), trips.end());
  unsigned int pos=0;
  unsigned int cur_h1 = trips[pos].hit1;
  unsigned int cur_h2 = trips[pos].hit2;
  while(pos < trips.size())
  {
    unsigned int next_pos = pos+1;
    if(next_pos >= trips.size()){break;}
    while( trips[pos] == trips[next_pos] )
    {
      next_pos+=1;
      if(next_pos >= trips.size()){break;}
    }
    if((next_pos - pos) > 1)
    {
      float best_chi2 = trips[pos].chi2;
      float next_chi2 = trips[pos+1].chi2;
      unsigned int best_pos = pos;
      for(unsigned int i=(pos+1);i<next_pos;++i)
      {
        if(input[trips[i].track].hits.size() < input[trips[best_pos].track].hits.size())
        {
          continue;
        }
        else if( (input[trips[i].track].hits.size() > input[trips[best_pos].track].hits.size()) || ( input[trips[i].track].hits.back().layer > input[trips[best_pos].track].hits.back().layer ) )
        {
          next_chi2 = best_chi2;
          best_chi2 = trips[i].chi2;
          best_pos = i;
          continue;
        }
        if((trips[i].chi2 < best_chi2) || ( usetrack[trips[best_pos].track]==false ))
        {
          next_chi2 = best_chi2;
          best_chi2 = trips[i].chi2;
          best_pos = i;
        }
        else if(trips[i].chi2 < next_chi2)
        {
          next_chi2 = trips[i].chi2;
        }
      }
      for(unsigned int i=pos;i<next_pos;++i)
      {
        if(i != best_pos)
        {
          usetrack[trips[i].track] = false;
        }
        else
        {
          next_best_chi2[trips[i].track] = next_chi2;
        }
      }
    }
    pos = next_pos;
    cur_h1 = trips[pos].hit1;
    cur_h2 = trips[pos].hit2;
  }
}


sPHENIXTracker::sPHENIXTracker(unsigned int n_phi,
			       unsigned int n_d,
			       unsigned int n_k,
			       unsigned int n_dzdl,
			       unsigned int n_z0,
			       HelixResolution& min_resolution,
			       HelixResolution& max_resolution,
			       HelixRange& range,
			       vector<float>& material,
			       vector<float>& radius,
			       float Bfield)
  : HelixHough(n_phi, n_d, n_k, n_dzdl, n_z0, min_resolution, max_resolution, range),
    fast_chi2_cut_par0(12.),
    fast_chi2_cut_par1(0.),
    fast_chi2_cut_max(FLT_MAX),
    chi2_cut(3.),
    chi2_removal_cut(1.),
    n_removal_hits(0),
    seeding(false),
    verbosity(0),
    cut_on_dca(false),
    dca_cut(0.01),
    vertex_x(0.),
    vertex_y(0.),
    vertex_z(0.),
    required_layers(0),
    reject_ghosts(false),
    nfits(0), findtracksiter(0),
    prev_max_k(0.),
    prev_max_dzdl(0.),
    prev_p_inv(0.),
    seed_layer(0),
    ca_chi2_cut(2.0),
    cosang_cut(0.985)
{
  vector<float> detector_material;
  
  for(unsigned int i=0;i<radius.size();++i)
  {
    detector_radii.push_back(radius[i]);
  }
  for(unsigned int i=0;i<material.size();++i)
  {
    detector_scatter.push_back(1.41421356237309515*0.0136*sqrt( 3.*material[i] ));
    detector_material.push_back(3.*material[i]);
  }
  
  detector_B_field = Bfield;
  
  integrated_scatter.assign(detector_scatter.size(),0.);
  float total_scatter_2 = 0.;
  for(unsigned int l=0;l<detector_scatter.size();++l)
  {
    total_scatter_2 += detector_scatter[l]*detector_scatter[l];
    integrated_scatter[l] = sqrt(total_scatter_2);
  }
  
  kalman = new CylinderKalman(detector_radii, detector_material, detector_B_field);
  
  vector<SimpleHit3D> one_layer;
  layer_sorted.assign(n_layers, one_layer);
  for(unsigned int i=0;i<4;++i){layer_sorted_1[i].assign(n_layers, one_layer);}
  temp_comb.assign(n_layers, 0);
}


sPHENIXTracker::sPHENIXTracker(vector<vector<unsigned int> >& zoom_profile,
			       unsigned int minzoom,
			       HelixRange& range,
			       vector<float>& material,
			       vector<float>& radius,
			       float Bfield,
			       bool parallel,
			       unsigned int num_threads)
  : HelixHough(zoom_profile, minzoom, range),
    fast_chi2_cut_par0(12.),
    fast_chi2_cut_par1(0.),
    fast_chi2_cut_max(FLT_MAX),
    chi2_cut(3.),
    chi2_removal_cut(1.),
    n_removal_hits(0),
    seeding(false),
    verbosity(0),
    cut_on_dca(false),
    dca_cut(0.01),
    vertex_x(0.),
    vertex_y(0.),
    vertex_z(0.),
    required_layers(0),
    reject_ghosts(false),
    nfits(0),
    findtracksiter(0),
    prev_max_k(0.),
    prev_max_dzdl(0.),
    prev_p_inv(0.),
    seed_layer(0),
    nthreads(num_threads),
    vssp(NULL),
    pins(NULL),
    is_parallel(parallel),
    is_thread(false),
    ca_chi2_cut(2.0),
    cosang_cut(0.985)
{
  vector<float> detector_material;
  
  for(unsigned int i=0;i<radius.size();++i)
  {
    detector_radii.push_back(radius[i]);
  }
  for(unsigned int i=0;i<material.size();++i)
  {
    detector_scatter.push_back(1.41421356237309515*0.0136*sqrt( 3.*material[i] ));
    detector_material.push_back(3.*material[i]);
  }
  
  detector_B_field = Bfield;
  
  integrated_scatter.assign(detector_scatter.size(),0.);
  float total_scatter_2 = 0.;
  for(unsigned int l=0;l<detector_scatter.size();++l)
  {
    total_scatter_2 += detector_scatter[l]*detector_scatter[l];
    integrated_scatter[l] = sqrt(total_scatter_2);
  }
  
  kalman = new CylinderKalman(detector_radii, detector_material, detector_B_field);
  
  vector<SimpleHit3D> one_layer;
  layer_sorted.assign(n_layers, one_layer);
  for(unsigned int i=0;i<4;++i){layer_sorted_1[i].assign(n_layers, one_layer);}
  temp_comb.assign(n_layers, 0);
  
  
  if(is_parallel == true)
  {
    Seamstress::init_vector(num_threads, vss);
    
    vssp = new vector<Seamstress*>();
    for(unsigned int i=0;i<vss.size();i++){vssp->push_back(&(vss[i]));}
    
    pins = new Pincushion<sPHENIXTracker>(this, vssp);
    
    vector<vector<unsigned int> > zoom_profile_new;
    for(unsigned int i=1;i<zoom_profile.size();++i)
    {
      zoom_profile_new.push_back(zoom_profile[i]);
    }
    
    for(unsigned int i=0;i<nthreads;++i)
    {
      thread_trackers.push_back(new sPHENIXTracker(zoom_profile, minzoom, range, material, radius, Bfield) );
      thread_trackers.back()->setThread();
      thread_trackers.back()->setStartZoom(1);
      thread_tracks.push_back(vector<SimpleTrack3D>());
      thread_ranges.push_back(HelixRange());
      thread_hits.push_back(vector<SimpleHit3D>());
      split_output_hits.push_back(new vector<vector<SimpleHit3D> >());
      split_ranges.push_back(new vector<HelixRange>());
      split_input_hits.push_back(vector<SimpleHit3D>());
    }
  }
}


sPHENIXTracker::~sPHENIXTracker()
{
  if ( kalman != NULL ) delete kalman;
  for(unsigned int i=0;i<vss.size();i++)
  {
    vss[i].stop();
  }
  for(unsigned int i=0;i<thread_trackers.size();++i)
  {
    delete thread_trackers[i];
    delete split_output_hits[i];
    delete split_ranges[i];
  }
  
  if ( pins != NULL ) delete pins;
  if ( vssp != NULL ) delete vssp;
}

float sPHENIXTracker::kappaToPt(float kappa) {  
  return detector_B_field / 333.6 / kappa;
}

float sPHENIXTracker::ptToKappa(float pt) {  
  return detector_B_field / 333.6 / pt;
}

// hel should be +- 1
static void xyTangent(SimpleHit3D& hit1, SimpleHit3D& hit2, float kappa, float hel, float& ux_out, float& uy_out, float& ux_in, float& uy_in)
{
  float x = hit2.x - hit1.x;
  float y = hit2.y - hit1.y;
  float D = sqrt(x*x + y*y);
  float ak = 0.5*kappa*D;
  float D_inv = 1./D;
  float hk = sqrt(1. - ak*ak);
  
  float kcx = (ak*x + hel*hk*y)*D_inv;
  float kcy = (ak*y - hel*hk*x)*D_inv;
  float ktx = -(kappa*y - kcy);
  float kty = kappa*x - kcx;
  float norm = 1./sqrt(ktx*ktx + kty*kty);
  ux_out = ktx*norm;
  uy_out = kty*norm;
  
  ktx = kcy;
  kty = -kcx;
  norm = 1./sqrt(ktx*ktx + kty*kty);
  ux_in = ktx*norm;
  uy_in = kty*norm;
}


// hel should be +- 1
static float cosScatter(SimpleHit3D& hit1, SimpleHit3D& hit2, SimpleHit3D& hit3, float kappa, float hel)
{
  float ux_in=0.;
  float uy_in=0.;
  float ux_out=0.;
  float uy_out=0.;
  
  float temp1=0.;
  float temp2=0.;
  
  xyTangent(hit1, hit2, kappa, hel, ux_in, uy_in, temp1, temp2);
  xyTangent(hit2, hit3, kappa, hel, temp1, temp2, ux_out, uy_out);
  
  return ux_in*ux_out + uy_in*uy_out;
}


static float dzdsSimple(SimpleHit3D& hit1, SimpleHit3D& hit2, float k)
{
  float x = hit2.x - hit1.x;
  float y = hit2.y - hit1.y;
  float D = sqrt(x*x + y*y);
  float s = 0.;
  float temp1 = k*D*0.5;temp1*=temp1;
  float temp2 = D*0.5;
  s += 2.*temp2;
  temp2*=temp1;
  s += temp2/3.;
  temp2*=temp1;
  s += (3./20.)*temp2;
  temp2*=temp1;
  s += (5./56.)*temp2;
  
  return (hit2.z - hit1.z)/s;
}


float sPHENIXTracker::dcaToVertexXY(SimpleTrack3D& track, float vx, float vy)
{
  float d_out = 0.;
  
  // find point at the dca to 0
  float x0 = track.d*cos(track.phi);
  float y0 = track.d*sin(track.phi);
  
  // change variables so x0,y0 -> 0,0
  float phi2 = atan2((1. + track.kappa*track.d)*sin(track.phi) - track.kappa*y0, (1. + track.kappa*track.d)*cos(track.phi) - track.kappa*x0);
  
  // translate so that (0,0) -> (x0 - vx , y0 - vy)
  float cosphi = cos(phi2);
  float sinphi = sin(phi2);
  float tx = cosphi + track.kappa*(x0-vx);
  float ty = sinphi + track.kappa*(y0-vy);
  float dk = sqrt( tx*tx + ty*ty ) - 1.;
  if(track.kappa == 0.){d_out = (x0-vx)*cosphi + (y0-vy)*sinphi;}
  else{d_out = dk/track.kappa;}
  return fabs(d_out);
}


void sPHENIXTracker::finalize(vector<SimpleTrack3D>& input, vector<SimpleTrack3D>& output)
{
  #ifdef AVXHOUGH
  if(findtracks_bin!=0)
  {
    findTracksBySegments_avx_run(input);
  }
  #endif
  
  if(is_thread == true)
  {
    for(unsigned int i=0;i<input.size();++i)
    {
      output.push_back(input[i]);
    }
    return;
  }
  
  unsigned int nt = input.size();
  vector<bool> usetrack;
  usetrack.assign(input.size(), true);
  vector<float> next_best_chi2;
  next_best_chi2.assign(input.size(), 99999.);
  
  if (reject_ghosts == true) {
    tripletRejection(input, output, usetrack, next_best_chi2);
  }
  
  vector<HelixKalmanState> states_new;
  
  for(unsigned int i=0;i<nt;++i)
  {
    if(usetrack[i] == true)
    {
      if( !(track_states[i].chi2 == track_states[i].chi2) ){continue;}
      
      output.push_back(input[i]);
      output.back().index = (output.size() - 1);
      states_new.push_back(track_states[i]);
      isolation_variable.push_back(next_best_chi2[i]);
    }
  }
  track_states = states_new;
  if(smooth_back == true)
  {
    for(unsigned int i=0;i<output.size();++i)
    {
      
      if(n_layers < 20)
      {
        HelixKalmanState state = track_states[i];
        
        track_states[i].C *= 3.;
        track_states[i].chi2 = 0.;
        track_states[i].x_int = 0.;
        track_states[i].y_int = 0.;
        track_states[i].z_int = 0.;
        //       track_states[i].position = output[i].hits.size();
        track_states[i].position = 0;
        for(int h=(output[i].hits.size() - 1);h>=0;--h)
        {
          SimpleHit3D hit = output[i].hits[h];
          float err_scale = 1.;
          int layer = hit.layer;
          if( (layer >= 0) && (layer < (int)(hit_error_scale.size()) ) ){err_scale = hit_error_scale[layer];}
          err_scale *= 3.0;//fudge factor, like needed due to non-gaussian errors
          hit.dx *= err_scale;hit.dy *= err_scale;hit.dz *= err_scale;
          kalman->addHit(hit, track_states[i]);
          track_states[i].position = h;
        }
        
        //       track_states[i].C = state.C;
        track_states[i].chi2 = state.chi2;
        track_states[i].C *= 2./3.;
        
        output[i].phi = track_states[i].phi;
        output[i].d = track_states[i].d;
        output[i].kappa = track_states[i].kappa;
        output[i].z0 = track_states[i].z0;
        output[i].dzdl = track_states[i].dzdl;
      }
      else
      {
        HelixKalmanState state = track_states[i];
        
        track_states[i].C *= 30;
        
        track_states[i].chi2 = 0.;
        track_states[i].x_int = 0.;
        track_states[i].y_int = 0.;
        track_states[i].z_int = 0.;
        track_states[i].position = output[i].hits.size();
        for(int h=(output[i].hits.size() - 1);h>=0;--h)
        {
          SimpleHit3D hit = output[i].hits[h];
          float err_scale = 0.66;
          hit.dx *= err_scale;hit.dy *= err_scale;hit.dz *= err_scale;
          kalman->addHit(hit, track_states[i]);
        }

        SimpleTrack3D temp_track = output[i];
        fitTrack(temp_track);
        track_states[i].kappa = temp_track.kappa;
        track_states[i].nu = sqrt(temp_track.kappa);
        
        if(!(track_states[i].kappa == track_states[i].kappa)){track_states[i] = state;}
        
        if(output[i].phi < 0.){output[i].phi += 2.*M_PI;}
        output[i].phi = track_states[i].phi;
        output[i].d = track_states[i].d;
        output[i].kappa = track_states[i].kappa;
        output[i].z0 = track_states[i].z0;
        output[i].dzdl = track_states[i].dzdl;
      }
    }
  }
  
  
  if(verbosity > 0)
  {
    cout<<"# fits = "<<nfits<<endl;
    cout<<"findTracks called "<<findtracksiter<<" times"<<endl;
    cout<<"CAtime = "<<CAtime<<endl;
    cout<<"KALime = "<<KALtime<<endl;
  }
  
}


void sPHENIXTracker::findTracks(vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks, const HelixRange& range)
{
  findtracksiter += 1;
  #ifdef AVXHOUGH
  findTracksBySegments_avx(hits,tracks,range);
  #else
  findTracksBySegments(hits,tracks,range);
  #endif
  
//   findTracksBySegments(hits,tracks,range);
}


void sPHENIXTracker::findSeededTracks(vector<SimpleTrack3D>& seeds, vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks, const HelixRange& range)
{
  findtracksiter += 1;
//   findSeededTracksByProjection(seeds, hits, tracks, range);
  findSeededTracksbySegments(seeds, hits, tracks, range);
}


bool sPHENIXTracker::breakRecursion(const vector<SimpleHit3D>& hits, const HelixRange& range)
{
  if(seeding == true){return false;}
  unsigned int layer_mask[4] = { 0, 0, 0, 0 };
  for(unsigned int i=0;i<hits.size();++i)
  {
    if(  hits[i].layer < 32  ) { layer_mask[0] = layer_mask[0] | (1<<hits[i].layer); }
    else if(  hits[i].layer < 64  ) { layer_mask[1] = layer_mask[1] | (1 << (hits[i].layer-32) ); }
    else if(  hits[i].layer < 96  ) { layer_mask[2] = layer_mask[2] | (1 << (hits[i].layer-64) ); }
    else if(  hits[i].layer < 128  ) { layer_mask[3] = layer_mask[3] | (1 << (hits[i].layer-96) ); }
  }
  unsigned int nlayers = __builtin_popcount( layer_mask[0] ) + __builtin_popcount( layer_mask[1] ) + __builtin_popcount( layer_mask[2] ) + __builtin_popcount( layer_mask[3] ) ;
  
  return (nlayers < required_layers);
  
}


float sPHENIXTracker::phiError(SimpleHit3D& hit, float min_k, float max_k, float min_d, float max_d, float min_z0, float max_z0, float min_dzdl, float max_dzdl, bool pairvoting)
{
  float Bfield_inv = 1./detector_B_field;
  float p_inv=0.;
  
  if((prev_max_k==max_k) && (prev_max_dzdl==max_dzdl))
  {
    p_inv=prev_p_inv;
  }
  else
  {
    prev_max_k=max_k;
    prev_max_dzdl=max_dzdl;
    prev_p_inv = 3.33333333333333314e+02*max_k*Bfield_inv*sqrt(1. - max_dzdl*max_dzdl);
    p_inv=prev_p_inv;
  }
  float total_scatter_2 = 0.;
  for(int i=seed_layer+1;i<=(hit.layer);++i)
  {
    float this_scatter = detector_scatter[i-1]*(detector_radii[i]-detector_radii[i-1])/detector_radii[i];
    total_scatter_2 += this_scatter*this_scatter;
  }
  float angle = p_inv*sqrt(total_scatter_2)*1.0;
  float dsize = 0.5*(max_d-min_d);
  float angle_from_d = dsize/detector_radii[hit.layer];
  float returnval = 0.;
  if(pairvoting==false)
  {
    if(angle_from_d > angle){returnval=0.;}
    else{returnval = (angle - angle_from_d);}
  }
  else
  {
    returnval = angle;
  }
  
  
  return returnval;
}


float sPHENIXTracker::dzdlError(SimpleHit3D& hit, float min_k, float max_k, float min_d, float max_d, float min_z0, float max_z0, float min_dzdl, float max_dzdl, bool pairvoting)
{
  float Bfield_inv = 1./detector_B_field;
  float p_inv=0.;
  
  if((prev_max_k==max_k) && (prev_max_dzdl==max_dzdl))
  {
    p_inv=prev_p_inv;
  }
  else
  {
    prev_max_k=max_k;
    prev_max_dzdl=max_dzdl;
    prev_p_inv = 3.33333333333333314e+02*max_k*Bfield_inv*sqrt(1. - max_dzdl*max_dzdl);
    p_inv=prev_p_inv;
  }
  float total_scatter_2 = 0.;
  for(int i=seed_layer+1;i<=(hit.layer);++i)
  {
    float this_scatter = detector_scatter[i-1]*(detector_radii[i]-detector_radii[i-1])/detector_radii[i];
    total_scatter_2 += this_scatter*this_scatter;
  }
  float angle = p_inv*sqrt(total_scatter_2)*1.0;
  float z0size = 0.5*(max_z0-min_z0);
  float angle_from_z0 = z0size/detector_radii[hit.layer];
  float returnval = 0.;
  if(pairvoting==false)
  {
    if(angle_from_z0 > angle){returnval=0.;}
    else{returnval = (angle - angle_from_z0);}
  }
  else
  {
    returnval = angle;
  }
  
  
  return returnval;
}


void sPHENIXTracker::setRangeFromSeed(HelixRange& range, SimpleTrack3D& seed)
{
  HelixKalmanState* state = &(seed_states[seed.index]);
  
  float dphi = 2.*sqrt(state->C(0,0));
  float dd = 2.*sqrt(state->C(1,1));
  float dk = 2.*state->C(2,2);
  float dz0 = 2.*sqrt(state->C(3,3));
  float ddzdl = 2.*sqrt(state->C(4,4));
  
  range.min_phi = seed.phi - dphi;
  range.max_phi = seed.phi + dphi;
  if(range.min_phi < 0.){range.min_phi = 0.;}
  if(range.max_phi > 2.*M_PI){range.max_phi = 2.*M_PI;}
  range.min_d = seed.d - dd;
  range.max_d = seed.d + dd;
  range.min_k = seed.kappa - dk;
  range.max_k = seed.kappa + dk;
  if(range.min_k < 0.){range.min_k = 0.;}
  
  range.min_k = range.min_k * range.min_k;
  range.max_k = range.max_k * range.max_k;
  
  range.min_dzdl = seed.dzdl - ddzdl;
  range.max_dzdl = seed.dzdl + ddzdl;
  range.min_z0 = seed.z0 - dz0;
  range.max_z0 = seed.z0 + dz0;
}


