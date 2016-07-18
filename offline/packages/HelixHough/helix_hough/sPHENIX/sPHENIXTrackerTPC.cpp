#include "sPHENIXTrackerTPC.h"
#include <float.h>
#include <sys/time.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "vector_math_inline.h"

using namespace std;
using namespace Eigen;
using namespace SeamStress;

class hitTriplet {
public:
  hitTriplet(unsigned int h1, unsigned int h2, unsigned int h3, unsigned int t,
             float c)
    : hit1(h1), hit2(h2), hit3(h3), track(t), chi2(c) {}
  ~hitTriplet() {}

  bool operator<(const hitTriplet& other) const {
    return (hit1 < other.hit1) ||
      ((hit2 < other.hit2) && (hit1 == other.hit1)) ||
      ((hit3 < other.hit3) && (hit1 == other.hit1) &&
       (hit2 == other.hit2));
  }

  bool operator==(const hitTriplet& other) const {
    return ((hit1 == other.hit1) && (hit2 == other.hit2) &&
            (hit3 == other.hit3));
  }

  unsigned int hit1, hit2, hit3, track;
  float chi2;
};

void sPHENIXTrackerTPC::tripletRejection(vector<SimpleTrack3D>& input,
                                         vector<SimpleTrack3D>& output,
                                         vector<bool>& usetrack,
                                         vector<float>& next_best_chi2) {
  vector<hitTriplet> trips;
  for (unsigned int i = 0; i < input.size(); ++i) {
    for (unsigned int h1 = 0; h1 < input[i].hits.size(); ++h1) {
      for (unsigned int h2 = (h1 + 1); h2 < input[i].hits.size(); ++h2) {
        for (unsigned int h3 = (h2 + 1); h3 < input[i].hits.size(); ++h3) {
          if (cut_on_dca == false) {
            trips.push_back(
			    hitTriplet(input[i].hits[h1].get_id(), input[i].hits[h2].get_id(),
				       input[i].hits[h3].get_id(), i, track_states[i].chi2));
          }
        }
      }
    }
  }
  if (trips.size() == 0) {
    return;
  }
  sort(trips.begin(), trips.end());
  unsigned int pos = 0;
  unsigned int cur_h1 = trips[pos].hit1;
  unsigned int cur_h2 = trips[pos].hit2;
  while (pos < trips.size()) {
    unsigned int next_pos = pos + 1;
    if (next_pos >= trips.size()) {
      break;
    }
    while (trips[pos] == trips[next_pos]) {
      next_pos += 1;
      if (next_pos >= trips.size()) {
        break;
      }
    }
    if ((next_pos - pos) > 1) {
      float best_chi2 = trips[pos].chi2;
      float next_chi2 = trips[pos + 1].chi2;
      unsigned int best_pos = pos;
      for (unsigned int i = (pos + 1); i < next_pos; ++i) {
        if (input[trips[i].track].hits.size() <
            input[trips[best_pos].track].hits.size()) {
          continue;
        } else if ((input[trips[i].track].hits.size() >
                    input[trips[best_pos].track].hits.size()) ||
                   (input[trips[i].track].hits.back().get_layer() >
                    input[trips[best_pos].track].hits.back().get_layer())) {
          next_chi2 = best_chi2;
          best_chi2 = trips[i].chi2;
          best_pos = i;
          continue;
        }
        if ((trips[i].chi2 < best_chi2) ||
            (usetrack[trips[best_pos].track] == false)) {
          next_chi2 = best_chi2;
          best_chi2 = trips[i].chi2;
          best_pos = i;
        } else if (trips[i].chi2 < next_chi2) {
          next_chi2 = trips[i].chi2;
        }
      }
      for (unsigned int i = pos; i < next_pos; ++i) {
        if (i != best_pos) {
          usetrack[trips[i].track] = false;
        } else {
          next_best_chi2[trips[i].track] = next_chi2;
        }
      }
    }
    pos = next_pos;
    cur_h1 = trips[pos].hit1;
    cur_h2 = trips[pos].hit2;
  }
}

sPHENIXTrackerTPC::sPHENIXTrackerTPC(unsigned int n_phi, unsigned int n_d,
                                     unsigned int n_k, unsigned int n_dzdl,
                                     unsigned int n_z0,
                                     HelixResolution& min_resolution,
                                     HelixResolution& max_resolution,
                                     HelixRange& range, vector<float>& material,
                                     vector<float>& radius, float Bfield)
  : HelixHough(n_phi, n_d, n_k, n_dzdl, n_z0, min_resolution, max_resolution,
	       range),
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
    ca_chi2_cut(2.0),
    cosang_cut(0.985),
    require_pixels(false) {
  vector<float> detector_material;

  for (unsigned int i = 0; i < radius.size(); ++i) {
    detector_radii.push_back(radius[i]);
  }
  for (unsigned int i = 0; i < material.size(); ++i) {
    detector_scatter.push_back(1.41421356237309515 * 0.0136 *
                               sqrt(3. * material[i]));
    detector_material.push_back(3. * material[i]);
  }

  detector_B_field = Bfield;

  integrated_scatter.assign(detector_scatter.size(), 0.);
  float total_scatter_2 = 0.;
  for (unsigned int l = 0; l < detector_scatter.size(); ++l) {
    total_scatter_2 += detector_scatter[l] * detector_scatter[l];
    integrated_scatter[l] = sqrt(total_scatter_2);
  }

  kalman =
    new CylinderKalman(detector_radii, detector_material, detector_B_field);

  vector<SimpleHit3D> one_layer;
  layer_sorted.assign(n_layers, one_layer);
  for (unsigned int i = 0; i < 4; ++i) {
    layer_sorted_1[i].assign(n_layers, one_layer);
  }
  temp_comb.assign(n_layers, 0);
}

sPHENIXTrackerTPC::sPHENIXTrackerTPC(
				     vector<vector<unsigned int> >& zoom_profile, unsigned int minzoom,
				     HelixRange& range, vector<float>& material, vector<float>& radius,
				     float Bfield, bool parallel, unsigned int num_threads)
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
    cosang_cut(0.985),
    require_pixels(false) {
  vector<float> detector_material;

  for (unsigned int i = 0; i < radius.size(); ++i) {
    detector_radii.push_back(radius[i]);
  }
  for (unsigned int i = 0; i < material.size(); ++i) {
    detector_scatter.push_back(1.41421356237309515 * 0.0136 *
                               sqrt(3. * material[i]));
    detector_material.push_back(3. * material[i]);
  }

  detector_B_field = Bfield;

  integrated_scatter.assign(detector_scatter.size(), 0.);
  float total_scatter_2 = 0.;
  for (unsigned int l = 0; l < detector_scatter.size(); ++l) {
    total_scatter_2 += detector_scatter[l] * detector_scatter[l];
    integrated_scatter[l] = sqrt(total_scatter_2);
  }

  kalman =
    new CylinderKalman(detector_radii, detector_material, detector_B_field);

  vector<SimpleHit3D> one_layer;
  layer_sorted.assign(n_layers, one_layer);
  for (unsigned int i = 0; i < 4; ++i) {
    layer_sorted_1[i].assign(n_layers, one_layer);
  }
  temp_comb.assign(n_layers, 0);

  if (is_parallel == true) {
    Seamstress::init_vector(num_threads, vss);

    vssp = new vector<Seamstress*>();
    for (unsigned int i = 0; i < vss.size(); i++) {
      vssp->push_back(&(vss[i]));
    }

    pins = new Pincushion<sPHENIXTrackerTPC>(this, vssp);

    vector<vector<unsigned int> > zoom_profile_new;
    for (unsigned int i = 1; i < zoom_profile.size(); ++i) {
      zoom_profile_new.push_back(zoom_profile[i]);
    }

    for (unsigned int i = 0; i < nthreads; ++i) {
      thread_trackers.push_back(new sPHENIXTrackerTPC(
						      zoom_profile, minzoom, range, material, radius, Bfield));
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

sPHENIXTrackerTPC::~sPHENIXTrackerTPC() {
  if (kalman != NULL) delete kalman;
  for (unsigned int i = 0; i < vss.size(); i++) {
    vss[i].stop();
  }
  for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
    delete thread_trackers[i];
    delete split_output_hits[i];
    delete split_ranges[i];
  }

  if (pins != NULL) delete pins;
  if (vssp != NULL) delete vssp;
}

float sPHENIXTrackerTPC::kappaToPt(float kappa) {
  return detector_B_field / 333.6 / kappa;
}

float sPHENIXTrackerTPC::ptToKappa(float pt) {
  return detector_B_field / 333.6 / pt;
}

// hel should be +- 1
static void xyTangent(SimpleHit3D& hit1, SimpleHit3D& hit2, float kappa,
                      float hel, float& ux_out, float& uy_out, float& ux_in,
                      float& uy_in) {
  float x = hit2.get_x() - hit1.get_x();
  float y = hit2.get_y() - hit1.get_y();
  float D = sqrt(x * x + y * y);
  float ak = 0.5 * kappa * D;
  float D_inv = 1. / D;
  float hk = sqrt(1. - ak * ak);

  float kcx = (ak * x + hel * hk * y) * D_inv;
  float kcy = (ak * y - hel * hk * x) * D_inv;
  float ktx = -(kappa * y - kcy);
  float kty = kappa * x - kcx;
  float norm = 1. / sqrt(ktx * ktx + kty * kty);
  ux_out = ktx * norm;
  uy_out = kty * norm;

  ktx = kcy;
  kty = -kcx;
  norm = 1. / sqrt(ktx * ktx + kty * kty);
  ux_in = ktx * norm;
  uy_in = kty * norm;
}

// hel should be +- 1
static float cosScatter(SimpleHit3D& hit1, SimpleHit3D& hit2, SimpleHit3D& hit3,
                        float kappa, float hel) {
  float ux_in = 0.;
  float uy_in = 0.;
  float ux_out = 0.;
  float uy_out = 0.;

  float temp1 = 0.;
  float temp2 = 0.;

  xyTangent(hit1, hit2, kappa, hel, ux_in, uy_in, temp1, temp2);
  xyTangent(hit2, hit3, kappa, hel, temp1, temp2, ux_out, uy_out);

  return ux_in * ux_out + uy_in * uy_out;
}

static float dzdsSimple(SimpleHit3D& hit1, SimpleHit3D& hit2, float k) {
  float x = hit2.get_x() - hit1.get_x();
  float y = hit2.get_y() - hit1.get_y();
  float D = sqrt(x * x + y * y);
  float s = 0.;
  float temp1 = k * D * 0.5;
  temp1 *= temp1;
  float temp2 = D * 0.5;
  s += 2. * temp2;
  temp2 *= temp1;
  s += temp2 / 3.;
  temp2 *= temp1;
  s += (3. / 20.) * temp2;
  temp2 *= temp1;
  s += (5. / 56.) * temp2;

  return (hit2.get_z() - hit1.get_z()) / s;
}

float sPHENIXTrackerTPC::dcaToVertexXY(SimpleTrack3D& track, float vx,
                                       float vy) {
  float d_out = 0.;

  // find point at the dca to 0
  float x0 = track.d * cos(track.phi);
  float y0 = track.d * sin(track.phi);

  // change variables so x0,y0 -> 0,0
  float phi2 =
    atan2((1. + track.kappa * track.d) * sin(track.phi) - track.kappa * y0,
	  (1. + track.kappa * track.d) * cos(track.phi) - track.kappa * x0);

  // translate so that (0,0) -> (x0 - vx , y0 - vy)
  float cosphi = cos(phi2);
  float sinphi = sin(phi2);
  float tx = cosphi + track.kappa * (x0 - vx);
  float ty = sinphi + track.kappa * (y0 - vy);
  float dk = sqrt(tx * tx + ty * ty) - 1.;
  if (track.kappa == 0.) {
    d_out = (x0 - vx) * cosphi + (y0 - vy) * sinphi;
  } else {
    d_out = dk / track.kappa;
  }
  return fabs(d_out);
}

void sPHENIXTrackerTPC::finalize(vector<SimpleTrack3D>& input,
                                 vector<SimpleTrack3D>& output) {

  if (is_thread == true) {
    for (unsigned int i = 0; i < input.size(); ++i) {
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

  for (unsigned int i = 0; i < nt; ++i) {
    if (usetrack[i] == true) {
      if (!(track_states[i].chi2 == track_states[i].chi2)) {
        continue;
      }

      output.push_back(input[i]);
      output.back().index = (output.size() - 1);
      states_new.push_back(track_states[i]);
      isolation_variable.push_back(next_best_chi2[i]);
    }
  }

  track_states = states_new;
  if (smooth_back == true) {
    for (unsigned int i = 0; i < output.size(); ++i) {
  
      // HelixKalmanState state = track_states[i];

      // track_states[i].C *= 30;

      // track_states[i].chi2 = 0.;
      // track_states[i].x_int = 0.;
      // track_states[i].y_int = 0.;
      // track_states[i].z_int = 0.;
      // track_states[i].position = output[i].hits.size();
      // for(int h=(output[i].hits.size() - 1);h>=0;--h)
      // {
      //   SimpleHit3D hit = output[i].hits[h];
      //   float err_scale = 1.0;
      //   hit.dx *= err_scale;hit.dy *= err_scale;hit.dz *= err_scale;
      //   kalman->addHit(hit, track_states[i]);
      // }

      // SimpleTrack3D temp_track = output[i];
      // fitTrack(temp_track);
      // if( temp_track.kappa == temp_track.kappa )
      // {
      //   track_states[i].kappa = temp_track.kappa;
      //   track_states[i].nu = sqrt(temp_track.kappa);
      // }

      SimpleTrack3D temp_track = output[i];
      vector<SimpleHit3D> temp_hits;
      vector<float> chi2_hit;
      fitTrack(temp_track, chi2_hit);
      for (unsigned int i = 0; i < chi2_hit.size(); ++i) {
	if (chi2_hit[i] < 10. || temp_track.hits[i].get_layer() < 2) {
	  temp_hits.push_back(temp_track.hits[i]);
	}
      }

      temp_track.hits = temp_hits;

      fitTrack(temp_track, chi2_hit);

      if (temp_track.kappa == temp_track.kappa) {
	track_states[i].kappa = temp_track.kappa;
	track_states[i].nu = sqrt(temp_track.kappa);
      }

      HelixKalmanState state = track_states[i];

      track_states[i].C *= 10.;

      track_states[i].chi2 = 0.;
      track_states[i].x_int = 0.;
      track_states[i].y_int = 0.;
      track_states[i].z_int = 0.;
      track_states[i].position = output[i].hits.size();
      for (int h = (temp_track.hits.size() - 1); h >= 0; --h) {
	SimpleHit3D hit = temp_track.hits[h];
	float err_scale = 1.0;

	/// \todo fudge factor location
	
	hit.set_ex( 2.0*sqrt(hit.get_size(0,0)) * err_scale);
	hit.set_ey( 2.0*sqrt(hit.get_size(1,1)) * err_scale);
	hit.set_ez( 2.0*sqrt(hit.get_size(2,2)) * err_scale);
	kalman->addHit(hit, track_states[i]);
      }

      if (fabs(track_states[i].d) < 0.01) {
	SimpleHit3D vertex_hit;
	vertex_hit.set_x(0.0);
	vertex_hit.set_y(0.0);
	vertex_hit.set_z(0.0);
	vertex_hit.set_ex(0.0001);
	vertex_hit.set_ey(0.0001);
	vertex_hit.set_ez(0.0001);
	temp_track.hits.push_back(vertex_hit);

	fitTrack(temp_track, chi2_hit);

	temp_track.hits.pop_back();
      }

      if (temp_track.kappa == temp_track.kappa) {
	track_states[i].kappa = temp_track.kappa;
	track_states[i].nu = sqrt(temp_track.kappa);
      }

      if (!(track_states[i].kappa == track_states[i].kappa)) {
	track_states[i] = state;
	if (temp_track.kappa == temp_track.kappa) {
	  track_states[i].kappa = temp_track.kappa;
	  track_states[i].nu = sqrt(temp_track.kappa);
	}
      }

      if (output[i].phi < 0.) {
	output[i].phi += 2. * M_PI;
      }
      output[i].phi = track_states[i].phi;
      output[i].d = track_states[i].d;
      output[i].kappa = track_states[i].kappa;
      output[i].z0 = track_states[i].z0;
      output[i].dzdl = track_states[i].dzdl;
    }
    
  }

  if (verbosity > 0) {
    cout << "# fits = " << nfits << endl;
    cout << "findTracks called " << findtracksiter << " times" << endl;
    cout << "CAtime = " << CAtime << endl;
    cout << "KALime = " << KALtime << endl;
  }
}

void sPHENIXTrackerTPC::findTracks(vector<SimpleHit3D>& hits,
                                   vector<SimpleTrack3D>& tracks,
                                   const HelixRange& range) {
  findtracksiter += 1;
  findTracksBySegments(hits, tracks, range);

  //   findTracksBySegments(hits,tracks,range);
}

void sPHENIXTrackerTPC::findSeededTracks(vector<SimpleTrack3D>& seeds,
                                         vector<SimpleHit3D>& hits,
                                         vector<SimpleTrack3D>& tracks,
                                         const HelixRange& range) {
  findtracksiter += 1;
  //   findSeededTracksByProjection(seeds, hits, tracks, range);
  findSeededTracksbySegments(seeds, hits, tracks, range);
}

bool sPHENIXTrackerTPC::breakRecursion(const vector<SimpleHit3D>& hits,
                                       const HelixRange& range) {
  if (seeding == true) {
    return false;
  }
  unsigned int layer_mask[4] = {0, 0, 0, 0};
  for (unsigned int i = 0; i < hits.size(); ++i) {
    if (hits[i].get_layer() < 32) {
      layer_mask[0] = layer_mask[0] | (1 << hits[i].get_layer());
    } else if (hits[i].get_layer() < 64) {
      layer_mask[1] = layer_mask[1] | (1 << (hits[i].get_layer() - 32));
    } else if (hits[i].get_layer() < 96) {
      layer_mask[2] = layer_mask[2] | (1 << (hits[i].get_layer() - 64));
    } else if (hits[i].get_layer() < 128) {
      layer_mask[3] = layer_mask[3] | (1 << (hits[i].get_layer() - 96));
    }
  }
  unsigned int nlayers =
    __builtin_popcount(layer_mask[0]) + __builtin_popcount(layer_mask[1]) +
    __builtin_popcount(layer_mask[2]) + __builtin_popcount(layer_mask[3]);

  if (require_pixels == true) {
    if (((layer_mask[0] & 1) == 0) || ((layer_mask[0] & 2) == 0)) {
      return true;
    }
  }

  return (nlayers < required_layers);
}

float sPHENIXTrackerTPC::phiError(SimpleHit3D& hit, float min_k, float max_k,
                                  float min_d, float max_d, float min_z0,
                                  float max_z0, float min_dzdl, float max_dzdl,
                                  bool pairvoting) {
  float Bfield_inv = 1. / detector_B_field;
  float p_inv = 0.;

  if ((prev_max_k == max_k) && (prev_max_dzdl == max_dzdl)) {
    p_inv = prev_p_inv;
  } else {
    prev_max_k = max_k;
    prev_max_dzdl = max_dzdl;
    prev_p_inv = 3.33333333333333314e+02 * max_k * Bfield_inv *
      sqrt(1. - max_dzdl * max_dzdl);
    p_inv = prev_p_inv;
  }
  float total_scatter_2 = 0.;
  for (int i = seed_layer + 1; i <= (hit.get_layer()); ++i) {
    float this_scatter = detector_scatter[i - 1] *
      (detector_radii[i] - detector_radii[i - 1]) /
      detector_radii[i];
    total_scatter_2 += this_scatter * this_scatter;
  }
  float angle = p_inv * sqrt(total_scatter_2) * 1.0;
  float dsize = 0.5 * (max_d - min_d);
  float angle_from_d = dsize / detector_radii[hit.get_layer()];
  float returnval = 0.;
  if (pairvoting == false) {
    if (angle_from_d > angle) {
      returnval = 0.;
    } else {
      returnval = (angle - angle_from_d);
    }
  } else {
    returnval = angle;
  }

  return returnval;
}

float sPHENIXTrackerTPC::dzdlError(SimpleHit3D& hit, float min_k, float max_k,
                                   float min_d, float max_d, float min_z0,
                                   float max_z0, float min_dzdl, float max_dzdl,
                                   bool pairvoting) {
  float Bfield_inv = 1. / detector_B_field;
  float p_inv = 0.;

  if ((prev_max_k == max_k) && (prev_max_dzdl == max_dzdl)) {
    p_inv = prev_p_inv;
  } else {
    prev_max_k = max_k;
    prev_max_dzdl = max_dzdl;
    prev_p_inv = 3.33333333333333314e+02 * max_k * Bfield_inv *
      sqrt(1. - max_dzdl * max_dzdl);
    p_inv = prev_p_inv;
  }
  float total_scatter_2 = 0.;
  for (int i = seed_layer + 1; i <= (hit.get_layer()); ++i) {
    float this_scatter = detector_scatter[i - 1] *
      (detector_radii[i] - detector_radii[i - 1]) /
      detector_radii[i];
    total_scatter_2 += this_scatter * this_scatter;
  }
  float angle = p_inv * sqrt(total_scatter_2) * 1.0;
  float z0size = 0.5 * (max_z0 - min_z0);
  float angle_from_z0 = z0size / detector_radii[hit.get_layer()];
  float returnval = 0.;
  if (pairvoting == false) {
    if (angle_from_z0 > angle) {
      returnval = 0.;
    } else {
      returnval = (angle - angle_from_z0);
    }
  } else {
    returnval = angle;
  }

  return returnval;
}

void sPHENIXTrackerTPC::setRangeFromSeed(HelixRange& range,
                                         SimpleTrack3D& seed) {
  HelixKalmanState* state = &(seed_states[seed.index]);

  float dphi = 2. * sqrt(state->C(0, 0));
  float dd = 2. * sqrt(state->C(1, 1));
  float dk = 2. * state->C(2, 2);
  float dz0 = 2. * sqrt(state->C(3, 3));
  float ddzdl = 2. * sqrt(state->C(4, 4));

  range.min_phi = seed.phi - dphi;
  range.max_phi = seed.phi + dphi;
  if (range.min_phi < 0.) {
    range.min_phi = 0.;
  }
  if (range.max_phi > 2. * M_PI) {
    range.max_phi = 2. * M_PI;
  }
  range.min_d = seed.d - dd;
  range.max_d = seed.d + dd;
  range.min_k = seed.kappa - dk;
  range.max_k = seed.kappa + dk;
  if (range.min_k < 0.) {
    range.min_k = 0.;
  }

  range.min_k = range.min_k * range.min_k;
  range.max_k = range.max_k * range.max_k;

  range.min_dzdl = seed.dzdl - ddzdl;
  range.max_dzdl = seed.dzdl + ddzdl;
  range.min_z0 = seed.z0 - dz0;
  range.max_z0 = seed.z0 + dz0;
}

float sPHENIXTrackerTPC::fitTrack(SimpleTrack3D& track,
				  float scale) {
  vector<float> chi2_hit;
  return sPHENIXTrackerTPC::fitTrack(track, chi2_hit, scale);
}

float sPHENIXTrackerTPC::fitTrack(SimpleTrack3D& track,
                                  vector<float>& chi2_hit,
				  float scale) {
  
  chi2_hit.clear();
  vector<float> xyres;
  vector<float> xyres_inv;
  vector<float> zres;
  vector<float> zres_inv;
  for (unsigned int i = 0; i < track.hits.size(); i++) {

    /// \todo rescale fudge factor
    
    float ex = (2.0*sqrt(track.hits[i].get_size(0,0))) * scale;
    float ey = (2.0*sqrt(track.hits[i].get_size(1,1))) * scale;
    float ez = (2.0*sqrt(track.hits[i].get_size(2,2))) * scale;
        
    xyres.push_back(sqrt( ex * ex + ey * ey ));
    xyres_inv.push_back(1. / xyres.back());
    zres.push_back( ez );
    zres_inv.push_back(1. / zres.back());
  }

  chi2_hit.resize(track.hits.size(), 0.);

  MatrixXf y = MatrixXf::Zero(track.hits.size(), 1);
  for (unsigned int i = 0; i < track.hits.size(); i++) {
    y(i, 0) = (pow(track.hits[i].get_x(), 2) + pow(track.hits[i].get_y(), 2));
    y(i, 0) *= xyres_inv[i];
  }

  MatrixXf X = MatrixXf::Zero(track.hits.size(), 3);
  for (unsigned int i = 0; i < track.hits.size(); i++) {
    X(i, 0) = track.hits[i].get_x();
    X(i, 1) = track.hits[i].get_y();
    X(i, 2) = -1.;
    X(i, 0) *= xyres_inv[i];
    X(i, 1) *= xyres_inv[i];
    X(i, 2) *= xyres_inv[i];
  }

  MatrixXf Xt = X.transpose();

  MatrixXf prod = Xt * X;

  MatrixXf Xty = Xt * y;
  MatrixXf beta = prod.ldlt().solve(Xty);

  float cx = beta(0, 0) * 0.5;
  float cy = beta(1, 0) * 0.5;
  float r = sqrt(cx * cx + cy * cy - beta(2, 0));

  float phi = atan2(cy, cx);
  float d = sqrt(cx * cx + cy * cy) - r;
  float k = 1. / r;

  MatrixXf diff = y - (X * beta);
  MatrixXf chi2 = (diff.transpose()) * diff;

  float dx = d * cos(phi);
  float dy = d * sin(phi);

  MatrixXf y2 = MatrixXf::Zero(track.hits.size(), 1);
  for (unsigned int i = 0; i < track.hits.size(); i++) {
    y2(i, 0) = track.hits[i].get_z();
    y2(i, 0) *= zres_inv[i];
  }

  MatrixXf X2 = MatrixXf::Zero(track.hits.size(), 2);
  for (unsigned int i = 0; i < track.hits.size(); i++) {
    float D = sqrt(pow(dx - track.hits[i].get_x(), 2) + pow(dy - track.hits[i].get_y(), 2));
    float s = 0.0;

    if (0.5 * k * D > 0.1) {
      float v = 0.5 * k * D;
      if (v >= 0.999999) {
        v = 0.999999;
      }
      s = 2. * asin(v) / k;
    } else {
      float temp1 = k * D * 0.5;
      temp1 *= temp1;
      float temp2 = D * 0.5;
      s += 2. * temp2;
      temp2 *= temp1;
      s += temp2 / 3.;
      temp2 *= temp1;
      s += (3. / 20.) * temp2;
      temp2 *= temp1;
      s += (5. / 56.) * temp2;
    }

    X2(i, 0) = s;
    X2(i, 1) = 1.0;

    X2(i, 0) *= zres_inv[i];
    X2(i, 1) *= zres_inv[i];
  }

  MatrixXf Xt2 = X2.transpose();
  MatrixXf prod2 = Xt2 * X2;

  MatrixXf Xty2 = Xt2 * y2;
  MatrixXf beta2 = prod2.ldlt().solve(Xty2);

  MatrixXf diff2 = y2 - (X2 * beta2);
  MatrixXf chi2_z = (diff2.transpose()) * diff2;

  float z0 = beta2(1, 0);
  float dzdl = beta2(0, 0) / sqrt(1. + beta2(0, 0) * beta2(0, 0));

  track.phi = phi;
  track.d = d;
  track.kappa = k;
  track.dzdl = dzdl;
  track.z0 = z0;

  if (track.kappa != 0.) {
    r = 1. / track.kappa;
  } else {
    r = 1.0e10;
  }

  cx = (track.d + r) * cos(track.phi);
  cy = (track.d + r) * sin(track.phi);

  float chi2_tot = 0.;
  for (unsigned int h = 0; h < track.hits.size(); h++) {
    float dx1 = track.hits[h].get_x() - cx;
    float dy1 = track.hits[h].get_y() - cy;

    float dx2 = track.hits[h].get_x() + cx;
    float dy2 = track.hits[h].get_y() + cy;

    float xydiff1 = sqrt(dx1 * dx1 + dy1 * dy1) - r;
    float xydiff2 = sqrt(dx2 * dx2 + dy2 * dy2) - r;
    float xydiff = xydiff2;
    if (fabs(xydiff1) < fabs(xydiff2)) {
      xydiff = xydiff1;
    }

    float ls_xy = xyres[h];

    chi2_hit[h] = 0.;
    chi2_hit[h] += xydiff * xydiff / (ls_xy * ls_xy);
    chi2_hit[h] += diff2(h, 0) * diff2(h, 0);

    chi2_tot += chi2_hit[h];
  }

  unsigned int deg_of_freedom = 2 * track.hits.size() - 5;

  return (chi2_tot) / ((float)(deg_of_freedom));
}


static inline double sign(double x) {
  return ((double)(x > 0.)) - ((double)(x < 0.));
}

void sPHENIXTrackerTPC::projectToLayer(SimpleTrack3D& seed, unsigned int layer,
                                       float& x, float& y, float& z) {
  float phi = seed.phi;
  float d = seed.d;
  float k = seed.kappa;
  float z0 = seed.z0;
  float dzdl = seed.dzdl;

  float hitx = seed.hits.back().get_x();
  float hity = seed.hits.back().get_y();

  float rad_det = detector_radii[layer];

  float cosphi = cos(phi);
  float sinphi = sin(phi);

  k = fabs(k);

  float kd = (d * k + 1.);
  float kcx = kd * cosphi;
  float kcy = kd * sinphi;
  float kd_inv = 1. / kd;
  float R2 = rad_det * rad_det;
  float a = 0.5 * (k * R2 + (d * d * k + 2. * d)) * kd_inv;
  float tmp1 = a * kd_inv;
  float P2x = kcx * tmp1;
  float P2y = kcy * tmp1;

  float h = sqrt(R2 - a * a);

  float ux = -kcy * kd_inv;
  float uy = kcx * kd_inv;

  float x1 = P2x + ux * h;
  float y1 = P2y + uy * h;
  float x2 = P2x - ux * h;
  float y2 = P2y - uy * h;
  float diff1 = (x1 - hitx) * (x1 - hitx) + (y1 - hity) * (y1 - hity);
  float diff2 = (x2 - hitx) * (x2 - hitx) + (y2 - hity) * (y2 - hity);
  float signk = 0.;
  if (diff1 < diff2) {
    signk = 1.;
  } else {
    signk = -1.;
  }
  x = P2x + signk * ux * h;
  y = P2y + signk * uy * h;

  double sign_dzdl = sign(dzdl);
  double onedzdl2_inv = 1. / (1. - dzdl * dzdl);
  double startx = d * cosphi;
  double starty = d * sinphi;
  double D = sqrt((startx - x) * (startx - x) + (starty - y) * (starty - y));
  double D_inv = 1. / D;
  double v = 0.5 * k * D;
  z = 0.;
  if (v > 0.1) {
    if (v >= 0.999999) {
      v = 0.999999;
    }
    double s = 2. * asin(v) / k;
    double s_inv = 1. / s;
    double sqrtvv = sqrt(1 - v * v);
    double dz = sqrt(s * s * dzdl * dzdl / (1. - dzdl * dzdl));
    z = z0 + sign_dzdl * dz;
  } else {
    double s = 0.;
    double temp1 = k * D * 0.5;
    temp1 *= temp1;
    double temp2 = D * 0.5;
    s += 2. * temp2;
    temp2 *= temp1;
    s += temp2 / 3.;
    temp2 *= temp1;
    s += (3. / 20.) * temp2;
    temp2 *= temp1;
    s += (5. / 56.) * temp2;
    double s_inv = 1. / s;
    double dz = sqrt(s * s * dzdl * dzdl / (1. - dzdl * dzdl));
    z = z0 + sign_dzdl * dz;
  }
}


void sPHENIXTrackerTPC::initSplitting(vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits)
{
  initEvent(hits, min_hits);
  (*(hits_vec[0])) = hits;
  zoomranges.clear();
  for(unsigned int z=0;z<=max_zoom;z++)
    {
      zoomranges.push_back(top_range);
    }
}


void sPHENIXTrackerTPC::findHelicesParallelOneHelicity(vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits, vector<SimpleTrack3D>& tracks)
{
  unsigned int hits_per_thread = (hits.size() + 2*nthreads)/nthreads;
  unsigned int pos=0;
  while(pos < hits.size())
    {
      for(unsigned int i=0;i<nthreads;++i)
	{
	  if(pos>=hits.size()){break;}
	  for(unsigned int j=0;j<hits_per_thread;++j)
	    {
	      if(pos>=hits.size()){break;}
	      split_input_hits[i].push_back(hits[pos]);
	      pos+=1;
	    }
	}
    }
  for(unsigned int i=0;i<nthreads;++i)
    {
      thread_trackers[i]->setTopRange(top_range);
      thread_trackers[i]->initSplitting(split_input_hits[i], thread_min_hits, thread_max_hits);
    }
  pins->sewStraight(&sPHENIXTrackerTPC::splitHitsParallelThread, nthreads);
  thread_ranges.clear();
  thread_hits.clear();
  
  unsigned int nbins = split_output_hits[0]->size();
  for(unsigned int b=0;b<nbins;++b)
    {
      thread_ranges.push_back( (*(split_ranges[0]))[b] );
      thread_hits.push_back(vector<SimpleHit3D>());
      for(unsigned int i=0;i<nthreads;++i)
	{
	  for(unsigned int j=0;j<(*(split_output_hits[i]))[b].size();++j)
	    {
	      thread_hits.back().push_back( (*(split_output_hits[i]))[b][j] );
	    }
	}
    }
  
  pins->sewStraight(&sPHENIXTrackerTPC::findHelicesParallelThread, nthreads);
}


void sPHENIXTrackerTPC::findHelicesParallel(vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits, vector<SimpleTrack3D>& tracks)
{
  thread_min_hits = min_hits;
  thread_max_hits = max_hits;
  
  for(unsigned int i=0;i<nthreads;++i)
    {
      thread_tracks[i].clear();
      thread_trackers[i]->clear();
      if(cluster_start_bin!=0){thread_trackers[i]->setClusterStartBin(cluster_start_bin-1);}
      else{thread_trackers[i]->setClusterStartBin(0);}
    }
  
  initSplitting(hits, min_hits, max_hits);
  
  if(separate_by_helicity==true)
    {
      for(unsigned int i=0;i<nthreads;++i)
	{
	  thread_trackers[i]->setSeparateByHelicity(true);
	  thread_trackers[i]->setOnlyOneHelicity(true);
	  thread_trackers[i]->setHelicity(true);
	  split_output_hits[i]->clear();
	  split_input_hits[i].clear();
	}
      findHelicesParallelOneHelicity(hits, min_hits, max_hits, tracks);
    
      for(unsigned int i=0;i<nthreads;++i)
	{
	  thread_trackers[i]->setSeparateByHelicity(true);
	  thread_trackers[i]->setOnlyOneHelicity(true);
	  thread_trackers[i]->setHelicity(false);
	  split_output_hits[i]->clear();
	  split_input_hits[i].clear();
	}
      findHelicesParallelOneHelicity(hits, min_hits, max_hits, tracks);
    }
  else
    {
      for(unsigned int i=0;i<nthreads;++i)
	{
	  thread_trackers[i]->setSeparateByHelicity(false);
	  thread_trackers[i]->setOnlyOneHelicity(false);
	  split_output_hits[i]->clear();
	  split_input_hits[i].clear();
	}
    
      findHelicesParallelOneHelicity(hits, min_hits, max_hits, tracks);
    }
  
  vector<SimpleTrack3D> temp_tracks;
  for(unsigned int i=0;i<nthreads;++i)
    {
      vector<HelixKalmanState>* states = &(thread_trackers[i]->getKalmanStates());
      for(unsigned int j=0;j<thread_tracks[i].size();++j)
	{
	  track_states.push_back((*states)[j]);
	  temp_tracks.push_back(thread_tracks[i][j]);
	}
    }
  finalize(temp_tracks, tracks);
}


void sPHENIXTrackerTPC::splitHitsParallelThread(void* arg)
{
  unsigned long int w = (*((unsigned long int *)arg));
  thread_trackers[w]->splitIntoBins(thread_min_hits, thread_max_hits, *(split_ranges[w]), *(split_output_hits[w]), 0);
}


void sPHENIXTrackerTPC::findHelicesParallelThread(void* arg)
{
  unsigned long int w = (*((unsigned long int *)arg));
  
  for(unsigned int i=w;i<thread_ranges.size();i+=nthreads)
    {
      if(thread_hits[i].size() == 0){continue;}
      thread_trackers[w]->setTopRange(thread_ranges[i]);
      thread_trackers[w]->findHelices(thread_hits[i], thread_min_hits, thread_max_hits, thread_tracks[w]);
    }
}



void sPHENIXTrackerTPC::calculateKappaTangents(
					       float* x1_a, float* y1_a, float* z1_a, float* x2_a, float* y2_a,
					       float* z2_a, float* x3_a, float* y3_a, float* z3_a, float* dx1_a,
					       float* dy1_a, float* dz1_a, float* dx2_a, float* dy2_a, float* dz2_a,
					       float* dx3_a, float* dy3_a, float* dz3_a, float* kappa_a, float* dkappa_a,
					       float* ux_mid_a, float* uy_mid_a, float* ux_end_a, float* uy_end_a,
					       float* dzdl_1_a, float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a) {
  static const __m128 two = {2., 2., 2., 2.};

  __m128 x1 = _mm_load_ps(x1_a);
  __m128 x2 = _mm_load_ps(x2_a);
  __m128 x3 = _mm_load_ps(x3_a);
  __m128 y1 = _mm_load_ps(y1_a);
  __m128 y2 = _mm_load_ps(y2_a);
  __m128 y3 = _mm_load_ps(y3_a);
  __m128 z1 = _mm_load_ps(z1_a);
  __m128 z2 = _mm_load_ps(z2_a);
  __m128 z3 = _mm_load_ps(z3_a);

  __m128 dx1 = _mm_load_ps(dx1_a);
  __m128 dx2 = _mm_load_ps(dx2_a);
  __m128 dx3 = _mm_load_ps(dx3_a);
  __m128 dy1 = _mm_load_ps(dy1_a);
  __m128 dy2 = _mm_load_ps(dy2_a);
  __m128 dy3 = _mm_load_ps(dy3_a);
  __m128 dz1 = _mm_load_ps(dz1_a);
  __m128 dz2 = _mm_load_ps(dz2_a);
  __m128 dz3 = _mm_load_ps(dz3_a);

  __m128 D12 = _mm_sub_ps(x2, x1);
  D12 = _mm_mul_ps(D12, D12);
  __m128 tmp1 = _mm_sub_ps(y2, y1);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D12 = _mm_add_ps(D12, tmp1);
  D12 = _vec_sqrt_ps(D12);

  __m128 D23 = _mm_sub_ps(x3, x2);
  D23 = _mm_mul_ps(D23, D23);
  tmp1 = _mm_sub_ps(y3, y2);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D23 = _mm_add_ps(D23, tmp1);
  D23 = _vec_sqrt_ps(D23);

  __m128 D31 = _mm_sub_ps(x1, x3);
  D31 = _mm_mul_ps(D31, D31);
  tmp1 = _mm_sub_ps(y1, y3);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D31 = _mm_add_ps(D31, tmp1);
  D31 = _vec_sqrt_ps(D31);

  __m128 k = _mm_mul_ps(D12, D23);
  k = _mm_mul_ps(k, D31);
  k = _vec_rec_ps(k);
  tmp1 = (D12 + D23 + D31) * (D23 + D31 - D12) * (D12 + D31 - D23) *
    (D12 + D23 - D31);
  tmp1 = _vec_sqrt_ps(tmp1);
  k *= tmp1;

  __m128 tmp2 = _mm_cmpgt_ps(tmp1, zero);
  tmp1 = _mm_and_ps(tmp2, k);
  tmp2 = _mm_andnot_ps(tmp2, zero);
  k = _mm_xor_ps(tmp1, tmp2);

  _mm_store_ps(kappa_a, k);
  __m128 k_inv = _vec_rec_ps(k);

  __m128 D12_inv = _vec_rec_ps(D12);
  __m128 D23_inv = _vec_rec_ps(D23);
  __m128 D31_inv = _vec_rec_ps(D31);

  __m128 dr1 = dx1 * dx1 + dy1 * dy1;
  dr1 = _vec_sqrt_ps(dr1);
  __m128 dr2 = dx2 * dx2 + dy2 * dy2;
  dr2 = _vec_sqrt_ps(dr2);
  __m128 dr3 = dx3 * dx3 + dy3 * dy3;
  dr3 = _vec_sqrt_ps(dr3);

  __m128 dk1 = (dr1 + dr2) * D12_inv * D12_inv;
  __m128 dk2 = (dr2 + dr3) * D23_inv * D23_inv;
  __m128 dk = dk1 + dk2;
  _mm_store_ps(dkappa_a, dk);

  __m128 ux12 = (x2 - x1) * D12_inv;
  __m128 uy12 = (y2 - y1) * D12_inv;
  __m128 ux23 = (x3 - x2) * D23_inv;
  __m128 uy23 = (y3 - y2) * D23_inv;
  __m128 ux13 = (x3 - x1) * D31_inv;
  __m128 uy13 = (y3 - y1) * D31_inv;

  __m128 cosalpha = ux12 * ux13 + uy12 * uy13;
  __m128 sinalpha = ux13 * uy12 - ux12 * uy13;

  __m128 ux_mid = ux23 * cosalpha - uy23 * sinalpha;
  __m128 uy_mid = ux23 * sinalpha + uy23 * cosalpha;
  _mm_store_ps(ux_mid_a, ux_mid);
  _mm_store_ps(uy_mid_a, uy_mid);

  __m128 ux_end = ux23 * cosalpha + uy23 * sinalpha;
  __m128 uy_end = uy23 * cosalpha - ux23 * sinalpha;

  _mm_store_ps(ux_end_a, ux_end);
  _mm_store_ps(uy_end_a, uy_end);

  // asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  __m128 v = one - sinalpha * sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s2 = _vec_atan_ps(v);
  s2 *= two;
  s2 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s2);
  tmp1 = _mm_andnot_ps(tmp1, D23);
  s2 = _mm_xor_ps(tmp1, tmp2);

  // dz/dl = (dz/ds)/sqrt(1 + (dz/ds)^2)
  // = dz/sqrt(s^2 + dz^2)
  __m128 del_z_2 = z3 - z2;
  __m128 dzdl_2 = s2 * s2 + del_z_2 * del_z_2;
  dzdl_2 = _vec_rsqrt_ps(dzdl_2);
  dzdl_2 *= del_z_2;
  __m128 ddzdl_2 = (dz2 + dz3) * D23_inv;
  _mm_store_ps(dzdl_2_a, dzdl_2);
  _mm_store_ps(ddzdl_2_a, ddzdl_2);

  sinalpha = ux13 * uy23 - ux23 * uy13;
  v = one - sinalpha * sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s1 = _vec_atan_ps(v);
  s1 *= two;
  s1 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s1);
  tmp1 = _mm_andnot_ps(tmp1, D12);
  s1 = _mm_xor_ps(tmp1, tmp2);

  __m128 del_z_1 = z2 - z1;
  __m128 dzdl_1 = s1 * s1 + del_z_1 * del_z_1;
  dzdl_1 = _vec_rsqrt_ps(dzdl_1);
  dzdl_1 *= del_z_1;
  __m128 ddzdl_1 = (dz1 + dz2) * D12_inv;
  _mm_store_ps(dzdl_1_a, dzdl_1);
  _mm_store_ps(ddzdl_1_a, ddzdl_1);
}

void sPHENIXTrackerTPC::calculateKappaTangents(
					       float* x1_a, float* y1_a, float* z1_a, float* x2_a, float* y2_a,
					       float* z2_a, float* x3_a, float* y3_a, float* z3_a, float* dx1_a,
					       float* dy1_a, float* dz1_a, float* dx2_a, float* dy2_a, float* dz2_a,
					       float* dx3_a, float* dy3_a, float* dz3_a, float* kappa_a, float* dkappa_a,
					       float* ux_mid_a, float* uy_mid_a, float* ux_end_a, float* uy_end_a,
					       float* dzdl_1_a, float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a,
					       float sinang_cut, float cosang_diff_inv, float* cur_kappa_a,
					       float* cur_dkappa_a, float* cur_ux_a, float* cur_uy_a, float* cur_chi2_a,
					       float* chi2_a) {
  static const __m128 two = {2., 2., 2., 2.};

  __m128 x1 = _mm_load_ps(x1_a);
  __m128 x2 = _mm_load_ps(x2_a);
  __m128 x3 = _mm_load_ps(x3_a);
  __m128 y1 = _mm_load_ps(y1_a);
  __m128 y2 = _mm_load_ps(y2_a);
  __m128 y3 = _mm_load_ps(y3_a);
  __m128 z1 = _mm_load_ps(z1_a);
  __m128 z2 = _mm_load_ps(z2_a);
  __m128 z3 = _mm_load_ps(z3_a);

  __m128 dx1 = _mm_load_ps(dx1_a);
  __m128 dx2 = _mm_load_ps(dx2_a);
  __m128 dx3 = _mm_load_ps(dx3_a);
  __m128 dy1 = _mm_load_ps(dy1_a);
  __m128 dy2 = _mm_load_ps(dy2_a);
  __m128 dy3 = _mm_load_ps(dy3_a);
  __m128 dz1 = _mm_load_ps(dz1_a);
  __m128 dz2 = _mm_load_ps(dz2_a);
  __m128 dz3 = _mm_load_ps(dz3_a);

  __m128 D12 = _mm_sub_ps(x2, x1);
  D12 = _mm_mul_ps(D12, D12);
  __m128 tmp1 = _mm_sub_ps(y2, y1);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D12 = _mm_add_ps(D12, tmp1);
  D12 = _vec_sqrt_ps(D12);

  __m128 D23 = _mm_sub_ps(x3, x2);
  D23 = _mm_mul_ps(D23, D23);
  tmp1 = _mm_sub_ps(y3, y2);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D23 = _mm_add_ps(D23, tmp1);
  D23 = _vec_sqrt_ps(D23);

  __m128 D31 = _mm_sub_ps(x1, x3);
  D31 = _mm_mul_ps(D31, D31);
  tmp1 = _mm_sub_ps(y1, y3);
  tmp1 = _mm_mul_ps(tmp1, tmp1);
  D31 = _mm_add_ps(D31, tmp1);
  D31 = _vec_sqrt_ps(D31);

  __m128 k = _mm_mul_ps(D12, D23);
  k = _mm_mul_ps(k, D31);
  k = _vec_rec_ps(k);
  tmp1 = (D12 + D23 + D31) * (D23 + D31 - D12) * (D12 + D31 - D23) *
    (D12 + D23 - D31);
  tmp1 = _vec_sqrt_ps(tmp1);
  k *= tmp1;

  __m128 tmp2 = _mm_cmpgt_ps(tmp1, zero);
  tmp1 = _mm_and_ps(tmp2, k);
  tmp2 = _mm_andnot_ps(tmp2, zero);
  k = _mm_xor_ps(tmp1, tmp2);

  _mm_store_ps(kappa_a, k);
  __m128 k_inv = _vec_rec_ps(k);

  __m128 D12_inv = _vec_rec_ps(D12);
  __m128 D23_inv = _vec_rec_ps(D23);
  __m128 D31_inv = _vec_rec_ps(D31);

  __m128 dr1 = dx1 * dx1 + dy1 * dy1;
  dr1 = _vec_sqrt_ps(dr1);
  __m128 dr2 = dx2 * dx2 + dy2 * dy2;
  dr2 = _vec_sqrt_ps(dr2);
  __m128 dr3 = dx3 * dx3 + dy3 * dy3;
  dr3 = _vec_sqrt_ps(dr3);

  __m128 dk1 = (dr1 + dr2) * D12_inv * D12_inv;
  __m128 dk2 = (dr2 + dr3) * D23_inv * D23_inv;
  __m128 dk = dk1 + dk2;
  _mm_store_ps(dkappa_a, dk);

  __m128 ux12 = (x2 - x1) * D12_inv;
  __m128 uy12 = (y2 - y1) * D12_inv;
  __m128 ux23 = (x3 - x2) * D23_inv;
  __m128 uy23 = (y3 - y2) * D23_inv;
  __m128 ux13 = (x3 - x1) * D31_inv;
  __m128 uy13 = (y3 - y1) * D31_inv;

  __m128 cosalpha = ux12 * ux13 + uy12 * uy13;
  __m128 sinalpha = ux13 * uy12 - ux12 * uy13;

  __m128 ux_mid = ux23 * cosalpha - uy23 * sinalpha;
  __m128 uy_mid = ux23 * sinalpha + uy23 * cosalpha;
  _mm_store_ps(ux_mid_a, ux_mid);
  _mm_store_ps(uy_mid_a, uy_mid);

  __m128 ux_end = ux23 * cosalpha + uy23 * sinalpha;
  __m128 uy_end = uy23 * cosalpha - ux23 * sinalpha;

  _mm_store_ps(ux_end_a, ux_end);
  _mm_store_ps(uy_end_a, uy_end);

  // asin(x) = 2*atan( x/( 1 + sqrt( 1 - x*x ) ) )
  __m128 v = one - sinalpha * sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s2 = _vec_atan_ps(v);
  s2 *= two;
  s2 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s2);
  tmp1 = _mm_andnot_ps(tmp1, D23);
  s2 = _mm_xor_ps(tmp1, tmp2);

  // dz/dl = (dz/ds)/sqrt(1 + (dz/ds)^2)
  // = dz/sqrt(s^2 + dz^2)
  __m128 del_z_2 = z3 - z2;
  __m128 dzdl_2 = s2 * s2 + del_z_2 * del_z_2;
  dzdl_2 = _vec_rsqrt_ps(dzdl_2);
  dzdl_2 *= del_z_2;
  __m128 ddzdl_2 = (dz2 + dz3) * D23_inv;
  _mm_store_ps(dzdl_2_a, dzdl_2);
  _mm_store_ps(ddzdl_2_a, ddzdl_2);

  sinalpha = ux13 * uy23 - ux23 * uy13;
  v = one - sinalpha * sinalpha;
  v = _vec_sqrt_ps(v);
  v += one;
  v = _vec_rec_ps(v);
  v *= sinalpha;
  __m128 s1 = _vec_atan_ps(v);
  s1 *= two;
  s1 *= k_inv;
  tmp1 = _mm_cmpgt_ps(k, zero);
  tmp2 = _mm_and_ps(tmp1, s1);
  tmp1 = _mm_andnot_ps(tmp1, D12);
  s1 = _mm_xor_ps(tmp1, tmp2);

  __m128 del_z_1 = z2 - z1;
  __m128 dzdl_1 = s1 * s1 + del_z_1 * del_z_1;
  dzdl_1 = _vec_rsqrt_ps(dzdl_1);
  dzdl_1 *= del_z_1;
  __m128 ddzdl_1 = (dz1 + dz2) * D12_inv;
  _mm_store_ps(dzdl_1_a, dzdl_1);
  _mm_store_ps(ddzdl_1_a, ddzdl_1);

  __m128 c_dk = _mm_load_ps(cur_dkappa_a);
  __m128 c_k = _mm_load_ps(cur_kappa_a);
  __m128 c_ux = _mm_load_ps(cur_ux_a);
  __m128 c_uy = _mm_load_ps(cur_uy_a);
  __m128 c_chi2 = _mm_load_ps(cur_chi2_a);
  __m128 sinang = _mm_load1_ps(&sinang_cut);
  __m128 cosdiff = _mm_load1_ps(&cosang_diff_inv);

  __m128 kdiff = c_k - k;
  __m128 n_dk = c_dk + dk + sinang * k;
  __m128 chi2_k = kdiff * kdiff / (n_dk * n_dk);
  __m128 cos_scatter = c_ux * ux_mid + c_uy * uy_mid;
  __m128 chi2_ang =
    (one - cos_scatter) * (one - cos_scatter) * cosdiff * cosdiff;
  tmp1 = dzdl_1 * sinang;
  _vec_fabs_ps(tmp1);
  __m128 chi2_dzdl = (dzdl_1 - dzdl_2) / (ddzdl_1 + ddzdl_2 + tmp1);
  chi2_dzdl *= chi2_dzdl;
  chi2_dzdl *= one_o_2;

  __m128 n_chi2 = c_chi2 + chi2_ang + chi2_k + chi2_dzdl;
  _mm_store_ps(chi2_a, n_chi2);
}

struct TempComb {
  TempComb() {}
  HelixKalmanState state;
  SimpleTrack3D track;

  inline bool operator<(const TempComb& other) const {
    if (track.hits.size() > other.track.hits.size()) {
      return true;
    } else if (track.hits.size() == other.track.hits.size()) {
      return state.chi2 < other.state.chi2;
    } else {
      return false;
    }
  }
};

void sPHENIXTrackerTPC::initDummyHits(vector<SimpleHit3D>& dummies,
                                      const HelixRange& range,
                                      HelixKalmanState& init_state) {
  SimpleTrack3D dummy_track;
  dummy_track.hits.push_back(SimpleHit3D());
  dummy_track.kappa = 0.5 * (range.min_k + range.max_k);
  dummy_track.phi = 0.5 * (range.min_phi + range.max_phi);
  dummy_track.d = 0.5 * (range.min_d + range.max_d);
  dummy_track.dzdl = 0.5 * (range.min_dzdl + range.max_dzdl);
  dummy_track.z0 = 0.5 * (range.min_z0 + range.max_z0);

  init_state.kappa = dummy_track.kappa;
  init_state.nu = sqrt(dummy_track.kappa);
  init_state.phi = dummy_track.phi;
  init_state.d = dummy_track.d;
  init_state.dzdl = dummy_track.dzdl;
  init_state.z0 = dummy_track.z0;

  init_state.C = Matrix<float, 5, 5>::Zero(5, 5);
  init_state.C(0, 0) = pow(range.max_phi - range.min_phi, 2.);
  init_state.C(1, 1) = pow(range.max_d - range.min_d, 2.);
  init_state.C(2, 2) = pow(10. * sqrt(range.max_k - range.min_k), 2.);
  init_state.C(3, 3) = pow(range.max_z0 - range.min_z0, 2.);
  init_state.C(4, 4) = pow(range.max_dzdl - range.min_dzdl, 2.);
  init_state.chi2 = 0.;
  init_state.position = 0;
  init_state.x_int = 0.;
  init_state.y_int = 0.;
  init_state.z_int = 0.;

  for (unsigned int i = 0; i < n_layers; ++i) {
    float x, y, z;
    projectToLayer(dummy_track, i, x, y, z);
    dummies[i].set_x(x);
    dummies[i].set_ex(5.);
    dummies[i].set_y(x);
    dummies[i].set_ey(5.);
    dummies[i].set_z(x);
    dummies[i].set_ez(5.);
    dummies[i].set_layer(i);
  }
}

static bool next_combo_n(vector<int> const& lsizes, vector<int>& comb_n) {
  unsigned int n = lsizes.size() / 2;
  for (int l = 0; l < n; ++l) {
    if (comb_n[l] == (lsizes[l] - 1)) {
      comb_n[l] = 0;
    } else {
      comb_n[l] += 1;
      return true;
    }
  }
  return false;
}

static SimpleHit3D& get_hit(vector<SimpleHit3D>& hits,
                            vector<SimpleHit3D>& dummies, int index) {
  if (index >= 0) {
    return hits[index];
  } else {
    return dummies[(-index) - 1];
  }
}

class hit_triplet {
public:
  hit_triplet(unsigned int h1, unsigned int h2, unsigned int h3, unsigned int t,
              float c)
    : hit1(h1), hit2(h2), hit3(h3), track(t), chi2(c) {}
  ~hit_triplet() {}

  bool operator<(const hit_triplet& other) const {
    return (hit1 < other.hit1) ||
      ((hit2 < other.hit2) && (hit1 == other.hit1)) ||
      ((hit3 < other.hit3) && (hit1 == other.hit1) &&
       (hit2 == other.hit2));
  }

  bool operator==(const hit_triplet& other) const {
    return ((hit1 == other.hit1) && (hit2 == other.hit2) &&
            (hit3 == other.hit3));
  }

  unsigned int hit1, hit2, hit3, track;
  float chi2;
};

static void triplet_rejection(vector<SimpleTrack3D>& input,
                              vector<float>& chi2s, vector<bool>& usetrack) {
  vector<hit_triplet> trips;
  for (unsigned int i = 0; i < input.size(); ++i) {
    for (unsigned int h1 = 0; h1 < input[i].hits.size(); ++h1) {
      for (unsigned int h2 = (h1 + 1); h2 < input[i].hits.size(); ++h2) {
        for (unsigned int h3 = (h2 + 1); h3 < input[i].hits.size(); ++h3) {
          trips.push_back(hit_triplet(input[i].hits[h1].get_id(),
                                      input[i].hits[h2].get_id(),
                                      input[i].hits[h3].get_id(), i, chi2s[i]));
        }
      }
    }
  }
  if (trips.size() == 0) {
    return;
  }
  sort(trips.begin(), trips.end());
  unsigned int pos = 0;
  unsigned int cur_h1 = trips[pos].hit1;
  unsigned int cur_h2 = trips[pos].hit2;
  while (pos < trips.size()) {
    unsigned int next_pos = pos + 1;
    if (next_pos >= trips.size()) {
      break;
    }
    while (trips[pos] == trips[next_pos]) {
      next_pos += 1;
      if (next_pos >= trips.size()) {
        break;
      }
    }
    if ((next_pos - pos) > 1) {
      float best_chi2 = trips[pos].chi2;
      float next_chi2 = trips[pos + 1].chi2;
      unsigned int best_pos = pos;
      for (unsigned int i = (pos + 1); i < next_pos; ++i) {
        if (input[trips[i].track].hits.size() <
            input[trips[best_pos].track].hits.size()) {
          continue;
        } else if ((input[trips[i].track].hits.size() >
                    input[trips[best_pos].track].hits.size()) ||
                   (input[trips[i].track].hits.back().get_layer() >
                    input[trips[best_pos].track].hits.back().get_layer())) {
          next_chi2 = best_chi2;
          best_chi2 = trips[i].chi2;
          best_pos = i;
          continue;
        }
        if ((trips[i].chi2 < best_chi2) ||
            (usetrack[trips[best_pos].track] == false)) {
          next_chi2 = best_chi2;
          best_chi2 = trips[i].chi2;
          best_pos = i;
        } else if (trips[i].chi2 < next_chi2) {
          next_chi2 = trips[i].chi2;
        }
      }
      for (unsigned int i = pos; i < next_pos; ++i) {
        if (i != best_pos) {
          usetrack[trips[i].track] = false;
        }
      }
    }
    pos = next_pos;
    cur_h1 = trips[pos].hit1;
    cur_h2 = trips[pos].hit2;
  }
}

static bool remove_bad_hits(SimpleTrack3D& track, float cut, float scale = 1.0) {
  SimpleTrack3D temp_track = track;
  float fit_chi2 = 0.;
  vector<float> chi2_hit;
  vector<float> temp_hits;
  while (true) {
    temp_track = track;
    fit_chi2 = sPHENIXTrackerTPC::fitTrack(temp_track, chi2_hit, scale);
    bool all_good = true;
    track.hits.clear();
    for (int h = 0; h < temp_track.hits.size(); h += 1) {
      if (chi2_hit[h] < cut) {
        track.hits.push_back(temp_track.hits[h]);
      } else {
        all_good = false;
      }
    }
    if (track.hits.size() < 3) {
      return false;
    }
    if (all_good == true) {
      return true;
    }
  }
}

static void initial_combos(int nhits, vector<SimpleHit3D>& hits,
                           vector<vector<int> >& layer_indexes,
                           vector<TempComb>& cur_comb, float CHI2_CUT,
                           int n_layers, CylinderKalman* kalman) {
  vector<int> lsizes;
  for (int i = 0; i < nhits; ++i) {
    lsizes.push_back(layer_indexes[i].size());
  }
  vector<int> comb_n;
  comb_n.assign(nhits, 0);
  vector<TempComb> comb;
  while (true) {
    SimpleTrack3D temp_track;
    TempComb tc;
    float sqrt12_inv = 1. / sqrt(12.);
    for (unsigned int h = 0; h < nhits; ++h) {
      if (layer_indexes[h][comb_n[h]] < 0) {
        continue;
      }
      SimpleHit3D& hit = hits[layer_indexes[h][comb_n[h]]];
      temp_track.hits.push_back(hit);
      temp_track.hits.back().set_ex( temp_track.hits.back().get_ex() * sqrt12_inv);
      temp_track.hits.back().set_ey( temp_track.hits.back().get_ey() * sqrt12_inv);
      temp_track.hits.back().set_ez( temp_track.hits.back().get_ez() * sqrt12_inv);
    }
    if (temp_track.hits.size() >= 3) {
      float init_chi2 = sPHENIXTrackerTPC::fitTrack(temp_track);
      comb.push_back(TempComb());
      TempComb& curcomb = comb.back();
      if ((init_chi2 / (2. * (temp_track.hits.size()) - 5.) > CHI2_CUT) ||
          (!(temp_track.hits[0].get_layer() == 0 && temp_track.hits[1].get_layer() == 1))) {
        comb.pop_back();
      } else {
        curcomb.state.phi = temp_track.phi;
        if (curcomb.state.phi < 0.) {
          curcomb.state.phi += 2. * M_PI;
        }
        curcomb.state.d = temp_track.d;
        curcomb.state.kappa = temp_track.kappa;
        curcomb.state.nu = sqrt(curcomb.state.kappa);
        curcomb.state.z0 = temp_track.z0;
        curcomb.state.dzdl = temp_track.dzdl;
        curcomb.state.C = Matrix<float, 5, 5>::Zero(5, 5);
        curcomb.state.C(0, 0) = pow(0.01, 2.);
        curcomb.state.C(1, 1) = pow(0.5, 2.);
        curcomb.state.C(2, 2) = pow(0.3 * curcomb.state.nu, 2.);
        curcomb.state.C(3, 3) = pow(0.5, 2.);
        curcomb.state.C(4, 4) = pow(0.05, 2.);
        curcomb.state.chi2 = 0.;
        curcomb.state.position = n_layers;
        curcomb.state.x_int = 0.;
        curcomb.state.y_int = 0.;
        curcomb.state.z_int = 0.;

        curcomb.track = temp_track;

        for (int h = 0; h < curcomb.track.hits.size(); h += 1) {
          kalman->addHit(curcomb.track.hits[h], curcomb.state);
        }
      }
    }
    if (next_combo_n(lsizes, comb_n) == false) {
      break;
    }
  }

  vector<SimpleTrack3D> comb_tracks;
  vector<TempComb> comb2;
  vector<float> chi2s;
  for (int t = 0; t < comb.size(); t += 1) {
    comb_tracks.push_back(comb[t].track);
    chi2s.push_back(comb[t].state.chi2);
    comb2.push_back(comb[t]);
  }
  vector<bool> usetrack(comb2.size(), true);

  triplet_rejection(comb_tracks, chi2s, usetrack);

  for (int t = 0; t < comb_tracks.size(); t += 1) {
    if (usetrack[t]) {
      if (comb2[t].state.chi2 / (2. * (comb2[t].track.hits.size()) - 5.) <
          CHI2_CUT) {
        cur_comb.push_back(comb2[t]);
      }
    }
  }

  sort(cur_comb.begin(), cur_comb.end());
  if (cur_comb.size() > 4) {
    cur_comb.resize(4, TempComb());
  }
}

static void extend_combos(vector<SimpleHit3D>& hits,
                          vector<vector<int> >& layer_indexes,
                          vector<TempComb>& cur_comb, float CHI2_CUT,
                          CylinderKalman* kalman, int layer_begin,
                          int layer_end) {
  vector<TempComb> comb;
  vector<TempComb> comb2;
  vector<TempComb> comb3;
  for (int c = 0; c < cur_comb.size(); c += 1) {
    comb2.clear();
    comb2.push_back(cur_comb[c]);
    for (int l = layer_begin; l <= layer_end; l += 1) {
      comb3.clear();
      for (int i = 0; i < comb2.size(); i += 1) {
        comb3.push_back(comb2[i]);
        for (int h = 0; h < (layer_indexes[l].size() - 1); h += 1) {
          SimpleHit3D hit = hits[layer_indexes[l][h]];
          hit.set_ex( hit.get_ex() / sqrt(12.));
          hit.set_ey( hit.get_ey() / sqrt(12.));
          hit.set_ez( hit.get_ez() / sqrt(12.));
          comb3.push_back(comb2[i]);
          kalman->addHit(hit, comb3.back().state);
          if (comb3.back().state.chi2 >
              CHI2_CUT * (2. * (comb3.back().track.hits.size()) - 5.)) {
            comb3.pop_back();
          } else {
            comb3.back().track.hits.push_back(hit);
          }
        }
      }
      comb2 = comb3;
    }
    for (int i = 0; i < comb2.size(); i += 1) {
      comb.push_back(comb2[i]);
    }
  }
  cur_comb.clear();

  vector<SimpleTrack3D> comb_tracks;
  comb2.clear();
  vector<float> chi2s;
  for (int t = 0; t < comb.size(); t += 1) {
    if (remove_bad_hits(comb[t].track, CHI2_CUT + 2.) == false) {
      continue;
    }
    comb_tracks.push_back(comb[t].track);
    chi2s.push_back(comb[t].state.chi2);
    comb2.push_back(comb[t]);
  }
  vector<bool> usetrack(comb2.size(), true);

  triplet_rejection(comb_tracks, chi2s, usetrack);

  for (int t = 0; t < comb_tracks.size(); t += 1) {
    if (usetrack[t]) {
      if (comb2[t].state.chi2 / (2. * (comb2[t].track.hits.size()) - 5.) <
          CHI2_CUT) {
        cur_comb.push_back(comb2[t]);
      }
    }
  }

  sort(cur_comb.begin(), cur_comb.end());
  if (cur_comb.size() > 4) {
    cur_comb.resize(4, TempComb());
  }
}

static bool fit_all_update(vector<vector<int> >& layer_indexes,
                           SimpleTrack3D& temp_track, vector<int>& best_ind,
                           vector<float>& best_chi2, SimpleTrack3D& track,
                           float& chi2, int iter = 0,
			   float tempscale = 1.0, float scale = 1.0) {
  if (iter != 0) {
    if (remove_bad_hits(track, 5., scale) == false) {
      return false;
    }
    sPHENIXTrackerTPC::fitTrack(track, scale);
  }

  vector<float> chi2_hit;
  float r = 1. / track.kappa;
  float cx = (track.d + r) * cos(track.phi);
  float cy = (track.d + r) * sin(track.phi);

  float phi = atan2(cy, cx);
  float d = sqrt(cx * cx + cy * cy) - r;
  float dx = d * cos(phi);
  float dy = d * sin(phi);

  float k = track.kappa;

  best_ind.assign(layer_indexes.size(), -1);
  best_chi2.assign(layer_indexes.size(), -1.);
  for (int i = 0; i < temp_track.hits.size(); i += 1) {
    float dx1 = temp_track.hits[i].get_x() - cx;
    float dy1 = temp_track.hits[i].get_y() - cy;
    float dx2 = temp_track.hits[i].get_x() + cx;
    float dy2 = temp_track.hits[i].get_y() + cy;
    float xydiff1 = sqrt(dx1 * dx1 + dy1 * dy1) - r;
    float xydiff2 = sqrt(dx2 * dx2 + dy2 * dy2) - r;
    float xydiff = xydiff2;
    if (fabs(xydiff1) < fabs(xydiff2)) {
      xydiff = xydiff1;
    }

    float tex = temp_track.hits[i].get_ex() * tempscale;
    float tey = temp_track.hits[i].get_ey() * tempscale;
    float tez = temp_track.hits[i].get_ez() * tempscale;
    
    float dr2 = tex * tex + tey * tey;
    float chi2_xy = xydiff * xydiff / dr2;

    float D = sqrt(pow(dx - temp_track.hits[i].get_x(), 2) +
                   pow(dy - temp_track.hits[i].get_y(), 2));
    float s = 0.0;

    if (0.5 * k * D > 0.1) {
      float v = 0.5 * k * D;
      if (v >= 0.999999) {
        v = 0.999999;
      }
      s = 2. * asin(v) / k;
    } else {
      float temp1 = k * D * 0.5;
      temp1 *= temp1;
      float temp2 = D * 0.5;
      s += 2. * temp2;
      temp2 *= temp1;
      s += temp2 / 3.;
      temp2 *= temp1;
      s += (3. / 20.) * temp2;
      temp2 *= temp1;
      s += (5. / 56.) * temp2;
    }

    float x2beta =
      s * track.dzdl / (sqrt(1. - track.dzdl * track.dzdl)) + track.z0;
    float diff = temp_track.hits[i].get_z() - x2beta;

    float chi2_z =
      diff * diff / (tez * tez);

    chi2_xy += chi2_z;

    if (chi2_xy > 5.0) {
      continue;
    }

    if ((best_chi2[temp_track.hits[i].get_layer()] < -0.) ||
        (best_chi2[temp_track.hits[i].get_layer()] > chi2_xy)) {
      best_chi2[temp_track.hits[i].get_layer()] = chi2_xy;
      best_ind[temp_track.hits[i].get_layer()] = i;
    }
  }

  track.hits.clear();
  for (int i = 0; i < layer_indexes.size(); i += 1) {
    if (best_ind[i] < 0) {
      continue;
    }
    track.hits.push_back(temp_track.hits[best_ind[i]]);
  }
  // chi2 = sPHENIXTrackerTPC::fitTrack(track, chi2_hit, tempscale);
  // if(remove_bad_hits(track,6.0)==false)
  // {
  //   return false;
  // }
  if (track.hits.size() < 8) {
    return false;
  }
  return true;
}

static bool fit_all(vector<SimpleHit3D>& hits,
                    vector<vector<int> >& layer_indexes, SimpleTrack3D& track,
                    float& chi2) {

  float scale1 = 32.;
  float scale2 = sqrt(sqrt(0.5));

  SimpleTrack3D temp_track;
  vector<float> chi2_hit;
  for (int i = 0; i < hits.size(); i += 1) {
    if (hits[i].get_layer() < 3) {
      continue;
    }
    temp_track.hits.push_back(hits[i]);
    // temp_track.hits.back().get_ex() *= 1./sqrt(12.);
    // temp_track.hits.back().get_ey() *= 1./sqrt(12.);
    // temp_track.hits.back().get_ez() *= 1./sqrt(12.);

    // temp_track.hits.back().set_ex( temp_track.hits.back().get_ex() * scale1 );
    // temp_track.hits.back().set_ey( temp_track.hits.back().get_ey() * scale1 );
    // temp_track.hits.back().set_ez( temp_track.hits.back().get_ez() * scale1 );
  }

  vector<int> best_ind;
  vector<float> best_chi2;
  track.hits.clear();
  for (int i = 0; i < hits.size(); i += 1) {
    if (hits[i].get_layer() < 3) {
      continue;
    }
    if (layer_indexes[hits[i].get_layer()].size() != 2) {
      continue;
    }
    track.hits.push_back(hits[i]);
    // track.hits.back().get_ex() *= 1./sqrt(12.);
    // track.hits.back().get_ey() *= 1./sqrt(12.);
    // track.hits.back().get_ez() *= 1./sqrt(12.);

    // track.hits.back().set_ex( track.hits.back().get_ex() * scale1);
    // track.hits.back().set_ey( track.hits.back().get_ey() * scale1);
    // track.hits.back().set_ez( track.hits.back().get_ez() * scale1);
  }
  sPHENIXTrackerTPC::fitTrack(track, chi2_hit, scale1);

  float tempscale = scale1;
  for (int i = 0; i < 20; ++i) {
    if (fit_all_update(layer_indexes, temp_track, best_ind, best_chi2, track,
                       chi2, i, tempscale, tempscale) == false) {
      return false;
    }
    tempscale *= scale2;
    // for (unsigned int h = 0; h < temp_track.hits.size(); ++h) {
    //   temp_track.hits[h].set_ex( temp_track.hits[h].get_ex() * scale2);
    //   temp_track.hits[h].set_ey( temp_track.hits[h].get_ey() * scale2);
    //   temp_track.hits[h].set_ez( temp_track.hits[h].get_ez() * scale2);
    // }
    sPHENIXTrackerTPC::fitTrack(track, chi2_hit, tempscale);
  }

  if (fit_all_update(layer_indexes, temp_track, best_ind, best_chi2, track,
                     chi2, 1.0, 1.0) == false) {
    return false;
  }

  // for (unsigned int h = 0; h < temp_track.hits.size(); ++h) {
  //   temp_track.hits[h].set_ex( temp_track.hits[h].get_ex() * 1./sqrt(12.) );
  //   temp_track.hits[h].set_ex( temp_track.hits[h].get_ey() * 1./sqrt(12.) );
  //   temp_track.hits[h].set_ex( temp_track.hits[h].get_ez() * 1./sqrt(12.) );
  // }

  chi2 = sPHENIXTrackerTPC::fitTrack(track, chi2_hit, 1.0);//scale1); // maybe this should be one
  if ((chi2 < 10.0) && (track.hits.size() > ((layer_indexes.size() * 1) / 2))) {
    return true;
  }
  return false;
}

void sPHENIXTrackerTPC::findTracksByCombinatorialKalman(
    vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks,
    const HelixRange& range) {
  // cout<<"findTracksByCombinatorialKalman "<<hits.size()<<"
  // "<<tracks.size()<<endl;

  timeval t1, t2;
  double time1 = 0.;
  double time2 = 0.;

  gettimeofday(&t1, NULL);
  float CHI2_CUT = chi2_cut + 2.;

  vector<vector<int> > layer_indexes;
  layer_indexes.assign(n_layers, vector<int>());
  for (unsigned int i = 0; i < hits.size(); ++i) {
    layer_indexes[hits[i].get_layer()].push_back(i);
  }
  for (int i = 0; i < (int)n_layers; ++i) {
    layer_indexes[i].push_back(-(i + 1));
  }

  vector<TempComb> cur_comb;

  SimpleTrack3D fit_all_track;
  float fit_all_chi2 = -1.;
  bool fit_all_success =
    fit_all(hits, layer_indexes, fit_all_track, fit_all_chi2);
  if (fit_all_success == true) {
    HelixKalmanState state;
    state.phi = fit_all_track.phi;
    if (state.phi < 0.) {
      state.phi += 2. * M_PI;
    }
    state.d = fit_all_track.d;
    state.kappa = fit_all_track.kappa;
    state.nu = sqrt(state.kappa);
    state.z0 = fit_all_track.z0;
    state.dzdl = fit_all_track.dzdl;
    state.C = Matrix<float, 5, 5>::Zero(5, 5);
    state.C(0, 0) = pow(0.01, 2.);
    state.C(1, 1) = pow(0.01, 2.);
    state.C(2, 2) = pow(0.01 * state.nu, 2.);
    state.C(3, 3) = pow(0.01, 2.);
    state.C(4, 4) = pow(0.01, 2.);
    state.chi2 = 0.;
    state.position = n_layers;
    state.x_int = 0.;
    state.y_int = 0.;
    state.z_int = 0.;

    state.C *= 3.;

    SimpleTrack3D temp_track;
    vector<SimpleHit3D> temp_hits;

    for (int h = ((int)(fit_all_track.hits.size() - 1)); h >= 0; h -= 1) {
      HelixKalmanState temp_state(state);
      kalman->addHit(fit_all_track.hits[h], temp_state);
      if ((temp_state.chi2 - state.chi2) < 4.) {
        state = temp_state;
        temp_hits.push_back(fit_all_track.hits[h]);
      }
    }

    state.chi2 = 0.;
    state.position = n_layers;
    state.x_int = 0.;
    state.y_int = 0.;
    state.z_int = 0.;

    state.C *= 3.;
    temp_hits.clear();
    for (int h = ((int)(fit_all_track.hits.size() - 1)); h >= 0; h -= 1) {
      HelixKalmanState temp_state(state);
      kalman->addHit(fit_all_track.hits[h], temp_state);
      if ((temp_state.chi2 - state.chi2) < 10.) {
        state = temp_state;
        temp_hits.push_back(fit_all_track.hits[h]);
      }
    }
    for (int h = ((int)(temp_hits.size() - 1)); h >= 0; h -= 1) {
      temp_track.hits.push_back(temp_hits[h]);
    }

    state.C *= 10.;

    vector<vector<int> > inner_combos;
    for (int l = 0; l < 2; ++l) {
      for (int i0 = 0; i0 < layer_indexes[0].size(); ++i0) {
        for (int i1 = 0; i1 < layer_indexes[1].size(); ++i1) {
          for (int i2 = 0; i2 < layer_indexes[2].size(); ++i2) {
            inner_combos.push_back(vector<int>(3, 0));
            inner_combos.back()[0] = layer_indexes[0][i0];
            inner_combos.back()[1] = layer_indexes[1][i1];
            inner_combos.back()[2] = layer_indexes[2][i2];
          }
        }
      }
    }

    float best_chi2_p[4] = {-1., -1., -1., -1.};
    int best_ind_p[4] = {-1, -1, -1, -1};

    for (int c = 0; c < inner_combos.size(); ++c) {
      HelixKalmanState temp_state(state);
      int n_p = 0;
      for (int l = 2; l >= 0; l -= 1) {
        if (inner_combos[c][l] >= 0) {
          n_p += 1;
          kalman->addHit(hits[inner_combos[c][l]], temp_state);
        }
      }
      float chi2 = temp_state.chi2 - state.chi2;
      if ((best_chi2_p[n_p] == -1) || (chi2 < best_chi2_p[n_p])) {
        best_chi2_p[n_p] = chi2;
        best_ind_p[n_p] = c;
      }
    }

    int best_np = 0;
    for (int n = 3; n >= 0; n -= 1) {
      if ((best_chi2_p[n] > 0.) && (best_chi2_p[n] < 25.)) {
        best_np = n;
        break;
      }
    }

    vector<SimpleHit3D> inner_hits;
    if (best_np > 0) {
      for (int l = 2; l >= 0; l -= 1) {
        int c = best_ind_p[best_np];
        if (inner_combos[c][l] >= 0) {
          inner_hits.push_back(hits[inner_combos[c][l]]);
          kalman->addHit(hits[inner_combos[c][l]], state);
        }
      }
    }

    // vector<SimpleHit3D> inner_hits;
    // int npixels = 0;
    // for( int l=2;l>=0;l-=1 )
    // {
    //   float best_chi2 = -1.;
    //   int best_index = -1;
    //   for(unsigned int i=0;i<(layer_indexes[l].size()-1);++i)
    //   {
    //     if( layer_indexes[l][i]<0 ){continue;}
    //     HelixKalmanState temp_state(state);

    //     SimpleHit3D temp_hit = hits[layer_indexes[l][i]];
    //     temp_hit.get_ex() /= sqrt(12.);
    //     temp_hit.get_ey() /= sqrt(12.);
    //     temp_hit.get_ez() /= sqrt(12.);

    //     // kalman->addHit( temp_hit , temp_state );
    //     // best_chi2 = temp_state.chi2 - state.chi2;

    //     temp_track.hits.push_back( temp_hit );
    //     vector<float> chi2_hit;
    //     float chi2 = sPHENIXTrackerTPC::fitTrack(temp_track, chi2_hit);
    //     best_chi2 = chi2;
    //     temp_track.hits.pop_back();

    //     best_index = layer_indexes[l][i];
    //   }
    //   if( best_index== -1 ){continue;}
    //   if( best_chi2 < 4. )
    //   {
    //     npixels += 1;
    //     SimpleHit3D temp_hit = hits[best_index];
    //     temp_hit.get_ex() /= sqrt(12.);
    //     temp_hit.get_ey() /= sqrt(12.);
    //     temp_hit.get_ez() /= sqrt(12.);
    //     inner_hits.push_back(temp_hit);
    //     // kalman->addHit( temp_hit , state );
    //   }
    // }

    vector<SimpleHit3D> new_hits;
    for (int i = (((int)(inner_hits.size())) - 1); i >= 0; i -= 1) {
      new_hits.push_back(inner_hits[i]);
    }
    for (unsigned int i = 0; i < temp_track.hits.size(); ++i) {
      new_hits.push_back(temp_track.hits[i]);
    }
    temp_track.hits = new_hits;

    // float chi2 = sPHENIXTrackerTPC::fitTrack(temp_track);
    // state.chi2 = chi2*( 2.*( temp_track.hits.size() ) - 5. );

    // if( (state.chi2 < chi2_cut*( 2.*( temp_track.hits.size() ) - 5. )) &&
    // (temp_track.hits.size()>(n_layers/2)) && npixels>1 )
    if ((state.chi2 < chi2_cut * (2. * (temp_track.hits.size()) - 5.)) &&
        (temp_track.hits.size() > (n_layers / 2)) && best_np > 0) {
      tracks.push_back(temp_track);
      track_states.push_back(state);
      if (remove_hits == true) {
        for (unsigned int i = 0; i < tracks.back().hits.size(); ++i) {
          (*hit_used)[tracks.back().hits[i].get_id()] = true;
        }
      }
      gettimeofday(&t2, NULL);
      time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
      time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
      CAtime += (time2 - time1);
      return;
    } else {
      // cout<<best_np<<" "<<state.chi2/( 2.*( temp_track.hits.size() ) - 5.
      // )<<" "<<temp_track.hits.size()<<endl;
      // cout<<"no track 0"<<endl;
      // if(remove_hits == true)
      // {
      //   for(unsigned int i=0;i<hits.size();++i)
      //   {
      //     (*hit_used)[hits[i].get_id()] = true;
      //   }
      // }

      gettimeofday(&t2, NULL);
      time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
      time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
      CAtime += (time2 - time1);
      return;
    }
  } else {
    // cout<<"no track 1"<<endl;
    // if(remove_hits == true)
    // {
    //   for(unsigned int i=0;i<hits.size();++i)
    //   {
    //     (*hit_used)[hits[i].get_id()] = true;
    //   }
    // }

    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
    CAtime += (time2 - time1);
    return;
  }
}

void sPHENIXTrackerTPC::findTracksBySegments(vector<SimpleHit3D>& hits,
                                             vector<SimpleTrack3D>& tracks,
                                             const HelixRange& range) {
  if (n_layers > 20) {
    findTracksByCombinatorialKalman(hits, tracks, range);
    return;
  }

  vector<TrackSegment>* cur_seg = &segments1;
  vector<TrackSegment>* next_seg = &segments2;
  unsigned int curseg_size = 0;
  unsigned int nextseg_size = 0;

  vector<TrackSegment> complete_segments;

  unsigned int allowed_missing = n_layers - req_layers;

  for (unsigned int l = 0; l < n_layers; ++l) {
    layer_sorted[l].clear();
  }
  for (unsigned int i = 0; i < hits.size(); ++i) {
    unsigned int min = (hits[i].get_layer() - allowed_missing);
    if (allowed_missing > hits[i].get_layer()) {
      min = 0;
    }
    for (unsigned int l = min; l <= hits[i].get_layer(); l += 1) {
      layer_sorted[l].push_back(hits[i]);
      SimpleHit3D& hit = layer_sorted[l].back();
      float err_scale = 1.;
      int layer = hit.get_layer();
      if ((layer >= 0) && (layer < (int)(hit_error_scale.size()))) {
        err_scale *= 3. * hit_error_scale[layer];
      }
      hit.set_ex( hit.get_ex() * err_scale);
      hit.set_ey( hit.get_ey() * err_scale);
      hit.set_ez( hit.get_ez() * err_scale);
    }
  }
  for (unsigned int l = 0; l < n_layers; ++l) {
    if (layer_sorted[l].size() == 0) {
      return;
    }
  }

  timeval t1, t2;
  double time1 = 0.;
  double time2 = 0.;

  gettimeofday(&t1, NULL);

  float cosang_diff = 1. - cosang_cut;
  float cosang_diff_inv = 1. / cosang_diff;
  float sinang_cut = sqrt(1. - cosang_cut * cosang_cut);
  float easy_chi2_cut = ca_chi2_cut;

  vector<float> inv_layer;
  inv_layer.assign(n_layers, 1.);
  for (unsigned int l = 3; l < n_layers; ++l) {
    inv_layer[l] = 1. / (((float)l) - 2.);
  }

  unsigned int hit_counter = 0;
  float x1_a[4] __attribute__((aligned(16)));
  float x2_a[4] __attribute__((aligned(16)));
  float x3_a[4] __attribute__((aligned(16)));
  float y1_a[4] __attribute__((aligned(16)));
  float y2_a[4] __attribute__((aligned(16)));
  float y3_a[4] __attribute__((aligned(16)));
  float z1_a[4] __attribute__((aligned(16)));
  float z2_a[4] __attribute__((aligned(16)));
  float z3_a[4] __attribute__((aligned(16)));

  float dx1_a[4] __attribute__((aligned(16)));
  float dx2_a[4] __attribute__((aligned(16)));
  float dx3_a[4] __attribute__((aligned(16)));
  float dy1_a[4] __attribute__((aligned(16)));
  float dy2_a[4] __attribute__((aligned(16)));
  float dy3_a[4] __attribute__((aligned(16)));
  float dz1_a[4] __attribute__((aligned(16)));
  float dz2_a[4] __attribute__((aligned(16)));
  float dz3_a[4] __attribute__((aligned(16)));

  float kappa_a[4] __attribute__((aligned(16)));
  float dkappa_a[4] __attribute__((aligned(16)));

  float ux_mid_a[4] __attribute__((aligned(16)));
  float uy_mid_a[4] __attribute__((aligned(16)));
  float ux_end_a[4] __attribute__((aligned(16)));
  float uy_end_a[4] __attribute__((aligned(16)));

  float dzdl_1_a[4] __attribute__((aligned(16)));
  float dzdl_2_a[4] __attribute__((aligned(16)));
  float ddzdl_1_a[4] __attribute__((aligned(16)));
  float ddzdl_2_a[4] __attribute__((aligned(16)));

  float cur_kappa_a[4] __attribute__((aligned(16)));
  float cur_dkappa_a[4] __attribute__((aligned(16)));
  float cur_ux_a[4] __attribute__((aligned(16)));
  float cur_uy_a[4] __attribute__((aligned(16)));
  float cur_chi2_a[4] __attribute__((aligned(16)));
  float chi2_a[4] __attribute__((aligned(16)));

  unsigned int hit1[4];
  unsigned int hit2[4];
  unsigned int hit3[4];

  TrackSegment temp_segment;
  temp_segment.hits.assign(n_layers, 0);
  // make segments out of first 3 layers
  for (unsigned int i = 0, sizei = layer_sorted[0].size(); i < sizei; ++i) {
    for (unsigned int j = 0, sizej = layer_sorted[1].size(); j < sizej; ++j) {
      for (unsigned int k = 0, sizek = layer_sorted[2].size(); k < sizek; ++k) {
        if ((layer_sorted[0][i].get_layer() >= layer_sorted[1][j].get_layer()) ||
            (layer_sorted[1][j].get_layer() >= layer_sorted[2][k].get_layer())) {
          continue;
        }

        x1_a[hit_counter] = layer_sorted[0][i].get_x();
        y1_a[hit_counter] = layer_sorted[0][i].get_y();
        z1_a[hit_counter] = layer_sorted[0][i].get_z();
        dx1_a[hit_counter] = layer_sorted[0][i].get_ex();
        dy1_a[hit_counter] = layer_sorted[0][i].get_ey();
        dz1_a[hit_counter] = layer_sorted[0][i].get_ez();

        x2_a[hit_counter] = layer_sorted[1][j].get_x();
        y2_a[hit_counter] = layer_sorted[1][j].get_y();
        z2_a[hit_counter] = layer_sorted[1][j].get_z();
        dx2_a[hit_counter] = layer_sorted[1][j].get_ex();
        dy2_a[hit_counter] = layer_sorted[1][j].get_ey();
        dz2_a[hit_counter] = layer_sorted[1][j].get_ez();

        x3_a[hit_counter] = layer_sorted[2][k].get_x();
        y3_a[hit_counter] = layer_sorted[2][k].get_y();
        z3_a[hit_counter] = layer_sorted[2][k].get_z();
        dx3_a[hit_counter] = layer_sorted[2][k].get_ex();
        dy3_a[hit_counter] = layer_sorted[2][k].get_ey();
        dz3_a[hit_counter] = layer_sorted[2][k].get_ez();

        hit1[hit_counter] = i;
        hit2[hit_counter] = j;
        hit3[hit_counter] = k;

        hit_counter += 1;

        if (hit_counter == 4) {
          calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a,
                                 z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a,
                                 dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a,
                                 ux_mid_a, uy_mid_a, ux_end_a, uy_end_a,
                                 dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);

          for (unsigned int h = 0; h < hit_counter; ++h) {
            temp_segment.chi2 =
	      (dzdl_1_a[h] - dzdl_2_a[h]) /
	      (ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h] * sinang_cut));
            temp_segment.chi2 *= temp_segment.chi2;
            if (temp_segment.chi2 > 2.0) {
              continue;
            }
            temp_segment.ux = ux_end_a[h];
            temp_segment.uy = uy_end_a[h];
            temp_segment.kappa = kappa_a[h];
            if (temp_segment.kappa > top_range.max_k) {
              continue;
            }
            temp_segment.dkappa = dkappa_a[h];
            temp_segment.hits[0] = hit1[h];
            temp_segment.hits[1] = hit2[h];
            temp_segment.hits[2] = hit3[h];
            temp_segment.n_hits = 3;
            unsigned int outer_layer =
	      layer_sorted[2][temp_segment.hits[2]].get_layer();
            if ((outer_layer - 2) > allowed_missing) {
              continue;
            }
            if ((n_layers - 3) <= allowed_missing) {
              complete_segments.push_back(temp_segment);
            }
            if (next_seg->size() == nextseg_size) {
              next_seg->push_back(temp_segment);
              nextseg_size += 1;
            } else {
              (*next_seg)[nextseg_size] = temp_segment;
              nextseg_size += 1;
            }
          }

          hit_counter = 0;
        }
      }
    }
  }
  if (hit_counter != 0) {
    calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a,
                           dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a,
                           dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a,
                           ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a,
                           ddzdl_2_a);

    for (unsigned int h = 0; h < hit_counter; ++h) {
      temp_segment.chi2 =
	(dzdl_1_a[h] - dzdl_2_a[h]) /
	(ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h] * sinang_cut));
      temp_segment.chi2 *= temp_segment.chi2;
      if (temp_segment.chi2 > 2.0) {
        continue;
      }
      temp_segment.ux = ux_end_a[h];
      temp_segment.uy = uy_end_a[h];
      temp_segment.kappa = kappa_a[h];
      if (temp_segment.kappa > top_range.max_k) {
        continue;
      }
      temp_segment.dkappa = dkappa_a[h];
      temp_segment.hits[0] = hit1[h];
      temp_segment.hits[1] = hit2[h];
      temp_segment.hits[2] = hit3[h];
      temp_segment.n_hits = 3;
      unsigned int outer_layer = layer_sorted[2][temp_segment.hits[2]].get_layer();
      if ((outer_layer - 2) > allowed_missing) {
        continue;
      }
      if ((n_layers - 3) <= allowed_missing) {
        complete_segments.push_back(temp_segment);
      }
      if (next_seg->size() == nextseg_size) {
        next_seg->push_back(temp_segment);
        nextseg_size += 1;
      } else {
        (*next_seg)[nextseg_size] = temp_segment;
        nextseg_size += 1;
      }
    }

    hit_counter = 0;
  }

  swap(cur_seg, next_seg);
  swap(curseg_size, nextseg_size);

  // add hits to segments layer-by-layer, cutting out bad segments
  unsigned int whichseg[4];
  for (unsigned int l = 3; l < n_layers; ++l) {
    if (l == (n_layers - 1)) {
      easy_chi2_cut *= 0.25;
    }
    nextseg_size = 0;
    for (unsigned int i = 0, sizei = curseg_size; i < sizei; ++i) {
      for (unsigned int j = 0, sizej = layer_sorted[l].size(); j < sizej; ++j) {
        if ((layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_layer() >=
             layer_sorted[l][j].get_layer())) {
          continue;
        }

        x1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_x();
        y1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_y();
        z1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_z();
        x2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_x();
        y2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_y();
        z2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_z();
        x3_a[hit_counter] = layer_sorted[l][j].get_x();
        y3_a[hit_counter] = layer_sorted[l][j].get_y();
        z3_a[hit_counter] = layer_sorted[l][j].get_z();

        dx1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_ex();
        dy1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_ey();
        dz1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_ez();
        dx2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_ex();
        dy2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_ey();
        dz2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_ez();
        dx3_a[hit_counter] = layer_sorted[l][j].get_ex();
        dy3_a[hit_counter] = layer_sorted[l][j].get_ey();
        dz3_a[hit_counter] = layer_sorted[l][j].get_ez();

        cur_kappa_a[hit_counter] = (*cur_seg)[i].kappa;
        cur_dkappa_a[hit_counter] = (*cur_seg)[i].dkappa;
        cur_ux_a[hit_counter] = (*cur_seg)[i].ux;
        cur_uy_a[hit_counter] = (*cur_seg)[i].uy;
        cur_chi2_a[hit_counter] = (*cur_seg)[i].chi2;

        whichseg[hit_counter] = i;
        hit1[hit_counter] = j;

        hit_counter += 1;
        if (hit_counter == 4) {
          calculateKappaTangents(
				 x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a,
				 dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a,
				 dkappa_a, ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a,
				 dzdl_2_a, ddzdl_1_a, ddzdl_2_a, sinang_cut, cosang_diff_inv,
				 cur_kappa_a, cur_dkappa_a, cur_ux_a, cur_uy_a, cur_chi2_a,
				 chi2_a);

          for (unsigned int h = 0; h < hit_counter; ++h) {
            if ((chi2_a[h]) * inv_layer[l] < easy_chi2_cut) {
              temp_segment.chi2 = chi2_a[h];
              temp_segment.ux = ux_end_a[h];
              temp_segment.uy = uy_end_a[h];
              temp_segment.kappa = kappa_a[h];
              if (temp_segment.kappa > top_range.max_k) {
                continue;
              }
              temp_segment.dkappa = dkappa_a[h];
              for (unsigned int ll = 0; ll < l; ++ll) {
                temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];
              }
              temp_segment.hits[l] = hit1[h];
              unsigned int outer_layer =
		layer_sorted[l][temp_segment.hits[l]].get_layer();
              temp_segment.n_hits = l + 1;
              if ((n_layers - (l + 1)) <= allowed_missing) {
                complete_segments.push_back(temp_segment);
              }
              if ((outer_layer - l) > allowed_missing) {
                continue;
              }
              if (next_seg->size() == nextseg_size) {
                next_seg->push_back(temp_segment);
                nextseg_size += 1;
              } else {
                (*next_seg)[nextseg_size] = temp_segment;
                nextseg_size += 1;
              }
            }
          }
          hit_counter = 0;
        }
      }
    }
    if (hit_counter != 0) {
      calculateKappaTangents(
			     x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a, dx1_a, dy1_a,
			     dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a,
			     ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a,
			     ddzdl_2_a, sinang_cut, cosang_diff_inv, cur_kappa_a, cur_dkappa_a,
			     cur_ux_a, cur_uy_a, cur_chi2_a, chi2_a);

      for (unsigned int h = 0; h < hit_counter; ++h) {
        if ((chi2_a[h]) * inv_layer[l] < easy_chi2_cut) {
          temp_segment.chi2 = chi2_a[h];
          temp_segment.ux = ux_end_a[h];
          temp_segment.uy = uy_end_a[h];
          temp_segment.kappa = kappa_a[h];
          if (temp_segment.kappa > top_range.max_k) {
            continue;
          }

          temp_segment.dkappa = dkappa_a[h];
          for (unsigned int ll = 0; ll < l; ++ll) {
            temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];
          }
          temp_segment.hits[l] = hit1[h];
          unsigned int outer_layer =
	    layer_sorted[l][temp_segment.hits[l]].get_layer();
          temp_segment.n_hits = l + 1;
          if ((n_layers - (l + 1)) <= allowed_missing) {
            complete_segments.push_back(temp_segment);
          }
          if ((outer_layer - l) > allowed_missing) {
            continue;
          }
          if (next_seg->size() == nextseg_size) {
            next_seg->push_back(temp_segment);
            nextseg_size += 1;
          } else {
            (*next_seg)[nextseg_size] = temp_segment;
            nextseg_size += 1;
          }
        }
      }
      hit_counter = 0;
    }
    swap(cur_seg, next_seg);
    swap(curseg_size, nextseg_size);
  }

  for (unsigned int i = 0; i < complete_segments.size(); ++i) {
    if (cur_seg->size() == curseg_size) {
      cur_seg->push_back(complete_segments[i]);
      curseg_size += 1;
    } else {
      (*cur_seg)[curseg_size] = complete_segments[i];
      curseg_size += 1;
    }
  }

  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
  CAtime += (time2 - time1);
  SimpleTrack3D temp_track;
  temp_track.hits.assign(n_layers, SimpleHit3D());
  vector<SimpleHit3D> temp_hits;
  for (unsigned int i = 0, sizei = curseg_size; i < sizei; ++i) {
    temp_track.hits.assign((*cur_seg)[i].n_hits, SimpleHit3D());

    temp_comb.assign((*cur_seg)[i].n_hits, 0);
    for (unsigned int l = 0; l < (*cur_seg)[i].n_hits; ++l) {
      temp_comb[l] = layer_sorted[l][(*cur_seg)[i].hits[l]].get_id();
    }
    sort(temp_comb.begin(), temp_comb.end());
    set<vector<unsigned int> >::iterator it = combos.find(temp_comb);
    if (it != combos.end()) {
      continue;
    }
    if (combos.size() > 10000) {
      combos.clear();
    }
    combos.insert(temp_comb);

    for (unsigned int l = 0; l < (*cur_seg)[i].n_hits; ++l) {
      temp_track.hits[l] = layer_sorted[l][(*cur_seg)[i].hits[l]];
    }

    gettimeofday(&t1, NULL);

    // unsigned int layer_out = temp_track.hits.size()-1;
    // unsigned int layer_mid = layer_out/2;
    // temp_track_3hits.hits[0] = temp_track.hits[0];
    // temp_track_3hits.hits[1] = temp_track.hits[layer_mid];
    // temp_track_3hits.hits[2] = temp_track.hits[layer_out];

    float init_chi2 = fitTrack(temp_track);

    if (init_chi2 > fast_chi2_cut_max) {
      if (init_chi2 > fast_chi2_cut_par0 +
	  fast_chi2_cut_par1 / kappaToPt(temp_track.kappa)) {
        gettimeofday(&t2, NULL);
        time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
        time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
        KALtime += (time2 - time1);
        continue;
      }
    }
    HelixKalmanState state;
    state.phi = temp_track.phi;
    if (state.phi < 0.) {
      state.phi += 2. * M_PI;
    }
    state.d = temp_track.d;
    state.kappa = temp_track.kappa;
    state.nu = sqrt(state.kappa);
    state.z0 = temp_track.z0;
    state.dzdl = temp_track.dzdl;
    state.C = Matrix<float, 5, 5>::Zero(5, 5);
    state.C(0, 0) = pow(0.01, 2.);
    state.C(1, 1) = pow(0.01, 2.);
    state.C(2, 2) = pow(0.01 * state.nu, 2.);
    state.C(3, 3) = pow(0.05, 2.);
    state.C(4, 4) = pow(0.05, 2.);
    state.chi2 = 0.;
    state.position = 0;
    state.x_int = 0.;
    state.y_int = 0.;
    state.z_int = 0.;

    for (unsigned int h = 0; h < temp_track.hits.size(); ++h) {
      kalman->addHit(temp_track.hits[h], state);
      nfits += 1;
    }

    // fudge factor for non-gaussian hit sizes
    state.C *= 3.;
    state.chi2 *= 6.;

    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
    KALtime += (time2 - time1);

    if (!(temp_track.kappa == temp_track.kappa)) {
      continue;
    }
    if (temp_track.kappa > top_range.max_k) {
      continue;
    }
    if (!(state.chi2 == state.chi2)) {
      continue;
    }
    if (state.chi2 / (2. * ((float)(temp_track.hits.size())) - 5.) > chi2_cut) {
      continue;
    }

    if (cut_on_dca == true) {
      if (fabs(temp_track.d) > dca_cut) {
        continue;
      }
      if (fabs(temp_track.z0) > dca_cut) {
        continue;
      }
    }

    tracks.push_back(temp_track);
    track_states.push_back(state);
    if ((remove_hits == true) && (state.chi2 < chi2_removal_cut) &&
        (temp_track.hits.size() >= n_removal_hits)) {
      for (unsigned int i = 0; i < temp_track.hits.size(); ++i) {
        (*hit_used)[temp_track.hits[i].get_id()] = true;
      }
    }
  }
}

void sPHENIXTrackerTPC::findSeededTracksbySegments(
						   vector<SimpleTrack3D>& seeds, vector<SimpleHit3D>& hits,
						   vector<SimpleTrack3D>& tracks, const HelixRange& range) {
  unsigned int allowed_missing = n_layers - seed_layer - req_layers;

  for (unsigned int l = 0; l < n_layers; ++l) {
    layer_sorted[l].clear();
  }
  for (unsigned int i = 0; i < hits.size(); ++i) {
    unsigned int min = (hits[i].get_layer() - allowed_missing);
    if (allowed_missing > hits[i].get_layer()) {
      min = 0;
    }
    for (unsigned int l = min; l <= hits[i].get_layer(); l += 1) {
      layer_sorted[l].push_back(hits[i]);
    }
  }

  findSeededTracksbySegments_run(seeds, hits, tracks);
}

void sPHENIXTrackerTPC::findSeededTracksbySegments_run(
						       vector<SimpleTrack3D>& seeds, vector<SimpleHit3D>& hits,
						       vector<SimpleTrack3D>& tracks) {
  if (seeds.size() == 0) {
    return;
  }

  timeval t1, t2;
  double time1 = 0.;
  double time2 = 0.;

  gettimeofday(&t1, NULL);

  unsigned int first_new_layer = seed_layer;

  unsigned int allowed_missing = n_layers - seed_layer - req_layers;

  vector<TrackSegment>* cur_seg = &segments1;
  vector<TrackSegment>* next_seg = &segments2;
  unsigned int curseg_size = 0;
  unsigned int nextseg_size = 0;

  vector<TrackSegment> complete_segments;

  float cosang_diff = 1. - cosang_cut;
  float cosang_diff_inv = 1. / cosang_diff;
  float sinang_cut = sqrt(1. - cosang_cut * cosang_cut);
  float easy_chi2_cut = ca_chi2_cut;

  unsigned int hit_counter = 0;
  float x1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float x2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float x3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float y1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float y2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float y3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float z1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float z2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float z3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};

  float dx1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float dx2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float dx3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float dy1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float dy2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float dy3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float dz1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float dz2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float dz3_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};

  float kappa_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float dkappa_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};

  float ux_mid_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float uy_mid_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float ux_end_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float uy_end_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};

  float dzdl_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float dzdl_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float ddzdl_1_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};
  float ddzdl_2_a[4] __attribute__((aligned(16))) = {0., 0., 0., 0.};

  TrackSegment temp_segment;
  temp_segment.hits.assign(n_layers, 0);
  unsigned int whichhit[4];
  unsigned int whichseed[4];
  for (unsigned int seed = 0, seedsize = seeds.size(); seed < seedsize;
       ++seed) {
    unsigned int firsthit = seeds[seed].hits.size();

    x1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_x();
    y1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_y();
    z1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_z();
    x2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_x();
    y2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_y();
    z2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_z();
    x3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_x();
    y3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_y();
    z3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_z();

    dx1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_ex();
    dy1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_ey();
    dz1_a[hit_counter] = seeds[seed].hits[firsthit - 3].get_ez();
    dx2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_ex();
    dy2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_ey();
    dz2_a[hit_counter] = seeds[seed].hits[firsthit - 2].get_ez();
    dx3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_ex();
    dy3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_ey();
    dz3_a[hit_counter] = seeds[seed].hits[firsthit - 1].get_ez();

    whichseed[hit_counter] = seed;
    hit_counter += 1;
    if (hit_counter == 4) {
      calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a,
                             z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a,
                             dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a,
                             uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a,
                             ddzdl_1_a, ddzdl_2_a);

      for (unsigned int h = 0; h < hit_counter; ++h) {
        temp_segment.chi2 = 0.;
        temp_segment.ux = ux_end_a[h];
        temp_segment.uy = uy_end_a[h];
        temp_segment.kappa = kappa_a[h];
        temp_segment.dkappa = dkappa_a[h];
        temp_segment.seed = whichseed[h];
        if (next_seg->size() == nextseg_size) {
          next_seg->push_back(temp_segment);
          nextseg_size += 1;
        } else {
          (*next_seg)[nextseg_size] = temp_segment;
          nextseg_size += 1;
        }
      }

      hit_counter = 0;
    }
  }
  if (hit_counter != 0) {
    calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a, z3_a,
                           dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a,
                           dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a,
                           ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a,
                           ddzdl_2_a);

    for (unsigned int h = 0; h < hit_counter; ++h) {
      temp_segment.chi2 = 0.;
      temp_segment.ux = ux_end_a[h];
      temp_segment.uy = uy_end_a[h];
      temp_segment.kappa = kappa_a[h];
      temp_segment.dkappa = dkappa_a[h];
      temp_segment.seed = whichseed[h];
      if (next_seg->size() == nextseg_size) {
        next_seg->push_back(temp_segment);
        nextseg_size += 1;
      } else {
        (*next_seg)[nextseg_size] = temp_segment;
        nextseg_size += 1;
      }
    }

    hit_counter = 0;
  }
  swap(cur_seg, next_seg);
  swap(curseg_size, nextseg_size);
  unsigned int whichseg[4];
  for (unsigned int l = first_new_layer; l < n_layers; ++l) {
    nextseg_size = 0;
    for (unsigned int j = 0, sizej = curseg_size; j < sizej; ++j) {
      unsigned int firsthit = seeds[(*cur_seg)[j].seed].hits.size();
      for (unsigned int i = 0, sizei = layer_sorted[l].size(); i < sizei; ++i) {
        if ((l - 2) <= (first_new_layer - 1)) {
          x1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 2 - (first_new_layer - firsthit)]
	    .get_x();
          y1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 2 - (first_new_layer - firsthit)]
	    .get_y();
          z1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 2 - (first_new_layer - firsthit)]
	    .get_z();
          dx1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 2 - (first_new_layer - firsthit)]
	    .get_ex();
          dy1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 2 - (first_new_layer - firsthit)]
	    .get_ey();
          dz1_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 2 - (first_new_layer - firsthit)]
	    .get_ez();
        } else {
          x1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_x();
          y1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_y();
          z1_a[hit_counter] = layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_z();
          dx1_a[hit_counter] =
	    layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_ex();
          dy1_a[hit_counter] =
	    layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_ey();
          dz1_a[hit_counter] =
	    layer_sorted[l - 2][(*cur_seg)[j].hits[l - 2]].get_ez();
        }
        if ((l - 1) <= (first_new_layer - 1)) {
          x2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 1 - (first_new_layer - firsthit)]
	    .get_x();
          y2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 1 - (first_new_layer - firsthit)]
	    .get_y();
          z2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 1 - (first_new_layer - firsthit)]
	    .get_z();
          dx2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 1 - (first_new_layer - firsthit)]
	    .get_ex();
          dy2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 1 - (first_new_layer - firsthit)]
	    .get_ey();
          dz2_a[hit_counter] = seeds[(*cur_seg)[j].seed]
	    .hits[l - 1 - (first_new_layer - firsthit)]
	    .get_ez();
        } else {
          x2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_x();
          y2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_y();
          z2_a[hit_counter] = layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_z();
          dx2_a[hit_counter] =
	    layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_ex();
          dy2_a[hit_counter] =
	    layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_ey();
          dz2_a[hit_counter] =
	    layer_sorted[l - 1][(*cur_seg)[j].hits[l - 1]].get_ez();
        }
        x3_a[hit_counter] = layer_sorted[l][i].get_x();
        y3_a[hit_counter] = layer_sorted[l][i].get_y();
        z3_a[hit_counter] = layer_sorted[l][i].get_z();
        dx3_a[hit_counter] = layer_sorted[l][i].get_ex();
        dy3_a[hit_counter] = layer_sorted[l][i].get_ey();
        dz3_a[hit_counter] = layer_sorted[l][i].get_ez();

        whichhit[hit_counter] = i;
        whichseg[hit_counter] = j;
        hit_counter += 1;
        if (hit_counter == 4) {
          calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a,
                                 z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a,
                                 dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a,
                                 ux_mid_a, uy_mid_a, ux_end_a, uy_end_a,
                                 dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);

          for (unsigned int h = 0; h < hit_counter; ++h) {
            float kdiff = (*cur_seg)[whichseg[h]].kappa - kappa_a[h];
            float dk = (*cur_seg)[whichseg[h]].dkappa + dkappa_a[h];
            dk += sinang_cut * kappa_a[h];
            float chi2_k = kdiff * kdiff / (dk * dk);
            float cos_scatter = (*cur_seg)[whichseg[h]].ux * ux_mid_a[h] +
	      (*cur_seg)[whichseg[h]].uy * uy_mid_a[h];
            float chi2_ang = (1. - cos_scatter) * (1. - cos_scatter) *
	      cosang_diff_inv * cosang_diff_inv;
            float chi2_dzdl =
	      (dzdl_1_a[h] - dzdl_2_a[h]) /
	      (ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h] * sinang_cut));
            chi2_dzdl *= chi2_dzdl;
            if (((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl) /
		((float)l - 2.) <
                easy_chi2_cut) {
              temp_segment.chi2 =
		(*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl;
              temp_segment.ux = ux_end_a[h];
              temp_segment.uy = uy_end_a[h];
              temp_segment.kappa = kappa_a[h];
              temp_segment.dkappa = dkappa_a[h];
              for (unsigned int ll = 0; ll < l; ++ll) {
                temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];
              }
              temp_segment.hits[l] = whichhit[h];
              temp_segment.seed = (*cur_seg)[whichseg[h]].seed;
              unsigned int outer_layer =
		layer_sorted[l][temp_segment.hits[l]].get_layer();
              temp_segment.n_hits = l + 1;
              if ((n_layers - (l + 1)) <= allowed_missing) {
                complete_segments.push_back(temp_segment);
              }
              if ((outer_layer - l) > allowed_missing) {
                continue;
              }
              if (next_seg->size() == nextseg_size) {
                next_seg->push_back(temp_segment);
                nextseg_size += 1;
              } else {
                (*next_seg)[nextseg_size] = temp_segment;
                nextseg_size += 1;
              }
            }
          }
          hit_counter = 0;
        }
      }
    }
    if (hit_counter != 0) {
      calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a,
                             z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a,
                             dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a,
                             uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a,
                             ddzdl_1_a, ddzdl_2_a);

      for (unsigned int h = 0; h < hit_counter; ++h) {
        float kdiff = (*cur_seg)[whichseg[h]].kappa - kappa_a[h];
        float dk = (*cur_seg)[whichseg[h]].dkappa + dkappa_a[h];
        dk += sinang_cut * kappa_a[h];
        float chi2_k = kdiff * kdiff / (dk * dk);
        float cos_scatter = (*cur_seg)[whichseg[h]].ux * ux_mid_a[h] +
	  (*cur_seg)[whichseg[h]].uy * uy_mid_a[h];
        float chi2_ang = (1. - cos_scatter) * (1. - cos_scatter) *
	  cosang_diff_inv * cosang_diff_inv;
        float chi2_dzdl =
	  (dzdl_1_a[h] - dzdl_2_a[h]) /
	  (ddzdl_1_a[h] + ddzdl_2_a[h] + fabs(dzdl_1_a[h] * sinang_cut));
        chi2_dzdl *= chi2_dzdl;
        if (((*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl) /
	    ((float)l - 2.) <
            easy_chi2_cut) {
          temp_segment.chi2 =
	    (*cur_seg)[whichseg[h]].chi2 + chi2_ang + chi2_k + chi2_dzdl;
          temp_segment.ux = ux_end_a[h];
          temp_segment.uy = uy_end_a[h];
          temp_segment.kappa = kappa_a[h];
          temp_segment.dkappa = dkappa_a[h];
          for (unsigned int ll = 0; ll < l; ++ll) {
            temp_segment.hits[ll] = (*cur_seg)[whichseg[h]].hits[ll];
          }
          temp_segment.hits[l] = whichhit[h];
          temp_segment.seed = (*cur_seg)[whichseg[h]].seed;
          unsigned int outer_layer =
	    layer_sorted[l][temp_segment.hits[l]].get_layer();
          temp_segment.n_hits = l + 1;
          if ((n_layers - (l + 1)) <= allowed_missing) {
            complete_segments.push_back(temp_segment);
          }
          if ((outer_layer - l) > allowed_missing) {
            continue;
          }
          if (next_seg->size() == nextseg_size) {
            next_seg->push_back(temp_segment);
            nextseg_size += 1;
          } else {
            (*next_seg)[nextseg_size] = temp_segment;
            nextseg_size += 1;
          }
        }
      }
      hit_counter = 0;
    }

    swap(cur_seg, next_seg);
    swap(curseg_size, nextseg_size);
  }

  for (unsigned int i = 0; i < complete_segments.size(); ++i) {
    if (cur_seg->size() == curseg_size) {
      cur_seg->push_back(complete_segments[i]);
      curseg_size += 1;
    } else {
      (*cur_seg)[curseg_size] = complete_segments[i];
      curseg_size += 1;
    }
  }

  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
  CAtime += (time2 - time1);

  gettimeofday(&t1, NULL);

  SimpleTrack3D temp_track;
  temp_track.hits.assign(n_layers, SimpleHit3D());
  for (unsigned int i = 0, sizei = curseg_size; i < sizei; ++i) {
    temp_track.hits.clear();
    for (unsigned int l = 0; l < seeds[(*cur_seg)[i].seed].hits.size(); ++l) {
      temp_track.hits.push_back(seeds[(*cur_seg)[i].seed].hits[l]);
    }
    for (unsigned int l = seed_layer; l < (*cur_seg)[i].n_hits; ++l) {
      temp_track.hits.push_back(layer_sorted[l][(*cur_seg)[i].hits[l]]);
    }

    HelixKalmanState state;

    fitTrack(temp_track);
    state.C = Matrix<float, 5, 5>::Zero(5, 5);
    state.phi = temp_track.phi;
    if (state.phi < 0.) {
      state.phi += 2. * M_PI;
    }
    state.d = temp_track.d;
    state.kappa = temp_track.kappa;
    state.nu = sqrt(state.kappa);
    state.z0 = temp_track.z0;
    state.dzdl = temp_track.dzdl;
    state.C(0, 0) = pow(0.1 * state.phi, 2.);
    state.C(1, 1) = pow(0.1 * state.d, 2.);
    state.C(2, 2) = pow(0.1 * state.nu, 2.);
    state.C(3, 3) = pow(0.1 * state.z0, 2.);
    state.C(4, 4) = pow(0.1 * state.dzdl, 2.);

    state.chi2 = 0.;
    state.position = 0;
    state.x_int = 0.;
    state.y_int = 0.;
    state.z_int = 0.;
    for (unsigned int l = 0; l < temp_track.hits.size(); ++l) {
      kalman->addHit(temp_track.hits[l], state);
      nfits += 1;
    }

    if (state.chi2 / (2. * ((float)(temp_track.hits.size() + 1)) - 5.) >
        chi2_cut) {
      continue;
    }
    if (state.chi2 != state.chi2) {
      continue;
    }

    temp_track.phi = state.phi;
    if (temp_track.phi < 0.) {
      temp_track.phi += 2. * M_PI;
    }
    if (temp_track.phi > 2. * M_PI) {
      temp_track.phi -= 2. * M_PI;
    }
    temp_track.d = state.d;
    temp_track.kappa = state.kappa;
    temp_track.z0 = state.z0;
    temp_track.dzdl = state.dzdl;

    if ((remove_hits == true) && (state.chi2 < chi2_removal_cut) &&
        (temp_track.hits.size() >= n_removal_hits)) {
      for (unsigned int l = (seeds[(*cur_seg)[i].seed].hits.size());
           l < temp_track.hits.size(); ++l) {
        (*hit_used)[temp_track.hits[l].get_id()] = true;
        temp_track.hits[l].set_id(index_mapping[temp_track.hits[l].get_id()]);
      }
    } else {
      for (unsigned int l = (seeds[(*cur_seg)[i].seed].hits.size());
           l < temp_track.hits.size(); ++l) {
        temp_track.hits[l].set_id(index_mapping[temp_track.hits[l].get_id()]);
      }
    }
    tracks.push_back(temp_track);
    track_states.push_back(state);
    seed_used[seeds[(*cur_seg)[i].seed].index] = true;
  }

  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
  KALtime += (time2 - time1);
}
