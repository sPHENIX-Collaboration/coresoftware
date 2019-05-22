#include "sPHENIXTrackerTpc.h"
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

void sPHENIXTrackerTpc::tripletRejection(vector<SimpleTrack3D>& input,
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

sPHENIXTrackerTpc::sPHENIXTrackerTpc(unsigned int n_phi, unsigned int n_d,
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

sPHENIXTrackerTpc::sPHENIXTrackerTpc(
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
    require_pixels(false){
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

    pins = new Pincushion<sPHENIXTrackerTpc>(this, vssp);

    vector<vector<unsigned int> > zoom_profile_new;
    for (unsigned int i = 1; i < zoom_profile.size(); ++i) {
      zoom_profile_new.push_back(zoom_profile[i]);
    }

    for (unsigned int i = 0; i < nthreads; ++i) {
      thread_trackers.push_back(new sPHENIXTrackerTpc(
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

sPHENIXTrackerTpc::~sPHENIXTrackerTpc() {
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

float sPHENIXTrackerTpc::kappaToPt(float kappa) {
  return detector_B_field / 333.6 / kappa;
}

float sPHENIXTrackerTpc::ptToKappa(float pt) {
  return detector_B_field / 333.6 / pt;
}


void sPHENIXTrackerTpc::finalize(vector<SimpleTrack3D>& input,
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
	kalman->addHit(hit, track_states[i]);
      }

      if (fabs(track_states[i].d) < 0.01) {
	SimpleHit3D vertex_hit;
	vertex_hit.set_x(0.0);
	vertex_hit.set_y(0.0);
	vertex_hit.set_z(0.0);
	//vertex_hit.set_ex(0.0001);
	//vertex_hit.set_ey(0.0001);
	//vertex_hit.set_ez(0.0001);
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

bool sPHENIXTrackerTpc::breakRecursion(const vector<SimpleHit3D>& hits,
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

float sPHENIXTrackerTpc::phiError(SimpleHit3D& hit, float min_k, float max_k,
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

float sPHENIXTrackerTpc::dzdlError(SimpleHit3D& hit, float min_k, float max_k,
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

void sPHENIXTrackerTpc::findTracks(vector<SimpleHit3D>& hits,
                                   vector<SimpleTrack3D>& tracks,
                                   const HelixRange& range) {
  cout << "findTracks " << endl;
  findtracksiter += 1;
  if( n_layers < 10 ) {
    findTracksBySegments(hits, tracks, range);
  } else {
    findTracksByCombinatorialKalman(hits, tracks, range);
  }
}

float sPHENIXTrackerTpc::fitTrack(SimpleTrack3D& track,
				  float scale) {
  vector<float> chi2_hit;
  return sPHENIXTrackerTpc::fitTrack(track, chi2_hit, scale);
}

float sPHENIXTrackerTpc::fitTrack(SimpleTrack3D& track,
                                  vector<float>& chi2_hit,
				  float scale) {
  
  chi2_hit.clear();
  vector<float> xyres;
  vector<float> xyres_inv;
  vector<float> zres;
  vector<float> zres_inv;
  for (unsigned int i = 0; i < track.hits.size(); i++) {

    float ex = (2.0*sqrt(track.hits[i].get_size(0,0))) * scale;
    float ey = (2.0*sqrt(track.hits[i].get_size(1,1))) * scale;
    float ez = (2.0*sqrt(track.hits[i].get_size(2,2))) * scale;

    if (track.hits[i].get_layer() < 0) {
      ex = 0.0001 * scale;
      ey = 0.0001 * scale;
      ez = 0.0001 * scale;
    }
    
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

void sPHENIXTrackerTpc::initSplitting(vector<SimpleHit3D>& hits,
                                      unsigned int min_hits,
                                      unsigned int max_hits) {
  initEvent(hits, min_hits);
  (*(hits_vec[0])) = hits;
  zoomranges.clear();
  for (unsigned int z = 0; z <= max_zoom; z++) {
    zoomranges.push_back(top_range);
  }
}

static bool remove_bad_hits(SimpleTrack3D& track, float cut, float scale = 1.0) {
  SimpleTrack3D temp_track = track;
  float fit_chi2 = 0.;
  vector<float> chi2_hit;
  vector<float> temp_hits;
  while (true) {
    temp_track = track;
    fit_chi2 = sPHENIXTrackerTpc::fitTrack(temp_track, chi2_hit, scale);
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

static bool fit_all_update(vector<vector<int> >& layer_indexes,
                           SimpleTrack3D& temp_track, vector<int>& best_ind,
                           vector<float>& best_chi2, SimpleTrack3D& track,
                           float& chi2, int iter = 0,
			   float tempscale = 1.0, float scale = 1.0) {
  if (iter != 0) {
    if (remove_bad_hits(track, 5., scale) == false) {
      return false;
    }
    sPHENIXTrackerTpc::fitTrack(track, scale);
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
   
    float tex = (2.0*sqrt(temp_track.hits[i].get_size(0,0))) * tempscale;
    float tey = (2.0*sqrt(temp_track.hits[i].get_size(1,1))) * tempscale;
    float tez = (2.0*sqrt(temp_track.hits[i].get_size(2,2))) * tempscale;
    
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
  }
  sPHENIXTrackerTpc::fitTrack(track, chi2_hit, scale1);

  float tempscale = scale1;
  for (int i = 0; i < 20; ++i) {
    if (fit_all_update(layer_indexes, temp_track, best_ind, best_chi2, track,
                       chi2, i, tempscale, tempscale) == false) {
      return false;
    }
    tempscale *= scale2;

    sPHENIXTrackerTpc::fitTrack(track, chi2_hit, tempscale);
  }

  if (fit_all_update(layer_indexes, temp_track, best_ind, best_chi2, track,
                     chi2, 1.0, 1.0) == false) {
    return false;
  }

  chi2 = sPHENIXTrackerTpc::fitTrack(track, chi2_hit, 1.0);//scale1); // maybe this should be one
  if ((chi2 < 10.0) && (track.hits.size() > ((layer_indexes.size() * 1) / 2))) {
    return true;
  }
  return false;
}

void sPHENIXTrackerTpc::findTracksByCombinatorialKalman(
    vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks,
    const HelixRange& range) {

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

    vector<SimpleHit3D> new_hits;
    for (int i = (((int)(inner_hits.size())) - 1); i >= 0; i -= 1) {
      new_hits.push_back(inner_hits[i]);
    }
    for (unsigned int i = 0; i < temp_track.hits.size(); ++i) {
      new_hits.push_back(temp_track.hits[i]);
    }
    temp_track.hits = new_hits;

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

      gettimeofday(&t2, NULL);
      time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
      time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
      CAtime += (time2 - time1);
      return;
    }
  } else {

    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
    CAtime += (time2 - time1);
    return;
  }
}

void sPHENIXTrackerTpc::findTracksBySegments(vector<SimpleHit3D>& hits,
                                             vector<SimpleTrack3D>& tracks,
                                             const HelixRange& range) {

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

        /// \todo location of a fudge scale factor

        dx1_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(0,0));
        dy1_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(1,1));
        dz1_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(2,2));

        x2_a[hit_counter] = layer_sorted[1][j].get_x();
        y2_a[hit_counter] = layer_sorted[1][j].get_y();
        z2_a[hit_counter] = layer_sorted[1][j].get_z();

        dx2_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(0,0));
        dy2_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(1,1));
        dz2_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(2,2));

        x3_a[hit_counter] = layer_sorted[2][k].get_x();
        y3_a[hit_counter] = layer_sorted[2][k].get_y();
        z3_a[hit_counter] = layer_sorted[2][k].get_z();
        dx3_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(0,0));
        dy3_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(1,1));
        dz3_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(2,2));

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

        dx1_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(0,0));
        dy1_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(1,1));
        dz1_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(2,2));
        dx2_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(0,0));
        dy2_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(1,1));
        dz2_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(2,2));
        dx3_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(0,0));
        dy3_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(1,1));
        dz3_a[hit_counter] = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(2,2));

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

void sPHENIXTrackerTpc::calculateKappaTangents(
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


void sPHENIXTrackerTpc::calculateKappaTangents(
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
