#include "sPHENIXSeedFinder.h"

#include "CylinderKalman.h"
#include "Pincushion.h"
#include "SimpleTrack3D.h"
#include "vector_math_inline.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <sys/time.h>
#include <utility>
#include <xmmintrin.h>

class HelixResolution;

#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp
#define LogWarning(exp)	std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp

using namespace std;
using namespace Eigen;
using namespace SeamStress;

class hit_triplet {
public:
	hit_triplet(unsigned int h1, unsigned int h2, unsigned int h3,
			unsigned int t, float c) :
			hit1(h1), hit2(h2), hit3(h3), track(t), chi2(c) {
	}
	~hit_triplet() {
	}

	bool operator<(const hit_triplet& other) const {
		return (hit1 < other.hit1)
				|| ((hit2 < other.hit2) && (hit1 == other.hit1))
				|| ((hit3 < other.hit3) && (hit1 == other.hit1)
						&& (hit2 == other.hit2));
	}

	bool operator==(const hit_triplet& other) const {
		return ((hit1 == other.hit1) && (hit2 == other.hit2)
				&& (hit3 == other.hit3));
	}

	unsigned int hit1, hit2, hit3, track;
	float chi2;
};

class hitTriplet {
public:
	hitTriplet(unsigned int h1, unsigned int h2, unsigned int h3,
			unsigned int t, float c) :
			hit1(h1), hit2(h2), hit3(h3), track(t), chi2(c) {
	}
	~hitTriplet() {
	}

	bool operator<(const hitTriplet& other) const {
		return (hit1 < other.hit1)
				|| ((hit2 < other.hit2) && (hit1 == other.hit1))
				|| ((hit3 < other.hit3) && (hit1 == other.hit1)
						&& (hit2 == other.hit2));
	}

	bool operator==(const hitTriplet& other) const {
		return ((hit1 == other.hit1) && (hit2 == other.hit2)
				&& (hit3 == other.hit3));
	}

	unsigned int hit1, hit2, hit3, track;
	float chi2;
};

void sPHENIXSeedFinder::tripletRejection(vector<SimpleTrack3D>& input,
		vector<SimpleTrack3D>& output, vector<bool>& usetrack,
		vector<float>& next_best_chi2) {

	vector<hitTriplet> trips;
	for (unsigned int i = 0; i < input.size(); ++i) {
		for (unsigned int h1 = 0; h1 < input[i].hits.size(); ++h1) {
			for (unsigned int h2 = (h1 + 1); h2 < input[i].hits.size(); ++h2) {
				for (unsigned int h3 = (h2 + 1); h3 < input[i].hits.size(); ++h3) {
					//if (cut_on_dca == false)
					{
						trips.push_back(hitTriplet(
										input[i].hits[h1].get_id(), input[i].hits[h2].get_id(),
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

sPHENIXSeedFinder::sPHENIXSeedFinder(unsigned int n_phi, unsigned int n_d,
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
      cosang_cut(0.985) {
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

sPHENIXSeedFinder::sPHENIXSeedFinder(vector<vector<unsigned int> >& zoom_profile,
                               unsigned int minzoom, HelixRange& range,
                               vector<float>& material, vector<float>& radius,
                               float Bfield, bool parallel,
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
      cosang_cut(0.985) {
	vector<float> detector_material;

	for (unsigned int i = 0; i < radius.size(); ++i) {
		detector_radii.push_back(radius[i]);
	}
	for (unsigned int i = 0; i < material.size(); ++i) {
		detector_scatter.push_back(
				1.41421356237309515 * 0.0136 * sqrt(3. * material[i]));
		detector_material.push_back(3. * material[i]);
	}

	detector_B_field = Bfield;

	integrated_scatter.assign(detector_scatter.size(), 0.);
	float total_scatter_2 = 0.;
	for (unsigned int l = 0; l < detector_scatter.size(); ++l) {
		total_scatter_2 += detector_scatter[l] * detector_scatter[l];
		integrated_scatter[l] = sqrt(total_scatter_2);
	}

	kalman = new CylinderKalman(detector_radii, detector_material,
			detector_B_field);

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

		pins = new Pincushion<sPHENIXSeedFinder>(this, vssp);

		vector<vector<unsigned int> > zoom_profile_new;
		for (unsigned int i = 1; i < zoom_profile.size(); ++i) {
			zoom_profile_new.push_back(zoom_profile[i]);
		}

		for (unsigned int i = 0; i < nthreads; ++i) {
			thread_trackers.push_back(
					new sPHENIXSeedFinder(zoom_profile, minzoom, range,
							material, radius, Bfield));
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

sPHENIXSeedFinder::~sPHENIXSeedFinder() {
	if (kalman != NULL)
		delete kalman;
	for (unsigned int i = 0; i < vss.size(); i++) {
		vss[i].stop();
	}
	for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
		delete thread_trackers[i];
		delete split_output_hits[i];
		delete split_ranges[i];
	}

	if (pins != NULL)
		delete pins;
	if (vssp != NULL)
		delete vssp;
}

float sPHENIXSeedFinder::kappaToPt(float kappa) {
	return detector_B_field / 333.6 / kappa;
}

float sPHENIXSeedFinder::ptToKappa(float pt) {
	return detector_B_field / 333.6 / pt;
}

// hel should be +- 1
static void xyTangent(SimpleHit3D& hit1, SimpleHit3D& hit2, float kappa,
		float hel, float& ux_out, float& uy_out, float& ux_in, float& uy_in) {
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

float sPHENIXSeedFinder::dcaToVertexXY(SimpleTrack3D& track, float vx,
		float vy) {
	float d_out = 0.;

	// find point at the dca to 0
	float x0 = track.d * cos(track.phi);
	float y0 = track.d * sin(track.phi);

	// change variables so x0,y0 -> 0,0
	float phi2 = atan2(
			(1. + track.kappa * track.d) * sin(track.phi) - track.kappa * y0,
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

void sPHENIXSeedFinder::finalize(vector<SimpleTrack3D>& input,
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

	//LogDebug("input size: ")<<input.size()<<endl;

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

	//LogDebug("output size: ")<<output.size()<<endl;

	track_states = states_new;
	if (smooth_back == true) {
		for (unsigned int i = 0; i < output.size(); ++i) {

			HelixKalmanState state = track_states[i];

			track_states[i].C *= 3.;
			track_states[i].chi2 = 0.;
			track_states[i].x_int = 0.;
			track_states[i].y_int = 0.;
			track_states[i].z_int = 0.;
			//       track_states[i].position = output[i].hits.size();
			track_states[i].position = 0;
			for (int h = (output[i].hits.size() - 1); h >= 0; --h) {
				SimpleHit3D hit = output[i].hits[h];
				float err_scale = 1.;
				int layer = hit.get_layer();
				if ((layer >= 0) && (layer < (int) (hit_error_scale.size()))) {
					err_scale = hit_error_scale[layer];
				}
				err_scale *= 3.0; // fudge factor, like needed due to non-gaussian errors

				// \todo location of a rescale fudge factor

				kalman->addHit(hit, track_states[i]);
				track_states[i].position = h;
			}

			//       track_states[i].C = state.C;
			track_states[i].chi2 = state.chi2;
			track_states[i].C *= 2. / 3.;

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

void sPHENIXSeedFinder::findTracks(vector<SimpleHit3D>& hits,
		vector<SimpleTrack3D>& tracks, const HelixRange& range) {

	findtracksiter += 1;
	findTracksBySegments(hits, tracks, range);
}

bool sPHENIXSeedFinder::breakRecursion(const vector<SimpleHit3D>& hits,
		const HelixRange& range) {
	if (seeding == true) {
		return false;
	}
	unsigned int layer_mask[4] = { 0, 0, 0, 0 };
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
	unsigned int nlayers = __builtin_popcount(layer_mask[0])
			+ __builtin_popcount(layer_mask[1])
			+ __builtin_popcount(layer_mask[2])
			+ __builtin_popcount(layer_mask[3]);

	return (nlayers < required_layers);
}

float sPHENIXSeedFinder::phiError(SimpleHit3D& hit, float min_k, float max_k,
		float min_d, float max_d, float min_z0, float max_z0, float min_dzdl,
		float max_dzdl, bool pairvoting) {
	float Bfield_inv = 1. / detector_B_field;
	float p_inv = 0.;

	if ((prev_max_k == max_k) && (prev_max_dzdl == max_dzdl)) {
		p_inv = prev_p_inv;
	} else {
		prev_max_k = max_k;
		prev_max_dzdl = max_dzdl;
		prev_p_inv = 3.33333333333333314e+02 * max_k * Bfield_inv
				* sqrt(1. - max_dzdl * max_dzdl);
		p_inv = prev_p_inv;
	}
	float total_scatter_2 = 0.;
	for (int i = seed_layer + 1; i <= (hit.get_layer()); ++i) {
		float this_scatter = detector_scatter[i - 1]
				* (detector_radii[i] - detector_radii[i - 1])
				/ detector_radii[i];
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

float sPHENIXSeedFinder::dzdlError(SimpleHit3D& hit, float min_k, float max_k,
		float min_d, float max_d, float min_z0, float max_z0, float min_dzdl,
		float max_dzdl, bool pairvoting) {
	float Bfield_inv = 1. / detector_B_field;
	float p_inv = 0.;

	if ((prev_max_k == max_k) && (prev_max_dzdl == max_dzdl)) {
		p_inv = prev_p_inv;
	} else {
		prev_max_k = max_k;
		prev_max_dzdl = max_dzdl;
		prev_p_inv = 3.33333333333333314e+02 * max_k * Bfield_inv
				* sqrt(1. - max_dzdl * max_dzdl);
		p_inv = prev_p_inv;
	}
	float total_scatter_2 = 0.;
	for (int i = seed_layer + 1; i <= (hit.get_layer()); ++i) {
		float this_scatter = detector_scatter[i - 1]
				* (detector_radii[i] - detector_radii[i - 1])
				/ detector_radii[i];
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

void sPHENIXSeedFinder::setRangeFromSeed(HelixRange& range,
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

float sPHENIXSeedFinder::fitTrack(SimpleTrack3D& track) {
	vector<float> chi2_hit;
	return sPHENIXSeedFinder::fitTrack(track, chi2_hit);
}

float sPHENIXSeedFinder::fitTrack(SimpleTrack3D& track,
		vector<float>& chi2_hit) {
	chi2_hit.clear();
	vector<float> xyres;
	vector<float> xyres_inv;
	vector<float> zres;
	vector<float> zres_inv;
	for (unsigned int i = 0; i < track.hits.size(); i++) {

		/// \todo fudge factor location

		xyres.push_back(
				sqrt(
						(0.5 * sqrt(12.) * sqrt(track.hits[i].get_size(0, 0)))
								* (0.5 * sqrt(12.)
										* sqrt(track.hits[i].get_size(0, 0)))
								+ (0.5 * sqrt(12.)
										* sqrt(track.hits[i].get_size(1, 1)))
										* (0.5 * sqrt(12.)
												* sqrt(
														track.hits[i].get_size(
																1, 1)))));
		xyres_inv.push_back(1. / xyres.back());
		zres.push_back((0.5 * sqrt(12.) * sqrt(track.hits[i].get_size(2, 2))));
		zres_inv.push_back(1. / zres.back());
	}

	chi2_hit.resize(track.hits.size(), 0.);

	MatrixXf y = MatrixXf::Zero(track.hits.size(), 1);
	for (unsigned int i = 0; i < track.hits.size(); i++) {
		y(i, 0) =
				(pow(track.hits[i].get_x(), 2) + pow(track.hits[i].get_y(), 2));
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
		float D = sqrt(
				pow(dx - track.hits[i].get_x(), 2)
						+ pow(dy - track.hits[i].get_y(), 2));
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

	return (chi2_tot) / ((float) (deg_of_freedom));
}

void sPHENIXSeedFinder::calculateKappaTangents(float* x1_a, float* y1_a,
		float* z1_a, float* x2_a, float* y2_a, float* z2_a, float* x3_a,
		float* y3_a, float* z3_a, float* dx1_a, float* dy1_a, float* dz1_a,
		float* dx2_a, float* dy2_a, float* dz2_a, float* dx3_a, float* dy3_a,
		float* dz3_a, float* kappa_a, float* dkappa_a, float* ux_mid_a,
		float* uy_mid_a, float* ux_end_a, float* uy_end_a, float* dzdl_1_a,
		float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a) {
	static const __m128 two = { 2., 2., 2., 2. };

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
	tmp1 = (D12 + D23 + D31) * (D23 + D31 - D12) * (D12 + D31 - D23)
			* (D12 + D23 - D31);
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

void sPHENIXSeedFinder::calculateKappaTangents(float* x1_a, float* y1_a,
		float* z1_a, float* x2_a, float* y2_a, float* z2_a, float* x3_a,
		float* y3_a, float* z3_a, float* dx1_a, float* dy1_a, float* dz1_a,
		float* dx2_a, float* dy2_a, float* dz2_a, float* dx3_a, float* dy3_a,
		float* dz3_a, float* kappa_a, float* dkappa_a, float* ux_mid_a,
		float* uy_mid_a, float* ux_end_a, float* uy_end_a, float* dzdl_1_a,
		float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a, float sinang_cut,
		float cosang_diff_inv, float* cur_kappa_a, float* cur_dkappa_a,
		float* cur_ux_a, float* cur_uy_a, float* cur_chi2_a, float* chi2_a) {
	static const __m128 two = { 2., 2., 2., 2. };

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
	tmp1 = (D12 + D23 + D31) * (D23 + D31 - D12) * (D12 + D31 - D23)
			* (D12 + D23 - D31);
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
	__m128 chi2_ang = (one - cos_scatter) * (one - cos_scatter) * cosdiff
			* cosdiff;
	tmp1 = dzdl_1 * sinang;
	_vec_fabs_ps(tmp1);
	__m128 chi2_dzdl = (dzdl_1 - dzdl_2) / (ddzdl_1 + ddzdl_2 + tmp1);
	chi2_dzdl *= chi2_dzdl;
	chi2_dzdl *= one_o_2;

	__m128 n_chi2 = c_chi2 + chi2_ang + chi2_k + chi2_dzdl;
	_mm_store_ps(chi2_a, n_chi2);
}

struct TempComb {
	TempComb() {
	}
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

void sPHENIXSeedFinder::initDummyHits(vector<SimpleHit3D>& dummies,
		const HelixRange& range, HelixKalmanState& init_state) {
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
		dummies[i].set_y(x);
		dummies[i].set_z(x);
		dummies[i].set_layer(i);
	}
}

static void triplet_rejection(vector<SimpleTrack3D>& input,
		vector<float>& chi2s, vector<bool>& usetrack) {
	vector<hit_triplet> trips;
	for (unsigned int i = 0; i < input.size(); ++i) {
		for (unsigned int h1 = 0; h1 < input[i].hits.size(); ++h1) {
			for (unsigned int h2 = (h1 + 1); h2 < input[i].hits.size(); ++h2) {
				for (unsigned int h3 = (h2 + 1); h3 < input[i].hits.size();
						++h3) {
					trips.push_back(
							hit_triplet(input[i].hits[h1].get_id(),
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
				if (input[trips[i].track].hits.size()
						< input[trips[best_pos].track].hits.size()) {
					continue;
				} else if ((input[trips[i].track].hits.size()
						> input[trips[best_pos].track].hits.size())
						|| (input[trips[i].track].hits.back().get_layer()
								> input[trips[best_pos].track].hits.back().get_layer())) {
					next_chi2 = best_chi2;
					best_chi2 = trips[i].chi2;
					best_pos = i;
					continue;
				}
				if ((trips[i].chi2 < best_chi2)
						|| (usetrack[trips[best_pos].track] == false)) {
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

void sPHENIXSeedFinder::findTracksBySegments(vector<SimpleHit3D>& hits,
		vector<SimpleTrack3D>& tracks, const HelixRange& range) {

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
		inv_layer[l] = 1. / (((float) l) - 2.);
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
		for (unsigned int j = 0, sizej = layer_sorted[1].size(); j < sizej;
				++j) {
			for (unsigned int k = 0, sizek = layer_sorted[2].size(); k < sizek;
					++k) {
				if ((layer_sorted[0][i].get_layer()
						>= layer_sorted[1][j].get_layer())
						|| (layer_sorted[1][j].get_layer()
								>= layer_sorted[2][k].get_layer())) {
					continue;
				}

				x1_a[hit_counter] = layer_sorted[0][i].get_x();
				y1_a[hit_counter] = layer_sorted[0][i].get_y();
				z1_a[hit_counter] = layer_sorted[0][i].get_z();

				/// \todo location of a fudge scale factor

				dx1_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[0][i].get_size(0, 0));
				dy1_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[0][i].get_size(1, 1));
				dz1_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[0][i].get_size(2, 2));

				x2_a[hit_counter] = layer_sorted[1][j].get_x();
				y2_a[hit_counter] = layer_sorted[1][j].get_y();
				z2_a[hit_counter] = layer_sorted[1][j].get_z();

				dx2_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[1][j].get_size(0, 0));
				dy2_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[1][j].get_size(1, 1));
				dz2_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[1][j].get_size(2, 2));

				x3_a[hit_counter] = layer_sorted[2][k].get_x();
				y3_a[hit_counter] = layer_sorted[2][k].get_y();
				z3_a[hit_counter] = layer_sorted[2][k].get_z();

				dx3_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[2][k].get_size(0, 0));
				dy3_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[2][k].get_size(1, 1));
				dz3_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[2][k].get_size(2, 2));

				hit1[hit_counter] = i;
				hit2[hit_counter] = j;
				hit3[hit_counter] = k;

				hit_counter += 1;

				if (hit_counter == 4) {
					calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a,
							x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a,
							dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a,
							ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a,
							dzdl_2_a, ddzdl_1_a, ddzdl_2_a);

					for (unsigned int h = 0; h < hit_counter; ++h) {
						temp_segment.chi2 = (dzdl_1_a[h] - dzdl_2_a[h])
								/ (ddzdl_1_a[h] + ddzdl_2_a[h]
										+ fabs(dzdl_1_a[h] * sinang_cut));
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
		calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a, y3_a,
				z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a, dy3_a,
				dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a, ux_end_a,
				uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a, ddzdl_2_a);

		for (unsigned int h = 0; h < hit_counter; ++h) {
			temp_segment.chi2 = (dzdl_1_a[h] - dzdl_2_a[h])
					/ (ddzdl_1_a[h] + ddzdl_2_a[h]
							+ fabs(dzdl_1_a[h] * sinang_cut));
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
			for (unsigned int j = 0, sizej = layer_sorted[l].size(); j < sizej;
					++j) {
				if ((layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_layer()
						>= layer_sorted[l][j].get_layer())) {
					continue;
				}

				x1_a[hit_counter] =
						layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_x();
				y1_a[hit_counter] =
						layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_y();
				z1_a[hit_counter] =
						layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_z();
				x2_a[hit_counter] =
						layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_x();
				y2_a[hit_counter] =
						layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_y();
				z2_a[hit_counter] =
						layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_z();
				x3_a[hit_counter] = layer_sorted[l][j].get_x();
				y3_a[hit_counter] = layer_sorted[l][j].get_y();
				z3_a[hit_counter] = layer_sorted[l][j].get_z();

				dx1_a[hit_counter] =
						0.5 * sqrt(12.0)
								* sqrt(
										layer_sorted[l - 2][(*cur_seg)[i].hits[l
												- 2]].get_size(0, 0));
				dy1_a[hit_counter] =
						0.5 * sqrt(12.0)
								* sqrt(
										layer_sorted[l - 2][(*cur_seg)[i].hits[l
												- 2]].get_size(1, 1));
				dz1_a[hit_counter] =
						0.5 * sqrt(12.0)
								* sqrt(
										layer_sorted[l - 2][(*cur_seg)[i].hits[l
												- 2]].get_size(2, 2));
				dx2_a[hit_counter] =
						0.5 * sqrt(12.0)
								* sqrt(
										layer_sorted[l - 1][(*cur_seg)[i].hits[l
												- 1]].get_size(0, 0));
				dy2_a[hit_counter] =
						0.5 * sqrt(12.0)
								* sqrt(
										layer_sorted[l - 1][(*cur_seg)[i].hits[l
												- 1]].get_size(1, 1));
				dz2_a[hit_counter] =
						0.5 * sqrt(12.0)
								* sqrt(
										layer_sorted[l - 1][(*cur_seg)[i].hits[l
												- 1]].get_size(2, 2));
				dx3_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[l][j].get_size(0, 0));
				dy3_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[l][j].get_size(1, 1));
				dz3_a[hit_counter] = 0.5 * sqrt(12.0)
						* sqrt(layer_sorted[l][j].get_size(2, 2));

				cur_kappa_a[hit_counter] = (*cur_seg)[i].kappa;
				cur_dkappa_a[hit_counter] = (*cur_seg)[i].dkappa;
				cur_ux_a[hit_counter] = (*cur_seg)[i].ux;
				cur_uy_a[hit_counter] = (*cur_seg)[i].uy;
				cur_chi2_a[hit_counter] = (*cur_seg)[i].chi2;

				whichseg[hit_counter] = i;
				hit1[hit_counter] = j;

				hit_counter += 1;
				if (hit_counter == 4) {
					calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a,
							x3_a, y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a,
							dz2_a, dx3_a, dy3_a, dz3_a, kappa_a, dkappa_a,
							ux_mid_a, uy_mid_a, ux_end_a, uy_end_a, dzdl_1_a,
							dzdl_2_a, ddzdl_1_a, ddzdl_2_a, sinang_cut,
							cosang_diff_inv, cur_kappa_a, cur_dkappa_a,
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
								temp_segment.hits[ll] =
										(*cur_seg)[whichseg[h]].hits[ll];
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
			calculateKappaTangents(x1_a, y1_a, z1_a, x2_a, y2_a, z2_a, x3_a,
					y3_a, z3_a, dx1_a, dy1_a, dz1_a, dx2_a, dy2_a, dz2_a, dx3_a,
					dy3_a, dz3_a, kappa_a, dkappa_a, ux_mid_a, uy_mid_a,
					ux_end_a, uy_end_a, dzdl_1_a, dzdl_2_a, ddzdl_1_a,
					ddzdl_2_a, sinang_cut, cosang_diff_inv, cur_kappa_a,
					cur_dkappa_a, cur_ux_a, cur_uy_a, cur_chi2_a, chi2_a);

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
						temp_segment.hits[ll] =
								(*cur_seg)[whichseg[h]].hits[ll];
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
	time1 = ((double) (t1.tv_sec) + (double) (t1.tv_usec) / 1000000.);
	time2 = ((double) (t2.tv_sec) + (double) (t2.tv_usec) / 1000000.);
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
			if (init_chi2
					> fast_chi2_cut_par0
							+ fast_chi2_cut_par1
									/ kappaToPt(temp_track.kappa)) {
				gettimeofday(&t2, NULL);
				time1 =
						((double) (t1.tv_sec) + (double) (t1.tv_usec) / 1000000.);
				time2 =
						((double) (t2.tv_sec) + (double) (t2.tv_usec) / 1000000.);
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
		time1 = ((double) (t1.tv_sec) + (double) (t1.tv_usec) / 1000000.);
		time2 = ((double) (t2.tv_sec) + (double) (t2.tv_usec) / 1000000.);
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
		if (state.chi2 / (2. * ((float) (temp_track.hits.size())) - 5.)
				> chi2_cut) {
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

		/*
		cout << "   add track with hits size " << temp_track.hits.size() 
		     << " z0 " << temp_track.z0
		     << endl;
		for (unsigned int i = 0; i < temp_track.hits.size(); ++i) 
		  {
		    cout << " temp_track.hits.get_id() " << temp_track.hits[i].get_id() << endl;
		  }
		*/

		tracks.push_back(temp_track);
		track_states.push_back(state);
		if ((remove_hits == true) && (state.chi2 < chi2_removal_cut)
				&& (temp_track.hits.size() >= n_removal_hits)) {
			for (unsigned int i = 0; i < temp_track.hits.size(); ++i) {
				(*hit_used)[temp_track.hits[i].get_id()] = true;
			}
		}
	}
}

static inline double sign(double x) {
	return ((double) (x > 0.)) - ((double) (x < 0.));
}

void sPHENIXSeedFinder::projectToLayer(SimpleTrack3D& seed, unsigned int layer,
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

void sPHENIXSeedFinder::initSplitting(vector<SimpleHit3D>& hits,
		unsigned int min_hits, unsigned int max_hits) {
	initEvent(hits, min_hits);
	(*(hits_vec[0])) = hits;
	zoomranges.clear();
	for (unsigned int z = 0; z <= max_zoom; z++) {
		zoomranges.push_back(top_range);
	}
}

void sPHENIXSeedFinder::findHelicesParallelOneHelicity(
		vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits,
		vector<SimpleTrack3D>& tracks) {
	unsigned int hits_per_thread = (hits.size() + 2 * nthreads) / nthreads;
	unsigned int pos = 0;
	while (pos < hits.size()) {
		for (unsigned int i = 0; i < nthreads; ++i) {
			if (pos >= hits.size()) {
				break;
			}
			for (unsigned int j = 0; j < hits_per_thread; ++j) {
				if (pos >= hits.size()) {
					break;
				}
				split_input_hits[i].push_back(hits[pos]);
				pos += 1;
			}
		}
	}
	for (unsigned int i = 0; i < nthreads; ++i) {
		thread_trackers[i]->setTopRange(top_range);
		thread_trackers[i]->initSplitting(split_input_hits[i], thread_min_hits,
				thread_max_hits);
	}
	pins->sewStraight(&sPHENIXSeedFinder::splitHitsParallelThread, nthreads);
	thread_ranges.clear();
	thread_hits.clear();

	unsigned int nbins = split_output_hits[0]->size();
	for (unsigned int b = 0; b < nbins; ++b) {
		thread_ranges.push_back((*(split_ranges[0]))[b]);
		thread_hits.push_back(vector<SimpleHit3D>());
		for (unsigned int i = 0; i < nthreads; ++i) {
			for (unsigned int j = 0; j < (*(split_output_hits[i]))[b].size();
					++j) {
				thread_hits.back().push_back((*(split_output_hits[i]))[b][j]);
			}
		}
	}

	pins->sewStraight(&sPHENIXSeedFinder::findHelicesParallelThread, nthreads);
}

void sPHENIXSeedFinder::findHelicesParallel(vector<SimpleHit3D>& hits,
		unsigned int min_hits, unsigned int max_hits,
		vector<SimpleTrack3D>& tracks) {
	thread_min_hits = min_hits;
	thread_max_hits = max_hits;

	for (unsigned int i = 0; i < nthreads; ++i) {
		thread_tracks[i].clear();
		thread_trackers[i]->clear();
		if (cluster_start_bin != 0) {
			thread_trackers[i]->setClusterStartBin(cluster_start_bin - 1);
		} else {
			thread_trackers[i]->setClusterStartBin(0);
		}
	}

	initSplitting(hits, min_hits, max_hits);

	if (separate_by_helicity == true) {
		for (unsigned int i = 0; i < nthreads; ++i) {
			thread_trackers[i]->setSeparateByHelicity(true);
			thread_trackers[i]->setOnlyOneHelicity(true);
			thread_trackers[i]->setHelicity(true);
			split_output_hits[i]->clear();
			split_input_hits[i].clear();
		}
		findHelicesParallelOneHelicity(hits, min_hits, max_hits, tracks);

		for (unsigned int i = 0; i < nthreads; ++i) {
			thread_trackers[i]->setSeparateByHelicity(true);
			thread_trackers[i]->setOnlyOneHelicity(true);
			thread_trackers[i]->setHelicity(false);
			split_output_hits[i]->clear();
			split_input_hits[i].clear();
		}
		findHelicesParallelOneHelicity(hits, min_hits, max_hits, tracks);
	} else {
		for (unsigned int i = 0; i < nthreads; ++i) {
			thread_trackers[i]->setSeparateByHelicity(false);
			thread_trackers[i]->setOnlyOneHelicity(false);
			split_output_hits[i]->clear();
			split_input_hits[i].clear();
		}

		findHelicesParallelOneHelicity(hits, min_hits, max_hits, tracks);
	}

	vector<SimpleTrack3D> temp_tracks;
	for (unsigned int i = 0; i < nthreads; ++i) {
		vector<HelixKalmanState>* states =
				&(thread_trackers[i]->getKalmanStates());
		for (unsigned int j = 0; j < thread_tracks[i].size(); ++j) {
			track_states.push_back((*states)[j]);
			temp_tracks.push_back(thread_tracks[i][j]);
		}
	}
	finalize(temp_tracks, tracks);
}

void sPHENIXSeedFinder::splitHitsParallelThread(void* arg) {
	unsigned long int w = (*((unsigned long int*) arg));
	thread_trackers[w]->splitIntoBins(thread_min_hits, thread_max_hits,
			*(split_ranges[w]), *(split_output_hits[w]), 0);
}

void sPHENIXSeedFinder::findHelicesParallelThread(void* arg) {
	unsigned long int w = (*((unsigned long int*) arg));

	for (unsigned int i = w; i < thread_ranges.size(); i += nthreads) {
		if (thread_hits[i].size() == 0) {
			continue;
		}
		thread_trackers[w]->setTopRange(thread_ranges[i]);
		thread_trackers[w]->findHelices(thread_hits[i], thread_min_hits,
				thread_max_hits, thread_tracks[w]);
	}
}
