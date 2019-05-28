#ifndef SPHENIX_SPHENIXSEEDFINDER_H
#define SPHENIX_SPHENIXSEEDFINDER_H

#include "HelixHough.h"
#include "HelixKalmanState.h"
#include "SimpleHit3D.h"
#include "Seamstress.h"


#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <vector>

class CylinderKalman;
class HelixResolution;
class SimpleTrack3D;

namespace SeamStress { template <class TClass> class Pincushion; }

class AngleIndexPair {
public:
	AngleIndexPair(float ang, unsigned int idx) :
			angle(ang), index(idx) {
		float twopi = 2. * M_PI;
		int a = (int) (angle / twopi);
		angle -= a * twopi;
		while (angle < 0.) {
			angle += twopi;
		}
		while (angle >= twopi) {
			angle -= twopi;
		}
	}
	~AngleIndexPair() {
	}

	bool operator<(const AngleIndexPair& other) const {
		return angle < other.angle;
	}

	static float absDiff(float angle1, float angle2) {
		float diff = (angle1 - angle2);
		while (diff > M_PI) {
			diff -= 2. * M_PI;
		}
		while (diff < -M_PI) {
			diff += 2. * M_PI;
		}
		diff = fabs(diff);
		return diff;
	}

	float angle;
	unsigned int index;
};

class AngleIndexList {
public:
	AngleIndexList() :
			sorted(false) {
	}
	~AngleIndexList() {
	}

	void addPair(AngleIndexPair& angind) {
		sorted = false;
		vec.push_back(angind);
	}

	void getRangeListSimple(float angle, float error,
			std::vector<AngleIndexPair*>& result) {
		result.clear();

		for (unsigned int i = 0; i < vec.size(); i++) {
			if (AngleIndexPair::absDiff(angle, vec[i].angle) <= error) {
				result.push_back(&(vec[i]));
			}
		}
	}

	void getRangeList(float angle, float error,
			std::vector<AngleIndexPair*>& result) {
		float twopi = 2. * M_PI;
		int a = (int) (angle / twopi);
		angle -= a * twopi;
		while (angle < 0.) {
			angle += twopi;
		}
		while (angle >= twopi) {
			angle -= twopi;
		}

		if (vec.size() <= 4) {
			return getRangeListSimple(angle, error, result);
		}

		result.clear();

		unsigned int closest = findClosest(angle);
		// first, traverse upward
		unsigned int current = closest;
		unsigned int lowest = 0;
		unsigned int highest = vec.size() - 1;
		while (true) {
			if (AngleIndexPair::absDiff(angle, vec[current].angle) <= error) {
				result.push_back(&(vec[current]));
				current = (current + 1) % (vec.size());
				if (current == closest) {
					return;
				}
			} else {
				break;
			}
		}

		// now, traverse downward
		if (current <= closest) {
			lowest = current;
		} else {
			highest = current;
		}
		current = ((closest + vec.size()) - 1) % (vec.size());
		while (true) {
			if ((current == lowest) || (current == highest)) {
				break;
			}
			if (AngleIndexPair::absDiff(angle, vec[current].angle) <= error) {
				result.push_back(&(vec[current]));
				current = ((current + vec.size()) - 1) % (vec.size());
			} else {
				break;
			}
		}
	}

	unsigned int findClosestSimple(float angle, unsigned int lower,
			unsigned int upper) {
		unsigned int closest = lower;
		float diff = AngleIndexPair::absDiff(vec[closest].angle, angle);
		for (unsigned int i = (lower + 1); i <= upper; i++) {
			float tempdiff = AngleIndexPair::absDiff(vec[i].angle, angle);
			if (tempdiff < diff) {
				closest = i;
				diff = tempdiff;
			}
		}

		return closest;
	}

	unsigned int findClosest(float angle) {
		if (vec.size() <= 4) {
			return findClosestSimple(angle, 0, vec.size() - 1);
		}

		if (sorted == false) {
			std::sort(vec.begin(), vec.end());
			sorted = true;
		}

		unsigned int lower = 0;
		unsigned int upper = vec.size() - 1;
		unsigned int middle = vec.size() / 2;
		while (true) {
			if ((upper - lower) <= 4) {
				return findClosestSimple(angle, lower, upper);
			}

			if (angle <= vec[middle].angle) {
				upper = middle;
				middle = (lower + upper) / 2;
			} else {
				lower = middle;
				middle = (lower + upper) / 2;
			}
		}
	}

	std::vector<AngleIndexPair> vec;
	bool sorted;
};

class TrackSegment {
public:
	TrackSegment() :
			chi2(0.), ux(0.), uy(0.), kappa(0.), dkappa(0.), seed(0), n_hits(0) {
	}
	~TrackSegment() {
	}

	float chi2;
	float ux, uy;
	float kappa;
	float dkappa;
	std::vector<unsigned int> hits;
	unsigned int seed;
	unsigned int bin;
	unsigned int n_hits;
};

class sPHENIXSeedFinder: public HelixHough {
public:
	sPHENIXSeedFinder(unsigned int n_phi, unsigned int n_d, unsigned int n_k,
			unsigned int n_dzdl, unsigned int n_z0,
			HelixResolution& min_resolution, HelixResolution& max_resolution,
			HelixRange& range, std::vector<float>& material,
			std::vector<float>& radius, float Bfield);
	sPHENIXSeedFinder(std::vector<std::vector<unsigned int> >& zoom_profile,
			unsigned int minzoom, HelixRange& range,
			std::vector<float>& material, std::vector<float>& radius,
			float Bfield, bool parallel = false, unsigned int num_threads = 1);
	virtual ~sPHENIXSeedFinder();

	void finalize(std::vector<SimpleTrack3D>& input,
			std::vector<SimpleTrack3D>& output);

	void findTracks(std::vector<SimpleHit3D>& hits,
			std::vector<SimpleTrack3D>& tracks, const HelixRange& range);

	void findTracksBySegments(std::vector<SimpleHit3D>& hits,
			std::vector<SimpleTrack3D>& tracks, const HelixRange& range);

	void initEvent(std::vector<SimpleHit3D>& hits, unsigned int min_hits) {
		int min_layer = 999999;
		int max_layer = 0;
		for (unsigned int i = 0; i < hits.size(); ++i) {
			if (hits[i].get_layer() < min_layer) {
				min_layer = hits[i].get_layer();
			}
			if (hits[i].get_layer() > max_layer) {
				max_layer = hits[i].get_layer();
			}
		}
		setSeedLayer(min_layer);
		setNLayers(max_layer + 1);

		nfits = 0;
		combos.clear();
		fail_combos.clear();
		pass_combos.clear();
		seeding = false;
		//     required_layers = min_hits;
		required_layers = req_layers;
		findtracksiter = 0;
		CAtime = 0.;
		KALtime = 0.;
		findtracks_bin = 0;
	}

	void initSeeding(std::vector<SimpleTrack3D>& seeds) {
		seeding = true;
		seed_used.clear();
		seed_used.assign(seeds.size(), false);
		AngleIndexList templist;
		angle_list.clear();
		angle_list.assign(n_layers, templist);
	}

	void clear() {
		track_states.clear();
		isolation_variable.clear();
	}

	void findHelicesParallel(std::vector<SimpleHit3D>& hits,
			unsigned int min_hits, unsigned int max_hits,
			std::vector<SimpleTrack3D>& tracks);
	void findHelicesParallelOneHelicity(std::vector<SimpleHit3D>& hits,
			unsigned int min_hits, unsigned int max_hits,
			std::vector<SimpleTrack3D>& tracks);

	float dcaToVertexXY(SimpleTrack3D& track, float vx, float vy);

	bool breakRecursion(const std::vector<SimpleHit3D>& hits,
			const HelixRange& range);

	float phiError(SimpleHit3D& hit, float min_k, float max_k, float min_d,
			float max_d, float min_z0, float max_z0, float min_dzdl,
			float max_dzdl, bool pairvoting = false);
	float dzdlError(SimpleHit3D& hit, float min_k, float max_k, float min_d,
			float max_d, float min_z0, float max_z0, float min_dzdl,
			float max_dzdl, bool pairvoting = false);

	static float fitTrack(SimpleTrack3D& track);
	static float fitTrack(SimpleTrack3D& track, std::vector<float>& chi2_hit);

	void setVerbosity(int v) {
		verbosity = v;
	}

	void setCutOnDca(bool dcut) {
		cut_on_dca = dcut;
	}
	void setDcaCut(float dcut) {
		dca_cut = dcut;
	}
	void setVertex(float vx, float vy, float vz) {
		vertex_x = vx;
		vertex_y = vy;
		vertex_z = vz;
	}

	void setRejectGhosts(bool rg) {
		reject_ghosts = rg;
	}

	std::vector<float>& getIsolation() {
		return isolation_variable;
	}

	static void calculateKappaTangents(float* x1_a, float* y1_a, float* z1_a,
			float* x2_a, float* y2_a, float* z2_a, float* x3_a, float* y3_a,
			float* z3_a, float* dx1_a, float* dy1_a, float* dz1_a, float* dx2_a,
			float* dy2_a, float* dz2_a, float* dx3_a, float* dy3_a,
			float* dz3_a, float* kappa_a, float* dkappa_a, float* ux_mid_a,
			float* uy_mid_a, float* ux_end_a, float* uy_end_a, float* dzdl_1_a,
			float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a);
	void calculateKappaTangents(float* x1_a, float* y1_a, float* z1_a,
			float* x2_a, float* y2_a, float* z2_a, float* x3_a, float* y3_a,
			float* z3_a, float* dx1_a, float* dy1_a, float* dz1_a, float* dx2_a,
			float* dy2_a, float* dz2_a, float* dx3_a, float* dy3_a,
			float* dz3_a, float* kappa_a, float* dkappa_a, float* ux_mid_a,
			float* uy_mid_a, float* ux_end_a, float* uy_end_a, float* dzdl_1_a,
			float* dzdl_2_a, float* ddzdl_1_a, float* ddzdl_2_a,
			float sinang_cut, float cosang_diff_inv, float* cur_kappa_a,
			float* cur_dkappa_a, float* cur_ux_a, float* cur_uy_a,
			float* cur_chi2_a, float* chi2_a);
	void projectToLayer(SimpleTrack3D& seed, unsigned int layer, float& x,
			float& y, float& z);

	void setRangeFromSeed(HelixRange& range, SimpleTrack3D& seed);

	void pairRejection(std::vector<SimpleTrack3D>& input,
			std::vector<SimpleTrack3D>& output, std::vector<bool>& usetrack,
			std::vector<float>& next_best_chi2);
	void tripletRejection(std::vector<SimpleTrack3D>& input,
			std::vector<SimpleTrack3D>& output, std::vector<bool>& usetrack,
			std::vector<float>& next_best_chi2);

	bool seedWasUsed(unsigned int seed_index) {
		return seed_used.at(seed_index);
	}

	void setSeparateByHelicity(bool sbh) {
		separate_by_helicity = sbh;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setSeparateByHelicity(sbh);
			}
		}
	}

	void requireLayers(unsigned int nl) {
		check_layers = true;
		req_layers = nl;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->requireLayers(nl);
			}
		}
	}

	virtual void setMaxHitsPairs(unsigned int mhp) {
		max_hits_pairs = mhp;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setMaxHitsPairs(mhp);
			}
		}
	}

	void setBinScale(float b_scl) {
		bin_scale = b_scl;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setBinScale(b_scl);
			}
		}
	}
	void setZBinScale(float b_scl) {
		z_bin_scale = b_scl;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setZBinScale(b_scl);
			}
		}
	}
	void setRemoveHits(bool rh) {
		remove_hits = rh;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setRemoveHits(rh);
			}
		}
	}

	void setSeedLayer(int sl) {
		seed_layer = sl;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setSeedLayer(sl);
			}
		}
	}

	void setNLayers(unsigned int n) {
		n_layers = n;
		std::vector<SimpleHit3D> one_layer;
		layer_sorted.clear();
		layer_sorted.assign(n_layers, one_layer);
		temp_comb.assign(n_layers, 0);
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setNLayers(n);
			}
		}
	}

	void setFastChi2Cut(float par0, float par1, float max) {
		fast_chi2_cut_par0 = par0;
		fast_chi2_cut_par1 = par1;
		fast_chi2_cut_max = max;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setFastChi2Cut(par0, par1, max);
			}
		}
	}
	void setChi2Cut(float c) {
		chi2_cut = c;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setChi2Cut(c);
			}
		}
	}
	void setChi2RemovalCut(float c) {
		chi2_removal_cut = c;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setChi2RemovalCut(c);
			}
		}
	}

	void setNRemovalHits(unsigned int n) {
		n_removal_hits = n;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setNRemovalHits(n);
			}
		}
	}

	void setCosAngleCut(float cut) {
		cosang_cut = cut;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setCosAngleCut(cut);
			}
		}
	}

	void setCellularAutomatonChi2Cut(float cut) {
		ca_chi2_cut = cut;
		if (is_parallel == true) {
			for (unsigned int i = 0; i < thread_trackers.size(); ++i) {
				thread_trackers[i]->setCellularAutomatonChi2Cut(cut);
			}
		}
	}

	void setThread() {
		is_thread = true;
	}

	void initSplitting(std::vector<SimpleHit3D>& hits, unsigned int min_hits,
			unsigned int max_hits);

	void setHitErrorScale(unsigned int layer, float scale) {
		if (layer >= hit_error_scale.size()) {
			hit_error_scale.resize(layer + 1, 1.);
		}
		hit_error_scale[layer] = scale;
	}

private:
	float kappaToPt(float kappa);
	float ptToKappa(float pt);

	void findHelicesParallelThread(void* arg);
	void splitHitsParallelThread(void* arg);

	void initDummyHits(std::vector<SimpleHit3D>& dummies,
			const HelixRange& range, HelixKalmanState& init_state);

	float fast_chi2_cut_par0;
	float fast_chi2_cut_par1;
	float fast_chi2_cut_max;
	float chi2_cut;
	float chi2_removal_cut;
	unsigned int n_removal_hits;
	std::set<std::vector<unsigned int> > combos;
	std::map<std::vector<unsigned int>, HelixKalmanState> pass_combos;
	std::set<std::vector<unsigned int> > fail_combos;
	std::vector<float> detector_scatter;
	std::vector<float> integrated_scatter;
	std::vector<float> detector_radii;
	bool seeding;
	int verbosity;
	bool cut_on_dca;
	float dca_cut;
	float vertex_x, vertex_y, vertex_z;
	unsigned int required_layers;
	bool reject_ghosts;
	unsigned int nfits;

	CylinderKalman* kalman;

	std::vector<float> isolation_variable;

	unsigned int findtracksiter;

	float prev_max_k;
	float prev_max_dzdl;
	float prev_p_inv;

	float detector_B_field;

	std::vector<TrackSegment> segments1;
	std::vector<TrackSegment> segments2;

	std::vector<std::vector<SimpleHit3D> > layer_sorted;
	std::vector<unsigned int> temp_comb;

	int seed_layer;

	std::vector<bool> seed_used;

	unsigned int nthreads;
	std::vector<SeamStress::Seamstress*> *vssp;
	std::vector<SeamStress::Seamstress> vss;
	SeamStress::Pincushion<sPHENIXSeedFinder> *pins;
	std::vector<sPHENIXSeedFinder*> thread_trackers;
	std::vector<std::vector<SimpleTrack3D> > thread_tracks;
	std::vector<HelixRange> thread_ranges;
	std::vector<std::vector<SimpleHit3D> > thread_hits;
	std::vector<std::vector<SimpleHit3D> > split_input_hits;
	std::vector<std::vector<std::vector<SimpleHit3D> >*> split_output_hits;
	std::vector<std::vector<HelixRange>*> split_ranges;
	bool is_parallel;
	unsigned int thread_min_hits;
	unsigned int thread_max_hits;
	bool is_thread;

	std::vector<AngleIndexList> angle_list;

	double CAtime;
	double KALtime;

	std::vector<std::vector<SimpleHit3D> > layer_sorted_1[4];
	unsigned int findtracks_bin;

	float ca_chi2_cut;
	float cosang_cut;

	std::vector<float> hit_error_scale;
};

#endif
