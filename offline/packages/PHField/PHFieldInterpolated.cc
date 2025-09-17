#include "PHFieldInterpolated.h"

#include <phool/phool.h> // PHWHERE

#include <TFile.h>
#include <TTree.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <iostream>

#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <boost/format.hpp>

void
PHFieldInterpolated::load_fieldmap (
	std::string const& filename,
	float const& magfield_rescale
) {
	// Necessary as calling c_str() on an empty string returns NULL,
	// but TFile::Open doesn't have a guard clause for a NULL argument
	if (filename.empty()) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " Empty filename";
		throw std::runtime_error(what.str());
	}

	TFile* file = TFile::Open(filename.c_str(), "READ");
	if (!file) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " Could not open file " << filename;
		throw std::runtime_error(what.str());
	}

	TTree* tree{};
	file->GetObject("fieldmap", tree);
	if (!tree) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " Could not get tree " << "fieldmap"
			<< " from file " << filename;
		throw std::runtime_error(what.str());
	}

	Point_t point;
	Field_t field;
	std::map<std::string, float&> branches = {
		{"x", point(0)}, {"y", point(1)}, {"z", point(2)},
		{"bx", field(0)}, {"by", field(1)}, {"bz", field(2)},
	};

	for (auto& [name, reference] : branches) {
		if (tree->SetBranchAddress(name.c_str(), &reference) != TTree::kMatch) {
			std::stringstream what;
			what
				<< PHWHERE
				<< " Could not get branch " << name
				<< " from tree " << "fieldmap"
				<< " from file " << filename;
			throw std::runtime_error(what.str());
		}
	}

	if (Verbosity()) {
		std::cout
			<< PHWHERE << "\n"
			<< " Loading fieldmap from file " << filename
			<< std::endl;
	}

	// The unique values taken on coordinate axes
	// used to check that the tree is a perfect grid
	std::array<std::set<float>, 3> points;

	for (std::size_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
		tree->GetEntry(n);

		// Unit conversions
		point *= cm;
		for (int i = 0; i < 3; ++i) {
			points[i].insert(point(i));
		}
	}

	for (int i = 0; i < 3; ++i) {
		if (points[i].size() < 2) {
			std::stringstream what;
			what
				<< PHWHERE
				<< " not enough points";
			throw std::runtime_error(what.str());
		}
		m_N(i) = points[i].size();
		m_min(i) = *points[i].begin();
		m_max(i) = *points[i].rbegin();
		m_D(i) = (m_max(i) - m_min(i)) / (m_N(i) - 1);
	}

	if (static_cast<Long64_t>(m_N(0)) * static_cast<Long64_t>(m_N(1)) * static_cast<Long64_t>(m_N(2)) != tree->GetEntriesFast()) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " tree is not a grid";
		throw std::runtime_error(what.str());
	}

	std::lock_guard lock(m_mutex);
	m_field.resize(tree->GetEntriesFast());
	for (std::size_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
		tree->GetEntry(n);

		// Unit conversions
		point *= cm;
		field *= tesla * magfield_rescale;

		// Check that the deterministic computation
		// actually matches what we're reading in
		Indices_t indices = get_indices(point);
		if (!point.isApprox(get_point(indices))) {
			std::stringstream what;
			what
				<< PHWHERE
				<< " point does not round trip"
				<< " (with " << point.transpose() << ")";
			throw std::runtime_error(what.str());
		}

		m_field[get_index(indices)] = field;
	}

	file->Close();
}

void
PHFieldInterpolated::GetFieldValue (
	double const* point_as_arr, // pointer to (immutable) double[4]
	double* field_as_arr // pointer to (mutable) double[3]
) const {
	Point_t point {
		(float)point_as_arr[0],
		(float)point_as_arr[1],
		(float)point_as_arr[2],
	};

	for (int i = 0; i < 3; ++i) {
		field_as_arr[i] = 0;
	}

	// Catches points out of bounds and early returns leaving field as 0
	try {
		validate_point(point);
	} catch (std::exception const&) {
		if (1 < Verbosity()) {
			std::cout
				<< PHWHERE
				<< " Returning 0 for point out of fieldmap bounds"
				<< " (at " << point.transpose() << ")"
				<< std::endl;
		}
		return;
	}

	// This should not throw if point is in bounds, which we just checked
	// If this throws, there is a bug in class logic
	Field_t field = get_interpolated(point);
	for (int i = 0; i < 3; ++i) {
		field_as_arr[i] = field(i);
	}
}

Eigen::VectorXf
PHFieldInterpolated::get_design_vector (
	Point_t const& point,
	InterpolationCache const& cache
) {
	float x = point(0) - cache.m_center(0);
	float y = point(1) - cache.m_center(1);
	float z = point(2) - cache.m_center(2);

	Eigen::VectorXf design_vector(20);
	design_vector <<
		// O(0)
		1,

		// O(1)
		x, y, z,

		// O(2)
		x*x, x*y, x*z,
		y*y, y*z,
		z*z,

		// O(3)
		x*x*x, x*x*y, x*x*z, x*y*y, x*y*z, x*z*z,
		y*y*y, y*y*z, y*z*z,
		z*z*z;

	return design_vector;
}

void
PHFieldInterpolated::cache_interpolation (
	Point_t const& point,
	InterpolationCache& cache
) const {
	Indices_t indices = get_indices(point);

	// We get neighbors offset by [-1, +2] relative to the buffered point
	// If we're close enough to an edge, this range isn't valid
	// buffer about a shifted voxel further in instead
	for (int i = 0; i < 3; ++i ) {
		while (indices(i) - 1 < 0) { ++indices(i); }
		while (m_N(i) < indices(i) + 3) { --indices(i); }
	}

	// Our coefficients have already been computed for the voxel we want to evaluate in
	if ((indices - cache.m_buffered_indices).norm() == 0) { return; }

	cache.m_buffered_indices = indices;

	// The point at the the center of the voxel
	// get_point gets the left-down-back corner
	cache.m_center = get_point(indices);
	for (int i = 0; i < 3; ++i) {
		cache.m_center(i) += m_D(i) / 2;
	}

	// Effectively we are solving three scalar interpolation problems
	// However, all systems of equations will share the same design matrix
	Eigen::MatrixXf M(64, 20);
	std::array<Eigen::VectorXf, 3> solution_vectors {
		Eigen::VectorXf(64),
		Eigen::VectorXf(64),
		Eigen::VectorXf(64),
	};

	// Get the 64 (nearest) neighboring points about the cell containing point and its neighboring cells
	// Note that the indices are of the left-down-back corner of this cell
	for (int row = 0; row < 64; ++row) {
		indices = {
			cache.m_buffered_indices(0) + (row / 16) - 1,
			cache.m_buffered_indices(1) + ((row / 4) % 4) - 1,
			cache.m_buffered_indices(2) + (row % 4) - 1,
		};

		M.row(row) = get_design_vector(get_point(indices), cache);
		for (int i = 0; i < 3; ++i) {
			solution_vectors[i](row) = get_field(indices)(i);
		}
	}

	// Use Eigen to solve the least squares problem we've set up
	for (int i = 0; i < 3; ++i) {
		cache.m_coefficients[i] = M.bdcSvd (
			Eigen::ComputeThinU | Eigen::ComputeThinV
		).solve (solution_vectors[i]);
	}
}

PHFieldInterpolated::InterpolationCache&
PHFieldInterpolated::get_cache (
) const {
	// Update the access counts
	InterpolationCache& this_cache = m_caches[std::this_thread::get_id()];
	for (auto& [thread_id, cache] : m_caches) {
		if (cache.m_queue_index < this_cache.m_queue_index) { ++cache.m_queue_index; }
	}
	this_cache.m_queue_index = 0;

	return this_cache;
}

PHFieldInterpolated::Field_t
PHFieldInterpolated::get_interpolated (
	Point_t const& point
) const {

	InterpolationCache this_cache;

	// Get a copy of the interpolation information this thread is using
	{
		std::lock_guard lock(m_mutex);
		this_cache = get_cache();
	}

	// Update the cache to be about the point
	cache_interpolation(point, this_cache);

	// Update its place in the member map
	{
		std::lock_guard lock(m_mutex);
		get_cache() = this_cache;

		// Prune map entries which haven't been used in a while
		std::erase_if (m_caches, [](auto const& key_val_pair) {
			return MAX_THREADS <= key_val_pair.second.m_queue_index;
		});
	}

	return {
		get_design_vector(point, this_cache).dot(this_cache.m_coefficients[0]),
		get_design_vector(point, this_cache).dot(this_cache.m_coefficients[1]),
		get_design_vector(point, this_cache).dot(this_cache.m_coefficients[2]),
	};
}

PHFieldInterpolated::Indices_t
PHFieldInterpolated::get_indices (
	std::size_t const& index
) const {
	return {
		((int)index / (m_N(1) * m_N(2))) % m_N(0),
		((int)index / m_N(2)) % m_N(1),
		(int)index % m_N(2),
	};
}

std::size_t
PHFieldInterpolated::get_index (
	Indices_t const& indices
) const {
	validate_indices(indices);
	return
		(indices(0) * m_N(1) * m_N(2)) +
		(indices(1) * m_N(2)) +
		indices(2);
}

PHFieldInterpolated::Indices_t
PHFieldInterpolated::get_indices (
	Point_t const& point
) const {
	validate_point(point);
	return {
		(int)std::floor(point(0) / m_D(0)) + (m_N(0) / 2),
		(int)std::floor(point(1) / m_D(1)) + (m_N(1) / 2),
		(int)std::floor(point(2) / m_D(2)) + (m_N(2) / 2),
	};
}

PHFieldInterpolated::Point_t
PHFieldInterpolated::get_point (
	Indices_t const& indices
) const {
	validate_indices(indices);
	return {
		static_cast<float>(indices(0) - ((m_N(0) - 1.0) / 2.0)) * m_D(0),
		static_cast<float>(indices(1) - ((m_N(1) - 1.0) / 2.0)) * m_D(1),
		static_cast<float>(indices(2) - ((m_N(2) - 1.0) / 2.0)) * m_D(2),
	};
}

void
PHFieldInterpolated::validate_indices (
	Indices_t const& indices
) const {
	for (int i = 0; i < 3; ++i) {
		if (indices(i) < 0 || indices(i) >= m_N(i)) {
			std::stringstream what;
			what
				<< PHWHERE
				<< " Component out of range"
				<< " (at " << indices.transpose() << ")";
			throw std::runtime_error(what.str());
		}
	}
}

void
PHFieldInterpolated::validate_point (
	Point_t const& point 
) const {
	for (int i = 0; i < 3; ++i) {
		if (m_max(i) < point(i) || point(i) < m_min(i)) {
			std::stringstream what;
			what
				<< PHWHERE
				<< " Component out of range"
				<< " (at " << point.transpose() << ")";
			throw std::runtime_error(what.str());
		}
	}
}

void
PHFieldInterpolated::print (
	std::ostream& stream
) const {
	std::lock_guard lock(m_mutex);

	stream
		<< PHWHERE << "\n"
		<< " size: " << m_field.size()
		<< " num caches: " << m_caches.size()
		<< std::endl;

	if (Verbosity() < 1) { return; }

	for (int i = 0; i < 3; ++i) {
		stream
			<< " component: " << i
			<< " min: " << m_min(i)
			<< " max: " << m_max(i)
			<< " step: " << m_D(i)
			<< " count: " << m_N(i)
			<< std::endl;
	}

	for (auto const& [thread_id, cache] : m_caches) {
		stream
			<< " thread id: " << thread_id
			<< " queue index: " << cache.m_queue_index
			<< std::endl;
	}

	if (Verbosity() < 2) { return; }

	for (std::size_t index = 0; index < m_field.size(); ++index) {
		Indices_t indices = get_indices(index);
		stream
			<< " indices: " << indices.transpose()
			<< " point: " << get_point(indices).transpose()
			<< " field: " << get_field(indices).transpose()
			<< std::endl;

		if (Verbosity() < 3 && index == 10) {
			stream << " ..." << std::endl;
			return;
		}
	}
}

void
PHFieldInterpolated::print_coefficients (
	std::ostream& stream
) const {
	std::lock_guard lock(m_mutex);
	InterpolationCache const& cache = m_caches[std::this_thread::get_id()];
	for (int i = 0; i < 3; ++i) {
		stream
			<< cache.m_coefficients[i].transpose()
			<< std::endl;
	}
}

