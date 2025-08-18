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
				<< " point does not round trip,"
				<< " tree is not a grid"
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
	try {
		Field_t field = get_interpolated ({
			(float)point_as_arr[0],
			(float)point_as_arr[1],
			(float)point_as_arr[2],
		});
		for (int i = 0; i < 3; ++i) {
			field_as_arr[i] = field(i);
		}
	} catch (std::exception const& e) {
		std::cout
			<< PHWHERE << "\n"
			<< "\t" << e.what() << "\n"
			<< std::endl;
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
	if ((indices - cache.m_buffered_indices).norm() == 0) { return; }
	cache.m_buffered_indices = indices;

	// The point at the the center of the cell
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

	// Get the 64 neighboring points about the cell containing point and its neighboring cells
	// Note that the indices are of the left-down-back corner of this cell
	for (int row = 0; row < 64; ++row) {
		indices = {
			cache.m_buffered_indices(0) + (row / 16) - 1,
			cache.m_buffered_indices(1) + ((row / 4) % 4) - 1,
			cache.m_buffered_indices(2) + (row % 4) - 1,
		};

		// Possible to throw while searching neighbors
		// Re-throw an error, but with a different message
		try {
			validate_indices(indices);
		} catch (std::exception const&) {
			cache.m_buffered_indices = {-1, -1, -1};
			std::stringstream what;
			what
				<< PHWHERE
				<< " Point too close to edge of fieldmap "
				<< " (at " << point.transpose() << ")";
			throw std::runtime_error(what.str());
		}

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

PHFieldInterpolated::Field_t
PHFieldInterpolated::get_interpolated (
	Point_t const& point
) const {
	// Because a separate access can change the map,
	// this must be locked from the top
	std::lock_guard lock(m_mutex);

	// Update the access counts
	InterpolationCache& this_cache = m_caches[std::this_thread::get_id()];
	for (auto& [thread_id, cache] : m_caches) {
		if (cache.m_queue_index < this_cache.m_queue_index) { ++cache.m_queue_index; }
	}
	this_cache.m_queue_index = 0;

	// // Since c++20
	// std::erase_if (m_caches, [](auto const& key_val_pair) {
	// 	return MAX_THREADS <= key_val_pair.second.m_queue_index;
	// });

	// Manual loop for c++17 compatability
	for (auto itr = m_caches.begin(); itr != m_caches.end();) {
		if (MAX_THREADS <= itr->second.m_queue_index) {
			itr = m_caches.erase(itr);
		} else {
			++itr;
		}
	}

	// Reference could be dangling after erase, erase_if
	this_cache = m_caches[std::this_thread::get_id()];
	cache_interpolation(point, this_cache);
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
				<< " (at " << indices << ")";
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
				<< " (at " << point << ")";
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

