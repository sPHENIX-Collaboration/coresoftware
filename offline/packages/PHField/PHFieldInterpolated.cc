#include "PHFieldInterpolated.h"

#include <phool/phool.h> // PHWHERE

#include <TFile.h>
#include <TTree.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <iostream>
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

	for (std::size_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
		tree->GetEntry(n);

		// Unit conversions
		point *= cm;
		field *= tesla * magfield_rescale;

		// Check if what we read in is consistent with our deterministic computations of the domain
		Point_t expected = get_point(get_indices(n));
		if (1.0E-4 < (point - expected).norm()) {
			std::stringstream what;
			what
				<< PHWHERE
				<< " Read point from file which is dissimilar to calculated value"
				<< " entry: " << n
				<< " expected: " << expected.transpose()
				<< " read: " << point.transpose();
			throw std::runtime_error(what.str());
		};

		m_field.push_back(field);

		// Print information when < 5 and an ellipsis when == 5
		int print_index{6};
		switch (Verbosity()) {
			case 0:
			case 1:
				// Left at 6 and will print nothing
				break;
			case 2:
				// Prints the next 5 terms after every x coordinate roll over
				print_index = (int)n % (GRID_COUNT * GRID_COUNT * GRID_COUNT);
				break;
			case 3:
				// Prints the next 5 terms after every y coordinate roll over
				print_index = (int)n % (GRID_COUNT * GRID_COUNT);
				break;
			case 4:
				// Prints the next 5 terms after every z coordinate roll over
				print_index = (int)n % (GRID_COUNT);
				break;
			default:
				// Prints the whole map
				print_index = 0;
				break;
		}

		if (print_index < 5) {
			std::cout
				<< " index: " << n
				<< " point: " << point.transpose()
				<< " field: " << field.transpose()
				<< std::endl;
		}

		if (print_index == 5) {
			std::cout
				<< " ... "
				<< std::endl;
		}
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
		// If you're from a printout,
		// you'll need to look around this file
		// for the literal "e.what()"
		std::cout
			<< PHWHERE << "\n"
			<< "\t" << e.what() << "\n"
			<< std::endl;
	}
}

Eigen::VectorXf
PHFieldInterpolated::get_design_vector (
	Point_t const& point
) const {
	float x = point(0) - m_center(0);
	float y = point(1) - m_center(1);
	float z = point(2) - m_center(2);

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
	Point_t const& point
) const {
	Indices_t indices = get_indices(point);
	if ((indices - m_buffered_indices).norm() == 0) { return; }
	m_buffered_indices = indices;

	// The point at the the center of the cell
	// get_point gets the left-down-back corner
	m_center = get_point(indices);
	for (int i = 0; i < 3; ++i) {
		m_center(i) += GRID_STEP / 2;
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
			m_buffered_indices(0) + (row / 16) - 1,
			m_buffered_indices(1) + ((row / 4) % 4) - 1,
			m_buffered_indices(2) + (row % 4) - 1,
		};

		// Possible to throw while searching neighbors
		// Re-throw an error, but with a different message
		try {
			validate_indices(indices);
		} catch (std::exception const&) {
			m_buffered_indices = {-1, -1, -1};
			std::stringstream what;
			what
				<< PHWHERE
				<< " Point too close to edge of fieldmap "
				<< " (at " << point << ")";
			throw std::runtime_error(what.str());
		}

		M.row(row) = get_design_vector(get_point(indices));
		for (int i = 0; i < 3; ++i) {
			solution_vectors[i](row) = get_field(indices)(i);
		}
	}

	for (int i = 0; i < 3; ++i) {
		m_coefficients[i] = M.bdcSvd (
			Eigen::ComputeThinU | Eigen::ComputeThinV
		).solve (solution_vectors[i]);
	}
}

PHFieldInterpolated::Field_t
PHFieldInterpolated::get_interpolated (
	Point_t const& point
) const {
	cache_interpolation(point);
	return {
		get_design_vector(point).dot(m_coefficients[0]),
		get_design_vector(point).dot(m_coefficients[1]),
		get_design_vector(point).dot(m_coefficients[2]),
	};
}

PHFieldInterpolated::Indices_t
PHFieldInterpolated::get_indices (
	std::size_t const& index
) {
	return {
		((int)index / (GRID_COUNT * GRID_COUNT)),
		((int)index / GRID_COUNT) % GRID_COUNT,
		((int)index) % GRID_COUNT,
	};
}

std::size_t
PHFieldInterpolated::get_index (
	Indices_t const& indices
) {
	validate_indices(indices);
	return
		(indices(0) * GRID_COUNT * GRID_COUNT) +
		(indices(1) * GRID_COUNT) +
		indices(2);
}

PHFieldInterpolated::Indices_t
PHFieldInterpolated::get_indices (
	Point_t const& point
) {
	validate_point(point);
	return {
		(int)std::floor(point(0) / GRID_STEP) + (GRID_COUNT / 2),
		(int)std::floor(point(1) / GRID_STEP) + (GRID_COUNT / 2),
		(int)std::floor(point(2) / GRID_STEP) + (GRID_COUNT / 2),
	};
}

PHFieldInterpolated::Point_t
PHFieldInterpolated::get_point (
	Indices_t const& indices
) {
	validate_indices(indices);
	return {
		// NOLINTNEXTLINE(bugprone-integer-division)
		(indices(0) - (GRID_COUNT / 2)) * GRID_STEP,
		// NOLINTNEXTLINE(bugprone-integer-division)
		(indices(1) - (GRID_COUNT / 2)) * GRID_STEP,
		// NOLINTNEXTLINE(bugprone-integer-division)
		(indices(2) - (GRID_COUNT / 2)) * GRID_STEP,
	};
}

void
PHFieldInterpolated::validate_indices (
	Indices_t const& indices
) {
	for (int i = 0; i < 3; ++i) {
		if (indices(i) < 0 || indices(i) >= GRID_COUNT) {
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
) {
	for (int i = 0; i < 3; ++i) {
		if (GRID_MAX < point(i) || point(i) < GRID_MIN) {
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
PHFieldInterpolated::print_map (
) const {
	for (std::size_t index = 0; index < m_field.size(); ++index) {
		Indices_t indices = get_indices(index);
		std::cout
			<< " indices: " << indices.transpose()
			<< " point: " << get_point(indices).transpose()
			<< " field: " << get_field(indices).transpose()
			<< std::endl;
	}
}

void
PHFieldInterpolated::print_coefficients (
) const {
	for (int i = 0; i < 3; ++i) {
		std::cout
			<< m_coefficients[i].transpose()
			<< std::endl;
	}
}

