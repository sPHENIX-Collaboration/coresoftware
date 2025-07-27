#pragma once

#include "PHField.h"

#if defined(__GNUC__) && !defined(__clang__)
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
	#include <Eigen/Dense>
	#pragma GCC diagnostic pop
#else
	#include <Eigen/Dense>
#endif

#include <array>
#include <deque>
#include <map>
#include <optional>
#include <string>

//! Provides a best-fit cubic approximation of the Field
//! In 3D cartesian coordinates
class PHFieldInterpolated : public PHField
{
public:
	typedef Eigen::Vector3i Indices_t;
	typedef Eigen::Vector3f Point_t;
	typedef Eigen::Vector3f Field_t;

	//! There is only one fieldmap, and all points lie on a perfect grid
	static int const GRID_COUNT = 111;
	static float constexpr GRID_STEP = 20;
	static float constexpr GRID_MIN = - GRID_STEP * (GRID_COUNT / 2);
	static float constexpr GRID_MAX = + GRID_STEP * (GRID_COUNT / 2);

	//! constructor
	//! @param[in] path to calibration file
	//! @param[in] magnetic field scaling factor (default 1.0)
	explicit PHFieldInterpolated () = default;

	//! destructor
	~PHFieldInterpolated () override = default;

	//! access field value
	//! Follow the convention of G4ElectroMagneticField
	//! @param[in]	pointer to double[4] representing (x, y, z, t), where field is to be evaluated (in Geant4/CLHEP units) (cm)
	//! @param[out] pointer to double[3] representing (b_x, b_y, b_z), mutated to contain the value of the field (in Geant4/CLHEP units) (T)
	void GetFieldValue (double const*, double*) const override;

	//! Loads the fieldmap from a ROOT file, expects specific contents
	void load_fieldmap (
		std::string const& = "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.499/share/calibrations/Field/Map/sphenix3dtrackingmapxyz.root",
		const float& = 1.0
	);

	//! Returns the 20 element row of a design matrix for position point about the stored value m_center
	Eigen::VectorXf get_design_vector (Point_t const&) const;
	//! Caches the interpolation for the cell containing given point
	void cache_interpolation (Point_t const&) const;
	//! Returns an O(3) best fit interpolation of the field at the point, updating the cached interpolated parameters as needed
	Field_t get_interpolated (Point_t const&) const;

	//! Gets the 3D grid indices from a 1D deque index
	static Indices_t get_indices (std::size_t const&);
	//! Gets the 1D deque index from 3D grid indices
	static std::size_t get_index (Indices_t const&);

	//! return the 3D grid indices of the left-down-back corner of the grid cell containing point
	static Indices_t get_indices (Point_t const&);
	//! return point at the left-down-back corner of the grid cell at the given 3D grid indices
	static Point_t get_point (Indices_t const&);

	//! return field at given 1D deque index
	Field_t get_field (std::size_t const& index) const { return m_field.at(index); }
	//! return field at the left-down-back corner of the grid cell at given 3D grid indices
	Field_t get_field (Indices_t const& indices) const { return get_field(get_index(indices)); }
	//! return the field at the left-down-back corner of the grid cell containing point
	Field_t get_field (Point_t const& point) const { return get_field(get_indices(point)); }

	//! throw if the 3D grid indices are not in the grid and would give an invalid index
	static void validate_indices (Indices_t const&);
	//! throw if the point lies outside the grid and would give an invalid index
	static void validate_point (Point_t const&);

	//! prints the contents of the loaded map to std::cout
	void print_map() const;

	//! prints the cached Taylor coefficients used for interpolation to std::cout
	void print_coefficients() const;

private:
	//! The entire fieldmap, loaded into memory (yikes)
	//! deque for O(1) access, without requiring contiguous allocation
	std::deque<Field_t> m_field;

	// the following should be wrapped in a class
	// a set of these should be maintained to buffer more cells

	//! The indices of the left-down-back corner of the cell that the interpolation is cached for
	mutable Indices_t m_buffered_indices = {-1, -1, -1};
	//! The point about which the field is currently interpolated
	mutable Point_t m_center;
	//! The 20 coefficients for an O(3) Taylor expansion of a function in 3D
	//! (one per component of the magnetic field)
	mutable std::array<Eigen::VectorXf, 3>	m_coefficients;
};
