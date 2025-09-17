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
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <thread>

//! Provides a best-fit cubic approximation of the Field
//! In 3D cartesian coordinates
class PHFieldInterpolated : public PHField
{
public:
	typedef Eigen::Vector3i Indices_t;
	typedef Eigen::Vector3f Point_t;
	typedef Eigen::Vector3f Field_t;

	//! constructor
	explicit PHFieldInterpolated () = default;

	//! destructor
	~PHFieldInterpolated () override = default;

	//! access field value
	//! Follow the convention of G4ElectroMagneticField
	//! @param[in]	pointer to double[4] representing (x, y, z, t), where field is to be evaluated (in Geant4/CLHEP units) (cm)
	//! @param[out] pointer to double[3] representing (b_x, b_y, b_z), mutated to contain the value of the field (in Geant4/CLHEP units) (T)
	void GetFieldValue (double const*, double*) const override;

	//! Returns an O(3) best fit interpolation of the field at the point, updating the cached interpolated parameters as needed
	//! Thread-safe top-level access
	Field_t get_interpolated (Point_t const&) const;

	//! Loads the fieldmap from a ROOT file, expects specific contents
	void load_fieldmap (
		std::string const& = "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.499/share/calibrations/Field/Map/sphenix3dtrackingmapxyz.root",
		const float& = 1.0
	);

	//! prints class info based on verbosity level
	//! thread-safe
	void print(std::ostream& = std::cout) const;

	//! prints the cached Taylor coefficients
	//! thread-safe
	void print_coefficients(std::ostream& = std::cout) const;

private:

	// Threads using the class get their own cache
	struct InterpolationCache {

		//! Tracks how "important" this cache is
		//! When a thread accesses this class, its corresponding cache is moved to the front of the queue
		//! Instances at the end of the "queue" are culled
		//! Maximal value allows idential logic to shuffle indexes for existing or new thread accesses
		std::size_t m_queue_index{std::numeric_limits<std::size_t>::max()};

		//! The indices of the left-down-back corner of the cell that the interpolation is cached for
		//! Initialized out-of-bounds so the coefficients are always computed on first call
		Indices_t m_buffered_indices = {-1, -1, -1};

		//! The point about which the field is currently interpolated
		Point_t m_center;

		//! The 20 coefficients for an O(3) Taylor expansion of a function in 3D
		//! (one per component of the magnetic field)
		std::array<Eigen::VectorXf, 3>	m_coefficients;
	};

	//! Returns the 20 element row of a design matrix for position point about the stored value m_center
	static Eigen::VectorXf get_design_vector (Point_t const&, InterpolationCache const&);

	//! Returns a mutable reference to the interpolation information cached for the calling thread
	InterpolationCache& get_cache () const;

	//! Caches the interpolation for the cell containing given point
	//! Should really be a member funciton of InterpolationCache
	void cache_interpolation (Point_t const&, InterpolationCache&) const;

	//! Gets the 3D grid indices from a 1D deque index
	Indices_t get_indices (std::size_t const&) const;
	//! Gets the 1D deque index from 3D grid indices
	std::size_t get_index (Indices_t const&) const;

	//! return the 3D grid indices of the left-down-back corner of the grid cell containing point
	Indices_t get_indices (Point_t const&) const;
	//! return point at the left-down-back corner of the grid cell at the given 3D grid indices
	Point_t get_point (Indices_t const&) const;

	//! return field at given 1D deque index
	Field_t get_field (std::size_t const& index) const { return m_field.at(index); }
	//! return field at the left-down-back corner of the grid cell at given 3D grid indices
	Field_t get_field (Indices_t const& indices) const { return get_field(get_index(indices)); }
	//! return the field at the left-down-back corner of the grid cell containing point
	Field_t get_field (Point_t const& point) const { return get_field(get_indices(point)); }

	//! throw if the 3D grid indices are not in the grid and would give an invalid index
	void validate_indices (Indices_t const&) const;
	//! throw if the point lies outside the grid and would give an invalid index
	void validate_point (Point_t const&) const;

	//! The most threads will we accomodate before pruning the cache map between calls
	static std::size_t const MAX_THREADS = 8;

	Indices_t m_N; // Number of points along an axis
	Point_t m_D; // Grid spacing of an axis
	Point_t m_min; // Minimum value of an axis
	Point_t m_max; // Maximum value of an axis

	//! The entire fieldmap, loaded into memory (yikes)
	//! deque for O(1) access, without requiring contiguous allocation
	std::deque<Field_t> m_field;

	// Map of thread ids to caches
	typedef std::map<std::thread::id, InterpolationCache> map_t;
	mutable map_t m_caches;

	// mutex for access
	mutable std::mutex m_mutex;
};
