#ifndef WEIGHTED_TRACK_H
#define WEIGHTED_TRACK_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>
#include <Eigen/Dense>

#ifdef __clang__
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wundefined-internal"
	#include <Eigen/Geometry>
	#pragma GCC diagnostic pop
#else
	#include <Eigen/Geometry>
#endif

#include <Math/Minimizer.h>

#include <iostream>
#include <utility>
#include <vector>

class ClusterFitPoint : public PHObject {
public:
	ClusterFitPoint() = default;
	virtual ~ClusterFitPoint() = default;

	PHObject* CloneMe() const override { return new ClusterFitPoint(*this); }

	/// The key of the associated cluster
	TrkrDefs::cluskey cluster_key{};
	/// The position of the cluster in global coordinates
	Eigen::Vector3d cluster_position = Eigen::Vector3d::Zero();
	/// The size of the cluster along the local coordinate axes
	Eigen::Vector2d cluster_errors = Eigen::Vector2d::Ones();

	/// The augmented affine transform matrix that acts on vectors in the sensor coordinates to convert them to global coordinates
	Eigen::Affine3d sensor_local_to_global_transform = Eigen::Affine3d::Identity();

	/// Returns the residuals w.r.t. a given intersection point (in global coordinates)
	Eigen::Vector2d get_residuals (Eigen::Vector3d const& /*intersection*/) const;
	Eigen::Vector2d get_weighted_residuals (Eigen::Vector3d const& /*intersection*/) const;

private:
	ClassDefOverride(ClusterFitPoint, 1)
};

class WeightedTrack : public PHObject, private std::vector<ClusterFitPoint> {
public:
	virtual ~WeightedTrack() = default;

	// Force derived classes to re-implement this method, since the class is pure virtual
	virtual PHObject* CloneMe() const override = 0;

	/// Sets member parameters (of derived classes which actually implement the track) using an array of double
	/// Unfortunately, helper members of derived classes must be declared mutable so this can be declared const
	/// All the following methods assume the parameters have been set using this method
	virtual void set_parameters (double const* /*params*/) = 0;
	/// Mutates the members of the given array such that calling set_parameters would result in the same track parameterization (inverse of set_parameters)
	virtual void get_parameters (double* /*params*/) const = 0;
	/// Returns the number of parameters needed for the fit, e.g. the minimum length of the array passed to set_parameters, get_parameters
	virtual std::size_t get_n_parameters () const = 0;
	/// Modifies the minimizer so it is ready to begin the fit optimization, including setting initial parameters and ranges
	virtual void configure_minimizer (ROOT::Math::Minimizer& /*minimizer*/) = 0;

	/// Returns the partial derivative of the track position w.r.t. the parameter given by index evaluated at the given path length
	/// Used as a helper to get_residual_derivatives, which accounts for the implicit dependence of the intersection on the path length
	virtual Eigen::Vector3d get_partial_derivative (int /*index*/, double /*path_length*/) const = 0;

	/// Evaluates the track position at the given path length
	virtual Eigen::Vector3d get_position_at_path_length (double /*path_length*/) const = 0;
	/// Evaluates the track slope at the given path length
	virtual Eigen::Vector3d get_slope_at_path_length (double /*path_length*/) const = 0;

	/// Finds the path length that gives the intersection with the plane described by the given transform
	/// For use with the sensor_transform member of ClusterFitPoint to find track states
	virtual double get_path_length_of_intersection (Eigen::Affine3d const& /*plane*/) const = 0;
	/// Finds the path length that gives the point of closest approach to the given point
	virtual double get_path_length_of_pca (Eigen::Vector3d const& /*point*/) const = 0;


	/// Returns the 3D position (in global coordinates) of the intersection of this track, with its current parameters, and the plane defined by the given Affine3d
	/// For use with the sensor_transform member of ClusterFitPoint to find track states
	Eigen::Vector3d get_intersection (Eigen::Affine3d const& plane) const { return get_position_at_path_length(get_path_length_of_intersection(plane)); }
	/// Returns the 3D position (in global coordinates) of the point on the track, with its current parameters, closest to the given point
	Eigen::Vector3d get_pca (Eigen::Vector3d const& point) const { return get_position_at_path_length(get_path_length_of_pca(point)); }

	/// Returns the sum of the weighted errors, for use as a functor for a minimizer instance
	double get_squared_error (double const*);
	/// Returns the projection vectors associated with a state as a matrix
	/// See Eqn. 31 of "Global $\Chi^{2}$ approach to the Alignment of the ATLAS Silicon Tracking Detectors"
	Eigen::Matrix<double, 2, 3> get_projection (Eigen::Affine3d const& /*plane*/) const;
	/// Returns the derivatives of the residuals w.r.t. the parameter given by index for the track intersection with the given plane
	Eigen::Vector2d get_residual_derivative (int /*index*/, Eigen::Affine3d const& /*plane*/) const;

	/// Sets whether to use the member vertex for the fit
	void use_vertex (bool use = true) { m_use_vertex = use; }
	/// Returns the vertex position
	Eigen::Vector3d get_vertex_position () { return m_vertex_position; }
	/// Returns the vertex covariance
	Eigen::Matrix3d get_vertex_covariance () { return m_vertex_covariance; }
	/// Sets the vertex position
	void get_vertex_position (Eigen::Vector3d&& vertex) { m_vertex_position = std::forward<Eigen::Vector3d>(vertex); }
	/// Sets the vertex covariance
	void set_vertex_covariance (Eigen::Matrix3d&& covariance) { m_vertex_covariance = std::forward<Eigen::Matrix3d>(covariance); }
	/// Component-wise access to the vertex position
	double& vertex_position (std::size_t i) { return m_vertex_position(i); }
	/// Component-wise access to the vertex covariance
	double& vertex_covariance (std::size_t i, std::size_t j) { return m_vertex_covariance(i, j); }


	// Types
	using std::vector<ClusterFitPoint>::value_type;
	using std::vector<ClusterFitPoint>::allocator_type;
	using std::vector<ClusterFitPoint>::size_type;
	using std::vector<ClusterFitPoint>::difference_type;
	using std::vector<ClusterFitPoint>::reference;
	using std::vector<ClusterFitPoint>::const_reference;
	using std::vector<ClusterFitPoint>::pointer;
	using std::vector<ClusterFitPoint>::const_pointer;
	using std::vector<ClusterFitPoint>::iterator;
	using std::vector<ClusterFitPoint>::const_iterator;
	using std::vector<ClusterFitPoint>::reverse_iterator;
	using std::vector<ClusterFitPoint>::const_reverse_iterator;

	// Element access
	using std::vector<ClusterFitPoint>::at;
	using std::vector<ClusterFitPoint>::operator[];
	using std::vector<ClusterFitPoint>::front;
	using std::vector<ClusterFitPoint>::back;
	using std::vector<ClusterFitPoint>::data;

	// Iterators
	using std::vector<ClusterFitPoint>::begin;
	using std::vector<ClusterFitPoint>::cbegin;
	using std::vector<ClusterFitPoint>::end;
	using std::vector<ClusterFitPoint>::cend;
	using std::vector<ClusterFitPoint>::rbegin;
	using std::vector<ClusterFitPoint>::crbegin;
	using std::vector<ClusterFitPoint>::rend;
	using std::vector<ClusterFitPoint>::crend;

	// Capacity
	using std::vector<ClusterFitPoint>::empty;
	using std::vector<ClusterFitPoint>::size;
	using std::vector<ClusterFitPoint>::max_size;
	using std::vector<ClusterFitPoint>::reserve;
	using std::vector<ClusterFitPoint>::capacity;
	using std::vector<ClusterFitPoint>::shrink_to_fit;

	// Modifiers
	using std::vector<ClusterFitPoint>::clear;
	using std::vector<ClusterFitPoint>::insert;
	// using std::vector<ClusterFitPoint>::insert_range; c++23
	using std::vector<ClusterFitPoint>::emplace;
	using std::vector<ClusterFitPoint>::erase;
	using std::vector<ClusterFitPoint>::push_back;
	using std::vector<ClusterFitPoint>::emplace_back;
	// using std::vector<ClusterFitPoint>::append_range; c++23
	using std::vector<ClusterFitPoint>::pop_back;
	using std::vector<ClusterFitPoint>::resize;
	using std::vector<ClusterFitPoint>::swap;

protected:
	WeightedTrack() = default;

private:
	bool m_use_vertex{false};
	Eigen::Vector3d m_vertex_position = Eigen::Vector3d::Zero();
	Eigen::Matrix3d m_vertex_covariance = Eigen::Matrix3d::Identity();

	ClassDefOverride(WeightedTrack, 1)
};

#endif//WEIGHTED_TRACK_H
