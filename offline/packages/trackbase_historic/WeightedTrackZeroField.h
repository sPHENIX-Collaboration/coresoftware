#ifndef WEIGHTED_TRACK_ZERO_FIELD_H
#define WEIGHTED_TRACK_ZERO_FIELD_H

#include "WeightedTrack.h"

class WeightedTrackZeroField : public WeightedTrack {
public:
	WeightedTrackZeroField () = default;
	virtual ~WeightedTrackZeroField () = default;

	PHObject* CloneMe() const override { return new WeightedTrackZeroField(*this); }

	/// Sets member parameters (of derived classes which actually implement the track) using an array of double
	/// Unfortunately, helper members of derived classes must be declared mutable so this can be declared const
	/// All the following methods assume the parameters have been set using this method
	void set_parameters (double const* /*params*/) override;
	/// Mutates the members of the given array such that calling set_parameters would result in the same track parameterization (inverse of set_parameters)
	void get_parameters (double* /*params*/) const override;
	/// Returns the number of parameters needed for the fit, e.g. the minimum length of the array passed to set_parameters, get_parameters
	std::size_t get_n_parameters () const override { return 4; }
	/// Modifies the minimizer so it is ready to begin the fit optimization, including setting initial parameters and ranges
	void configure_minimizer (ROOT::Math::Minimizer& /*minimizer*/) override;

	/// Returns the partial derivative of the track position w.r.t. the parameter given by index evaluated at the given path length
	/// Used as a helper to get_residual_derivatives, which accounts for the implicit dependence of the intersection on the path length
	Eigen::Vector3d get_partial_derivative (int /*index*/, double /*path_length*/) const override;

	/// Evaluates the track position at the given path length
	Eigen::Vector3d get_position_at_path_length (double /*path_length*/) const override;
	/// Evaluates the track slope at the given path length
	Eigen::Vector3d get_slope_at_path_length (double /*path_length*/) const override;

	/// Finds the path length that gives the intersection with the plane described by the given transform
	/// For use with the sensor_transform member of ClusterFitPoint to find track states
	double get_path_length_of_intersection (Eigen::Affine3d const& /*plane*/) const override;
	/// Finds the path length that gives the point of closest approach to the given point
	double get_path_length_of_pca (Eigen::Vector3d const& /*point*/) const override;

private:
	Eigen::Vector3d m_slope = Eigen::Vector3d::Zero();
	Eigen::Vector3d m_intercept = Eigen::Vector3d::Zero();

	double m_cos_theta{1.0};
	double m_sin_theta{0.0};
	double m_cos_phi{1.0};
	double m_sin_phi{0.0};

	ClassDefOverride(WeightedTrackZeroField, 1);
};

#endif//WEIGHTED_TRACK_ZERO_FIELD_H
