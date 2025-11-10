#include "WeightedTrack.h"

double
WeightedTrack::get_squared_error (
	double const* parameters
) {
	set_parameters(parameters);

	double error{0};
	for (auto const& point : *this) {
		Eigen::Vector3d intersection = get_intersection(point.sensor_local_to_global_transform);
		Eigen::Vector2d weighted_residuals = point.get_weighted_residuals(intersection);
		error += weighted_residuals.squaredNorm();
	}

	if (m_use_vertex) {
		error += m_vertex_position.transpose() * m_vertex_covariance * m_vertex_position;
	}

	return error;
}

Eigen::Matrix<double, 2, 3>
WeightedTrack::get_projection (
	Eigen::Affine3d const& plane
) const {
	Eigen::Matrix<double, 2, 3> projection = Eigen::Matrix<double, 2, 3>::Zero();
	double path_length = get_path_length_of_intersection(plane);
	Eigen::Vector3d slope = get_slope_at_path_length(path_length);

	Eigen::Vector3d z = plane.rotation().col(2);
	for (int i = 0; i < 2; ++i) {
		Eigen::Vector3d alpha = plane.rotation().col(i);
		/// See Eqn. 31 of "Global $\Chi^{2}$ approach to the Alignment of the ATLAS Silicon Tracking Detectors"
		projection.row(i) = alpha - (slope.dot(alpha) / slope.dot(z)) * z;
	}

	return projection;
}

Eigen::Vector2d
WeightedTrack::get_residual_derivative (
	int index,
	Eigen::Affine3d const& plane
) const {
	double path_length = get_path_length_of_intersection(plane);
	return get_projection(plane) * get_partial_derivative(index, path_length);
}

Eigen::Vector2d
ClusterFitPoint::get_residuals (
	Eigen::Vector3d const& intersection 
) const {
	Eigen::Vector3d displacement = intersection - cluster_position; // Defines the sign of the residual
	return Eigen::Vector2d {
		sensor_local_to_global_transform.rotation().col(0).dot(displacement),
		sensor_local_to_global_transform.rotation().col(1).dot(displacement),
	};
}

Eigen::Vector2d
ClusterFitPoint::get_weighted_residuals (
	Eigen::Vector3d const& intersection 
) const {
	Eigen::Vector2d residuals = get_residuals(intersection);
	residuals(0) /= cluster_errors(0);
	residuals(1) /= cluster_errors(1);
	return residuals;
}
