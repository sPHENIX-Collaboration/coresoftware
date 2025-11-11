#include "WeightedTrackZeroField.h"

#include <phool/phool.h>

#include <stdexcept>

void
WeightedTrackZeroField::set_parameters (
	double const* params
) {
	m_cos_theta = std::cos(params[2]);
	m_sin_theta = std::sin(params[2]);
	m_cos_phi = std::cos(params[3]);
	m_sin_phi = std::sin(params[3]);

	m_intercept = Eigen::Vector3d {
		params[0], // x-component of xy plane intercept
		params[1], // y-component of xy plane intercept
		0,
	};

	m_slope = Eigen::Vector3d {
		m_sin_theta * m_cos_phi,
		m_sin_theta * m_sin_phi,
		m_cos_theta
	};
}

void
WeightedTrackZeroField::get_parameters (
	double* params
) const {
	params[0] = m_intercept(0);
	params[1] = m_intercept(1);
	params[2] = std::acos(m_slope(2));
	params[3] = std::atan2(m_slope(1), m_slope(0));
}

void
WeightedTrackZeroField::configure_minimizer (
	ROOT::Math::Minimizer& minimizer
) {
	if (size() < 2) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " too few points to perform fit"
			<< std::endl;
		throw std::runtime_error(what.str());
	}

	// Initial slope
	m_slope = (operator[](1).cluster_position - operator[](0).cluster_position).normalized();
	double phi = std::atan2(m_slope(1), m_slope(0));
	double theta = std::acos(m_slope(2));
	m_cos_theta = m_slope(2);
	m_sin_theta = std::sqrt(m_slope(0) * m_slope(0) + m_slope(1) * m_slope(1));
	m_cos_phi = m_slope(0) / m_sin_theta;
	m_sin_phi = m_slope(1) / m_sin_theta;

	// Initial intercept
	// Solve for s:
	// 0 = (m_slope * s + operator[](0).cluster_position).dot({0, 0, 1})
	double s = -operator[](0).cluster_position(2) / m_slope(2);
	m_intercept = m_slope * s + operator[](0).cluster_position;

	minimizer.SetVariable(0, "x",     m_intercept(0), 1E-6);
	minimizer.SetVariable(1, "y",     m_intercept(1), 1E-6);
	minimizer.SetVariable(2, "theta", theta,          1E-6);
	minimizer.SetVariable(3, "phi",   phi,            1E-6);

	// Bounds for angular parameters
	minimizer.SetVariableLimits(2,          0,     3.1416);
	minimizer.SetVariableLimits(3, phi - 1.58, phi + 1.58);
}

Eigen::Vector3d
WeightedTrackZeroField::get_partial_derivative (
	int index,
	double path_length
) const {
	switch (index) {
		case 0: return Eigen::Vector3d {1.0, 0.0, 0.0};
		case 1: return Eigen::Vector3d {0.0, 1.0, 0.0};
		case 2: return Eigen::Vector3d {  m_cos_theta * m_cos_phi, m_cos_theta * m_sin_phi, -m_sin_theta} * path_length;
		case 3: return Eigen::Vector3d { -m_sin_theta * m_sin_phi, m_sin_theta * m_cos_phi, 0.0} * path_length;
		default: return Eigen::Vector3d::Zero();
	}
}

Eigen::Vector3d
WeightedTrackZeroField::get_position_at_path_length (
	double path_length
) const {
	return m_slope * path_length + m_intercept;
}

Eigen::Vector3d
WeightedTrackZeroField::get_slope_at_path_length (
	double /*path_length*/
) const {
	return m_slope;
}

double
WeightedTrackZeroField::get_path_length_of_intersection (
	Eigen::Affine3d const& plane
) const {
	// Solve s: 
	// 0 = (m_slope * s + m_intercept - plane.translation()).dot(plane.rotation().col(2))
	return (plane.translation() - m_intercept).dot(plane.rotation().col(2)) / m_slope.dot(plane.rotation().col(2));
}

double
WeightedTrackZeroField::get_path_length_of_pca (
	Eigen::Vector3d const& point
) const {
	// Solve s:
	// 0 = (m_slope * s + m_intercept - point).dot(m_slope)
	return (point - m_intercept).dot(m_slope);
}

