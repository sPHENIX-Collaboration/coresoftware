#ifndef INTT_SURVEY_MAP_H
#define INTT_SURVEY_MAP_H

#include "InttMapping.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/SVD>

namespace InttSurvey
{
	Eigen::Affine3d GetTransform(struct Intt::Offline_s const&);
}

#endif//INTT_SURVEY_MAP_H
