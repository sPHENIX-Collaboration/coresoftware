#ifndef FITNEWTON_GAUSSIANINTEGRALGRADHESSIAN_H
#define FITNEWTON_GAUSSIANINTEGRALGRADHESSIAN_H

#include "FunctionGradHessian.h"

#include <Eigen/Core>

namespace FitNewton
{
  class GaussianIntegralGradHessian : public FunctionGradHessian
  {
    public:
      GaussianIntegralGradHessian();
      ~GaussianIntegralGradHessian();
      
      bool calcValGradHessian(const Eigen::VectorXd& x, double& val, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian);
      
      FunctionGradHessian* Clone() const;
  };
}




#endif

