#ifndef FITNEWTON_GAUSSIANGRADHESSIAN_H
#define FITNEWTON_GAUSSIANGRADHESSIAN_H

#include "FunctionGradHessian.h"

#include <Eigen/Core>

namespace FitNewton
{
  class GaussianGradHessian : public FunctionGradHessian
  {
    public:
      GaussianGradHessian();
      ~GaussianGradHessian();
      
      bool calcValGradHessian(const Eigen::VectorXd& x, double& val, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian);
      
      FunctionGradHessian* Clone() const;
  };
}




#endif

