#ifndef __GAUSSIANINTEGRALGRADHESSIAN__
#define __GAUSSIANINTEGRALGRADHESSIAN__

#include "FunctionGradHessian.h"

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

