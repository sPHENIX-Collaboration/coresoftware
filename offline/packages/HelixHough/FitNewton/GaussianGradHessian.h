#ifndef __GAUSSIANGRADHESSIAN__
#define __GAUSSIANGRADHESSIAN__

#include "FunctionGradHessian.h"

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

