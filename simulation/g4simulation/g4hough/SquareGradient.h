#ifndef __SQUAREGRADIENT__
#define __SQUAREGRADIENT__

#include "FunctionGradHessian.h"

namespace NewtonMinimizer
{
  class SquareGradient : public FunctionGradHessian
  {
    public:
      SquareGradient(FunctionGradHessian* func);
      ~SquareGradient();
      
      // this class is intended to be used only with line-search, and so this function doesn't fill the Hessian, and really should not be called
      bool calcValGradHessian(const Eigen::VectorXd& x, double& val, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian);
      
      // fill the gradient of this function based on "function"
      // this class represents the square of the magnitude of the gradient of "function",
      // so we can calculate the val and grad from the val,grad,hessian of function
      bool calcValGrad(const Eigen::VectorXd& x, double& val, Eigen::VectorXd& grad);
      
      FunctionGradHessian* Clone() const;
      
      bool verbose;
      
    protected:
      FunctionGradHessian* function;
  };
}




#endif
