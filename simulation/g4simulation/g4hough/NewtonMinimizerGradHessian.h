#ifndef G4HOUGH_NEWTONMINIMIZERGRADHESSIAN_H
#define G4HOUGH_NEWTONMINIMIZERGRADHESSIAN_H

#include <Eigen/Core>  // for VectorXd

#include <vector>      // for vector

namespace NewtonMinimizer
{
  class FunctionGradHessian;
  class NewtonMinimizerGradHessian
  {
    public:
      NewtonMinimizerGradHessian();
      ~NewtonMinimizerGradHessian();
      
      void setFunction(FunctionGradHessian* func);
      
      bool minimize(const Eigen::VectorXd& start_point, Eigen::VectorXd& min_point, double tol=0x1.0p-30, unsigned int max_iter=1024, double abs_tol=0.);
      bool findSaddlePoint(const Eigen::VectorXd& start_point, Eigen::VectorXd& min_point, double tol=0x1.0p-30, unsigned int max_iter=1024, double abs_tol=0.);
      
      void fixParameter(unsigned int par);
      void unfixParameter(unsigned int par);
      
    private:
      FunctionGradHessian* function;
      bool zoom(const double& wolfe1, const double& wolfe2, double& lo, double& hi, Eigen::VectorXd& try_grad, Eigen::VectorXd& direction, double& grad0_dir, double& val0, double& val_lo, Eigen::VectorXd& init_params, Eigen::VectorXd& try_params, unsigned int max_iter, double& result);
      bool lineSearch(double& alpha, const double& wolfe1, const double& wolfe2, Eigen::VectorXd& try_grad, Eigen::VectorXd& direction, double& grad0_dir, double& val0, Eigen::VectorXd& init_params, Eigen::VectorXd& try_params, const double& precision, const double& accuracy, unsigned int max_iter, double& result);
      
      std::vector<int> fixparameter;
      
  };
  
}

#endif

