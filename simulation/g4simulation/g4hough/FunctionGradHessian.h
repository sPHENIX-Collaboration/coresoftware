#ifndef __FUNCTIONGRADHESSIAN__
#define __FUNCTIONGRADHESSIAN__

#include <vector>
#include <Eigen/Core>

namespace NewtonMinimizer
{
  class FunctionGradHessian
  {
    public:
      FunctionGradHessian(unsigned int nparams=1, unsigned int nfixedparams=1) : npars(nparams), nfixedpars(nfixedparams)
      {
        fixedpars.assign(nfixedpars, 0.);
      }
      
      virtual ~FunctionGradHessian(){}
      
      //false if parameters x are not in allowed region
      virtual bool calcValGradHessian(const Eigen::VectorXd& x, double& val, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian) = 0;
      virtual bool calcValGrad(const Eigen::VectorXd& x, double& val, Eigen::VectorXd& grad)
      {
        Eigen::MatrixXd hess = Eigen::MatrixXd::Zero(grad.size(), grad.size());
        return calcValGradHessian(x, val, grad, hess);
      }
      
      unsigned int nPars(){return npars;}
      unsigned int nFixedPars(){return nfixedpars;}
      std::vector<double> getFixedPars(){return fixedpars;}
      
      void setFixedPar(unsigned int coor, double val){fixedpars[coor]=val;}
      
      //create a new instance of this class from the current instance
      virtual FunctionGradHessian* Clone() const = 0;
      
      virtual void computeCovariance(const double& val, const Eigen::MatrixXd& hessian){}
      
      virtual void rescaleMove(const Eigen::VectorXd& pars, Eigen::VectorXd& move){}
      
      
      
    protected:
      unsigned int npars;
      unsigned int nfixedpars;
      
      std::vector<double> fixedpars;
  };
}




#endif

