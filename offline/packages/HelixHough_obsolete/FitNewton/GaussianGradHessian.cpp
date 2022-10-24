#include "GaussianGradHessian.h"

#include <cmath>
#include <vector>

using namespace std;
using namespace Eigen;

namespace FitNewton
{
  GaussianGradHessian::GaussianGradHessian() : FunctionGradHessian(3, 1)
  {
    
  }
  
  
  GaussianGradHessian::~GaussianGradHessian()
  {
    
  }
  
  
  FunctionGradHessian* GaussianGradHessian::Clone() const
  {
    GaussianGradHessian* clone = new GaussianGradHessian();
    clone->setFixedPar(0, fixedpars[0]);
    return clone;
  }
  
  
  bool GaussianGradHessian::calcValGradHessian(const VectorXd& x, double& val, VectorXd& grad, MatrixXd& hessian)
  {
    double temp0 = (fixedpars[0] - x(1));
    double temp1 = temp0*temp0;
    double temp4 = 1./(x(2));
    double temp2 = temp4*temp4;
    temp4*=temp2;
    double temp3 = 1./x(0);
    
    val = (x(0))*exp(-temp1*0.5*temp2);
    
    grad(0) = temp3*val;
    grad(1) = val*temp0*temp2;
    grad(2) = val*temp1*temp4;
    
    hessian(0,0) = 0.;
    hessian(0,1) = grad(1)*temp3;
    hessian(0,2) = grad(2)*temp3;
    hessian(1,1) = temp2*(grad(1)*temp0 - val);
    hessian(1,2) = temp0*(grad(2)*temp2 - 2.*temp4*val);
    hessian(2,2) = temp1*(grad(2)*temp4 - 3.*temp2*temp2*val);
    hessian(1,0) = hessian(0,1);
    hessian(2,0) = hessian(0,2);
    hessian(2,1) = hessian(1,2);
    return true;
  }
}


// f = x0*exp(-(t-x1)^2/(2*x2^2))
// 
// df/dx0 = exp(-(t-x1)^2/(2*x2^2)) = f/x0
// df/dx1 = x0*exp(-(t-x1)^2/(2*x2^2))*(t1-x1)/(x2^2) = f*((t-x1)/x2^2)
// df/dx2 = x0*exp(-(t-x1)^2/(2*x2^2))*(t1-x1)^2/(x2^3) = f*(t-x1)^2/(x2^3)
// 
// d^2f/dx0^2 = 0
// d^2f/dx0dx1 = (df/dx1)/x0
// d^2f/dx0dx2 = (df/dx2)/x0
// d^2f/dx1^2 = (df/dx1)*((t-x1)/x2^2) - f/x2^2 = (1/x2^2)*((df/dx1)*(t-x1) - f)
// d^2f/dx1dx2 = (df/dx2)*((t-x1)/x2^2) - 2*f*((t-x1)/x2^3)
// d^2f/dx2^2 = (df/dx2)*(t-x1)^2/(x2^3) - 3*f*(t-x1)^2/(x2^4)

