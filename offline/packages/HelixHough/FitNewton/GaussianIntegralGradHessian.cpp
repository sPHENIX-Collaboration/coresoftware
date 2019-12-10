#include "GaussianIntegralGradHessian.h"

#include <cmath>
#include <vector>

using namespace std;
using namespace Eigen;

namespace FitNewton
{
  GaussianIntegralGradHessian::GaussianIntegralGradHessian() : FunctionGradHessian(3, 2)
  {
    
  }
  
  
  GaussianIntegralGradHessian::~GaussianIntegralGradHessian()
  {
    
  }
  
  
  FunctionGradHessian* GaussianIntegralGradHessian::Clone() const
  {
    GaussianIntegralGradHessian* clone = new GaussianIntegralGradHessian();
    clone->setFixedPar(0, fixedpars[0]);
    clone->setFixedPar(1, fixedpars[1]);
    return clone;
  }
  
  
  bool GaussianIntegralGradHessian::calcValGradHessian(const VectorXd& x, double& val, VectorXd& grad, MatrixXd& hessian)
  {
    double temp0 = (fixedpars[0] - x(1));
    double temp0_2 = temp0*temp0;
    double temp1 = (fixedpars[1] - x(1));
    double temp1_2 = temp1*temp1;
    double temp2 = 1./x(2);
    double sqrt_pi_o_2 = 0x1.40d931ff627059657ca41fae722cep0;
    double sqrt_1_o_2 = 0xb.504f333f9de6484597d89b3754abep-4;
    
    val = sqrt_pi_o_2*x(0)*x(2)*( erf(sqrt_1_o_2*temp1*temp2) - erf(sqrt_1_o_2*temp0*temp2) );
    
    double temp2_2 = temp2*temp2;
    double etemp1 = exp(-temp1_2*0.5*temp2_2);
    double etemp0 = exp(-temp0_2*0.5*temp2_2);
    
    double temp4 = 1./x(0);
    double temp5 = ( -etemp1*temp1 + etemp0*temp0 );
    
    grad(0) = val*temp4;
    grad(1) = x(0)*(-etemp1 + etemp0);
    grad(2) = temp2*( val + x(0)*temp5 );
    
    hessian(0,0) = 0.;
    hessian(0,1) = temp4*grad(1);
    hessian(0,2) = temp4*grad(2);
    hessian(1,1) = x(0)*temp2_2*temp5;
    hessian(1,2) = x(0)*temp2*temp2_2*( -etemp1*temp1_2 + etemp0*temp0_2 );
    hessian(2,2) = x(0)*temp2_2*temp2_2*( -etemp1*temp1_2*temp1 + -etemp0*temp0_2*temp0 );
    return true;
  }
}


// f = sqrt(pi/2.)*x0*x2*(erf((1/sqrt(2))*(t1-x1)/x2) - erf((1/sqrt(2))*(t0-x1)/x2))
// 
// df/dx0 = f/x0
// df/dx1 = x0*[ -exp(-(t1-x1)^2/(2*x2^2)) + exp(-(t0-x1)^2/(2*x2^2)) ]
// df/dx2 = (1/x2)*(f + x0*[ -exp(-(t1-x1)^2/(2*x2^2))*(t1-x1) + exp(-(t0-x1)^2/(2*x2^2))*(t0-x1) ])
// 
// d^2f/dx0^2 = 0
// d^2f/dx0dx1 = (df/dx1)/x0
// d^2f/dx0dx2 = (df/dx2)/x0
// d^2f/dx1^2 = (x0/x2^2)*[ -exp(-(t1-x1)^2/(2*x2^2))*(t1-x1) + exp(-(t0-x1)^2/(2*x2^2))*(t0-x1) ]
// d^2f/dx1dx2 = (x0/x2^3)*[ -exp(-(t1-x1)^2/(2*x2^2))*(t1-x1)^2 + exp(-(t0-x1)^2/(2*x2^2))*(t0-x1)^2 ]
// d^2f/dx2^2 = (x0/x2^4)*[ -exp(-(t1-x1)^2/(2*x2^2))*(t1-x1)^3 + exp(-(t0-x1)^2/(2*x2^2))*(t0-x1)^3  ]

