#include "SquareGradient.h"

#include <vector>

using namespace std;
using namespace Eigen;


namespace FitNewton
{
  SquareGradient::SquareGradient(FunctionGradHessian* func) : FunctionGradHessian(func->nPars(), func->nFixedPars()), verbose(true), function(func)
  {
    fixedpars = func->getFixedPars();
  }
  
  
  SquareGradient::~SquareGradient()
  {
    
  }
  
  
  FunctionGradHessian* SquareGradient::Clone() const
  {
    SquareGradient* clone = new SquareGradient(function);
    return clone;
  }
  
  
  bool SquareGradient::calcValGrad(const VectorXd& x, double& val, VectorXd& grad)
  {
    // this class is the square of the gradient of function
    // val = sum_i fgrad(i)^2
    // grad(j) = sum_i 2*fgrad(i)*fhessian(i,j)
    
    double fval=0.;
    VectorXd fgrad = VectorXd::Zero(npars);
    MatrixXd fhessian = MatrixXd::Zero(npars,npars);
    function->calcValGradHessian(x, fval, fgrad, fhessian);
    
    val = 0.;
    grad = VectorXd::Zero(npars);
    
    for(unsigned int i=0;i<npars;++i)
    {
      val += fgrad(i)*fgrad(i);
    }
    
    for(unsigned int j=0;j<npars;++j)
    {
      for(unsigned int i=0;i<npars;++i)
      {
        grad(j) += fgrad(i)*fhessian(i,j);
      }
      grad(j) *= 2.;
    }
    
//     if(verbose == true)
//     {
//       for(unsigned int i=0;i<npars;++i)
//       {
//         VectorXd xtemp = x;
//         xtemp(i) += 1.0e-6;
//         
//         double fval=0.;
//         VectorXd fgrad = VectorXd::Zero(npars);
//         double tempval = 0.;
//         function->calcValGrad(xtemp, fval, fgrad);
//         
//         for(unsigned int j=0;j<npars;++j)
//         {
//           tempval += fgrad(j)*fgrad(j);
//         }
//         
//         cout<<"grad test "<<i<<" "<<(tempval-val)*1.0e6<<" "<<grad(i)<<endl;
//       }
//       cout<<endl;
//     }
    
    
    
    
    
    
    
    
    
    
    return true;
  }
  
  
  bool SquareGradient::calcValGradHessian(const VectorXd& x, double& val, VectorXd& grad, MatrixXd& hessian)
  {
//     verbose = false;
//     
//     this->calcValGrad(x, val, grad);
//     
//     hessian = MatrixXd::Zero(npars,npars);
//     for(unsigned int j=0;j<npars;++j)
//     {
//       VectorXd xtemp = x;
//       xtemp(j) += 1.0e-8;
//       VectorXd gradtemp = VectorXd::Zero(npars);
//       double valtemp=0.;
//       
//       this->calcValGrad(xtemp, valtemp, gradtemp);
//       
//       for(unsigned int i=0;i<npars;++i)
//       {
//         hessian(i,j) = (gradtemp(i) - grad(i))/(1.0e-8);
//       }
//     }
//     
//     verbose = true;
//     
//     return true;
    
    
    
    
    
    
    
    
    
    
    
    return calcValGrad(x,val,grad);
  }
  
  
  
}
