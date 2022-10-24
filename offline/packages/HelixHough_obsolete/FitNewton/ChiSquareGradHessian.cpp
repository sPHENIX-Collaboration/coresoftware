#include "ChiSquareGradHessian.h"

#include "Pincushion.h"

#include <Eigen/LU>

#include <memory>

using namespace std;
using namespace Eigen;
using namespace SeamStress;

namespace FitNewton
{
  ChiSquareGradHessian::ChiSquareGradHessian(FunctionGradHessian* func_instance, unsigned long int numthreads) : FunctionGradHessian(func_instance->nPars(), 0), func(func_instance), data_has_errors(false), errors_are_weights(false)
  {
    Seamstress::init_vector(numthreads, vss);
    
    vssp = new vector<Seamstress*>();
    for(unsigned long int i=0;i<vss.size();i++){vssp->push_back(&(vss[i]));}
    
    pins = new Pincushion<ChiSquareGradHessian>(this, vssp);
    
    nthreads = numthreads;
    
    pthread_mutexattr_init(&mattr);
    pthread_mutexattr_settype(&mattr, PTHREAD_MUTEX_NORMAL);
    pthread_mutex_init(&mutex, &mattr);
  }


  ChiSquareGradHessian::~ChiSquareGradHessian()
  {
    for(unsigned int i=0;i<nthreads;i++)
    {
      vss[i].stop();
    }
    delete pins;
    delete vssp;
  }
  
  
  FunctionGradHessian* ChiSquareGradHessian::Clone() const
  {
    ChiSquareGradHessian* clone = new ChiSquareGradHessian(func, nthreads);
    clone->setData(*data);
    clone->setPoints(*points);
    if(data_has_errors==true){clone->setErrors(*data_errors);}
    return clone;
  }
  
  
  void ChiSquareGradHessian::computeCovariance(const double& val, const MatrixXd& hessian)
  {
    covariance = hessian.inverse();
//     hessian.computeInverse(&covariance);
    if(data_has_errors==false){covariance *= (val/( (double)(points->size() - npars) ) );}
    else if(errors_are_weights==true)
    {
      double scale = 0.;
      for(unsigned int i=0;i<data_errors->size();i++)
      {
        double temp = 1./((*data_errors)[i]);temp*=temp;
        scale += temp;
      }
      covariance *= val/(scale - (double)npars);
    }
  }
  
  
  void ChiSquareGradHessian::getCovariance(MatrixXd& cov)
  {
    cov = covariance;
  }


  bool ChiSquareGradHessian::calcValGradHessian(const VectorXd& x, double& val, VectorXd& grad, MatrixXd& hessian)
  {
    current_eval = &x;
    bool bounds = true;
    FunctionGradHessian* func_instance = func->Clone();
    for(unsigned int j=0;j<(*points)[0].size();j++){func_instance->setFixedPar(j, (*points)[0][j]);}
    bounds = func_instance->calcValGrad((*current_eval), val, grad);
    delete func_instance;
    
    
    val=0.;
    grad = VectorXd::Zero(npars);
    hessian = MatrixXd::Zero(npars, npars);
    
    val_output = &val;
    grad_output = &grad;
    hessian_output = &hessian;
    
    covariance = MatrixXd::Zero(npars, npars);
    
    niter = points->size();
    if(nthreads>niter){thread_tot=niter;}
    else{thread_tot=nthreads;}
    pins->sewStraight(&ChiSquareGradHessian::calcValGradHessianThread1, thread_tot);
    return bounds;
  }


  void ChiSquareGradHessian::calcValGradHessianThread1(void* arg)
  {
    unsigned long int w = (*((unsigned long int *)arg));
    unsigned long int part = niter/thread_tot;
    unsigned long int rem = niter - part*thread_tot;
    unsigned long int start = 0;
    unsigned long int end = 0;
    if(w<rem)
    {
      start = (part+1)*w;
      end = start + part;
    }
    else
    {
      start = (part+1)*rem;
      start += (part*(w-rem));
      end = start + part - 1;
    }
    
    //initialize the val, grad, and hessian to be added to those of the main thread
    double temp_val=0.;
    VectorXd temp_grad = VectorXd::Zero(npars);
    MatrixXd temp_hessian = MatrixXd::Zero(npars, npars);
    
    //initialize the val, grad, and hessian of the function with fixed parameters at the given points vector
    double temp_val2=0.;
    VectorXd temp_grad2 = VectorXd::Zero(npars);
    MatrixXd temp_hessian2 = MatrixXd::Zero(npars, npars);
    
    VectorXd cur_grad = VectorXd::Zero(npars);
    
    //f += (t - v)^2
    //df/dx += 2*[ t - v ]*(dt/dx)
    //d^2f/dxdy += 2*[(dt/dy)*(dt/dx) + [t-v]*(d^2t/dxdy) ]
    
    FunctionGradHessian* func_instance = func->Clone();
    
    for(unsigned long int i=(start);i<=(end);i++)
    {
      for(unsigned int j=0;j<(*points)[i].size();j++){func_instance->setFixedPar(j, (*points)[i][j]);}
      func_instance->calcValGradHessian((*current_eval), temp_val2, temp_grad2, temp_hessian2);
      double inv_variance = 1.;
      if(data_has_errors==true){inv_variance = 1./(*data_errors)[i];inv_variance*=inv_variance;}
      
      double temp1 = (temp_val2 - (*data)[i]);
      temp_val += temp1*temp1*inv_variance;
      
      for(unsigned int j=0;j<npars;j++)
      {
        cur_grad(j) = 2.*temp1*temp_grad2(j)*inv_variance;
        temp_grad(j) += cur_grad(j);
      }
      
      for(unsigned int k=0;k<npars;k++)
      {
        for(unsigned int j=0;j<npars;j++)
        {
          temp_hessian(j,k) += 2.*inv_variance*(temp_grad2(j)*temp_grad2(k) + temp1*temp_hessian2(j,k));
        }
      }
    }
    
    
    pthread_mutex_lock(&mutex);
    (*val_output) += temp_val;
    (*grad_output) += temp_grad;
    (*hessian_output) += temp_hessian;
    pthread_mutex_unlock(&mutex);
    
    
    
    delete func_instance;
  }

  
  bool ChiSquareGradHessian::calcValGrad(const VectorXd& x, double& val, VectorXd& grad)
  {
    current_eval = &x;
    bool bounds = true;
    FunctionGradHessian* func_instance = func->Clone();
    for(unsigned int j=0;j<(*points)[0].size();j++){func_instance->setFixedPar(j, (*points)[0][j]);}
    bounds = func_instance->calcValGrad((*current_eval), val, grad);
    delete func_instance;
    
    val=0.;
    grad = VectorXd::Zero(npars);
    
    val_output = &val;
    grad_output = &grad;
    
    niter = points->size();
    if(nthreads>niter){thread_tot=niter;}
    else{thread_tot=nthreads;}
    pins->sewStraight(&ChiSquareGradHessian::calcValGradThread1, thread_tot);
    return bounds;
  }
  
  
  void ChiSquareGradHessian::calcValGradThread1(void* arg)
  {
    unsigned long int w = (*((unsigned long int *)arg));
    unsigned long int part = niter/thread_tot;
    unsigned long int rem = niter - part*thread_tot;
    unsigned long int start = 0;
    unsigned long int end = 0;
    if(w<rem)
    {
      start = (part+1)*w;
      end = start + part;
    }
    else
    {
      start = (part+1)*rem;
      start += (part*(w-rem));
      end = start + part - 1;
    }
    
    //initialize the val, grad, and hessian to be added to those of the main thread
    double temp_val=0.;
    VectorXd temp_grad = VectorXd::Zero(npars);
    
    //initialize the val, grad, and hessian of the function with fixed parameters at the given points vector
    double temp_val2=0.;
    VectorXd temp_grad2 = VectorXd::Zero(npars);
    
    VectorXd cur_grad = VectorXd::Zero(npars);
    
    //f += (t - v)^2
    //df/dx += 2*[ t - v ]*(dt/dx)
    //d^2f/dxdy += 2*[(dt/dy)*(dt/dx) + [t-v]*(d^2t/dxdy) ]
    
    FunctionGradHessian* func_instance = func->Clone();
    
    for(unsigned long int i=(start);i<=(end);i++)
    {
      for(unsigned int j=0;j<(*points)[i].size();j++){func_instance->setFixedPar(j, (*points)[i][j]);}
      func_instance->calcValGrad((*current_eval), temp_val2, temp_grad2);
      double inv_variance = 1.;
      if(data_has_errors==true){inv_variance = 1./(*data_errors)[i];inv_variance*=inv_variance;}
      
      double temp1 = (temp_val2 - (*data)[i]);
      temp_val += temp1*temp1*inv_variance;
      
      for(unsigned int j=0;j<npars;j++)
      {
        cur_grad(j) = 2.*temp1*temp_grad2(j)*inv_variance;
        temp_grad(j) += cur_grad(j);
      }
    }
    
    
    pthread_mutex_lock(&mutex);
    (*val_output) += temp_val;
    (*grad_output) += temp_grad;
    pthread_mutex_unlock(&mutex);
    
    
    
    delete func_instance;
  }
  
}


