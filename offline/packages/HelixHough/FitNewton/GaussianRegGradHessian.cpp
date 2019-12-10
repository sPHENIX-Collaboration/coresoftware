#include "GaussianRegGradHessian.h"
#include "Pincushion.h"

#include <Eigen/LU>

#include <cmath>
#include <memory>

using namespace std;
using namespace Eigen;
using namespace SeamStress;

namespace FitNewton
{
  GaussianRegGradHessian::GaussianRegGradHessian(FunctionGradHessian* func_instance, double sig, unsigned long int numthreads) : FunctionGradHessian(func_instance->nPars(), 0), func(func_instance), data_has_errors(false), errors_are_weights(false), sigma(sig)
  {
    Seamstress::init_vector(numthreads, vss);
    
    vssp = new vector<Seamstress*>();
    for(unsigned long int i=0;i<vss.size();i++){vssp->push_back(&(vss[i]));}
    
    pins = new Pincushion<GaussianRegGradHessian>(this, vssp);
    
    nthreads = numthreads;
    
    pthread_mutexattr_init(&mattr);
    pthread_mutexattr_settype(&mattr, PTHREAD_MUTEX_NORMAL);
    pthread_mutex_init(&mutex, &mattr);
    
    vector<vector<double> >* tempvec = NULL;
    thread_points.assign(nthreads, tempvec);
    vector<double>* tempvec2 = NULL;
    thread_data_errors.assign(nthreads, tempvec2);
    thread_data.assign(nthreads, tempvec2);
    
    resetPointsDataErrors();
  }
  
  
  GaussianRegGradHessian::~GaussianRegGradHessian()
  {
    for(unsigned int i=0;i<nthreads;i++)
    {
      vss[i].stop();
    }
    delete pins;
    delete vssp;
    
    for(unsigned long int i=0;i<nthreads;i++)
    {
      if(thread_points[i] != NULL){delete thread_points[i];}
      if(thread_data_errors[i] != NULL){delete thread_data_errors[i];}
      if(thread_data[i] != NULL){delete thread_data[i];}
    }
  }
  
  
  void GaussianRegGradHessian::setPoints(const vector<vector<double> >& POINTS)
  {
    if(POINTS.size() == 0){return;}
    points = &POINTS;
    niter = points->size();
    if(nthreads>niter){thread_tot=niter;}
    else{thread_tot=nthreads;}
    pins->sewStraight(&GaussianRegGradHessian::setPointsThread1, nthreads);
  }
  
  
  void GaussianRegGradHessian::setPointsThread1(void* arg)
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
    if(thread_points[w] != NULL){delete thread_points[w];}
    thread_points[w] = new vector<vector<double> >();
    vector<double> tempvec;
    tempvec.assign((*points)[0].size(), 0.);
    for(unsigned long int i=(start);i<=(end);i++)
    {
      thread_points[w]->push_back(tempvec);
      for(unsigned long int j=0;j<thread_points[w]->back().size();j++)
      {
        thread_points[w]->back()[j] = (*points)[i][j];
      }
    }
  }
  
  
  void GaussianRegGradHessian::setErrors(const vector<double>& ERROR)
  {
    if(ERROR.size() == 0){return;}
    data_errors = &ERROR;data_has_errors=true;
    niter = data_errors->size();
    if(nthreads>niter){thread_tot=niter;}
    else{thread_tot=nthreads;}
    pins->sewStraight(&GaussianRegGradHessian::setErrorsThread1, nthreads);
  }
  
  
  void GaussianRegGradHessian::setErrorsThread1(void* arg)
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
    
    if(thread_data_errors[w] != NULL){delete thread_data_errors[w];}
    thread_data_errors[w] = new vector<double>();
    for(unsigned long int i=(start);i<=(end);i++)
    {
      thread_data_errors[w]->push_back((*data_errors)[i]);
    }
  }
  
  
  void GaussianRegGradHessian::setData(const vector<double>& DATA)
  {
    if(DATA.size() == 0){return;}
    data = &DATA;
    niter = data->size();
    if(nthreads>niter){thread_tot=niter;}
    else{thread_tot=nthreads;}
    pins->sewStraight(&GaussianRegGradHessian::setDataThread1, nthreads);
  }
  
  
  void GaussianRegGradHessian::setDataThread1(void* arg)
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
    
    if(thread_data[w] != NULL){delete thread_data[w];}
    thread_data[w] = new vector<double>();
    for(unsigned long int i=(start);i<=(end);i++)
    {
      thread_data[w]->push_back((*data)[i]);
    }
  }
  
  
  FunctionGradHessian* GaussianRegGradHessian::Clone() const
  {
    GaussianRegGradHessian* clone = new GaussianRegGradHessian(func, sigma, nthreads);
    clone->setData(*data);
    clone->setPoints(*points);
    if(data_has_errors==true){clone->setErrors(*data_errors);}
    return clone;
  }
  
  
  void GaussianRegGradHessian::computeCovariance(const double& val, const MatrixXd& hessian)
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
  
  
  void GaussianRegGradHessian::getCovariance(MatrixXd& cov)
  {
    cov = covariance;
  }
  
  
  bool GaussianRegGradHessian::calcValGradHessian(const VectorXd& x, double& val, VectorXd& grad, MatrixXd& hessian)
  {
    current_eval = &x;
    
    val=0.;
    grad = VectorXd::Zero(npars);
    hessian = MatrixXd::Zero(npars, npars);
    
    val_output = &val;
    grad_output = &grad;
    hessian_output = &hessian;
    covariance = MatrixXd::Zero(npars, npars);
    if(points == NULL){return true;}
    niter = points->size();
    if(nthreads>niter){thread_tot=niter;}
    else{thread_tot=nthreads;}
    pins->sewStraight(&GaussianRegGradHessian::calcValGradHessianThread1, thread_tot);
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
//       scale/=( (double)(points->size() - npars) );
//       covariance *= scale;
      covariance *= val/(scale - (double)npars);
    }
    return true;
  }
  
  
  void GaussianRegGradHessian::calcValGradHessianThread1(void* arg)
  {
    unsigned long int w = (*((unsigned long int *)arg));
    //initialize the val, grad, and hessian to be added to those of the main thread
    double temp_val=0.;
    VectorXd temp_grad = VectorXd::Zero(npars);
    MatrixXd temp_hessian = MatrixXd::Zero(npars, npars);
    
    //initialize the val, grad, and hessian of the function with fixed parameters at the given points vector
    double temp_val2=0.;
    VectorXd temp_grad2 = VectorXd::Zero(npars);
    MatrixXd temp_hessian2 = MatrixXd::Zero(npars, npars);
    
    VectorXd cur_grad = VectorXd::Zero(npars);
    
    // f = -exp(-(t-v)^2/(2*s^2))
    // df/dx = -f*(1/s^2)*(t-v)*dt/dx
    // d^2f/dxdy = (-1/s^2)*[ (df/dy)*(t-v)*(dt/dx) + f*(dt/dy)*(dt/dx) + f*(t-v)*(d^2t/dxdy) ]
    
    FunctionGradHessian* func_instance = func->Clone();
    
    double inv_sig2 = 1./(sigma*sigma);
    
    
    for(unsigned long int i=0;i<thread_points[w]->size();i++)
    {
      for(unsigned int j=0;j<(*(thread_points[w]))[i].size();j++){func_instance->setFixedPar(j, (*(thread_points[w]))[i][j]);}
      func_instance->calcValGradHessian((*current_eval), temp_val2, temp_grad2, temp_hessian2);
      double inv_variance = 1.;
      if(data_has_errors==true){inv_variance = 1./(*(thread_data_errors[w]))[i];inv_variance*=inv_variance;}
      
      double temp1 = (temp_val2 - (*(thread_data[w]))[i]);
      double temp2 = -exp(-temp1*temp1*0.5*inv_sig2);
      
      temp_val += temp2*inv_variance;
      
      for(unsigned int j=0;j<npars;j++)
      {
        cur_grad(j) = -temp2*inv_sig2*temp1*temp_grad2(j)*inv_variance;
        temp_grad(j) += cur_grad(j);
      }
      
      for(unsigned int k=0;k<npars;k++)
      {
        for(unsigned int j=0;j<npars;j++)
        {
          temp_hessian(j, k) += -inv_sig2*( cur_grad(k)*temp1*temp_grad2(j) + temp2*temp_grad2(k)*temp_grad2(j) + temp2*temp1*temp_hessian2(j, k) )*inv_variance;
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
}






