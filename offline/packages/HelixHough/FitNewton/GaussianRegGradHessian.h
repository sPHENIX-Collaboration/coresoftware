#ifndef FITNEWTON_GAUSSIANREGGRADHESSIAN_H
#define FITNEWTON_GAUSSIANREGGRADHESSIAN_H

#include "FunctionGradHessian.h"
#include "Seamstress.h"

#include <Eigen/Core>

#include <cstddef>
#include <pthread.h>
#include <vector>

namespace SeamStress { template <class TClass> class Pincushion; }

namespace FitNewton
{
  class GaussianRegGradHessian : public FunctionGradHessian
  {
    public:
      GaussianRegGradHessian(FunctionGradHessian* func_instance, double sig, unsigned long int numthreads=1);
      
      
      ~GaussianRegGradHessian();
      
//       void setPoints(const std::vector<std::vector<double> >& POINTS){points = &POINTS;}
//       void setData(const std::vector<double>& DATA){data = &DATA;}
//       void setErrors(const std::vector<double>& ERROR){data_errors = &ERROR;data_has_errors=true;}
      
      void setPoints(const std::vector<std::vector<double> >& POINTS);
      void setData(const std::vector<double>& DATA);
      void setErrors(const std::vector<double>& ERROR);
      
      void setPointsThread1(void* arg);
      void setErrorsThread1(void* arg);
      void setDataThread1(void* arg);
      
      
      void setFunction(FunctionGradHessian* f){func = f;}
      
      void resetPointsDataErrors(){points=NULL;data=NULL;data_errors=NULL;data_has_errors=false;}
      
      bool calcValGradHessian(const Eigen::VectorXd& x, double& val, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian);
      void calcValGradHessianThread1(void* arg);
      
      void getCovariance(Eigen::MatrixXd& cov);
      
      FunctionGradHessian* Clone() const;
      
      void setWidth(double sig){sigma = sig;}
      
      //set if errors are really just inverse weights of the data points
      void setErrorsAreWeights(bool e_a_w){errors_are_weights = e_a_w;}
      
      void computeCovariance(const double& val, const Eigen::MatrixXd& hessian);
      
      
    private:
      FunctionGradHessian* func;
      const std::vector<std::vector<double> >* points;
      const std::vector<double>* data;
      const std::vector<double>* data_errors;
      
      std::vector<std::vector<std::vector<double> >* > thread_points;
      std::vector<std::vector<double>* > thread_data;
      std::vector<std::vector<double>* > thread_data_errors;
      
      const Eigen::VectorXd* current_eval;
      
      double* val_output;
      Eigen::VectorXd* grad_output;
      Eigen::MatrixXd* hessian_output;
      
      Eigen::MatrixXd covariance;
      
      bool data_has_errors;
      bool errors_are_weights;
      
      std::vector<SeamStress::Seamstress*> *vssp;
      std::vector<SeamStress::Seamstress> vss;
      SeamStress::Pincushion<GaussianRegGradHessian> *pins;
      unsigned long int nthreads;
      
      unsigned long int thread_tot;
      unsigned long int niter;
      pthread_mutex_t mutex;
      pthread_mutexattr_t mattr;
      
      double sigma;
  };
}



#endif

