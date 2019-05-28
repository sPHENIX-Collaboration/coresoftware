#ifndef FITNEWTON_CHISQUAREGRADHESSIAN_H
#define FITNEWTON_CHISQUAREGRADHESSIAN_H

#include "FunctionGradHessian.h"
#include "Seamstress.h"

#include <Eigen/Core>

#include <cstddef>
#include <pthread.h>
#include <vector>

namespace SeamStress { template <class TClass> class Pincushion; }

namespace FitNewton
{
  class ChiSquareGradHessian : public FunctionGradHessian
  {
    public:
      ChiSquareGradHessian(FunctionGradHessian* func_instance, unsigned long int numthreads=1);
      
      
      ~ChiSquareGradHessian();
      
      void setPoints(const std::vector<std::vector<double> >& POINTS){points = &POINTS;}
      void setData(const std::vector<double>& DATA){data = &DATA;}
      void setErrors(const std::vector<double>& ERROR){data_errors = &ERROR;data_has_errors=true;}
      
      void setFunction(FunctionGradHessian* f){func = f;}
      
      void resetPointsDataErrors(){points=NULL;data=NULL;data_errors=NULL;data_has_errors=false;}
      
      bool calcValGradHessian(const Eigen::VectorXd& x, double& val, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian);
      void calcValGradHessianThread1(void* arg);
      
      bool calcValGrad(const Eigen::VectorXd& x, double& val, Eigen::VectorXd& grad);
      void calcValGradThread1(void* arg);
      
      void getCovariance(Eigen::MatrixXd& cov);
      
      FunctionGradHessian* Clone() const;
      
      //set if errors are really just inverse weights of the data points
      void setErrorsAreWeights(bool e_a_w){errors_are_weights = e_a_w;}
      
      void computeCovariance(const double& val, const Eigen::MatrixXd& hessian);
      
      
    private:
      FunctionGradHessian* func;
      const std::vector<std::vector<double> >* points;
      const std::vector<double>* data;
      const std::vector<double>* data_errors;
      
      const Eigen::VectorXd* current_eval;
      
      double* val_output;
      Eigen::VectorXd* grad_output;
      Eigen::MatrixXd* hessian_output;
      
      Eigen::MatrixXd covariance;
      
      bool data_has_errors;
      bool errors_are_weights;
      
      std::vector<SeamStress::Seamstress*> *vssp;
      std::vector<SeamStress::Seamstress> vss;
      SeamStress::Pincushion<ChiSquareGradHessian> *pins;
      unsigned long int nthreads;
      
      unsigned long int thread_tot;
      unsigned long int niter;
      pthread_mutex_t mutex;
      pthread_mutexattr_t mattr;
  };
}



#endif

