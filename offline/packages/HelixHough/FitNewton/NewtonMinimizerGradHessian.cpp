#include "NewtonMinimizerGradHessian.h"

#include "FunctionGradHessian.h"
#include "SquareGradient.h"

#include <Eigen/LU>

#include <algorithm>                 // for fill
#include <cmath>
#include <iostream>
#include <utility>

using namespace std;
using namespace Eigen;

namespace FitNewton
{
  NewtonMinimizerGradHessian::NewtonMinimizerGradHessian()
  {
    
  }


  NewtonMinimizerGradHessian::~NewtonMinimizerGradHessian()
  {
    
  }
  
  
  void NewtonMinimizerGradHessian::setFunction(FunctionGradHessian* func)
  {
    function=func;
    fixparameter.clear();
    fixparameter.assign(func->nPars(), 0);
  }
  
  
  void NewtonMinimizerGradHessian::fixParameter(unsigned int par)
  {
    if(par < fixparameter.size())
    {
      fixparameter[par] = 1;
    }
  }
  
  
  void NewtonMinimizerGradHessian::unfixParameter(unsigned int par)
  {
    if(par < fixparameter.size())
    {
      fixparameter[par] = 0;
    }
  }
  
  
  bool NewtonMinimizerGradHessian::zoom(const double& wolfe1, const double& wolfe2, double& lo, double& hi, VectorXd& try_grad, VectorXd& direction, double& grad0_dir, double& val0, double& val_lo, VectorXd& init_params, VectorXd& try_params, unsigned int max_iter, double& result)
  {
    double tryval = val0;
    double alpha = tryval;
    double temp1 = tryval;
    double temp2 = tryval;
    double temp3 = tryval;
    bool bounds = true;
    
    unsigned int counter = 1;
    while(true)
    {
      alpha = lo + hi;
      alpha*=0.5;
      try_params = direction;
      try_params *= alpha;
      try_params += init_params;
      bounds = function->calcValGrad(try_params, tryval, try_grad);
      for(unsigned int i=0;i<fixparameter.size();++i)
      {
        if(fixparameter[i] != 0)
        {
          try_grad[i] = 0.;
        }
      }
      temp1 = wolfe1*alpha;
      temp1 *= grad0_dir;
      temp1 += val0;
      
      if( ( tryval > temp1 ) || ( tryval >= val_lo ) || (bounds == false) )
      {
//         if( (fabs((tryval - val_lo)/(tryval)) < 1.0e-4) && (tryval < val_lo) ){result = alpha;return true;}
        if( (fabs((tryval - val_lo)/(tryval)) < 1.0e-4) ){result = alpha;return true;}
        hi = alpha;
      }
      else
      {
        temp1 = try_grad.dot(direction);
        temp2 = -wolfe2*grad0_dir;
        temp3 = fabs(temp1);
        
        if( temp3 <= fabs(temp2) )
        {
          result = alpha;
          return bounds;
        }
        temp3 = hi - lo;
        temp1 *= temp3;
        if( temp1 >= 0.)
        {
          hi = lo;
        }
        lo = alpha;
        val_lo = tryval;
      }
      counter++;
      if(counter > max_iter){return false;}
    }
  }
  
  
  bool NewtonMinimizerGradHessian::lineSearch(double& alpha, const double& wolfe1, const double& wolfe2, VectorXd& try_grad, VectorXd& direction, double& grad0_dir, double& val0, VectorXd& init_params, VectorXd& try_params, const double& precision, const double& accuracy, unsigned int max_iter, double& result)
  {
    double tryval = val0;
    double prev_val = tryval;
    double prev_alpha = tryval;
    double lo = tryval;
    double hi = tryval;
    double temp1 = tryval;
    double temp2 = tryval;
    double temp3 = tryval;
    unsigned int i = 1;
    bool bounds = true;
    
    prev_alpha = 0.;
    while(true)
    {
      try_params = direction;
      try_params *= alpha;
      try_params += init_params;
      
      bounds = function->calcValGrad(try_params, tryval, try_grad);
      for(unsigned int j=0;j<fixparameter.size();++j)
      {
        if(fixparameter[j] != 0)
        {
          try_grad[j] = 0.;
        }
      }
      if(bounds == false)
      {
        while(true)
        {
          alpha *= 0.5;
          try_params = direction;
          try_params *= alpha;
          try_params += init_params;
          
          bounds = function->calcValGrad(try_params, tryval, try_grad);
          for(unsigned int j=0;j<fixparameter.size();++j)
          {
            if(fixparameter[j] != 0)
            {
              try_grad[j] = 0.;
            }
          }
          if(bounds==true)
          {
            if(tryval < val0)
            {
              result=alpha;
              return true;
            }
          }
          if(i>max_iter){return false;}
          i++;
        }
        
        
        
        alpha = 0.5*(prev_alpha + alpha);
        i++;
        if(i > max_iter){return false;}
        continue;
      }
      
      temp1 = wolfe1*alpha;
      temp1 *= grad0_dir;
      temp1 += val0;
      if( ( tryval > temp1 ) || ( ( tryval > prev_val ) && (i>1) ))
      {
        lo = prev_alpha;
        hi = alpha;
        return zoom(wolfe1, wolfe2, lo, hi, try_grad, direction, grad0_dir, val0, prev_val, init_params, try_params, max_iter, result);
      }
      temp1 = try_grad.dot(direction);
      temp2 = -wolfe2*grad0_dir;
      temp3 = fabs(temp1);
      
      if( temp3 <= fabs(temp2) )
      {
        result = alpha;
        return bounds;
      }
      if( temp1 >= 0. )
      {
        lo = alpha;
        hi = prev_alpha;
        return zoom(wolfe1, wolfe2, lo, hi, try_grad, direction, grad0_dir, val0, tryval, init_params, try_params, max_iter, result);
      }
      prev_val = tryval;
      prev_alpha = alpha;
      alpha *= 2.;
      i++;
      if(i > max_iter){return false;}
    }
  }
  
  
  bool NewtonMinimizerGradHessian::findSaddlePoint(const VectorXd& start_point, VectorXd& min_point, double tol, unsigned int max_iter, double abs_tol)
  {
    try
    {
      if(!function){throw (char*)("minimize called but function has not been set");}
    }
    catch(char* str)
    {
      cout<<"Exception from NewtonMinimizerGradHessian: "<<str<<endl;
      throw;
      return false;
    }
    
    unsigned int npars = function->nPars();
    
    try
    {
      if(!(npars>0)){throw (char*)("function to be minimized has zero dimensions");}
      if(start_point.rows()!=(int)npars){throw (char*)("input to minimizer not of dimension required by function to be minimized");}
    }
    catch(char* str)
    {
      cout<<"Exception from NewtonMinimizerGradHessian: "<<str<<endl;
      throw;
      return false;
    }
    
    //parameters used for the Wolfe conditions
    double c1 = 1.0e-6;
    double c2 = 1.0e-1;
    
    VectorXd working_points[2];
    working_points[0] = VectorXd::Zero(npars);
    working_points[1] = VectorXd::Zero(npars);
    
    VectorXd* current_point = &(working_points[0]);
    VectorXd* try_point = &(working_points[1]);
    
    (*current_point) = start_point;
    
    double value=0.;
    double prev_value=0.;
    
    double dir_der=0.;
    double scale = 1.;
    double scale_temp = 1.;
    double grad0_dir = 0.;
    
    bool good_value = true;
    bool bounds = true;
    
    unsigned int search_iter = 64;
    
    VectorXd grad[2];
    grad[0] = VectorXd::Zero(npars);
    grad[1] = VectorXd::Zero(npars);
    VectorXd* current_grad = &(grad[0]);
    VectorXd* try_grad = &(grad[1]);
    
    VectorXd newgrad = VectorXd::Zero(npars);
    
    MatrixXd hessian = MatrixXd::Zero(npars, npars);
    
    VectorXd move = VectorXd::Zero(npars);
    VectorXd unit_move = VectorXd::Zero(npars);
    
    FunctionGradHessian* orig_func = function;
    SquareGradient gradsquared(orig_func);
    function = &gradsquared;
    
    //try a Newton iteration
    orig_func->calcValGradHessian((*current_point), value, (*current_grad), hessian);
    move = -hessian.fullPivLu().solve(*current_grad);
    gradsquared.calcValGrad((*current_point), value, newgrad);
    good_value=true;
    for(unsigned int i=0;i<npars;++i){if(!(move(i) == move(i))){good_value=false;break;}}
    if(good_value == false){move = -newgrad;}
    dir_der = newgrad.dot(move);
    if(dir_der>0.){move = -move;}
    gradsquared.rescaleMove((*current_point), move);
    grad0_dir = newgrad.dot(move);
    scale_temp = 1.;
    //find scale, such that move*scale satisfies the strong Wolfe conditions
    bounds = lineSearch(scale_temp, c1, c2, (*try_grad), move, grad0_dir, value, (*current_point), (*try_point), tol, abs_tol, search_iter, scale);
    if(bounds == false){min_point = start_point; function=orig_func;return false;}
    move *= scale;
    (*try_point) = ((*current_point) + move);
    orig_func->calcValGradHessian((*try_point), prev_value, (*try_grad), hessian);
    gradsquared.calcValGrad((*try_point), prev_value, newgrad);
    swap(current_point, try_point);
    swap(current_grad, try_grad);
    swap(value, prev_value);
    
    unsigned long int count = 1;
    bool converged=false;
    while(converged==false)
    {
      if((fabs((prev_value - value)/prev_value)<tol || fabs(prev_value - value)<abs_tol)){converged=true;break;}
      prev_value = value;
      //try a Newton iteration
      move = -hessian.fullPivLu().solve(*current_grad);
      gradsquared.calcValGrad((*current_point), value, newgrad);
      
      good_value=true;
      for(unsigned int i=0;i<npars;++i){if(!(move(i) == move(i))){good_value=false;break;}}
      scale_temp = 1.;
      scale_temp = fabs(value/newgrad.dot(move));
      if(good_value == false){move = -newgrad;scale_temp = fabs(value/move.dot(move));}
      dir_der = newgrad.dot(move);
      if(dir_der>0.){move = -move;}
      gradsquared.rescaleMove((*current_point), move);
      grad0_dir = newgrad.dot(move);
      //find scale, such that move*scale satisfies the strong Wolfe conditions
      bounds = lineSearch(scale_temp, c1, c2, (*try_grad), move, grad0_dir, value, (*current_point), (*try_point), tol, abs_tol, search_iter, scale);
      if(bounds == false){min_point = (*current_point); function=orig_func;return false;}
      move *= scale;
      (*try_point) = ((*current_point) + move);
      orig_func->calcValGradHessian((*try_point), value, (*try_grad), hessian);
      gradsquared.calcValGrad((*try_point), value, newgrad);
      swap(current_point, try_point);
      swap(current_grad, try_grad);
      count++;
      if(count > max_iter){break;}
    }
    orig_func->computeCovariance(value, hessian);
    min_point = (*current_point);
    function=orig_func;return converged;
  }
  
  
  bool NewtonMinimizerGradHessian::minimize(const VectorXd& start_point, VectorXd& min_point, double tol, unsigned int max_iter, double abs_tol)
  {
    try
    {
      if(!function){throw (char*)("minimize called but function has not been set");}
    }
    catch(char* str)
    {
      cout<<"Exception from NewtonMinimizerGradHessian: "<<str<<endl;
      throw;
      return false;
    }
    
    unsigned int npars = function->nPars();
    
    try
    {
      if(!(npars>0)){throw (char*)("function to be minimized has zero dimensions");}
      if(start_point.rows()!=(int)npars){throw (char*)("input to minimizer not of dimension required by function to be minimized");}
    }
    catch(char* str)
    {
      cout<<"Exception from NewtonMinimizerGradHessian: "<<str<<endl;
      throw;
      return false;
    }
    
    
    //parameters used for the Wolfe conditions
    double c1 = 1.0e-4;
    double c2 = 0.9;
    
    
    VectorXd working_points[2];
    working_points[0] = VectorXd::Zero(npars);
    working_points[1] = VectorXd::Zero(npars);
    
    VectorXd* current_point = &(working_points[0]);
    VectorXd* try_point = &(working_points[1]);
    
    (*current_point) = start_point;
    
    double value=0.;
    double prev_value=0.;
    
    double dir_der=0.;
    double scale = 1.;
    double scale_temp = 1.;
    double grad0_dir = 0.;
    
    bool good_value = true;
    bool bounds = true;
    
    unsigned int search_iter = 32;
    
    VectorXd grad[2];
    grad[0] = VectorXd::Zero(npars);
    grad[1] = VectorXd::Zero(npars);
    VectorXd* current_grad = &(grad[0]);
    VectorXd* try_grad = &(grad[1]);
    
    
    MatrixXd hessian = MatrixXd::Zero(npars, npars);
    
    VectorXd move = VectorXd::Zero(npars);
    VectorXd unit_move = VectorXd::Zero(npars);
    
    //try a Newton iteration
    function->calcValGradHessian((*current_point), value, (*current_grad), hessian);
    for(unsigned int i=0;i<fixparameter.size();++i)
    {
      if(fixparameter[i] != 0)
      {
        (*current_grad)[i] = 0.;
        for(unsigned int j=0;j<npars;++j)
        {
          hessian(j,i) = 0.;
          hessian(i,j) = 0.;
        }
        hessian(i,i) = 1.;
      }
    }
    
    move = -hessian.fullPivLu().solve(*current_grad);
    good_value=true;
    for(unsigned int i=0;i<npars;++i){if(!(move(i) == move(i))){good_value=false;break;}}
    if(good_value == false){move = - (*current_grad);}
    dir_der = (*current_grad).dot(move);
    //if the inverse hessian times the negative gradient isn't even a descent direction, negate the direction
    if(dir_der>0.)
    {
      move = -move;
    }
    function->rescaleMove((*current_point), move);
    grad0_dir = (*current_grad).dot(move);
    scale_temp = 1.;
    //find scale, such that move*scale satisfies the strong Wolfe conditions
    bounds = lineSearch(scale_temp, c1, c2, (*try_grad), move, grad0_dir, value, (*current_point), (*try_point), tol, abs_tol, search_iter, scale);
    if(bounds == false){min_point = start_point; return false;}
    move *= scale;
    (*try_point) = ((*current_point) + move);
    function->calcValGradHessian((*try_point), prev_value, (*try_grad), hessian);
    for(unsigned int i=0;i<fixparameter.size();++i)
    {
      if(fixparameter[i] != 0)
      {
        (*try_grad)[i] = 0.;
        for(unsigned int j=0;j<npars;++j)
        {
          hessian(j,i) = 0.;
          hessian(i,j) = 0.;
        }
        hessian(i,i) = 1.;
      }
    }
    swap(current_point, try_point);
    swap(current_grad, try_grad);
    swap(value, prev_value);
    
    unsigned long int count = 1;
    bool converged=false;
    while(converged==false)
    {
      if((fabs((prev_value - value)/prev_value)<tol || fabs(prev_value - value)<abs_tol)){converged=true;break;}
      prev_value = value;
      //try a Newton iteration
      
      move = -hessian.fullPivLu().solve(*current_grad);
      good_value=true;
      for(unsigned int i=0;i<npars;++i){if(!(move(i) == move(i))){good_value=false;break;}}
      if(good_value == false){move = - (*current_grad);}
      dir_der = (*current_grad).dot(move);
      //if the inverse hessian times the negative gradient isn't even a descent direction, negate the direction
      if(dir_der>0.)
      {
        move = -move;
      }
      function->rescaleMove((*current_point), move);
      grad0_dir = (*current_grad).dot(move);
      scale_temp = 1.;
      //find scale, such that move*scale satisfies the strong Wolfe conditions
//       scale_temp = fabs(value/grad0_dir);
      bounds = lineSearch(scale_temp, c1, c2, (*try_grad), move, grad0_dir, value, (*current_point), (*try_point), tol, abs_tol, search_iter, scale);
      if(bounds == false){min_point = (*current_point); return false;}
      move *= scale;
      (*try_point) = ((*current_point) + move);
      function->calcValGradHessian((*try_point), value, (*try_grad), hessian);
      for(unsigned int i=0;i<fixparameter.size();++i)
      {
        if(fixparameter[i] != 0)
        {
          (*try_grad)[i] = 0.;
          for(unsigned int j=0;j<npars;++j)
          {
            hessian(j,i) = 0.;
            hessian(i,j) = 0.;
          }
          hessian(i,i) = 1.;
        }
      }
      swap(current_point, try_point);
      swap(current_grad, try_grad);
      
      count++;
      if(count > max_iter){break;}
    }
    function->computeCovariance(value, hessian);
    
    min_point = (*current_point);
    return converged;
  }
}


