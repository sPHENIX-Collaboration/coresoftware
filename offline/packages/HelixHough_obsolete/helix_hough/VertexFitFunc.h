#ifndef HELIXHOUGH_VERTEXFITFUNC_H
#define HELIXHOUGH_VERTEXFITFUNC_H

#include "FunctionGradHessian.h"

#include <Eigen/Core>

#include <vector>

class SimpleTrack3D;

class HelixDCAFunc : public FitNewton::FunctionGradHessian
{
  public:
    HelixDCAFunc();
    ~HelixDCAFunc();
    
    bool calcValGradHessian(const Eigen::VectorXd& x, double& val, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian);
    
    FitNewton::FunctionGradHessian* Clone() const {return new HelixDCAFunc();}
    
    double getTangent(unsigned int coor){return tangent[coor];}
    double getPoint(unsigned int coor){return point[coor];}
    Eigen::Matrix<float,3,3> getPointCovariance(){return point_covariance;}
    void setCovariance(Eigen::Matrix<float,5,5> const& cov){covariance = cov;}
    
  private:
    Eigen::Matrix<float,5,5> covariance;
    std::vector<double> tangent;
    std::vector<double> point;
    Eigen::Matrix<float,3,3> point_covariance;
};

class VertexFitFunc : public FitNewton::FunctionGradHessian
{
  public:
    VertexFitFunc();
    ~VertexFitFunc();
    
    bool calcValGradHessian(const Eigen::VectorXd& x, double& val, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian);
    
    FitNewton::FunctionGradHessian* Clone() const {return new VertexFitFunc();}
    
    void setTracks(std::vector<SimpleTrack3D>* trks){tracks = trks;}
    void setCovariances(std::vector<Eigen::Matrix<float,5,5> >* covs){covariances = covs;}
    
  private:
    std::vector<Eigen::Matrix<float,5,5> >* covariances;
    std::vector<SimpleTrack3D> *tracks;
};


#endif
