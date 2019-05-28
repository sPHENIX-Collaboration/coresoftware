#include "VertexFitFunc.h"

#include "NewtonMinimizerGradHessian.h"
#include "SimpleTrack3D.h"

#include <cmath>
#include <memory>

using namespace std;
using namespace FitNewton;
using namespace Eigen;


HelixDCAFunc::HelixDCAFunc() : FunctionGradHessian(1, 8)
{
  tangent.assign(3,0.);
  point.assign(3,0.);
}


HelixDCAFunc::~HelixDCAFunc()
{
  
}


bool HelixDCAFunc::calcValGradHessian(const VectorXd& x, double& val, VectorXd& grad, MatrixXd& hessian)
{
//   the 5 helix parameters are k,d,phi,z0,dzdl
//   l will be used as the independent variable
//   l=0 corresponds to (d*cos(phi), d*sin(phi), z0)
//   fixedpars[5,6,7] correspond to vx,vy,vz of the point we are calculating the DCA to
//   
//   a change in l corresponds to what change in azimuth?
//   
//   l^2 = s^2 + z^2
//   
//   1 = s^2/l^2 + dzdl^2
//   
//   s^2 = (1 - dzdl^2)*l^2
//   
//   change in s gives what change in azimuth?
//   
//   s = theta*r = theta/k
//   
//   theta = s*k
//   
//   s = 2.*asin(0.5*k*D)/k;
//   
//   s*k*0.5 = asin(0.5*k*D)
//   
//   sin(s*k/2) = k*D/2
//   
//   D = (2/k)*sin(s*k/2)
//   
//   
//   in what direction do we move by D in the xy plane?
//   
//   theta/2 = asin(0.5*k*D)
//   
//   l == x
  
  
  point_covariance = Matrix<float,3,3>::Zero(3,3);
  
  
  float k = fixedpars[2];
  float d = fixedpars[1];
  float phi = fixedpars[0];
  float z0 = fixedpars[3];
  float cosphi = cos(phi);
  float sinphi = sin(phi);
  //determine starting point on helix xp,yp,zp
  float xp = d*cosphi;
  float yp = d*sinphi;
  float zp = z0;
  
  // pc = (df/dx)*C*(df/dx)^T
  Matrix<float,3,5> dpdx =  Matrix<float,3,5>::Zero(3,5);
  // dxp/dphi = -d*sinphi , dyp/dphi = d*cosphi
  dpdx(0,0) = -d*sinphi;dpdx(1,0) = d*cosphi;
  // dxp/dd = cosphi;dyp/dd = sinphi
  dpdx(0,1) = cosphi;dpdx(1,1) = sinphi;
  dpdx(2,3) = 1.;
  point_covariance = dpdx * covariance * (dpdx.transpose());
  
  
  //determine tangent direction ux,uy in xy plane
  float ux = sinphi;
  float uy = -cosphi;
  float duxdphi = -cosphi;
  float duydphi = sinphi;
  if(d <= 0.)
  {
    ux=-ux;
    uy=-uy;
    duxdphi=-duxdphi;
    duydphi=-duydphi;
  }
  
  float dzdl = fixedpars[4];
  //determine point on the helix with parameter l == x(0)
  float dsdl = sqrt(1. - dzdl*dzdl);
  float s = dsdl*x(0);
  float dsddzdl = -x(0) * (1./sqrt(1. - dzdl*dzdl)) * dzdl;
  
  float psi = 0.5*s*k;
  float D = 0.;
  Matrix<float,1,5> dDdx = Matrix<float,1,5>::Zero(1,5);
  
  if(psi > 0.1)
  {
    float sinpsi = sin(psi);
    float cospsi = cos(psi);
    D = (2./k)*sinpsi;
    // dsinpsi/dk = cos(psi)*0.5*s
    dDdx(0,2) = -2.*sinpsi/(k*k) + (cospsi/k)*s;
    // dsinpsi/ddzdl = cospsi * 0.5 * k * ds/ddzdl
    dDdx(0,4) = (2./k) * ( cospsi * 0.5 * k * dsddzdl );
  }
  else
  {
    //sin(del) = del - del^3/6 + del^5/120 - del^7/5040. + ...
    //(2/k)*sin(s*k/2) = s - s^3*(k/2)^2/6 + s^5*(k/2)^4/120 - s^7*(k/2)^6/5040
    
    float srun=s;
    float s2 = s*s;
    float khalf2 = 0.5*k;khalf2*=khalf2;
    float krun = khalf2;
    
    D = s;
    dDdx(0,4) = dsddzdl;
    srun *= s2;
    D -= srun*krun/6.;
    dDdx(0,4) -= ( (krun/6.) * 3.*s*s*dsddzdl );
    // krun = (1/4) * k^2
    dDdx(0,2) -= ( (srun/6.) * (2./4.) * k );
    srun *= s2;
    krun *= khalf2;
    D += srun*krun/120.;
    dDdx(0,4) += ( (krun/120.) * 5.*s2*s2*dsddzdl );
    dDdx(0,2) -= ( (srun/120.) * (4./16.)*k*k*k );
    srun *= s2;
    krun *= khalf2;
    D -= srun*krun/5040.;
    dDdx(0,4) -= ( (krun/5040.) * 7.*s2*s2*s2*dsddzdl );
    dDdx(0,2) -= ( (srun/5040.) * (6./64.)*k*k*k*k*k );
  }
  float dpsidk = 0.5*s;float dpsiddzdl = 0.5*k*dsddzdl;
  if( ((int)((floor(fabs(s*k/M_PI)))))%2 == 1  )
  {
    psi = M_PI - psi;
    dpsidk = -dpsidk;
    dpsiddzdl = -dpsiddzdl;
  }
  float cospsi = cos(psi);
  float sinpsi = sin(psi);
  Matrix<float,1,5> dux2dx = Matrix<float,1,5>::Zero(1,5);
  Matrix<float,1,5> duy2dx = Matrix<float,1,5>::Zero(1,5);
  float ux2 = cospsi*ux - sinpsi*uy;
  dux2dx(0,0) = cospsi*duxdphi - sinpsi*duydphi;
  dux2dx(0,2) = ux*-sinpsi*dpsidk - uy*cospsi*dpsidk;
  dux2dx(0,4) = ux*-sinpsi*dpsiddzdl - uy*cospsi*dpsiddzdl;
  float uy2 = sinpsi*ux + cospsi*uy;
  duy2dx(0,0) = sinpsi*duxdphi + cospsi*duydphi;
  duy2dx(0,2) = ux*cospsi*dpsidk - uy*sinpsi*dpsidk;
  duy2dx(0,4) = ux*cospsi*dpsiddzdl - uy*sinpsi*dpsiddzdl;
  
  xp += D*ux2;
  yp += D*uy2;
  zp += dzdl*x(0);
  
  Matrix<float,1,5> temp_1_5 = Matrix<float,1,5>::Zero(1,5);
  temp_1_5 = D*dux2dx + ux2*dDdx;
  dpdx = Matrix<float,3,5>::Zero(3,5);
  for(int i=0;i<5;++i){dpdx(0,i) = temp_1_5(0,i);}
  temp_1_5 = D*duy2dx + uy2*dDdx;
  for(int i=0;i<5;++i){dpdx(1,i) = temp_1_5(0,i);}
  dpdx(2,4) = x(0);
  
  point_covariance += dpdx * covariance * (dpdx.transpose());
  
  point[0] = xp;
  point[1] = yp;
  point[2] = zp;
  
  
  //calculate the new tangent angle
  if(psi>0. || psi < M_PI){psi = 2.*psi;}
  else if(psi < 0.)
  {
    psi = -psi;
    psi = M_PI - psi;
    psi = 2.*psi;
  }
  else
  {
    psi -= M_PI;
    psi = 2.*psi;
  }
  cospsi = cos(psi);
  sinpsi = sin(psi);
  ux2 = cospsi*ux - sinpsi*uy;
  uy2 = sinpsi*ux + cospsi*uy;
  
  
  
  
//   f = (vx - xp)^2 + (vy - yp)^2 + (vz - zp)^2
//   df/dl = 2*(vx - xp)*dxp/dl + 2*(vy - yp)*dyp/dl + 2*(vz - zp)*dzp/dl
//   
//   
//   dx/dt = r*ux2 = ux2/k
//   dy/dt = r*uy2 = uy2/k
//   
//   dx/ds = ux2
//   dy/ds = uy2
//   
//   dx/dl = (dx/ds)*(ds/dl)
//   
//   s^2 = (1 - dzdl^2)*l^2
//   
//   s = sqrt(1 - dzdl^2)*l
//   
//   ds/dl = sqrt(1 - dzdl^2)
//   
//   
//   x = cx + r*cos(t)
//   y = cy + r*sin(t)
//   
//   ux2 = -sin(t)
//   uy2 = cos(t)
//   
//   dx/dl = sqrt(1 - dzdl^2)*ux2
//   d^2x/dl^2 = sqrt(1 - dzdl^2)*d(ux2)/dl
//   
//   dux2/dt = -cos(t) = -uy2
//   dt/ds = k
//   ds/dl = sqrt(1 - dzdl^2)
//   
//   d^2x/dl^2 = -(1 - dzdl^2)*uy2*k
//   
//   d^2y/dl^2 = sqrt(1 - dzdl^2)*d(uy2)/dl
//   d^2y/dl^2 = (1 - dzdl^2)*ux2*k
//   
//   f = (vx - xp)^2 + (vy - yp)^2 + (vz - zp)^2
//   df/dl = -2*(vx - xp)*dxp/dl + -2*(vy - yp)*dyp/dl + -2*(vz - zp)*dzp/dl
//   d^2f/dl^2 = 2*(dxp/dl)^2 + -2*(vx-xp)*d^2xp/dl^2 + ditto for y and z
//   
//   dzp/dl = dzdl
  
  float dx = xp - fixedpars[5];
  float dy = yp - fixedpars[6];
  float dz = zp - fixedpars[7];
  
  val = dx*dx + dy*dy + dz*dz;
  
  grad(0) = 2.*dx*dsdl*ux2 + 2.*dy*dsdl*uy2 + 2.*dz*dzdl;
  
  hessian(0,0) = 2.*dsdl*ux2*dsdl*ux2;
  hessian(0,0) += 2.*dsdl*uy2*dsdl*uy2;
  hessian(0,0) += 2.*dzdl*dzdl;
  hessian(0,0) += -2.*dx*dsdl*dsdl*uy2*k;
  hessian(0,0) += 2.*dy*dsdl*dsdl*ux2*k;
  
  float xylen = sqrt(1-dzdl*dzdl);
  tangent[0] = ux2*xylen;
  tangent[1] = uy2*xylen;
  tangent[2] = dzdl;
  
  return true;
}





//fixedpars[0] is a gaussian regulator sigma
VertexFitFunc::VertexFitFunc() : FunctionGradHessian(3,1)
{
  
}


VertexFitFunc::~VertexFitFunc()
{
  
}


bool VertexFitFunc::calcValGradHessian(const VectorXd& x, double& val, VectorXd& grad, MatrixXd& hessian)
{
  NewtonMinimizerGradHessian minimizer;
  HelixDCAFunc helixfunc;
  minimizer.setFunction(&helixfunc);
  
  val = 0;
  for(int i=0;i<3;i++){grad(i)=0.;}
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      hessian(j,i)=0.;
    }
  }
  
  for(unsigned int i=0;i<tracks->size();i++)
  {
    if(covariances->size() > 0){helixfunc.setCovariance( covariances->at(i) );}
    
    helixfunc.setFixedPar(0, tracks->at(i).phi);
    helixfunc.setFixedPar(1, tracks->at(i).d);
    helixfunc.setFixedPar(2, tracks->at(i).kappa);
    helixfunc.setFixedPar(3, tracks->at(i).z0);
    helixfunc.setFixedPar(4, tracks->at(i).dzdl);
    helixfunc.setFixedPar(5, x(0));
    helixfunc.setFixedPar(6, x(1));
    helixfunc.setFixedPar(7, x(2));
    VectorXd start_point = VectorXd::Zero(1);
    VectorXd min_point = VectorXd::Zero(1);
    //find the point on the helix closest to the point x
    minimizer.minimize(start_point, min_point, 0x1.0p-30, 16, 0x1.0p-40);
    
    //now calculate the chi-square contribution from this track
    double tval=0.;
    VectorXd tgrad = VectorXd::Zero(1);
    MatrixXd thessian = MatrixXd::Zero(1,1);
    
    helixfunc.calcValGradHessian(min_point, tval, tgrad, thessian);
    
    float point[3];float tangent[3];
    point[0] = helixfunc.getPoint(0);
    point[1] = helixfunc.getPoint(1);
    point[2] = helixfunc.getPoint(2);
    tangent[0] = helixfunc.getTangent(0);
    tangent[1] = helixfunc.getTangent(1);
    tangent[2] = helixfunc.getTangent(2);
    Eigen::Matrix<float,3,3> point_covariance = helixfunc.getPointCovariance();
    
//     cout<<"point_covariance = "<<endl<<point_covariance<<endl<<endl;
    
    
    float d[3];
    d[0] = x(0) - point[0];
    d[1] = x(1) - point[1];
    d[2] = x(2) - point[2];
    
    float F = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    
//     Eigen::Matrix<float,1,3> dFdp = Eigen::Matrix<float,1,3>::Zero(1,3);
//     dFdp(0,0) = -2.*d[0];
//     dFdp(0,1) = -2.*d[1];
//     dFdp(0,2) = -2.*d[2];
    
//     float F_cov = dFdp * point_covariance * (dFdp.transpose());
    float F_cov = point_covariance(0,0) + point_covariance(1,1) + point_covariance(2,2);
    float F_err_inv = 1./(F_cov);
    if(covariances->size() == 0){F_err_inv = 1.;}
    
    float f = F;
    
    float g0 = -2.*tangent[0]*( d[0]*tangent[0] + d[1]*tangent[1] + d[2]*tangent[2] ) + 2*d[0];
    float g1 = -2.*tangent[1]*( d[0]*tangent[0] + d[1]*tangent[1] + d[2]*tangent[2] ) + 2*d[1];
    float g2 = -2.*tangent[2]*( d[0]*tangent[0] + d[1]*tangent[1] + d[2]*tangent[2] ) + 2*d[2];
    float h0 = 2.*(1. - tangent[0]*tangent[0]);
    float h1 = 2.*(1. - tangent[1]*tangent[1]);
    float h2 = 2.*(1. - tangent[2]*tangent[2]);
    
    float invs2 = 1./(2.*fixedpars[0]*fixedpars[0]);
    
    float g = -exp(-f*invs2);
    
    if( !(g == g) ){continue;}
    
    val += g*F_err_inv;
    
    float grd0 = -invs2*g*g0*F_err_inv;
    float grd1 = -invs2*g*g1*F_err_inv;
    float grd2 = -invs2*g*g2*F_err_inv;
    
    grad(0) += grd0;
    grad(1) += grd1;
    grad(2) += grd2;
    
    hessian(0,0) += -invs2*F_err_inv*( grd0*g0 + g*h0 );
    hessian(1,1) += -invs2*F_err_inv*( grd1*g1 + g*h1 );
    hessian(2,2) += -invs2*F_err_inv*( grd2*g2 + g*h2 );
    float htemp = -invs2*F_err_inv*( grd1*g0);
    hessian(0,1) += htemp;
    hessian(1,0) += htemp;
    htemp = -invs2*F_err_inv*( grd2*g0 );
    hessian(0,2) += htemp;
    hessian(2,0) += htemp;
    htemp = -invs2*F_err_inv*( grd2*g1 );
    hessian(1,2) += htemp;
    hessian(2,1) += htemp;
  }
  
  return true;
}


