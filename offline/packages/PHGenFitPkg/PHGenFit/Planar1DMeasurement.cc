/*!
 * \file		Planar1DMeasurement.cc
 * \brief		Handles the palnar type of measurements.
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "Planar1DMeasurement.h"

#include <GenFit/PlanarMeasurement.h>
#include <GenFit/DetPlane.h>
#include <GenFit/SharedPlanePtr.h>     // for SharedPlanePtr
#include <GenFit/StateOnPlane.h>

#include <TMatrixDSymfwd.h>            // for TMatrixDSym
#include <TMatrixTSym.h>               // for TMatrixTSym
#include <TVectorDfwd.h>               // for TVectorD
#include <TVectorT.h>                  // for TVectorT
#include <TVector3.h>

namespace 
{ 
  template<class T> 
    inline constexpr T square( const T& x ) { return x*x; }
}
  
namespace PHGenFit
{
  //___________________________________________________________
  void Planar1DMeasurement::init(const TVector3& pos, const TVector3& u, const TVector3& v, double error)
  {
    int nDim = 1;
    TVectorD hitCoords(nDim);
    TMatrixDSym hitCov(nDim);
    
    hitCoords(0) = 0;
    hitCov(0, 0) = square( error );
    
    genfit::SharedPlanePtr plane( new genfit::DetPlane(pos, u, v));
    
    static constexpr int measurementCounter = 0;
    auto measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, -1, measurementCounter, nullptr);
    measurement->setPlane(plane, measurementCounter);
    _measurement = measurement;
  }
  
  //___________________________________________________________
  Planar1DMeasurement::Planar1DMeasurement(const TVector3& pos, const TVector3& u, const TVector3& v, double error)
  { init(pos, u, v, error); }
  
  //___________________________________________________________
  Planar1DMeasurement::Planar1DMeasurement(const TVector3& pos, const TVector3& n, double error)
  {
    /*!
    *  Mainly for "n" in xy plane, or n.z() = 0;
    *  In this case, According to https://root.cern.ch/doc/master/TVector3_8h_source.html#l00301
    *  u = n.Orthogonal() is always in xy plane
    *  v = n x u, is along z direction.
    *  If z is not 0, but z <= x or y, then u will still be in xy plane, phi direction error, but v will not be along z axis.
    */
    TVector3 u = n.Orthogonal();
    TVector3 v = n.Cross(u);
    init(pos, u, v, error);
  }
  
  //___________________________________________________________
  void Planar1DMeasurement::setStripV(bool value) 
  { static_cast<genfit::PlanarMeasurement*>( _measurement )->setStripV( value ); }

}
