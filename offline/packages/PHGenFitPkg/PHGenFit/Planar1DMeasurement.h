#ifndef PHGENFIT_PLANAR1DMEASUREMENT_H
#define PHGENFIT_PLANAR1DMEASUREMENT_H

/*!
 * \file		Planar1DMeasurement.h
 * \brief		Handles one dimentional planar measurements.
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "Measurement.h"

class TVector3;

namespace PHGenFit
{
  class Planar1DMeasurement : public Measurement
  {
    public:
    //! constructor
    Planar1DMeasurement(const TVector3& pos, const TVector3& u, const TVector3& v, double error);
    
    //! constructor 
    /*! \param n the vector normal to the plane. u and v plane vectors are defined as such: 
     * u is orthogonal to n, and in the xy plane
     * v is perpendicular to n and u, and thus along z
     * This constructor assumes that the measurement is along the u direction as defined above
     * to make the measurement along the v direction, use the setStripV method
     */
    Planar1DMeasurement(const TVector3& pos, const TVector3& n, double error);
    
    /** @brief Use if the coordinate for 1D hits measured in V direction.
    *
    * Per default for 1D planar hits, the coordinate is measured in U direction.
    * With this function you can set it to be measured in V direction.
    * This affects the outcoe of constructHMatrix().
    */
    void setStripV(bool = true);
    
    private:
    
    //* actual setup of the measurement
    void init(const TVector3& pos, const TVector3& u, const TVector3& v, const double error);

    
  };
  
}  // namespace PHGenFit

#endif
