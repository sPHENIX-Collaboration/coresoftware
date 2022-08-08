
/***********************************************************************
* Copyright 1998-2020 CERN for the benefit of the EvtGen authors       *
*                                                                      *
* This file is part of EvtGen.                                         *
*                                                                      *
* EvtGen is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* EvtGen is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     *
***********************************************************************/

#ifndef EVTSPINDENSITY_HH
#define EVTSPINDENSITY_HH
#include "EvtGenBase/EvtComplex.hh"

// Description: This class holds a spin density matrix, it is
//              a complex nxn matrix.

class EvtSpinDensity {
  public:
    EvtSpinDensity( const EvtSpinDensity& density );
    EvtSpinDensity& operator=( const EvtSpinDensity& density );
    virtual ~EvtSpinDensity();

    EvtSpinDensity();
    void setDim( int n );
    int getDim() const;
    void set( int i, int j, const EvtComplex& rhoij );
    const EvtComplex& get( int i, int j ) const;
    double normalizedProb( const EvtSpinDensity& d );
    friend std::ostream& operator<<( std::ostream& s, const EvtSpinDensity& d );
    void setDiag( int n );

    int check();

  private:
    EvtComplexPtrPtr rho;
    int dim;
};

#endif
