
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

#include "EvtGenBase/EvtPartProp.hh"

#include "EvtGenBase/EvtAbsLineShape.hh"
#include "EvtGenBase/EvtFlatLineShape.hh"
#include "EvtGenBase/EvtManyDeltaFuncLineShape.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRelBreitWignerBarrierFact.hh"

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
using std::fstream;

EvtPartProp::EvtPartProp() :
    _id( -1, -1 ), _idchgconj( -1, -1 ), _chg3( 0 ), _stdhep( 0 ), _lundkc( 0 )
{
    _ctau = 0.0;
    _name = "*******";
    _spintype = EvtSpinType::SCALAR;
}

EvtPartProp::EvtPartProp( const EvtPartProp& x )
{
    _lineShape.reset( x._lineShape ? x._lineShape->clone() : nullptr );
    _ctau = x._ctau;
    _name = x._name;
    _spintype = x._spintype;
    _id = x._id;
    _idchgconj = x._idchgconj;
    _chg3 = x._chg3;
    _stdhep = x._stdhep;
    _lundkc = x._lundkc;
}

void EvtPartProp::setName( std::string pname )
{
    _name = pname;
}

EvtPartProp& EvtPartProp::operator=( const EvtPartProp& x )
{
    _lineShape.reset( x._lineShape ? x._lineShape->clone() : nullptr );

    _ctau = x._ctau;
    _name = x._name;
    _chg3 = x._chg3;
    _spintype = x._spintype;
    return *this;
}

void EvtPartProp::initLineShape( double mass, double width, double maxRange )
{
    _lineShape = std::make_unique<EvtRelBreitWignerBarrierFact>( mass, width,
                                                                 maxRange,
                                                                 _spintype );
}

void EvtPartProp::newLineShape( std::string type )
{
    double m = _lineShape->getMass();
    double w = _lineShape->getWidth();
    double mR = _lineShape->getMaxRange();
    EvtSpinType::spintype st = _lineShape->getSpinType();
    if ( type == "RELBW" ) {
        _lineShape = std::make_unique<EvtRelBreitWignerBarrierFact>( m, w, mR,
                                                                     st );
    } else if ( type == "NONRELBW" ) {
        _lineShape = std::make_unique<EvtAbsLineShape>( m, w, mR, st );
    } else if ( type == "FLAT" ) {
        _lineShape = std::make_unique<EvtFlatLineShape>( m, w, mR, st );
    } else if ( type == "MANYDELTAFUNC" ) {
        _lineShape = std::make_unique<EvtManyDeltaFuncLineShape>( m, w, mR, st );
    } else {
        _lineShape.reset();
    }
}

void EvtPartProp::reSetMass( double mass )
{
    if ( !_lineShape )
        ::abort();
    _lineShape->reSetMass( mass );
}
void EvtPartProp::reSetWidth( double width )
{
    if ( !_lineShape )
        ::abort();
    _lineShape->reSetWidth( width );
}

void EvtPartProp::setPWForDecay( int spin, EvtId d1, EvtId d2 )
{
    if ( !_lineShape )
        ::abort();
    _lineShape->setPWForDecay( spin, d1, d2 );
}

void EvtPartProp::setPWForBirthL( int spin, EvtId par, EvtId othD )
{
    if ( !_lineShape )
        ::abort();
    _lineShape->setPWForBirthL( spin, par, othD );
}

void EvtPartProp::reSetMassMin( double mass )
{
    if ( !_lineShape )
        ::abort();
    _lineShape->reSetMassMin( mass );
}
void EvtPartProp::reSetMassMax( double mass )
{
    if ( !_lineShape )
        ::abort();
    _lineShape->reSetMassMax( mass );
}
void EvtPartProp::reSetBlatt( double blatt )
{
    if ( !_lineShape )
        ::abort();
    _lineShape->reSetBlatt( blatt );
}
void EvtPartProp::reSetBlattBirth( double blatt )
{
    if ( !_lineShape )
        ::abort();
    _lineShape->reSetBlattBirth( blatt );
}
void EvtPartProp::includeBirthFactor( bool yesno )
{
    if ( !_lineShape )
        ::abort();
    _lineShape->includeBirthFactor( yesno );
}
void EvtPartProp::includeDecayFactor( bool yesno )
{
    if ( !_lineShape )
        ::abort();
    _lineShape->includeDecayFactor( yesno );
}
