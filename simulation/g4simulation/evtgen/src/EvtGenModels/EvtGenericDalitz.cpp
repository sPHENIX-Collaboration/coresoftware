
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

#include "EvtGenModels/EvtGenericDalitz.hh"

#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"

#include "EvtGenModels/EvtDalitzTable.hh"

std::string EvtGenericDalitz::getName()
{
    return "GENERIC_DALITZ";
}

EvtDecayBase* EvtGenericDalitz::clone()
{
    return new EvtGenericDalitz();
}

void EvtGenericDalitz::init()
{
    checkNArg( 1 );

    EvtId parnum = getParentId();
    EvtId d1 = getDaug( 0 );
    EvtId d2 = getDaug( 1 );
    EvtId d3 = getDaug( 2 );

    std::vector<EvtDalitzDecayInfo> decays =
        EvtDalitzTable::getInstance( getArgStr( 0 ) )->getDalitzTable( parnum );

    std::vector<EvtDalitzDecayInfo>::iterator i = decays.begin();
    for ( ; i != decays.end(); i++ ) {
        EvtId daughter1 = ( *i ).daughter1();
        EvtId daughter2 = ( *i ).daughter2();
        EvtId daughter3 = ( *i ).daughter3();

        if ( d1 == daughter1 && d2 == daughter2 && d3 == daughter3 ) {
            _d1 = 0;
            _d2 = 1;
            _d3 = 2;
        } else if ( d1 == daughter1 && d2 == daughter3 && d3 == daughter2 ) {
            _d1 = 0;
            _d2 = 2;
            _d3 = 1;
        } else if ( d1 == daughter2 && d2 == daughter1 && d3 == daughter3 ) {
            _d1 = 1;
            _d2 = 0;
            _d3 = 2;
        } else if ( d1 == daughter2 && d2 == daughter3 && d3 == daughter1 ) {
            _d1 = 1;
            _d2 = 2;
            _d3 = 0;
        } else if ( d1 == daughter3 && d2 == daughter1 && d3 == daughter2 ) {
            _d1 = 2;
            _d2 = 0;
            _d3 = 1;
        } else if ( d1 == daughter3 && d2 == daughter2 && d3 == daughter1 ) {
            _d1 = 2;
            _d2 = 1;
            _d3 = 0;
        } else {
            continue;
        }

        _resonances = ( *i ).getResonances();
        setProbMax( ( *i ).getProbMax() );
        return;
    }
}

void EvtGenericDalitz::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtVector4R p4_d1 = p->getDaug( _d1 )->getP4();
    EvtVector4R p4_d2 = p->getDaug( _d2 )->getP4();
    EvtVector4R p4_d3 = p->getDaug( _d3 )->getP4();

    double mA = p->getDaug( _d1 )->mass();
    double mB = p->getDaug( _d2 )->mass();
    double mC = p->getDaug( _d3 )->mass();

    double m2AB = ( p4_d1 + p4_d2 ).mass2();
    double m2CA = ( p4_d1 + p4_d3 ).mass2();
    double m2BC = ( p4_d2 + p4_d3 ).mass2();

    EvtDalitzPoint point( mA, mB, mC, m2AB, m2BC, m2CA );

    EvtComplex amp( 0, 0 );
    std::vector<std::pair<EvtComplex, EvtDalitzReso>>::iterator i =
        _resonances.begin();
    for ( ; i != _resonances.end(); i++ ) {
        std::pair<EvtComplex, EvtDalitzReso> res = ( *i );
        amp += res.first * res.second.evaluate( point );
    }

    vertex( amp );
    return;
}

std::string EvtGenericDalitz::getParamName( int i )
{
    switch ( i ) {
        case 0:
            return "xmlFile";
        default:
            return "";
    }
}
