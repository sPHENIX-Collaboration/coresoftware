
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

#include "EvtGenBase/EvtSecondary.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <iostream>
using std::endl;
using std::ostream;

void EvtSecondary::init()
{
    _npart = 0;
}

int EvtSecondary::getNPart()
{
    return _npart;
}

void EvtSecondary::createSecondary( int stdhepindex, EvtParticle* prnt )
{
    _stdhepindex[_npart] = stdhepindex;
    if ( prnt->getNDaug() == 0 ) {
        _id1[_npart] = 0;
        _id2[_npart] = 0;
        _id3[_npart] = 0;
        _npart++;
        return;
    }
    if ( prnt->getNDaug() == 1 ) {
        _id1[_npart] = EvtPDL::getStdHep( prnt->getDaug( 0 )->getId() );
        _id2[_npart] = 0;
        _id3[_npart] = 0;
        _npart++;
        return;
    }
    if ( prnt->getNDaug() == 2 ) {
        _id1[_npart] = EvtPDL::getStdHep( prnt->getDaug( 0 )->getId() );
        _id2[_npart] = EvtPDL::getStdHep( prnt->getDaug( 1 )->getId() );
        _id3[_npart] = 0;
        _npart++;
        return;
    }
    if ( prnt->getNDaug() == 3 ) {
        _id1[_npart] = EvtPDL::getStdHep( prnt->getDaug( 0 )->getId() );
        _id2[_npart] = EvtPDL::getStdHep( prnt->getDaug( 1 )->getId() );
        _id3[_npart] = EvtPDL::getStdHep( prnt->getDaug( 2 )->getId() );
        _npart++;
        return;
    }

    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "More than 3 decay products in a secondary particle!" << endl;
}

ostream& operator<<( ostream& s, const EvtSecondary& secondary )
{
    s << endl;
    s << "Secondary decays:" << endl;

    int i;
    for ( i = 0; i < secondary._npart; i++ ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << i << " " << secondary._stdhepindex[i] << " " << secondary._id1[i]
            << " " << secondary._id2[i] << " " << secondary._id3[i] << endl;
    }

    s << endl;

    return s;
}
