
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

#include "EvtGenBase/EvtParticleDecayList.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtStatus.hh"

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
using std::endl;
using std::fstream;

EvtParticleDecayList::EvtParticleDecayList( const EvtParticleDecayList& o )
{
    _nmode = o._nmode;
    _rawbrfrsum = o._rawbrfrsum;
    _decaylist = new EvtParticleDecayPtr[_nmode];

    int i;
    for ( i = 0; i < _nmode; i++ ) {
        _decaylist[i] = new EvtParticleDecay;

        EvtDecayBase* tModel = o._decaylist[i]->getDecayModel();

        EvtDecayBase* tModelNew = tModel->clone();
        if ( tModel->getPHOTOS() ) {
            tModelNew->setPHOTOS();
        }
        if ( tModel->verbose() ) {
            tModelNew->setVerbose();
        }
        if ( tModel->summary() ) {
            tModelNew->setSummary();
        }
        std::vector<std::string> args;
        int j;
        for ( j = 0; j < tModel->getNArg(); j++ ) {
            args.push_back( tModel->getArgStr( j ) );
        }
        tModelNew->saveDecayInfo( tModel->getParentId(), tModel->getNDaug(),
                                  tModel->getDaugs(), tModel->getNArg(), args,
                                  tModel->getModelName(),
                                  tModel->getBranchingFraction() );
        _decaylist[i]->setDecayModel( tModelNew );

        _decaylist[i]->setBrfrSum( o._decaylist[i]->getBrfrSum() );
        _decaylist[i]->setMassMin( o._decaylist[i]->getMassMin() );
    }
}

EvtParticleDecayList::~EvtParticleDecayList()
{
    int i;
    for ( i = 0; i < _nmode; i++ ) {
        delete _decaylist[i];
    }

    if ( _decaylist != 0 )
        delete[] _decaylist;
}

void EvtParticleDecayList::printSummary()
{
    int i;
    for ( i = 0; i < _nmode; i++ ) {
        _decaylist[i]->printSummary();
    }
}

void EvtParticleDecayList::removeDecay()
{
    int i;
    for ( i = 0; i < _nmode; i++ ) {
        delete _decaylist[i];
    }

    delete[] _decaylist;
    _decaylist = 0;
    _nmode = 0;
    _rawbrfrsum = 0.0;
}

EvtDecayBase* EvtParticleDecayList::getDecayModel( int imode )
{
    EvtDecayBase* theModel( 0 );
    if ( imode >= 0 && imode < _nmode ) {
        EvtParticleDecay* theDecay = _decaylist[imode];
        if ( theDecay != 0 ) {
            theModel = theDecay->getDecayModel();
        }
    }

    return theModel;
}

EvtDecayBase* EvtParticleDecayList::getDecayModel( EvtParticle* p )
{
    if ( p->getNDaug() != 0 ) {
        assert( p->getChannel() >= 0 );
        return getDecay( p->getChannel() ).getDecayModel();
    }
    if ( p->getChannel() > ( -1 ) ) {
        return getDecay( p->getChannel() ).getDecayModel();
    }

    if ( getNMode() == 0 ) {
        return 0;
    }
    if ( getRawBrfrSum() < 0.00000001 ) {
        return 0;
    }

    if ( getNMode() == 1 ) {
        p->setChannel( 0 );
        return getDecay( 0 ).getDecayModel();
    }

    if ( p->getChannel() > ( -1 ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Internal error!!!" << endl;
        ::abort();
    }

    int j;

    for ( j = 0; j < 10000000; j++ ) {
        double u = EvtRandom::Flat();

        int i;
        bool breakL = false;
        for ( i = 0; i < getNMode(); i++ ) {
            if ( breakL )
                continue;
            if ( u < getDecay( i ).getBrfrSum() ) {
                breakL = true;
                //special case for decay of on particel to another
                // e.g. K0->K0S

                if ( getDecay( i ).getDecayModel()->getNDaug() == 1 ) {
                    p->setChannel( i );
                    return getDecay( i ).getDecayModel();
                }

                if ( p->hasValidP4() ) {
                    if ( getDecay( i ).getMassMin() < p->mass() ) {
                        p->setChannel( i );
                        return getDecay( i ).getDecayModel();
                    }
                } else {
                    //Lange apr29-2002 - dont know the mass yet
                    p->setChannel( i );
                    return getDecay( i ).getDecayModel();
                }
            }
        }
    }

    //Ok, we tried 10000000 times above to pick a decay channel that is
    //kinematically allowed! Now we give up and search all channels!
    //if that fails, the particle will not be decayed!

    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Tried 10000000 times to generate decay of "
        << EvtPDL::name( p->getId() ) << " with mass=" << p->mass() << endl;
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Will take first kinematically allowed decay in the decay table"
        << endl;

    int i;

    //Need to check that we don't use modes with 0 branching fractions.
    double previousBrSum = 0.0;
    for ( i = 0; i < getNMode(); i++ ) {
        if ( getDecay( i ).getBrfrSum() != previousBrSum ) {
            if ( getDecay( i ).getMassMin() < p->mass() ) {
                p->setChannel( i );
                return getDecay( i ).getDecayModel();
            }
        }
        previousBrSum = getDecay( i ).getBrfrSum();
    }

    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Could not decay:" << EvtPDL::name( p->getId() ).c_str()
        << " with mass:" << p->mass() << " will throw event away! " << endl;

    EvtStatus::setRejectFlag();
    return 0;
}

void EvtParticleDecayList::setNMode( int nmode )
{
    EvtParticleDecayPtr* _decaylist_new = new EvtParticleDecayPtr[nmode];

    if ( _nmode != 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Error _nmode not equal to zero!!!" << endl;
        ::abort();
    }
    if ( _decaylist != 0 ) {
        delete[] _decaylist;
    }
    _decaylist = _decaylist_new;
    _nmode = nmode;
}

EvtParticleDecay& EvtParticleDecayList::getDecay( int nchannel ) const
{
    if ( nchannel >= _nmode ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Error getting channel:" << nchannel << " with only " << _nmode
            << " stored!" << endl;
        ::abort();
    }
    return *( _decaylist[nchannel] );
}

void EvtParticleDecayList::makeChargeConj( EvtParticleDecayList* conjDecayList )
{
    _rawbrfrsum = conjDecayList->_rawbrfrsum;

    setNMode( conjDecayList->_nmode );

    int i;

    for ( i = 0; i < _nmode; i++ ) {
        _decaylist[i] = new EvtParticleDecay;
        _decaylist[i]->chargeConj( conjDecayList->_decaylist[i] );
    }
}

void EvtParticleDecayList::addMode( EvtDecayBase* decay, double brfrsum,
                                    double massmin )
{
    EvtParticleDecayPtr* newlist = new EvtParticleDecayPtr[_nmode + 1];

    int i;
    for ( i = 0; i < _nmode; i++ ) {
        newlist[i] = _decaylist[i];
    }

    _rawbrfrsum = brfrsum;

    newlist[_nmode] = new EvtParticleDecay;

    newlist[_nmode]->setDecayModel( decay );
    newlist[_nmode]->setBrfrSum( brfrsum );
    newlist[_nmode]->setMassMin( massmin );

    EvtDecayBase* newDec = newlist[_nmode]->getDecayModel();
    for ( i = 0; i < _nmode; i++ ) {
        if ( newDec->matchingDecay( *( newlist[i]->getDecayModel() ) ) ) {
            //sometimes its ok..
            if ( newDec->getModelName() == "JETSET" ||
                 newDec->getModelName() == "PYTHIA" )
                continue;
            if ( newDec->getModelName() == "JSCONT" ||
                 newDec->getModelName() == "PYCONT" )
                continue;
            if ( newDec->getModelName() == "PYGAGA" )
                continue;
            if ( newDec->getModelName() == "LUNDAREALAW" )
                continue;
            if ( newDec->getModelName() == "TAUOLA" )
                continue;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Two matching decays with same parent in decay table\n";
            EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Please fix that\n";
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Parent " << EvtPDL::name( newDec->getParentId() ).c_str()
                << endl;
            for ( int j = 0; j < newDec->getNDaug(); j++ )
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Daughter "
                    << EvtPDL::name( newDec->getDaug( j ) ).c_str() << endl;
            assert( 0 );
        }
    }

    if ( _nmode != 0 ) {
        delete[] _decaylist;
    }

    if ( ( _nmode == 0 ) && ( _decaylist != 0 ) )
        delete[] _decaylist;

    _nmode++;

    _decaylist = newlist;
}

void EvtParticleDecayList::finalize()
{
    if ( _nmode > 0 ) {
        if ( _rawbrfrsum < 0.000001 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Please give me a "
                << "branching fraction sum greater than 0\n";
            assert( 0 );
        }
        if ( fabs( _rawbrfrsum - 1.0 ) > 0.0001 ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Warning, sum of branching fractions for "
                << EvtPDL::name( _decaylist[0]->getDecayModel()->getParentId() )
                       .c_str()
                << " is " << _rawbrfrsum << endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "rescaled to one! " << endl;
        }

        int i;

        for ( i = 0; i < _nmode; i++ ) {
            double brfrsum = _decaylist[i]->getBrfrSum() / _rawbrfrsum;
            _decaylist[i]->setBrfrSum( brfrsum );
        }
    }
}

EvtParticleDecayList& EvtParticleDecayList::operator=( const EvtParticleDecayList& o )
{
    if ( this != &o ) {
        removeDecay();
        _nmode = o._nmode;
        _rawbrfrsum = o._rawbrfrsum;
        _decaylist = new EvtParticleDecayPtr[_nmode];

        int i;
        for ( i = 0; i < _nmode; i++ ) {
            _decaylist[i] = new EvtParticleDecay;

            EvtDecayBase* tModel = o._decaylist[i]->getDecayModel();

            EvtDecayBase* tModelNew = tModel->clone();
            if ( tModel->getPHOTOS() ) {
                tModelNew->setPHOTOS();
            }
            if ( tModel->verbose() ) {
                tModelNew->setVerbose();
            }
            if ( tModel->summary() ) {
                tModelNew->setSummary();
            }
            std::vector<std::string> args;
            int j;
            for ( j = 0; j < tModel->getNArg(); j++ ) {
                args.push_back( tModel->getArgStr( j ) );
            }
            tModelNew->saveDecayInfo( tModel->getParentId(), tModel->getNDaug(),
                                      tModel->getDaugs(), tModel->getNArg(),
                                      args, tModel->getModelName(),
                                      tModel->getBranchingFraction() );
            _decaylist[i]->setDecayModel( tModelNew );

            //_decaylist[i]->setDecayModel(tModel);
            _decaylist[i]->setBrfrSum( o._decaylist[i]->getBrfrSum() );
            _decaylist[i]->setMassMin( o._decaylist[i]->getMassMin() );
        }
    }
    return *this;
}

void EvtParticleDecayList::removeMode( EvtDecayBase* decay )
{
    // here we will delete a decay with the same final state particles
    // and recalculate the branching fractions for the remaining modes
    int match = -1;
    int i;
    double match_bf;

    for ( i = 0; i < _nmode; i++ ) {
        if ( decay->matchingDecay( *( _decaylist[i]->getDecayModel() ) ) ) {
            match = i;
        }
    }

    if ( match < 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " Attempt to remove undefined mode for" << endl
            << "Parent " << EvtPDL::name( decay->getParentId() ).c_str() << endl
            << "Daughters: ";
        for ( int j = 0; j < decay->getNDaug(); j++ )
            EvtGenReport( EVTGEN_ERROR, "" )
                << EvtPDL::name( decay->getDaug( j ) ).c_str() << " ";
        EvtGenReport( EVTGEN_ERROR, "" ) << endl;
        ::abort();
    }

    if ( match == 0 ) {
        match_bf = _decaylist[match]->getBrfrSum();
    } else {
        match_bf = ( _decaylist[match]->getBrfrSum() -
                     _decaylist[match - 1]->getBrfrSum() );
    }

    double divisor = 1 - match_bf;
    if ( divisor < 0.000001 && _nmode > 1 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Removing requested mode leaves "
            << EvtPDL::name( decay->getParentId() ).c_str()
            << " with zero sum branching fraction," << endl
            << "but more than one decay mode remains. Aborting." << endl;
        ::abort();
    }

    EvtParticleDecayPtr* newlist = new EvtParticleDecayPtr[_nmode - 1];

    for ( i = 0; i < match; i++ ) {
        newlist[i] = _decaylist[i];
        newlist[i]->setBrfrSum( newlist[i]->getBrfrSum() / divisor );
    }
    for ( i = match + 1; i < _nmode; i++ ) {
        newlist[i - 1] = _decaylist[i];
        newlist[i - 1]->setBrfrSum(
            ( newlist[i - 1]->getBrfrSum() - match_bf ) / divisor );
    }

    delete[] _decaylist;

    _nmode--;

    _decaylist = newlist;

    if ( _nmode == 0 ) {
        delete[] _decaylist;
    }
}

bool EvtParticleDecayList::isJetSet() const
{
    int i;
    EvtDecayBase* decayer;

    for ( i = 0; i < getNMode(); i++ ) {
        decayer = getDecay( i ).getDecayModel();
        if ( decayer->getModelName() == "PYTHIA" )
            return true;
    }

    return false;
}
