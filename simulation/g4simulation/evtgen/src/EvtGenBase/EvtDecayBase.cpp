
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

#include "EvtGenBase/EvtDecayBase.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtStatus.hh"

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
using std::endl;
using std::fstream;
void EvtDecayBase::checkQ()
{
    int i;
    int q = 0;
    int qpar;

    //If there are no daughters (jetset etc) then we do not
    //want to do this test.  Why?  Because sometimes the parent
    //will have a nonzero charge.

    if ( _ndaug != 0 ) {
        for ( i = 0; i < _ndaug; i++ ) {
            q += EvtPDL::chg3( _daug[i] );
        }
        qpar = EvtPDL::chg3( _parent );

        if ( q != qpar ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << _modelname.c_str() << " generator expected "
                << " charge to be conserved, found:" << endl;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Parent charge of " << ( qpar / 3 ) << endl;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Sum of daughter charge of " << ( q / 3 ) << endl;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "The parent is " << EvtPDL::name( _parent ).c_str() << endl;
            for ( i = 0; i < _ndaug; i++ ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Daughter " << EvtPDL::name( _daug[i] ).c_str() << endl;
            }
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Will terminate execution!" << endl;

            ::abort();
        }
    }
}

double EvtDecayBase::getProbMax( double prob )
{
    int i;

    //diagnostics
    sum_prob += prob;
    if ( prob > max_prob )
        max_prob = prob;

    if ( defaultprobmax && ntimes_prob <= 500 ) {
        //We are building up probmax with this iteration
        ntimes_prob += 1;
        if ( prob > probmax ) {
            probmax = prob;
        }
        if ( ntimes_prob == 500 ) {
            probmax *= 1.2;
        }
        return 1000000.0 * prob;
    }

    if ( prob > probmax * 1.0001 ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "prob > probmax:(" << prob << ">" << probmax << ")";
        EvtGenReport( EVTGEN_INFO, "" ) << "(" << _modelname.c_str() << ") ";
        EvtGenReport( EVTGEN_INFO, "" )
            << EvtPDL::name( _parent ).c_str() << " -> ";
        for ( i = 0; i < _ndaug; i++ ) {
            EvtGenReport( EVTGEN_INFO, "" )
                << EvtPDL::name( _daug[i] ).c_str() << " ";
        }
        EvtGenReport( EVTGEN_INFO, "" ) << endl;

        if ( defaultprobmax )
            probmax = prob;
    }

    ntimes_prob += 1;

    return probmax;

}    //getProbMax

double EvtDecayBase::resetProbMax( double prob )
{
    EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Reseting prob max\n";
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "prob > probmax:(" << prob << ">" << probmax << ")";
    EvtGenReport( EVTGEN_INFO, "" ) << "(" << _modelname.c_str() << ")";
    EvtGenReport( EVTGEN_INFO, "" ) << EvtPDL::getStdHep( _parent ) << "->";

    for ( int i = 0; i < _ndaug; i++ ) {
        EvtGenReport( EVTGEN_INFO, "" ) << EvtPDL::getStdHep( _daug[i] ) << " ";
    }
    EvtGenReport( EVTGEN_INFO, "" ) << endl;

    probmax = 0.0;
    defaultprobmax = 0;
    ntimes_prob = 0;

    return prob;
}

std::string EvtDecayBase::commandName()
{
    return std::string( "" );
}

void EvtDecayBase::command( std::string )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Should never call EvtDecayBase::command" << endl;
    ::abort();
}

std::string EvtDecayBase::getParamName( int i )
{
    switch ( i ) {
        case 0:
            return "param00";
        case 1:
            return "param01";
        case 2:
            return "param02";
        case 3:
            return "param03";
        case 4:
            return "param04";
        case 5:
            return "param05";
        case 6:
            return "param06";
        case 7:
            return "param07";
        case 8:
            return "param08";
        case 9:
            return "param09";
        default:
            return "";
    }
}

std::string EvtDecayBase::getParamDefault( int /*i*/ )
{
    return "";
}

void EvtDecayBase::init()
{
    //This default version of init does nothing;
    //A specialized version of this function can be
    //supplied for each decay model to do initialization.

    return;
}

void EvtDecayBase::initProbMax()
{
    //This function is called if the decay does not have a
    //specialized initialization.
    //The default is to set the maximum
    //probability to 0 and the number of times called to 0
    //and defaultprobmax to 1 such that the decay will be
    //generated many many times
    //in order to generate a reasonable maximum probability
    //for the decay.

    defaultprobmax = 1;
    ntimes_prob = 0;
    probmax = 0.0;

}    //initProbMax

void EvtDecayBase::saveDecayInfo( EvtId ipar, int ndaug, EvtId* daug, int narg,
                                  std::vector<std::string>& args,
                                  std::string name, double brfr )
{
    int i;

    _brfr = brfr;
    _ndaug = ndaug;
    _narg = narg;
    _parent = ipar;

    _dsum = 0;

    if ( _ndaug > 0 ) {
        _daug.resize( _ndaug );
        for ( i = 0; i < _ndaug; i++ ) {
            _daug[i] = daug[i];
            _dsum += daug[i].getAlias();
        }
    } else {
        _daug.clear();
    }

    if ( _narg > 0 ) {
        _args.resize( _narg + 1 );
        for ( i = 0; i < _narg; i++ ) {
            _args[i] = args[i];
        }
    } else {
        _args.clear();
    }

    _modelname = name;

    this->init();
    this->initProbMax();

    if ( _chkCharge ) {
        this->checkQ();
    }

    if ( defaultprobmax ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "No default probmax for ";
        EvtGenReport( EVTGEN_INFO, "" ) << "(" << _modelname.c_str() << ") ";
        EvtGenReport( EVTGEN_INFO, "" )
            << EvtPDL::name( _parent ).c_str() << " -> ";
        for ( i = 0; i < _ndaug; i++ ) {
            EvtGenReport( EVTGEN_INFO, "" )
                << EvtPDL::name( _daug[i] ).c_str() << " ";
        }
        EvtGenReport( EVTGEN_INFO, "" ) << endl;
        EvtGenReport( EVTGEN_INFO, "" )
            << "This is fine for development, but must be provided for production."
            << endl;
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Never fear though - the decay will use the \n";
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "500 iterations to build up a good probmax \n";
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "before accepting a decay. " << endl;
    }
}

EvtDecayBase::EvtDecayBase()
{
    //the default is that the user module does _not_ set
    // any probmax.
    defaultprobmax = 1;
    ntimes_prob = 0;
    probmax = 0.0;

    _photos = 0;
    _verbose = 0;
    _summary = 0;
    _parent = EvtId( -1, -1 );
    _ndaug = 0;
    _narg = 0;
    _modelname = "**********";

    //Default is to check that charge is conserved
    _chkCharge = 1;

    //statistics collection!

    max_prob = 0.0;
    sum_prob = 0.0;
}

void EvtDecayBase::printSummary() const
{
    if ( ntimes_prob > 0 ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Calls = " << ntimes_prob
            << " eff: " << sum_prob / ( probmax * ntimes_prob )
            << " frac. max:" << max_prob / probmax;
        EvtGenReport( EVTGEN_INFO, "" )
            << " probmax:" << probmax << " max:" << max_prob << " : ";
    }

    printInfo();
}

void EvtDecayBase::printInfo() const
{
    EvtGenReport( EVTGEN_INFO, "" ) << EvtPDL::name( _parent ).c_str() << " -> ";
    for ( int i = 0; i < _ndaug; i++ ) {
        EvtGenReport( EVTGEN_INFO, "" )
            << EvtPDL::name( _daug[i] ).c_str() << " ";
    }
    EvtGenReport( EVTGEN_INFO, "" ) << " (" << _modelname.c_str() << ")" << endl;
}

void EvtDecayBase::setProbMax( double prbmx )
{
    defaultprobmax = 0;
    probmax = prbmx;
}

void EvtDecayBase::noProbMax()
{
    defaultprobmax = 0;
}

double EvtDecayBase::findMaxMass( EvtParticle* p )
{
    double maxOkMass = EvtPDL::getMaxMass( p->getId() );

    //protect against vphotons
    if ( maxOkMass < 0.0000000001 )
        return 10000000.;
    //and against already determined masses
    if ( p->hasValidP4() )
        maxOkMass = p->mass();

    EvtParticle* par = p->getParent();
    if ( par ) {
        double maxParMass = findMaxMass( par );
        size_t i;
        double minDaugMass = 0.;
        for ( i = 0; i < par->getNDaug(); i++ ) {
            EvtParticle* dau = par->getDaug( i );
            if ( dau != p ) {
                // it might already have a mass
                if ( dau->isInitialized() || dau->hasValidP4() )
                    minDaugMass += dau->mass();
                else
                    //give it a bit of phase space
                    minDaugMass += 1.000001 * EvtPDL::getMinMass( dau->getId() );
            }
        }
        if ( maxOkMass > ( maxParMass - minDaugMass ) )
            maxOkMass = maxParMass - minDaugMass;
    }
    return maxOkMass;
}

// given list of daughters ( by number ) returns a
// list of viable masses.

void EvtDecayBase::findMass( EvtParticle* p )
{
    //Need to also check that this mass does not screw
    //up the parent
    //This code assumes that for the ith daughter, 0..i-1
    //already have a mass
    double maxOkMass = findMaxMass( p );

    int count = 0;
    double mass;
    bool massOk = false;
    size_t i;
    while ( !massOk ) {
        count++;
        if ( count > 10000 ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Can not find a valid mass for: "
                << EvtPDL::name( p->getId() ).c_str() << endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Now printing parent and/or grandparent tree\n";
            if ( p->getParent() ) {
                if ( p->getParent()->getParent() ) {
                    p->getParent()->getParent()->printTree();
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << p->getParent()->getParent()->mass() << endl;
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << p->getParent()->mass() << endl;
                } else {
                    p->getParent()->printTree();
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << p->getParent()->mass() << endl;
                }
            } else
                p->printTree();
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "maxokmass=" << maxOkMass << " "
                << EvtPDL::getMinMass( p->getId() ) << " "
                << EvtPDL::getMaxMass( p->getId() ) << endl;
            if ( p->getNDaug() ) {
                for ( i = 0; i < p->getNDaug(); i++ ) {
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << p->getDaug( i )->mass() << " ";
                }
                EvtGenReport( EVTGEN_INFO, "EvtGen" ) << endl;
            }
            if ( maxOkMass >= EvtPDL::getMinMass( p->getId() ) ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" )
                    << "taking a default value\n";
                p->setMass( maxOkMass );
                return;
            }
            assert( 0 );
        }
        mass = EvtPDL::getMass( p->getId() );
        //Just need to check that this mass is > than
        //the mass of all daughters
        double massSum = 0.;
        if ( p->getNDaug() ) {
            for ( i = 0; i < p->getNDaug(); i++ ) {
                massSum += p->getDaug( i )->mass();
            }
        }
        //some special cases are handled with 0 (stable) or 1 (k0->ks/kl) daughters
        if ( p->getNDaug() < 2 )
            massOk = true;
        if ( p->getParent() ) {
            if ( p->getParent()->getNDaug() == 1 )
                massOk = true;
        }
        if ( !massOk ) {
            if ( massSum < mass )
                massOk = true;
            if ( mass > maxOkMass )
                massOk = false;
        }
    }

    p->setMass( mass );
}

void EvtDecayBase::findMasses( EvtParticle* p, int ndaugs, EvtId daugs[10],
                               double masses[10] )
{
    int i;
    double mass_sum;

    int count = 0;

    if ( !( p->firstornot() ) ) {
        for ( i = 0; i < ndaugs; i++ ) {
            masses[i] = p->getDaug( i )->mass();
        }    //for
    }        //if
    else {
        p->setFirstOrNot();
        // if only one daughter do it

        if ( ndaugs == 1 ) {
            masses[0] = p->mass();
            return;
        }

        //until we get a combo whose masses are less than _parent mass.
        do {
            mass_sum = 0.0;

            for ( i = 0; i < ndaugs; i++ ) {
                masses[i] = EvtPDL::getMass( daugs[i] );
                mass_sum = mass_sum + masses[i];
            }

            count++;

            if ( count == 10000 ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Decaying particle:" << EvtPDL::name( p->getId() ).c_str()
                    << " (m=" << p->mass() << ")" << endl;
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "To the following daugthers" << endl;
                for ( i = 0; i < ndaugs; i++ ) {
                    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                        << EvtPDL::name( daugs[i] ).c_str() << endl;
                }
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Has been rejected " << count
                    << " times, will now take minimal masses "
                    << " of daugthers" << endl;

                mass_sum = 0.;
                for ( i = 0; i < ndaugs; i++ ) {
                    masses[i] = EvtPDL::getMinMass( daugs[i] );
                    mass_sum = mass_sum + masses[i];
                }
                if ( mass_sum > p->mass() ) {
                    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                        << "Parent mass=" << p->mass()
                        << "to light for daugthers." << endl
                        << "Will throw the event away." << endl;
                    //dont terminate - start over on the event.
                    EvtStatus::setRejectFlag();
                    mass_sum = 0.;
                    //	  ::abort();
                }
            }
        } while ( mass_sum > p->mass() );
    }    //else

    return;
}

void EvtDecayBase::checkNArg( int a1, int a2, int a3, int a4 )
{
    if ( _narg != a1 && _narg != a2 && _narg != a3 && _narg != a4 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << _modelname.c_str() << " generator expected " << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << a1 << endl;
        ;
        if ( a2 > -1 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << " or " << a2 << endl;
        }
        if ( a3 > -1 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << " or " << a3 << endl;
        }
        if ( a4 > -1 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << " or " << a4 << endl;
        }
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " arguments but found:" << _narg << endl;
        printSummary();
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}
void EvtDecayBase::checkNDaug( int d1, int d2 )
{
    if ( _ndaug != d1 && _ndaug != d2 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << _modelname.c_str() << " generator expected ";
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << d1;
        if ( d2 > -1 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << " or " << d2;
        }
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " daughters but found:" << _ndaug << endl;
        printSummary();
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}

void EvtDecayBase::checkSpinParent( EvtSpinType::spintype sp )
{
    EvtSpinType::spintype parenttype = EvtPDL::getSpinType( getParentId() );
    if ( parenttype != sp ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << _modelname.c_str() << " did not get the correct parent spin\n";
        printSummary();
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}

void EvtDecayBase::checkSpinDaughter( int d1, EvtSpinType::spintype sp )
{
    EvtSpinType::spintype parenttype = EvtPDL::getSpinType( getDaug( d1 ) );
    if ( parenttype != sp ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << _modelname.c_str()
            << " did not get the correct daughter spin d=" << d1 << endl;
        printSummary();
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}

double* EvtDecayBase::getArgs()
{
    if ( !_argsD.empty() )
        return _argsD.data();
    //The user has asked for a list of doubles - the arguments
    //better all be doubles...
    if ( _narg == 0 )
        return _argsD.data();

    _argsD.resize( _narg );
    for ( int i = 0; i < _narg; i++ ) {
        char* tc;
        _argsD[i] = strtod( _args[i].c_str(), &tc );
    }
    return _argsD.data();
}

double EvtDecayBase::getArg( unsigned int j )
{
    // Verify string

    if ( getParentId().getId() == 25 ) {
        int i = 0;
        ++i;
    }

    const char* str = _args[j].c_str();
    int i = 0;
    while ( str[i] != 0 ) {
        if ( isalpha( str[i] ) && str[i] != 'e' ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "String " << str << " is not a number" << endl;
            assert( 0 );
        }
        i++;
    }

    char** tc = 0;
    double result = strtod( _args[j].c_str(), tc );

    if ( _storedArgs.size() < j + 1 ) {    // then store the argument's value
        _storedArgs.push_back( result );
    }

    return result;
}

bool EvtDecayBase::matchingDecay( const EvtDecayBase& other ) const
{
    if ( _ndaug != other._ndaug )
        return false;
    if ( _parent != other._parent )
        return false;

    std::vector<int> useDs;
    for ( int i = 0; i < _ndaug; i++ )
        useDs.push_back( 0 );

    for ( int i = 0; i < _ndaug; i++ ) {
        bool foundIt = false;
        for ( int j = 0; j < _ndaug; j++ ) {
            if ( useDs[j] == 1 )
                continue;
            if ( _daug[i] == other._daug[j] &&
                 _daug[i].getAlias() == other._daug[j].getAlias() ) {
                foundIt = true;
                useDs[j] = 1;
                break;
            }
        }
        if ( foundIt == false )
            return false;
    }
    for ( int i = 0; i < _ndaug; i++ )
        if ( useDs[i] == 0 )
            return false;

    return true;
}
