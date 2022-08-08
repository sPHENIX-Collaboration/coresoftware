
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

#ifdef EVTGEN_TAUOLA

#include "EvtGenExternal/EvtTauolaEngine.hh"

#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSymTable.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "Tauola/Log.h"
#include "Tauola/Tauola.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

using std::endl;

EvtTauolaEngine::EvtTauolaEngine( bool useEvtGenRandom )
{
    // PDG standard code integer ID for tau particle
    _tauPDG = 15;
    // Number of possible decay modes in Tauola
    _nTauolaModes = 22;

    EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Setting up TAUOLA." << endl;

    // These three lines are not really necessary since they are the default.
    // But they are here so that we know what the initial conditions are.
    Tauolapp::Tauola::setDecayingParticle( _tauPDG );    // tau PDG code
    Tauolapp::Tauola::setSameParticleDecayMode(
        Tauolapp::Tauola::All );    // all modes allowed
    Tauolapp::Tauola::setOppositeParticleDecayMode(
        Tauolapp::Tauola::All );    // all modes allowed

    // Limit the number of warnings printed out. Can't choose zero here, unfortunately
    Tauolapp::Log::SetWarningLimit( 1 );

    // Initial the Tauola external generator
    if ( useEvtGenRandom == true ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Using EvtGen random number engine also for Tauola++" << endl;

        Tauolapp::Tauola::setRandomGenerator( EvtRandom::Flat );
    }

    // Use the BaBar-tuned chiral current calculations by default. Can be changed using the
    // TauolaCurrentOption keyword in decay files
    Tauolapp::Tauola::setNewCurrents( 1 );

    Tauolapp::Tauola::initialize();

    // Initialise various default parameters
    // Neutral and charged spin propagator choices
    _neutPropType = 0;
    _posPropType = 0;
    _negPropType = 0;

    // Set-up possible decay modes _after_ we have read the (user) decay file
    _initialised = false;
}

void EvtTauolaEngine::initialise()
{
    // Set up all possible tau decay modes.
    // This should be done just before the first doDecay() call,
    // since we want to make sure that any decay.dec files are processed
    // first to get lists of particle modes and their alias definitions
    // (for creating EvtParticles with the right history information).

    if ( _initialised == false ) {
        this->setUpPossibleTauModes();
        this->setOtherParameters();

        _initialised = true;
    }
}

void EvtTauolaEngine::setUpPossibleTauModes()
{
    // Get the decay table list defined by the decay.dec files.
    // Only look for the first tau particle decay mode definitions with the Tauola name,
    // since that generator only allows the same BFs for both tau+ and tau- decays.
    // We can not choose a specific tau decay event-by-event, since this is
    // only possible before we call Tauola::initialize().
    // Otherwise, we could have selected a random mode ourselves for tau- and tau+
    // separately (via selecting a random number and comparing it to be less than
    // the cumulative BF) for each event.

    int nPDL = EvtPDL::entries();
    int iPDL( 0 );

    bool gotAnyTauolaModes( false );

    for ( iPDL = 0; iPDL < nPDL; iPDL++ ) {
        EvtId particleId = EvtPDL::getEntry( iPDL );
        int PDGId = EvtPDL::getStdHep( particleId );

        if ( abs( PDGId ) == _tauPDG && gotAnyTauolaModes == false ) {
            int aliasInt = particleId.getAlias();

            // Get the list of decay modes for this tau particle (alias)
            int nModes = EvtDecayTable::getInstance()->getNModes( aliasInt );
            int iMode( 0 ), iTauMode( 0 );

            // Vector to store tau mode branching fractions.
            // The size of this vector equals the total number of possible
            // Tauola decay modes. Initialise all BFs to zero.
            std::vector<double> tauolaModeBFs( _nTauolaModes );

            for ( iTauMode = 0; iTauMode < _nTauolaModes; iTauMode++ ) {
                tauolaModeBFs[iTauMode] = 0.0;
            }

            double totalTauModeBF( 0.0 );

            int nNonTauolaModes( 0 );

            // Loop through each decay mode
            for ( iMode = 0; iMode < nModes; iMode++ ) {
                EvtDecayBase* decayModel =
                    EvtDecayTable::getInstance()->findDecayModel( aliasInt,
                                                                  iMode );
                if ( decayModel ) {
                    // Check that the decay model name matches TAUOLA
                    std::string modelName = decayModel->getName();
                    if ( modelName == "TAUOLA" ) {
                        if ( gotAnyTauolaModes == false ) {
                            gotAnyTauolaModes = true;
                        }

                        // Extract the decay mode integer type and branching fraction
                        double BF = decayModel->getBranchingFraction();
                        int modeArrayInt = this->getModeInt( decayModel ) - 1;

                        if ( modeArrayInt >= 0 && modeArrayInt < _nTauolaModes ) {
                            tauolaModeBFs[modeArrayInt] = BF;
                            totalTauModeBF += BF;
                        }

                    } else {
                        nNonTauolaModes++;
                    }

                }    // Decay mode exists

            }    // Loop over decay models

            if ( gotAnyTauolaModes == true && nNonTauolaModes > 0 ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Please remove all non-TAUOLA decay modes for particle "
                    << EvtPDL::name( particleId ) << endl;
                ::abort();
            }

            // Normalise all (non-zero) tau mode BFs to sum up to 1.0, and
            // let Tauola know about these normalised branching fractions
            if ( totalTauModeBF > 0.0 ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" )
                    << "Setting TAUOLA BF modes using the definitions for the particle "
                    << EvtPDL::name( particleId ) << endl;

                for ( iTauMode = 0; iTauMode < _nTauolaModes; iTauMode++ ) {
                    tauolaModeBFs[iTauMode] /= totalTauModeBF;
                    double modeBF = tauolaModeBFs[iTauMode];
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << "Setting TAUOLA BF for mode " << iTauMode + 1
                        << " = " << modeBF << endl;
                    Tauolapp::Tauola::setTauBr( iTauMode + 1, modeBF );
                }

                EvtGenReport( EVTGEN_INFO, "EvtGen" )
                    << "Any other TAUOLA BF modes for other tau particle decay mode definitions will be ignored!"
                    << endl;
            }

        }    // Got tau particle and have yet to get a TAUOLA mode

    }    // Loop through PDL entries
}

int EvtTauolaEngine::getModeInt( EvtDecayBase* decayModel )
{
    int modeInt( 0 );

    if ( decayModel ) {
        int nVars = decayModel->getNArg();

        if ( nVars > 0 ) {
            modeInt = static_cast<int>( decayModel->getArg( 0 ) );
        }
    }

    return modeInt;
}

void EvtTauolaEngine::setOtherParameters()
{
    // Set other Tauola parameters using the "Defined" keyword in the decay file. If any of
    // these are not found in the decay file, then default values are assumed/kept

    // 1) TauolaNeutralProp: Specify the neutral propagator type used for spin matrix calculations
    // "Z" (default), "Gamma", "Higgs" (H0), "PseudoHiggs" (A0), "MixedHiggs" (A0/H0)
    int iErr( 0 );
    std::string neutPropName = EvtSymTable::get( "TauolaNeutralProp", iErr );
    if ( neutPropName == "Z0" || neutPropName == "Z" ) {
        _neutPropType = Tauolapp::TauolaParticle::Z0;
    } else if ( neutPropName == "Gamma" ) {
        _neutPropType = Tauolapp::TauolaParticle::GAMMA;
    } else if ( neutPropName == "Higgs" ) {
        _neutPropType = Tauolapp::TauolaParticle::HIGGS;
    } else if ( neutPropName == "PseudoHiggs" ) {
        _neutPropType = Tauolapp::TauolaParticle::HIGGS_A;
    } else if ( neutPropName == "MixedHiggs" ) {
        _neutPropType = Tauolapp::Tauola::getHiggsScalarPseudoscalarPDG();
    }

    if ( _neutPropType != 0 ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "TAUOLA neutral spin propagator PDG id set to " << _neutPropType
            << endl;
    }

    // 2) TauolaChargedProp: Specify the charged propagator type used for spin matrix calculations
    // "W" (default), "Higgs" (H+/H-)
    std::string chargedPropName = EvtSymTable::get( "TauolaChargedProp", iErr );
    if ( chargedPropName == "W" ) {
        _negPropType = Tauolapp::TauolaParticle::W_MINUS;
        _posPropType = Tauolapp::TauolaParticle::W_PLUS;
    } else if ( chargedPropName == "Higgs" ) {
        _negPropType = Tauolapp::TauolaParticle::HIGGS_MINUS;
        _posPropType = Tauolapp::TauolaParticle::HIGGS_PLUS;
    }

    if ( _negPropType != 0 ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "TAUOLA negative charge spin propagator PDG id set to "
            << _negPropType << endl;
    }

    if ( _posPropType != 0 ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "TAUOLA positive charge spin propagator PDG id set to "
            << _posPropType << endl;
    }

    // 3) TauolaHiggsMixingAngle: Specify the mixing angle between the neutral scalar & pseudoscalar Higgs
    // A0/H0; the default mixing angle is pi/4 radians
    std::string mixString = EvtSymTable::get( "TauolaHiggsMixingAngle", iErr );
    // If the definition name is not found, get() just returns the first argument string
    if ( mixString != "TauolaHiggsMixingAngle" ) {
        double mixAngle = std::atof( mixString.c_str() );
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "TAUOLA Higgs mixing angle set to " << mixAngle << " radians"
            << endl;
        Tauolapp::Tauola::setHiggsScalarPseudoscalarMixingAngle( mixAngle );
    }

    // 4) TauolaBRi, where i = 1,2,3,4: Redefine sub-channel branching fractions using the setTaukle
    // function, after initialized() has been called. Default values = 0.5, 0.5, 0.5 and 0.6667
    int j( 1 );
    std::vector<double> BRVect;
    BRVect.push_back( 0.5 );
    BRVect.push_back( 0.5 );
    BRVect.push_back( 0.5 );
    BRVect.push_back( 0.6667 );

    for ( j = 1; j < 5; j++ ) {
        std::ostringstream o;
        o << j;
        std::string BRName = "TauolaBR" + o.str();
        std::string stringBR = EvtSymTable::get( BRName, iErr );

        // If the definition name is not found, get() just returns the first argument string
        if ( stringBR != BRName ) {
            BRVect[j - 1] = std::atof( stringBR.c_str() );
        }
    }

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "TAUOLA::setTaukle values are " << BRVect[0] << ", " << BRVect[1]
        << ", " << BRVect[2] << ", " << BRVect[3] << endl;

    Tauolapp::Tauola::setTaukle( BRVect[0], BRVect[1], BRVect[2], BRVect[3] );

    // 5) Specify the hadronic current option, e.g. orig CLEO = 0, BaBar-tuned = 1 (default), ...
    // No check is made by EvtGen on valid integer options - its just passed to Tauola
    std::string currentOption = EvtSymTable::get( "TauolaCurrentOption", iErr );
    // If the definition name is not found, get() just returns the first argument string
    if ( currentOption != "TauolaCurrentOption" ) {
        int currentOpt = std::atoi( currentOption.c_str() );
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "TAUOLA current option = " << currentOpt << endl;

        Tauolapp::Tauola::setNewCurrents( currentOpt );
    }
}

bool EvtTauolaEngine::doDecay( EvtParticle* tauParticle )
{
    if ( _initialised == false ) {
        this->initialise();
    }

    if ( tauParticle == 0 ) {
        return false;
    }

    // Check that we have a tau particle.
    EvtId partId = tauParticle->getId();
    if ( abs( EvtPDL::getStdHep( partId ) ) != _tauPDG ) {
        return false;
    }

    int nTauDaug = tauParticle->getNDaug();

    // If the number of tau daughters is not zero, then we have already decayed
    // it using Tauola/another decay algorithm.
    if ( nTauDaug > 0 ) {
        return true;
    }

    this->decayTauEvent( tauParticle );

    return true;
}

void EvtTauolaEngine::decayTauEvent( EvtParticle* tauParticle )
{
    // Either we have a tau particle within a decay chain, or a single particle.
    // Create a dummy HepMC event & vertex for the parent particle, containing the tau as
    // one of the outgoing particles. If we have a decay chain, the parent will be the
    // incoming particle, while the daughters, including the tau, are outgoing particles.
    // For the single particle case, the incoming particle is null, while the single tau
    // is the only outgoing particle.
    // We can then pass this event to Tauola which should then decay the tau particle.
    // We also consider all other tau particles from the parent decay in the logic below.

    // Create the dummy event.
    auto theEvent = std::make_unique<GenEvent>( Units::GEV, Units::MM );

    // Create the decay "vertex".
    GenVertexPtr theVertex = newGenVertexPtr();
    theEvent->add_vertex( theVertex );

    // Get the parent of this tau particle
    EvtParticle* theParent = tauParticle->getParent();
    GenParticlePtr hepMCParent( 0 );

    // Assign the parent particle as the incoming particle to the vertex.
    if ( theParent ) {
        hepMCParent = this->createGenParticle( theParent );
        theVertex->add_particle_in( hepMCParent );
    } else {
        // The tau particle has no parent. Set "itself" as the incoming particle for the first vertex.
        // This is needed, otherwise Tauola warns of momentum non-conservation for this (1st) vertex.
        GenParticlePtr tauGenInit = this->createGenParticle( tauParticle );
        theVertex->add_particle_in( tauGenInit );
    }

    // Find all daughter particles and assign them as outgoing particles to the vertex.
    // This will include the tau particle we are currently processing.
    // If the parent decay has more than one tau particle, we need to include them as well.
    // This is important since Tauola needs the correct physics correlations: we do not
    // want Tauola to decay each particle separately if they are from tau pair combinations.
    // Tauola will process the event, and we will create EvtParticles from all tau decay
    // products, i.e. the tau particle we currently have and any other tau particles.
    // EvtGen will try to decay the other tau particle(s) by calling EvtTauola and therefore
    // this function. However, we check to see if the tau candidate has any daughters already.
    // If it does, then we have already set the tau decay products from Tauola.

    // Map to store (HepMC,EvtParticle) pairs for each tau candidate from the parent
    // decay. This is needed to find out what EvtParticle corresponds to a given tau HepMC
    // candidate: we do not want to recreate existing EvtParticle pointers.
    std::map<GenParticlePtr, EvtParticle*> tauMap;

    // Keep track of the original EvtId of the parent particle, since we may need to set
    // the equivalent HepMCParticle has a gauge boson to let Tauola calculate spin effects
    EvtId origParentId( -1, -1 );

    if ( theParent ) {
        // Original parent id
        origParentId = EvtPDL::getId( theParent->getName() );

        // Find all tau particles in the decay tree and store them in the map.
        // Keep track of how many tau daughters this parent particle has
        int nTaus( 0 );
        int nDaug( theParent->getNDaug() );
        int iDaug( 0 );

        for ( iDaug = 0; iDaug < nDaug; iDaug++ ) {
            EvtParticle* theDaughter = theParent->getDaug( iDaug );

            if ( theDaughter ) {
                GenParticlePtr hepMCDaughter = this->createGenParticle(
                    theDaughter );
                theVertex->add_particle_out( hepMCDaughter );

                EvtId theId = theDaughter->getId();
                int PDGInt = EvtPDL::getStdHep( theId );

                if ( abs( PDGInt ) == _tauPDG ) {
                    // Delete any siblings for the tau particle
                    if ( theDaughter->getNDaug() > 0 ) {
                        theDaughter->deleteDaughters( false );
                    }
                    tauMap[hepMCDaughter] = theDaughter;
                    nTaus++;
                } else {
                    // Treat all other particles as "stable"
                    hepMCDaughter->set_status( Tauolapp::TauolaParticle::STABLE );
                }

            }    // theDaughter != 0
        }        // Loop over daughters

        // For the parent particle, artifically set the PDG to a boson with the same 4-momentum
        // so that spin correlations are calculated inside Tauola.
        // This leaves the original parent _EvtParticle_ unchanged
        if ( nTaus > 0 && hepMCParent ) {
            int parCharge = EvtPDL::chg3( origParentId ) /
                            3;    // (3*particle charge)/3 = particle charge
            if ( parCharge == 0 && _neutPropType != 0 ) {
                hepMCParent->set_pdg_id( _neutPropType );
            } else if ( parCharge == -1 && _negPropType != 0 ) {
                hepMCParent->set_pdg_id( _negPropType );
            } else if ( parCharge == 1 && _posPropType != 0 ) {
                hepMCParent->set_pdg_id( _posPropType );
            }
        }

    } else {
        // We only have the one tau particle. Store only this in the map.
        GenParticlePtr singleTau = this->createGenParticle( tauParticle );
        theVertex->add_particle_out( singleTau );
        tauMap[singleTau] = tauParticle;
    }

    // Now pass the event to Tauola for processing
    // Create a Tauola event object
#ifdef EVTGEN_HEPMC3
    Tauolapp::TauolaHepMC3Event tauolaEvent( theEvent.get() );
#else
    Tauolapp::TauolaHepMCEvent tauolaEvent( theEvent.get() );
#endif

    // Run the Tauola algorithm
    tauolaEvent.decayTaus();

    // Loop over all tau particles in the HepMC event and create their EvtParticle daughters.
    // Store all final "stable" descendent particles as the tau daughters, i.e.
    // let Tauola decay any resonances such as a_1 or rho.
    // If there is more than one tau particle in the event, then also create the
    // corresponding EvtParticles for their daughters as well. They will not be
    // re-decayed since we check at the start of this function if the tau particle has
    // any daughters before running Tauola decayTaus().

#ifdef EVTGEN_HEPMC3
    for ( auto aParticle : theEvent->particles() ) {
#else
    HepMC::GenEvent::particle_iterator eventIter;
    for ( eventIter = theEvent->particles_begin();
          eventIter != theEvent->particles_end(); ++eventIter ) {
        // Check to see if we have a tau particle
        HepMC::GenParticle* aParticle = ( *eventIter );
#endif

        if ( aParticle && abs( aParticle->pdg_id() ) == _tauPDG ) {
            // Find out what EvtParticle corresponds to the HepMC particle.
            // We need this to create and attach EvtParticle daughters.
            EvtParticle* tauEvtParticle = tauMap[aParticle];

            if ( tauEvtParticle ) {
                // Get the tau 4-momentum in the lab (first mother) frame. We need to boost
                // all the tau daughters to this frame, such that daug.getP4() is in the tau restframe.
                EvtVector4R tauP4CM = tauEvtParticle->getP4Lab();
                tauP4CM.set( tauP4CM.get( 0 ), -tauP4CM.get( 1 ),
                             -tauP4CM.get( 2 ), -tauP4CM.get( 3 ) );

                // Get the decay vertex for the tau particle
                GenVertexPtr endVertex = aParticle->end_vertex();

                std::vector<EvtId> daugIdVect;
                std::vector<EvtVector4R> daugP4Vect;

                // Loop through all descendants
#ifdef EVTGEN_HEPMC3
                for ( auto tauDaug :
                      HepMC3::Relatives::DESCENDANTS( endVertex ) ) {
#else
                HepMC::GenVertex::particle_iterator tauIter;
                // Loop through all descendants
                for ( tauIter = endVertex->particles_begin( HepMC::descendants );
                      tauIter != endVertex->particles_end( HepMC::descendants );
                      ++tauIter ) {
                    HepMC::GenParticle* tauDaug = ( *tauIter );
#endif
                    // Check to see if this descendant has its own decay vertex, e.g. rho resonance.
                    // If so, skip this daughter and continue looping through the descendant list
                    // until we reach the final "stable" products (e.g. pi pi from rho -> pi pi).
                    GenVertexPtr daugDecayVtx = tauDaug->end_vertex();
                    if ( daugDecayVtx ) {
                        continue;
                    }

                    // Store the particle id and 4-momentum
                    int tauDaugPDG = tauDaug->pdg_id();
                    EvtId daugId = EvtPDL::evtIdFromStdHep( tauDaugPDG );
                    daugIdVect.push_back( daugId );

                    FourVector tauDaugP4 = tauDaug->momentum();
                    double tauDaug_px = tauDaugP4.px();
                    double tauDaug_py = tauDaugP4.py();
                    double tauDaug_pz = tauDaugP4.pz();
                    double tauDaug_E = tauDaugP4.e();

                    EvtVector4R daugP4( tauDaug_E, tauDaug_px, tauDaug_py,
                                        tauDaug_pz );
                    daugP4Vect.push_back( daugP4 );

                }    // Loop over HepMC tau daughters

                // Create the tau EvtParticle daughters and assign their ids and 4-mtm
                int nDaug = daugIdVect.size();

                tauEvtParticle->makeDaughters( nDaug, daugIdVect );

                int iDaug( 0 );
                for ( iDaug = 0; iDaug < nDaug; iDaug++ ) {
                    EvtParticle* theDaugPart = tauEvtParticle->getDaug( iDaug );

                    if ( theDaugPart ) {
                        EvtId theDaugId = daugIdVect[iDaug];
                        EvtVector4R theDaugP4 = daugP4Vect[iDaug];
                        theDaugP4.applyBoostTo(
                            tauP4CM );    // Boost the 4-mtm to the tau rest frame
                        theDaugPart->init( theDaugId, theDaugP4 );
                    }

                }    // Loop over tau daughters
            }

        }    // We have a tau HepMC particle in the event
    }

    theEvent->clear();
}

GenParticlePtr EvtTauolaEngine::createGenParticle( EvtParticle* theParticle )
{
    // Method to create an HepMC::GenParticle version of the given EvtParticle.
    if ( theParticle == 0 ) {
        return 0;
    }

    // Get the 4-momentum (E, px, py, pz) for the EvtParticle
    EvtVector4R p4 = theParticle->getP4Lab();

    // Convert this to the HepMC 4-momentum
    double E = p4.get( 0 );
    double px = p4.get( 1 );
    double py = p4.get( 2 );
    double pz = p4.get( 3 );

    FourVector hepMC_p4( px, py, pz, E );

    int PDGInt = EvtPDL::getStdHep( theParticle->getId() );

    // Set the status flag for the particle.
    int status = Tauolapp::TauolaParticle::HISTORY;

    GenParticlePtr genParticle = newGenParticlePtr( hepMC_p4, PDGInt, status );

    return genParticle;
}

#endif
