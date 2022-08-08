
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

#ifndef EVT_INTERVAL_DECAY_AMP
#define EVT_INTERVAL_DECAY_AMP

#define VERBOSE true
#include "EvtGenBase/EvtAmpFactory.hh"
#include "EvtGenBase/EvtAmpPdf.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtMacros.hh"
#include "EvtGenBase/EvtMultiChannelParser.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPdf.hh"
#include "EvtGenBase/EvtReport.hh"

#include <iostream>
#include <string>
#include <vector>

// Decay model that uses the "amplitude on an interval"
// templatization

template <class T>
class EvtIntervalDecayAmp : public EvtDecayAmp {
  public:
    EvtIntervalDecayAmp() : _probMax( 0. ), _nScan( 0 ), _fact( 0 ) {}

    EvtIntervalDecayAmp( const EvtIntervalDecayAmp<T>& other ) :
        _probMax( other._probMax ), _nScan( other._nScan ), COPY_PTR( _fact )
    {
    }

    virtual ~EvtIntervalDecayAmp() { delete _fact; }

    // Initialize model

    void init() override
    {
        // Collect model parameters and parse them

        vector<std::string> args;
        int i;
        for ( i = 0; i < getNArg(); i++ )
            args.push_back( getArgStr( i ) );
        EvtMultiChannelParser parser;
        parser.parse( args );

        // Create factory and interval

        if ( VERBOSE )
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Create factory and interval" << std::endl;
        _fact = createFactory( parser );

        // Maximum PDF value over the Dalitz plot can be specified, or a scan
        // can be performed.

        _probMax = parser.pdfMax();
        _nScan = parser.nScan();
        if ( VERBOSE )
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Pdf maximum " << _probMax << std::endl;
        if ( VERBOSE )
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Scan number " << _nScan << std::endl;
    }

    void initProbMax() override
    {
        if ( 0 == _nScan ) {
            if ( _probMax > 0 )
                setProbMax( _probMax );
            else
                assert( 0 );
        } else {
            double factor = 1.2;    // increase maximum probability by 20%
            EvtAmpPdf<T> pdf( *_fact->getAmp() );
            EvtPdfSum<T>* pc = _fact->getPC();
            EvtPdfDiv<T> pdfdiv( pdf, *pc );
            printf( "Sampling %d points to find maximum\n", _nScan );
            EvtPdfMax<T> x = pdfdiv.findMax( *pc, _nScan );
            _probMax = factor * x.value();
            printf( "Found maximum %f\n", x.value() );
            printf( "Increase to   %f\n", _probMax );
            setProbMax( _probMax );
        }
    }

    void decay( EvtParticle* p ) override
    {
        // Set things up in most general way

        static EvtId B0 = EvtPDL::getId( "B0" );
        static EvtId B0B = EvtPDL::getId( "anti-B0" );
        double t;
        EvtId other_b;
        EvtComplex ampl( 0., 0. );

        // Sample using pole-compensator pdf

        EvtPdfSum<T>* pc = getPC();
        _x = pc->randomPoint();

        if ( _fact->isCPModel() ) {
            // Time-dependent Dalitz plot changes
            // Dec 2005 (ddujmic@slac.stanford.edu)

            EvtComplex A = _fact->getAmp()->evaluate( _x );
            EvtComplex Abar = _fact->getAmpConj()->evaluate( _x );

            EvtCPUtil::getInstance()->OtherB( p, t, other_b );

            double dm = _fact->dm();
            double mixAmpli = _fact->mixAmpli();
            double mixPhase = _fact->mixPhase();
            EvtComplex qoverp( cos( mixPhase ) * mixAmpli,
                               sin( mixPhase ) * mixAmpli );
            EvtComplex poverq( cos( mixPhase ) / mixAmpli,
                               -sin( mixPhase ) / mixAmpli );

            if ( other_b == B0B )
                ampl = A * cos( dm * t / ( 2 * EvtConst::c ) ) +
                       EvtComplex( 0., 1. ) * Abar *
                           sin( dm * t / ( 2 * EvtConst::c ) ) * qoverp;
            if ( other_b == B0 )
                ampl = Abar * cos( dm * t / ( 2 * EvtConst::c ) ) +
                       EvtComplex( 0., 1. ) * A *
                           sin( dm * t / ( 2 * EvtConst::c ) ) * poverq;

        } else {
            ampl = amplNonCP( _x );
        }

        // Pole-compensate

        double comp = sqrt( pc->evaluate( _x ) );
        assert( comp > 0 );
        vertex( ampl / comp );

        // Now generate random angles, rotate and setup
        // the daughters

        std::vector<EvtVector4R> v = initDaughters( _x );

        size_t N = p->getNDaug();
        if ( v.size() != N ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Number of daughters " << N << std::endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Momentum vector size " << v.size() << std::endl;
            assert( 0 );
        }

        for ( size_t i = 0; i < N; i++ ) {
            p->getDaug( i )->init( getDaugs()[i], v[i] );
        }
    }

    virtual EvtAmpFactory<T>* createFactory(
        const EvtMultiChannelParser& parser ) = 0;
    virtual std::vector<EvtVector4R> initDaughters( const T& p ) const = 0;

    // provide access to the decay point and to the amplitude of any decay point.
    // this is used by EvtBtoKD3P:
    const T& x() const { return _x; }
    EvtComplex amplNonCP( const T& x )
    {
        return _fact->getAmp()->evaluate( x );
    }
    EvtPdfSum<T>* getPC() { return _fact->getPC(); }

  protected:
    double _probMax;    // Maximum probability
    int _nScan;         // Number of points for max prob DP scan
    T _x;               // Decay point

    EvtAmpFactory<T>* _fact;    // factory
};

#endif
