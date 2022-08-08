
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

#ifndef EVTCPUTIL_HH
#define EVTCPUTIL_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPatches.hh"
class EvtParticle;
class EvtId;

class EvtCPUtil {
  public:
    EvtCPUtil( int mixingType );
    ~EvtCPUtil();

    enum MixingType
    {
        Coherent = 0,
        Incoherent = 1
    };

    static EvtCPUtil* getInstance();

    void setMixingType( int mixingType ) { _mixingType = mixingType; }
    int getMixingType() { return _mixingType; }

    void fractB0CP( EvtComplex Af, EvtComplex Abarf, double deltam, double beta,
                    double& fract );

    void fractB0nonCP( EvtComplex Af, EvtComplex Abarf, EvtComplex Afbar,
                       EvtComplex Abarfbar, double deltam, double beta,
                       int flip, double& fract );

    // Mark Whitehead 7/12/2009
    // Add required lines from EvtIncoherentMixing.hh to fix CPV

    // Functions to check if a B has mixed (comes from a B)
    bool isB0Mixed( EvtParticle* );
    bool isBsMixed( EvtParticle* );

    bool flipIsEnabled();
    void enableFlip();
    void disableFlip();

    void OtherB( EvtParticle* p, double& t, EvtId& otherb );

    void OtherCoherentB( EvtParticle* p, double& t, EvtId& otherb, double probB0 );
    void OtherIncoherentB( EvtParticle* p, double& t, EvtId& otherb,
                           double probB0 );

    void OtherB( EvtParticle* p, double& t, EvtId& otherb, double probB0 );

    //id is the produced particle
    //t returns the lifetime of the particle
    //and mix will be 1 if it mixed otherwise 0
    void incoherentMix( const EvtId id, double& t, int& mix );

    double getDeltaGamma( const EvtId id );
    double getDeltaM( const EvtId id );

  private:
    bool _enableFlip;
    int _mixingType;
};

#endif
