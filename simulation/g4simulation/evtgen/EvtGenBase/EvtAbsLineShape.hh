
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

#ifndef EVTABSLINESHAPE_HH
#define EVTABSLINESHAPE_HH

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtSpinType.hh"

#include <vector>

class EvtId;

class EvtAbsLineShape {
  public:
    EvtAbsLineShape() = default;
    EvtAbsLineShape( double mass, double width, double maxRange,
                     EvtSpinType::spintype sp );
    virtual ~EvtAbsLineShape() = default;
    EvtAbsLineShape& operator=( const EvtAbsLineShape& x );
    EvtAbsLineShape( const EvtAbsLineShape& x );

    double getMass() { return _mass; }
    double getMassMin() { return _massMin; }
    double getMassMax() { return _massMax; }
    double getMaxRange() { return _maxRange; }
    double getWidth() { return _width; }
    EvtSpinType::spintype getSpinType() { return _spin; }
    virtual double rollMass();
    virtual EvtAbsLineShape* clone();

    void reSetMass( double mass ) { _mass = mass; }
    void reSetWidth( double width ) { _width = width; }
    void reSetMassMin( double mass ) { _massMin = mass; }
    void reSetMassMax( double mass ) { _massMax = mass; }
    virtual void reSetBlatt( double /*blatt*/ ){};
    virtual void reSetBlattBirth( double /*blatt*/ ){};
    void includeBirthFactor( bool yesno ) { _includeBirthFact = yesno; }
    void includeDecayFactor( bool yesno ) { _includeDecayFact = yesno; }
    void setPWForDecay( int spin, EvtId d1, EvtId d2 )
    {
        _userSetPW.push_back( spin );
        _userSetPWD1.push_back( d1 );
        _userSetPWD2.push_back( d2 );
    }
    void setPWForBirthL( int spin, EvtId par, EvtId othD )
    {
        _userSetBirthPW.push_back( spin );
        _userSetBirthOthD.push_back( othD );
        _userSetBirthPar.push_back( par );
    }

    virtual double getRandMass( EvtId* parId, int nDaug, EvtId* dauId,
                                EvtId* othDaugId, double maxMass,
                                double* dauMasses );
    virtual double getMassProb( double mass, double massPar, int nDaug,
                                double* massDau );

  protected:
    bool _includeDecayFact;
    bool _includeBirthFact;
    double _mass;
    double _massMin;
    double _massMax;
    double _width;
    double _maxRange;

    // allow for special cases where the default method of picking the
    //lowest allowed partial wave for a decay is not the right answer.
    // string is "<spin> <daughter1> <daughter2>"
    //new 9/12/2003 Lange
    std::vector<EvtId> _userSetPWD1, _userSetPWD2;
    std::vector<int> _userSetPW;

    // also do it for birth factors
    std::vector<EvtId> _userSetBirthPar, _userSetBirthOthD;
    std::vector<int> _userSetBirthPW;

    EvtSpinType::spintype _spin;
};

#endif
