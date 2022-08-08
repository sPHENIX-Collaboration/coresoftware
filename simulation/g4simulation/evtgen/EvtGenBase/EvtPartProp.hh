
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

#ifndef EVTPARTPROP_HH
#define EVTPARTPROP_HH

#include "EvtGenBase/EvtAbsLineShape.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtSpinType.hh"

#include <memory>
#include <string>

class EvtPartProp {
  public:
    EvtPartProp();
    EvtPartProp( const EvtPartProp& x );

    double getMass() { return _lineShape->getMass(); }
    double getMassMin() { return _lineShape->getMassMin(); }
    double getMassMax() { return _lineShape->getMassMax(); }
    double getMaxRange() { return _lineShape->getMaxRange(); }
    double getWidth() { return _lineShape->getWidth(); }

    double getRandMass( EvtId* parId, int nDaug, EvtId* dauId, EvtId* othDauId,
                        double maxMass, double* dauMasses )
    {
        return _lineShape->getRandMass( parId, nDaug, dauId, othDauId, maxMass,
                                        dauMasses );
    }
    double getMassProb( double mass, double massPar, int nDaug, double* massDau )
    {
        return _lineShape->getMassProb( mass, massPar, nDaug, massDau );
    }

    double getctau() { return _ctau; }
    void setctau( double tau ) { _ctau = tau; }

    int getChg3() { return _chg3; }
    void setChg3( int c3 ) { _chg3 = c3; }

    EvtSpinType::spintype getSpinType() { return _spintype; }
    void setSpinType( EvtSpinType::spintype stype ) { _spintype = stype; }

    const std::string& getName() { return _name; }
    void setName( std::string pname );

    EvtId getId() { return _id; }
    void setId( EvtId id ) { _id = id; }

    EvtId getIdChgConj() { return _idchgconj; }
    void setIdChgConj( EvtId idchgconj ) { _idchgconj = idchgconj; }

    int getStdHep() { return _stdhep; }
    void setStdHep( int stdhep ) { _stdhep = stdhep; }

    int getLundKC() { return _lundkc; }
    void setLundKC( int lundkc ) { _lundkc = lundkc; }

    EvtAbsLineShape* getLineShape() { return _lineShape.get(); }
    void initLineShape( double mass, double width, double maxRange );
    //  void initLineShape(double mass, double width, double maxRange, double mDaug1, double mDaug2, int l);

    // setLineShape takes ownership of l
    void setLineShape( EvtAbsLineShape* l ) { _lineShape.reset( l ); }
    double rollMass() { return _lineShape->rollMass(); }

    EvtPartProp& operator=( const EvtPartProp& x );

    void reSetMass( double mass );
    void reSetWidth( double width );

    void reSetMassMin( double mass );
    void reSetMassMax( double mass );
    void reSetBlatt( double blatt );
    void reSetBlattBirth( double blatt );
    void includeBirthFactor( bool yesno );
    void includeDecayFactor( bool yesno );
    void newLineShape( std::string type );
    void setPWForDecay( int spin, EvtId d1, EvtId d2 );
    void setPWForBirthL( int spin, EvtId par, EvtId othD );

  private:
    std::unique_ptr<EvtAbsLineShape> _lineShape;

    double _ctau;
    EvtId _id;
    EvtId _idchgconj;
    EvtSpinType::spintype _spintype;
    int _chg3;
    int _stdhep;
    int _lundkc;
    std::string _name;
};

#endif
