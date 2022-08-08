
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

#ifndef EVTFOURBODYPHSP_HH
#define EVTFOURBODYPHSP_HH

#include <vector>
#include <array>
#include <utility>

#include "EvtGenBase/EvtDecayProb.hh"

class EvtParticle;

class EvtFourBodyPhsp : public EvtDecayProb {
  public:
    enum Shape
    {
        rectangle = 1,
        trapezoid = 2,
        pentagon = 3,
        variable = 4
    };

    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;

    void decay( EvtParticle* parent ) override;

  private:
    std::array<double, 4> phspFactor( const double mM, const double m12,
                                      const double m34,
                                      std::array<double, 4>& daughters ) const;
    Shape determineBoundaryShape( const double m12Min, const double m12Max,
                                  const double m34Max,
                                  const double mMother ) const;

    std::pair<double, double> generatePairMasses(
        const double m12Min, const double m12Max, const double m34Min,
        const double m34Max, const double mMother,
        const EvtFourBodyPhsp::Shape shape ) const;
    std::pair<double, double> generateRectangle( const double m12Min,
                                                 const double m12Max,
                                                 const double m34Min,
                                                 const double m34Max ) const;
    std::pair<double, double> generateTrapezoid( const double m12Min,
                                                 const double m12Max,
                                                 const double m34Min,
                                                 const double mMother ) const;

    std::array<double, 4> m_daughterMasses{ -1, -1, -1, -1 };

    double m_m12Min;
    double m_m12Max;
    double m_m34Min;
    double m_m34Max;

    double m_trapNorm;
    double m_trapCoeff1;
    double m_trapCoeff2;

    double m_pentagonSplit;
    double m_pentagonFraction;

    Shape m_boundaryShape;

    bool m_stableMother{true};
    bool m_stableDaughters{true};
    bool m_fixedBoundary{true};
};

#endif
