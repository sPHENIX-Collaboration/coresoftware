
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

#include "EvtGenModels/EvtBtoXsgammaRootFinder.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtItgSimpsonIntegrator.hh"
#include "EvtGenModels/EvtItgTwoCoeffFcn.hh"

#include <math.h>
using std::endl;

//-------------
// C Headers --
//-------------
extern "C" {
}

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

#define EVTITGROOTFINDER_MAXIT 100
#define EVTITGROOTFINDER_RELATIVEPRECISION 1.0e-16

double EvtBtoXsgammaRootFinder::GetRootSingleFunc(
    const EvtItgAbsFunction* theFunc, double functionValue, double lowerValue,
    double upperValue, double precision )
{
    // Use the bisection to find the root.
    // Iterates until find root to the accuracy of precision

    double xLower = 0.0, xUpper = 0.0;
    double root = 0;

    double f1 = theFunc->value( lowerValue ) - functionValue;
    double f2 = theFunc->value( upperValue ) - functionValue;

    if ( f1 * f2 > 0.0 ) {
        EvtGenReport( EVTGEN_WARNING, "EvtGen" )
            << "EvtBtoXsgammaRootFinder: No root in specified range !" << endl;
        return 0;
    }

    // Already have root
    if ( fabs( f1 ) < precision ) {
        root = lowerValue;
        return root;
    }
    if ( fabs( f2 ) < precision ) {
        root = upperValue;
        return root;
    }

    // Orient search so that f(xLower) < 0
    if ( f1 < 0.0 ) {
        xLower = lowerValue;
        xUpper = upperValue;
    } else {
        xLower = upperValue;
        xUpper = lowerValue;
    }

    double rootGuess = 0.5 * ( lowerValue + upperValue );
    double dxold = fabs( upperValue - lowerValue );
    double dx = dxold;

    double f = theFunc->value( rootGuess ) - functionValue;

    for ( int j = 0; j < EVTITGROOTFINDER_MAXIT; j++ ) {
        dxold = dx;
        dx = 0.5 * ( xUpper - xLower );
        rootGuess = xLower + dx;

        // If change in root is negligible, take it as solution.
        if ( fabs( xLower - rootGuess ) < precision ) {
            root = rootGuess;
            return root;
        }

        f = theFunc->value( rootGuess ) - functionValue;

        if ( f < 0.0 ) {
            xLower = rootGuess;
        } else {
            xUpper = rootGuess;
        }
    }

    EvtGenReport( EVTGEN_WARNING, "EvtGen" )
        << "EvtBtoXsgammaRootFinder: Maximum number of iterations "
        << "in EvtBtoXsgammaRootFinder::foundRoot exceeded!"
        << " Returning false." << endl;
    return 0;
}

double EvtBtoXsgammaRootFinder::GetGaussIntegFcnRoot(
    EvtItgAbsFunction* theFunc1, EvtItgAbsFunction* theFunc2,
    double integ1Precision, double integ2Precision, int maxLoop1, int maxLoop2,
    double integLower, double integUpper, double lowerValue, double upperValue,
    double precision )
{
    // Use the bisection to find the root.
    // Iterates until find root to the accuracy of precision

    //Need to work with integrators
    auto func1Integ = EvtItgSimpsonIntegrator{*theFunc1, integ1Precision,
                                              maxLoop1};
    auto func2Integ = EvtItgSimpsonIntegrator{*theFunc2, integ2Precision,
                                              maxLoop2};

    //coefficient 1 of the integrators is the root to be found
    //need to set this to lower value to start off with
    theFunc1->setCoeff( 1, 0, lowerValue );
    theFunc2->setCoeff( 1, 0, lowerValue );

    double f1 = func1Integ.evaluate( integLower, integUpper ) -
                theFunc2->getCoeff( 1, 2 ) *
                    func2Integ.evaluate( integLower, integUpper );
    theFunc1->setCoeff( 1, 0, upperValue );
    theFunc2->setCoeff( 1, 0, upperValue );
    double f2 = func1Integ.evaluate( integLower, integUpper ) -
                theFunc2->getCoeff( 1, 2 ) *
                    func2Integ.evaluate( integLower, integUpper );

    double xLower = 0.0, xUpper = 0.0;
    double root = 0;

    if ( f1 * f2 > 0.0 ) {
        EvtGenReport( EVTGEN_WARNING, "EvtGen" )
            << "EvtBtoXsgammaRootFinder: No root in specified range !" << endl;
        return false;
    }

    // Already have root
    if ( fabs( f1 ) < precision ) {
        root = lowerValue;
        return root;
    }
    if ( fabs( f2 ) < precision ) {
        root = upperValue;
        return root;
    }

    // Orient search so that f(xLower) < 0
    if ( f1 < 0.0 ) {
        xLower = lowerValue;
        xUpper = upperValue;
    } else {
        xLower = upperValue;
        xUpper = lowerValue;
    }

    double rootGuess = 0.5 * ( lowerValue + upperValue );
    double dxold = fabs( upperValue - lowerValue );
    double dx = dxold;

    theFunc1->setCoeff( 1, 0, rootGuess );
    theFunc2->setCoeff( 1, 0, rootGuess );
    double f = func1Integ.evaluate( integLower, integUpper ) -
               theFunc2->getCoeff( 1, 2 ) *
                   func2Integ.evaluate( integLower, integUpper );

    for ( int j = 0; j < EVTITGROOTFINDER_MAXIT; j++ ) {
        dxold = dx;
        dx = 0.5 * ( xUpper - xLower );
        rootGuess = xLower + dx;

        // If change in root is negligible, take it as solution.
        if ( fabs( xLower - rootGuess ) < precision ) {
            root = rootGuess;
            return root;
        }

        theFunc1->setCoeff( 1, 0, rootGuess );
        theFunc2->setCoeff( 1, 0, rootGuess );
        f = func1Integ.evaluate( integLower, integUpper ) -
            theFunc2->getCoeff( 1, 2 ) *
                func2Integ.evaluate( integLower, integUpper );

        if ( f < 0.0 ) {
            xLower = rootGuess;
        } else {
            xUpper = rootGuess;
        }
    }

    EvtGenReport( EVTGEN_WARNING, "EvtGen" )
        << "EvtBtoXsgammaRootFinder: Maximum number of iterations "
        << "in EvtBtoXsgammaRootFinder::foundRoot exceeded!"
        << " Returning false." << endl;
    return 0;
}
