
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

#ifndef EVT_PRED_GEN_HH
#define EVT_PRED_GEN_HH

#include <stdio.h>

// A predicate is applied to a generator to get another generator.
// Accept-reject can be implemented in this way.
//
//           Predicate
// Generator    ->     Generator

template <class Generator, class Predicate>
class EvtPredGen {
  public:
    typedef typename Generator::result_type result_type;

    EvtPredGen() : itsTried( 0 ), itsPassed( 0 ) {}

    EvtPredGen( Generator gen, Predicate pred ) :
        itsGen( gen ), itsPred( pred ), itsTried( 0 ), itsPassed( 0 )
    {
    }

    EvtPredGen( const EvtPredGen& other ) :
        itsGen( other.itsGen ),
        itsPred( other.itsPred ),
        itsTried( other.itsTried ),
        itsPassed( other.itsPassed )
    {
    }

    ~EvtPredGen() {}

    result_type operator()()
    {
        int i = 0;
        int MAX = 10000;
        while ( i++ < MAX ) {
            itsTried++;
            result_type point = itsGen();
            if ( itsPred( point ) ) {
                itsPassed++;
                return point;
            }
        }

        printf( "No random point generated after %d attempts\n", MAX );
        printf( "Sharp peak? Consider using pole compensation.\n" );
        printf( "I will now pick a point at random to return.\n" );
        return itsGen();
    }

    inline int getTried() const { return itsTried; }
    inline int getPassed() const { return itsPassed; }

  protected:
    Generator itsGen;
    Predicate itsPred;
    int itsTried;
    int itsPassed;
};

#endif
