
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

#ifndef EVT_ABS_BINNING_HH
#define EVT_ABS_BINNING_HH
#define BIN_OUTSIDE -1

#include <stdio.h>

/*
 * Data point to bin value mapping
 */

template <class T>
class EvtAbsBinning {
  public:
    EvtAbsBinning() {}
    EvtAbsBinning( const EvtAbsBinning<T>& other ) {}
    virtual ~EvtAbsBinning() {}

    virtual EvtAbsBinning<T>* clone() const = 0;
    virtual int getBin( const T& point ) const = 0;
    virtual T getBinPoint( int bin ) const = 0;
    virtual double size( int bin ) const = 0;

    virtual int nTypes() const = 0;

    virtual char* typeLabel( int i ) const
    {
        char* a = new char[128];
        sprintf( a, "%d", i );
        return a;
    }
};

#endif
