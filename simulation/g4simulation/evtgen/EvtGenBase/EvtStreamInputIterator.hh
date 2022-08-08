
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

#ifndef EVT_STREAM_INPUT_ITERATOR_HH
#define EVT_STREAM_INPUT_ITERATOR_HH

#include "EvtGenBase/EvtStreamAdapter.hh"

#include <cstddef>
#include <iterator>

using std::input_iterator_tag;

// Adapters are used to convert various types of input streams
// into an iteratable interface.

template <class Point>
class EvtStreamInputIterator {
  public:
    typedef input_iterator_tag iterator_category;
    typedef Point value_type;
    typedef ptrdiff_t difference_type;
    typedef const Point* pointer;
    typedef const Point& reference;

    EvtStreamInputIterator() : _counter( 0 ) {}

    EvtStreamInputIterator( const EvtStreamInputIterator& other ) :
        _counter( other._counter ? other._counter->clone() : 0 ),
        _currentValue( other._currentValue )
    {
    }

    EvtStreamInputIterator( EvtStreamAdapter<Point>& counter ) :
        _counter( counter.clone() )
    {
        _currentValue = _counter->currentValue();
    }

    ~EvtStreamInputIterator()
    {
        if ( _counter )
            delete _counter;
    }

    reference operator*() const { return _currentValue; }

    EvtStreamInputIterator& operator++()
    {
        _read();
        return *this;
    }

    EvtStreamInputIterator operator++( int )
    {
        EvtStreamInputIterator tmp = *this;
        _read();
        return tmp;
    }

    bool operator==( const EvtStreamInputIterator& other ) const
    {
        // Equality is only defined for two past the end iterators
        return ( pastEnd() && other.pastEnd() );
    }

  protected:
    EvtStreamAdapter<Point>* _counter;
    value_type _currentValue;

    bool pastEnd() const
    {
        bool ret = true;
        if ( _counter )
            ret = _counter->pastEnd();
        return ret;
    }

    // Advances the iterator

    void _read()
    {
        _counter->advance();
        _currentValue = _counter->currentValue();
    }
};

// For adaptable generators these shorthand functions can be used
// to construct iterators.

template <class Generator>
EvtStreamInputIterator<typename Generator::result_type> iter( Generator gen,
                                                              int N = 0 )
{
    typedef typename Generator::result_type Point;
    EvtGenStreamAdapter<Point, Generator> counter( gen, N );
    return EvtStreamInputIterator<Point>( counter );
}

#endif
