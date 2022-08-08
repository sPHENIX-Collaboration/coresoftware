
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

#include "EvtGenBase/EvtDecayMode.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <assert.h>
#include <iostream>
using std::endl;
using std::ostream;

using std::string;
using std::vector;

EvtDecayMode::EvtDecayMode( std::string mother, vector<string> dau ) :
    _mother( mother ), _dau( dau )
{
}

EvtDecayMode::EvtDecayMode( const EvtDecayMode& other ) :
    _mother( other._mother ), _dau( other._dau )
{
}

EvtDecayMode::EvtDecayMode( const char* decay )
{
    // Parse the decay string, it should be in a standard
    // format, e.g. "B+ -> pi+ pi+ pi-" with all spaces

    string s( decay );

    // mother

    string::size_type i = s.find_first_not_of( " " );
    string::size_type j = s.find_first_of( " ", i );

    if ( i == string::npos ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "No non-space character found" << endl;
        assert( 0 );
    }

    if ( j == string::npos ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "No space before -> found" << endl;
        assert( 0 );
    }

    _mother = string( s, i, j - i );

    i = s.find_first_not_of( " ", j );
    j = s.find_first_of( "->", j );
    if ( i != j ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Multiple mothers?" << i << "," << j << endl;
        assert( 0 );
    }
    j += 2;

    while ( 1 ) {
        i = s.find_first_not_of( " ", j );
        j = s.find_first_of( " ", i );

        if ( i == string::npos )
            break;
        if ( j == string::npos ) {
            _dau.push_back( string( s, i, s.size() - i + 1 ) );
            break;
        } else {
            _dau.push_back( string( s, i, j - i ) );
        }
    }
}

const char* EvtDecayMode::mother() const
{
    return _mother.c_str();
}

int EvtDecayMode::nD() const
{
    return _dau.size();
}

const char* EvtDecayMode::dau( int i ) const
{
    assert( 0 <= i && i < (int)_dau.size() );
    return _dau[i].c_str();
}

std::string EvtDecayMode::mode() const
{
    string ret = _mother + string( " -> " );

    for ( size_t i = 0; i < _dau.size() - 1; i++ ) {
        ret += string( _dau[i] ) + string( " " );
    }
    ret += _dau[_dau.size() - 1];
    return ret;
}

ostream& EvtDecayMode::print( ostream& os ) const
{
    os << _mother.c_str() << " ->";
    for ( size_t i = 0; i < _dau.size(); i++ ) {
        os << " " << _dau[i].c_str();
    }
    return os;
}

std::string EvtDecayMode::m( EvtCyclic3::Pair i ) const
{
    string s( "m(" );
    s.append( dau( EvtCyclic3::first( i ) ) );
    s.append( "," );
    s.append( dau( EvtCyclic3::second( i ) ) );
    s.append( ")" );
    return s;
}

std::string EvtDecayMode::q( EvtCyclic3::Pair i ) const
{
    string s( "q(" );
    s.append( dau( EvtCyclic3::first( i ) ) );
    s.append( "," );
    s.append( dau( EvtCyclic3::second( i ) ) );
    s.append( ")" );
    return s;
}

std::string EvtDecayMode::dal( EvtCyclic3::Pair i, EvtCyclic3::Pair j ) const
{
    string s( q( i ) );
    s.append( ":" );
    s.append( q( j ) );
    return s;
}

ostream& operator<<( ostream& os, const EvtDecayMode& mode )
{
    mode.print( os );
    return os;
}
