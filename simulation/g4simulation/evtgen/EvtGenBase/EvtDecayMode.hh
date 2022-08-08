
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

#ifndef EVT_DECAY_MODE_HH
#define EVT_DECAY_MODE_HH

#include "EvtGenBase/EvtCyclic3.hh"

#include <iosfwd>
#include <string>
#include <vector>

class EvtDecayMode final {
  public:
    EvtDecayMode( const char* decay );
    EvtDecayMode( const EvtDecayMode& other );
    EvtDecayMode( std::string mother, std::vector<std::string> dau );

    const char* mother() const;
    int nD() const;
    const char* dau( int i ) const;

    std::ostream& print( std::ostream& ) const;

    // Frequent name combinations

    std::string m( EvtCyclic3::Pair i ) const;
    std::string q( EvtCyclic3::Pair i ) const;
    std::string dal( EvtCyclic3::Pair i, EvtCyclic3::Pair j ) const;
    std::string mode() const;

  private:
    std::string _mother;
    std::vector<std::string> _dau;
};

std::ostream& operator<<( std::ostream&, const EvtDecayMode& );

#endif
