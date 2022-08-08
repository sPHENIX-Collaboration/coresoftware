
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

#ifndef EVTID_HH
#define EVTID_HH

#include <iostream>
#include <string>

class EvtId {
  public:
    //need a default constructor
    EvtId() : _id( -1 ), _alias( -1 ) {}

    EvtId( int id, int alias ) : _id( id ), _alias( alias ) {}

    friend std::ostream& operator<<( std::ostream& s, const EvtId& v );

    int operator==( const EvtId& id ) const { return _id == id._id; }
    int operator!=( const EvtId& id ) const { return _id != id._id; }
    int operator<( const EvtId& id ) const { return _id < id._id; }

    int isConjugate( const EvtId& id ) const;

    int getId() const { return _id; }

    int getAlias() const { return _alias; }

    int isAlias() const { return _id != _alias; }

    std::string getName() const;

  private:
    //particle number 0..n. The order of particles are determined
    //by the order in pdt.table
    int _id;
    //if the particle is an alias to another particle alias!=id
    //The only place where the alias should be used is for looking
    //up decays in the decay table.
    int _alias;
};

#endif
