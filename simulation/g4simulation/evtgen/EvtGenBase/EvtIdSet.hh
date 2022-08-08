
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

#ifndef EVTIDSET_HH
#define EVTIDSET_HH

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <string>
class EvtId;

class EvtIdSet {
  public:
    //need a default constructor

    EvtIdSet( const EvtId name1 );
    EvtIdSet( const std::string name1 );

    EvtIdSet( const EvtId name1, const EvtId name2 );

    EvtIdSet( const std::string name1, const std::string name2 );

    EvtIdSet( const EvtId name1, const EvtId name2, const EvtId name3 );

    EvtIdSet( const std::string name1, const std::string name2,
              const std::string name3 );

    EvtIdSet( const EvtId name1, const EvtId name2, const EvtId name3,
              const EvtId name4 );

    EvtIdSet( const std::string name1, const std::string name2,
              const std::string name3, const std::string name4 );

    EvtIdSet( const EvtId name1, const EvtId name2, const EvtId name3,
              const EvtId name4, const EvtId name5 );

    EvtIdSet( const std::string name1, const std::string name2,
              const std::string name3, const std::string name4,
              const std::string name5 );

    EvtIdSet( const EvtId name1, const EvtId name2, const EvtId name3,
              const EvtId name4, const EvtId name5, const EvtId name6 );

    EvtIdSet( const std::string name1, const std::string name2,
              const std::string name3, const std::string name4,
              const std::string name5, const std::string name6 );

    EvtIdSet( const EvtId name1, const EvtId name2, const EvtId name3,
              const EvtId name4, const EvtId name5, const EvtId name6,
              const EvtId name7 );

    EvtIdSet( const std::string name1, const std::string name2,
              const std::string name3, const std::string name4,
              const std::string name5, const std::string name6,
              const std::string name7 );

    EvtIdSet( const EvtId name1, const EvtId name2, const EvtId name3,
              const EvtId name4, const EvtId name5, const EvtId name6,
              const EvtId name7, const EvtId name8 );

    EvtIdSet( const std::string name1, const std::string name2,
              const std::string name3, const std::string name4,
              const std::string name5, const std::string name6,
              const std::string name7, const std::string name8 );

    EvtIdSet( const EvtId name1, const EvtId name2, const EvtId name3,
              const EvtId name4, const EvtId name5, const EvtId name6,
              const EvtId name7, const EvtId name8, const EvtId name9 );

    EvtIdSet( const std::string name1, const std::string name2,
              const std::string name3, const std::string name4,
              const std::string name5, const std::string name6,
              const std::string name7, const std::string name8,
              const std::string name9 );

    EvtIdSet( const EvtId name1, const EvtId name2, const EvtId name3,
              const EvtId name4, const EvtId name5, const EvtId name6,
              const EvtId name7, const EvtId name8, const EvtId name9,
              const EvtId name10 );

    EvtIdSet( const std::string name1, const std::string name2,
              const std::string name3, const std::string name4,
              const std::string name5, const std::string name6,
              const std::string name7, const std::string name8,
              const std::string name9, const std::string name10 );

    EvtIdSet( const EvtId name1, const EvtId name2, const EvtId name3,
              const EvtId name4, const EvtId name5, const EvtId name6,
              const EvtId name7, const EvtId name8, const EvtId name9,
              const EvtId name10, const EvtId name11 );

    EvtIdSet( const std::string name1, const std::string name2,
              const std::string name3, const std::string name4,
              const std::string name5, const std::string name6,
              const std::string name7, const std::string name8,
              const std::string name9, const std::string name10,
              const std::string name11 );

    EvtIdSet( const EvtId name1, const EvtId name2, const EvtId name3,
              const EvtId name4, const EvtId name5, const EvtId name6,
              const EvtId name7, const EvtId name8, const EvtId name9,
              const EvtId name10, const EvtId name11, const EvtId name12 );

    EvtIdSet( const std::string name1, const std::string name2,
              const std::string name3, const std::string name4,
              const std::string name5, const std::string name6,
              const std::string name7, const std::string name8,
              const std::string name9, const std::string name10,
              const std::string name11, const std::string name12 );

    ~EvtIdSet() { delete[] _list; }

    EvtIdSet( const EvtIdSet& set1 );
    EvtIdSet( const EvtIdSet& set1, const EvtIdSet& set2 );

    int contains( const EvtId id );
    int contains( const std::string id );

    void append( const EvtIdSet set1 );
    int sizeOfSet() const;
    EvtId getElem( const int i ) const;

  private:
    int _numInList;
    EvtId* _list;
};

#endif
