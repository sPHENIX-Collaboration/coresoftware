
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

#ifndef EVTMODELALIAS_HH
#define EVTMODELALIAS_HH

#include <string>
#include <vector>

class EvtModelAlias {
  public:
    EvtModelAlias(){};
    EvtModelAlias( std::string alias, std::string model,
                   std::vector<std::string> args );
    ~EvtModelAlias(){};
    EvtModelAlias( const EvtModelAlias& copyMe );
    EvtModelAlias operator=( const EvtModelAlias& copyMe );
    bool matchAlias( const std::string& cand )
    {
        if ( cand == _aliasName )
            return true;
        return false;
    }
    std::string getName() { return _model; }
    std::vector<std::string> getArgList();

  private:
    std::string _aliasName;
    std::string _model;
    std::vector<std::string> _modelArgs;
};
#endif
