
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

#ifndef EVTMODEL_HH
#define EVTMODEL_HH

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtStringHash.hh"

#include <map>
//#include <fstream.h>

//Class to read in and handle the decays available
//to EvtGen for each particle, and the model to be
//used for each one.

class EvtModel {
  public:
    static EvtModel& instance();

    void registerModel( EvtDecayBase* prototype );

    int isModel( std::string name );

    EvtDecayBase* getFcn( std::string model_name );

    int isCommand( std::string cmd );
    void storeCommand( std::string cmd, std::string cnfgstr );

  private:
    EvtModel();

    static EvtModel* _instance;

    std::map<std::string, EvtDecayBase*> _modelNameHash;
    std::map<std::string, EvtDecayBase*> _commandNameHash;
};

inline EvtModel& EvtModel::instance()
{
    if ( _instance == 0 )
        _instance = new EvtModel;
    return *_instance;
}

#endif
