
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

#ifndef EVTDECAYPARM_HH
#define EVTDECAYPARM_HH

#include <string>

// Description: Class to keep the arguments and daughters of a decay

class EvtParticle;

typedef void ( *fcnPtr )( EvtParticle*, int, int*, double* );

class EvtDecayParm {
  public:
    EvtDecayParm();
    ~EvtDecayParm();

    void init( fcnPtr pfcn, int ndaug, int* daugs, int narg, double* args,
               std::string name );

    int getNDaug() { return itsndaug; }
    int getNArg() { return itsnarg; }
    int* getDaugs() { return itsdaugs; }
    double* getArgs() { return itsargs; }
    fcnPtr getfcnPtr() { return itsfcn; }
    std::string getModelName() { return modelname; }

  private:
    fcnPtr itsfcn;
    int itsndaug;
    int* itsdaugs;
    int itsnarg;
    double* itsargs;
    std::string modelname;
};

#endif
