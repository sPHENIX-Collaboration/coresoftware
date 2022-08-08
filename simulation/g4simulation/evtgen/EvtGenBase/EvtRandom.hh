
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

#ifndef EVTRANDOM_HH
#define EVTRANDOM_HH

class EvtRandomEngine;

class EvtRandom {
  public:
    static double Flat();
    static double Flat( double max );
    static double Flat( double min, double max );

    //generate unit Gaussian
    static double Gaussian();

    static double random();

    //This class does not take ownership of the random engine;
    //the caller needs to make sure that the engine is not
    //destroyed.
    static void setRandomEngine( EvtRandomEngine* randomEngine );

  private:
    static EvtRandomEngine* _randomEngine;
};

#endif
