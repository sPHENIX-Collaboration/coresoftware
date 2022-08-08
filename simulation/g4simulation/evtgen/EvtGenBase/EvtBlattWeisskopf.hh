
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

#ifndef EVT_BLATT_WEISSKOPF_HH
#define EVT_BLATT_WEISSKOPF_HH

// Blatt-Weisskopf penetration form factor for a resonance R->AB.
// Taken from CLEO preprint 00-23 (hep-ex/0011065)

class EvtBlattWeisskopf {
  public:
    EvtBlattWeisskopf( int LL, double R, double p0 );
    EvtBlattWeisskopf( const EvtBlattWeisskopf& );

    double operator()( double p ) const;

  private:
    int _LL;           // angular momentum of daughters
    double _radial;    // resonance radial parameter
    double _p0;

    double _F0;    // formula evaluated at _p0
    double compute( double p ) const;
};

#endif
