
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

#ifndef EVTGENBASE_EVTDILOG_HH
#define EVTGENBASE_EVTDILOG_HH 1

// Include files

/** @namespace EvtDiLog EvtDiLog.hh EvtGenBase/EvtDiLog.hh
 *  Dilogarithm function (replaces CERNLIB DDILOG)
 *
 *  @author Patrick Robbe
 *  @date   2007-01-23
 */
namespace EvtDiLog {
    double DiLog( double x );

    // constants for computation
    static const double Z1 = 1.;
    static const double HF = Z1 / 2.;
    static const double PI = 3.14159265358979324;
    static const double PI3 = PI * PI / 3.;
    static const double PI6 = PI * PI / 6.;
    static const double PI12 = PI * PI / 12.;
    static const double C[20] = {
        0.42996693560813697,  0.40975987533077105,  -0.01858843665014592,
        0.00145751084062268,  -0.00014304184442340, 0.00001588415541880,
        -0.00000190784959387, 0.00000024195180854,  -0.00000003193341274,
        0.00000000434545063,  -0.00000000060578480, 0.00000000008612098,
        -0.00000000001244332, 0.00000000000182256,  -0.00000000000027007,
        0.00000000000004042,  -0.00000000000000610, 0.00000000000000093,
        -0.00000000000000014, 0.00000000000000002};
}    // namespace EvtDiLog
#endif    // EVTGENBASE_EVTDILOG_HH
