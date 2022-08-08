
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

#ifndef EVTBTOSLLWILCNLO_HH
#define EVTBTOSLLWILCNLO_HH
class EvtComplex;

// Description: The calculation of the Wilson coefficients for
//              b -> (d,s) ell+ ell-  transitions in the SM at NLO
//              according to the paper:
//                    A.J.Buras, M.Munz, Phys.Rev.D52, p.189 (1995).

class EvtbTosllWilsCoeffNLO {
  public:
    double As( double mu, int Nf, int ias );
    double Li2( double w );

    double C1( double mu, double Mw, int Nf, int ias );
    double C2( double mu, double Mw, int Nf, int ias );
    double C3( double mu, double Mw, int Nf, int ias );
    double C4( double mu, double Mw, int Nf, int ias );
    double C5( double mu, double Mw, int Nf, int ias );
    double C6( double mu, double Mw, int Nf, int ias );

    double A( double z );
    double B( double z );
    double C_Bur( double z );
    double D_Bur( double z );
    double E( double z );
    double F_Bur( double z );
    double Y( double z );
    double Z( double z );

    double C7gamma( double mu, double Mw, double mt, int Nf, int ias );

    double Pe( double eta );
    double P0ndr( double asW, double eta );

    double C9v( double mu, double Mw, double mt, int Nf, int ias );
    double C10a( double mt, double Mw );

    double Reh( double mu, double mQ, double q2 );
    double Imh( double mQ, double q2 );

    double ReResonant( double q2, double GV, double GllV, double MV );
    double ImResonant( double q2, double GV, double GllV, double MV );

    double ReHtot( int qflavour, int res_swch, int ias, int Nf, double mu,
                   double mQ, double q2, double ml, double Mw );
    double ImHtot( int qflavour, int res_swch, int ias, int Nf, double mu,
                   double mQ, double q2, double ml, double Mw );

    double omega( double q2, double m2 );

    double ReC9eff( int decay_id, int res_swch, int ias, int Nf, double q2,
                    double m2, double md, double mc, double mu, double mt,
                    double Mw, double ml, double Relambda_qu,
                    double Imlambda_qu );
    double ImC9eff( int decay_id, int res_swch, int ias, int Nf, double q2,
                    double m2, double md, double mc, double mu, double Mw,
                    double ml, double Relambda_qu, double Imlambda_qu );

    EvtComplex GetC9Eff( int decay_id, int res_swch, int ias, int Nf, double q2,
                         double m2, double md, double mc, double mu, double mt,
                         double Mw, double ml, double Relambda_qu,
                         double Imlambda_qu );
    EvtComplex GetC10Eff( double mt, double Mw );
    EvtComplex GetC7Eff( double mu, double Mw, double mt, int Nf, int ias );
};

#endif
