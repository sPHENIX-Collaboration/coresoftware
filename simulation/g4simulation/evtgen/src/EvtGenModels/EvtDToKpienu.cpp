
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

#include "EvtGenModels/EvtDToKpienu.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

std::string EvtDToKpienu::getName()
{
    return "DToKpienu";
}

EvtDecayBase* EvtDToKpienu::clone()
{
    return new EvtDToKpienu;
}

void EvtDToKpienu::init()
{
    checkNArg( 0 );
    checkNDaug( 4 );
    checkSpinParent( EvtSpinType::SCALAR );
    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "EvtDToKpienu ==> Initialization !" << std::endl;
    nAmps = 2;

    rS = -11.57;    // S-wave
    rS1 = 0.08;
    a_delta = 1.94;
    b_delta = -0.81;
    m0_1430_S = 1.425;
    width0_1430_S = 0.270;
    type[0] = 0;

    mV = 1.81;
    mA = 2.61;
    V_0 = 1.411;
    A1_0 = 1;
    A2_0 = 0.788;

    m0 = 0.8946;    // P-wave K*
    width0 = 0.04642;
    rBW = 3.07;
    rho = 1;
    phi = 0;
    type[1] = 1;

    m0_1410 = 1.414;    // P-wave K*(1410)
    width0_1410 = 0.232;
    rho_1410 = 0.1;
    phi_1410 = 0.;
    type[2] = 2;

    TV_0 = 1;    // D-wave K*2(1430)
    T1_0 = 1;
    T2_0 = 1;
    m0_1430 = 1.4324;
    width0_1430 = 0.109;
    rho_1430 = 15;
    phi_1430 = 0;
    type[3] = 3;

    mD = 1.86962;
    mPi = 0.13957;
    mK = 0.49368;
    Pi = atan2( 0.0, -1.0 );
    root2 = sqrt( 2. );
    root2d3 = sqrt( 2. / 3 );
    root1d2 = sqrt( 0.5 );
    root3d2 = sqrt( 1.5 );
}

void EvtDToKpienu::initProbMax()
{
    /*
   * This piece of code could in principle be used to calculate maximum
   * probablity on fly. But as it uses high number of points and model
   * deals with single final state, we keep hardcoded number for now rather
   * than adapting code to work here.

     double maxprob = 0.0;
     for(int ir=0;ir<=60000000;ir++){
        p->initializePhaseSpace(getNDaug(),getDaugs());
        EvtVector4R _K  = p->getDaug(0)->getP4();
        EvtVector4R _pi = p->getDaug(1)->getP4();
        EvtVector4R _e  = p->getDaug(2)->getP4();
        EvtVector4R _nu = p->getDaug(3)->getP4();
        int pid = EvtPDL::getStdHep(p->getDaug(0)->getId());
        int charm;
        if(pid == -321) charm = 1;
        else charm = -1;
        double m2, q2, cosV, cosL, chi;
        KinVGen(_K, _pi, _e, _nu, charm, m2, q2, cosV, cosL, chi);
        double _prob = calPDF(m2, q2, cosV, cosL, chi);
        if(_prob>maxprob) {
            maxprob=_prob;
            EvtGenReport(EVTGEN_INFO,"EvtGen") << "Max PDF = " << ir << " charm= " << charm << " prob= "
     << _prob << std::endl;
        }
     }
     EvtGenReport(EVTGEN_INFO,"EvtGen") << "Max!!!!!!!!!!! " << maxprob<< std::endl;
  */
    setProbMax( 22750.0 );
}

void EvtDToKpienu::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    EvtVector4R K = p->getDaug( 0 )->getP4();
    EvtVector4R pi = p->getDaug( 1 )->getP4();
    EvtVector4R e = p->getDaug( 2 )->getP4();
    EvtVector4R nu = p->getDaug( 3 )->getP4();

    int pid = EvtPDL::getStdHep( p->getDaug( 0 )->getId() );
    int charm;
    if ( pid == -321 ) {
        charm = 1;
    } else {
        charm = -1;
    }
    double m2, q2, cosV, cosL, chi;
    KinVGen( K, pi, e, nu, charm, m2, q2, cosV, cosL, chi );
    double prob = calPDF( m2, q2, cosV, cosL, chi );
    setProb( prob );
    return;
}

void EvtDToKpienu::KinVGen( const EvtVector4R& vp4_K, const EvtVector4R& vp4_Pi,
                            const EvtVector4R& vp4_Lep, const EvtVector4R& vp4_Nu,
                            const int charm, double& m2, double& q2,
                            double& cosV, double& cosL, double& chi ) const
{
    EvtVector4R vp4_KPi = vp4_K + vp4_Pi;
    EvtVector4R vp4_W = vp4_Lep + vp4_Nu;

    m2 = vp4_KPi.mass2();
    q2 = vp4_W.mass2();

    EvtVector4R boost;
    boost.set( vp4_KPi.get( 0 ), -vp4_KPi.get( 1 ), -vp4_KPi.get( 2 ),
               -vp4_KPi.get( 3 ) );
    EvtVector4R vp4_Kp = boostTo( vp4_K, boost );
    cosV = vp4_Kp.dot( vp4_KPi ) / ( vp4_Kp.d3mag() * vp4_KPi.d3mag() );

    boost.set( vp4_W.get( 0 ), -vp4_W.get( 1 ), -vp4_W.get( 2 ), -vp4_W.get( 3 ) );
    EvtVector4R vp4_Lepp = boostTo( vp4_Lep, boost );
    cosL = vp4_Lepp.dot( vp4_W ) / ( vp4_Lepp.d3mag() * vp4_W.d3mag() );

    EvtVector4R V = vp4_KPi / vp4_KPi.d3mag();
    EvtVector4R C = vp4_Kp.cross( V );
    C /= C.d3mag();
    EvtVector4R D = vp4_Lepp.cross( V );
    D /= D.d3mag();
    double sinx = C.cross( V ).dot( D );
    double cosx = C.dot( D );
    chi = sinx > 0 ? acos( cosx ) : -acos( cosx );
    if ( charm == -1 )
        chi = -chi;
}

double EvtDToKpienu::calPDF( const double m2, const double q2, const double cosV,
                             const double cosL, const double chi ) const
{
    double m = sqrt( m2 );
    double q = sqrt( q2 );

    // begin to calculate form factor
    EvtComplex F10( 0.0, 0.0 );
    EvtComplex F11( 0.0, 0.0 );
    EvtComplex F21( 0.0, 0.0 );
    EvtComplex F31( 0.0, 0.0 );
    EvtComplex F12( 0.0, 0.0 );
    EvtComplex F22( 0.0, 0.0 );
    EvtComplex F32( 0.0, 0.0 );

    EvtComplex f10( 0.0, 0.0 );
    EvtComplex f11( 0.0, 0.0 );
    EvtComplex f21( 0.0, 0.0 );
    EvtComplex f31( 0.0, 0.0 );
    EvtComplex f12( 0.0, 0.0 );
    EvtComplex f22( 0.0, 0.0 );
    EvtComplex f32( 0.0, 0.0 );
    EvtComplex coef( 0.0, 0.0 );
    double amplitude_temp, delta_temp;

    for ( int index = 0; index < nAmps; index++ ) {
        switch ( type[index] ) {
            case 0:    // calculate form factor of S wave
            {
                NRS( m, q, rS, rS1, a_delta, b_delta, mA, m0_1430_S,
                     width0_1430_S, amplitude_temp, delta_temp, f10 );
                F10 = F10 + f10;
                break;
            }
            case 1:    // calculate form factor of P wave (K*)
            {
                ResonanceP( m, q, mV, mA, V_0, A1_0, A2_0, m0, width0, rBW,
                            amplitude_temp, delta_temp, f11, f21, f31 );
                coef = getCoef( rho, phi );
                F11 = F11 + coef * f11;
                F21 = F21 + coef * f21;
                F31 = F31 + coef * f31;
                break;
            }
            case 2:    // calculate form factor of P wave (K*(1410))
            {
                ResonanceP( m, q, mV, mA, V_0, A1_0, A2_0, m0_1410, width0_1410,
                            rBW, amplitude_temp, delta_temp, f11, f21, f31 );
                coef = getCoef( rho_1410, phi_1410 );
                F11 = F11 + coef * f11;
                F21 = F21 + coef * f21;
                F31 = F31 + coef * f31;
                break;
            }
            case 3:    // calculate form factor of D wave
            {
                ResonanceD( m, q, mV, mA, TV_0, T1_0, T2_0, m0_1430, width0_1430,
                            rBW, amplitude_temp, delta_temp, f12, f22, f32 );
                coef = getCoef( rho_1430, phi_1430 );
                F12 = F12 + coef * f12;
                F22 = F22 + coef * f22;
                F32 = F32 + coef * f32;
                break;
            }
            default: {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "No such form factor type!!!" << std::endl;
                break;
            }
        }
    }

    // begin to calculate pdf value
    double I, I1, I2, I3, I4, I5, I6, I7, I8, I9;

    double cosV2 = cosV * cosV;
    double sinV = sqrt( 1.0 - cosV2 );
    double sinV2 = sinV * sinV;

    EvtComplex F1 = F10 + F11 * cosV + F12 * ( 1.5 * cosV2 - 0.5 );
    EvtComplex F2 = F21 * root1d2 + F22 * cosV * root3d2;
    EvtComplex F3 = F31 * root1d2 + F32 * cosV * root3d2;

    I1 = 0.25 * ( abs2( F1 ) + 1.5 * sinV2 * ( abs2( F2 ) + abs2( F3 ) ) );
    I2 = -0.25 * ( abs2( F1 ) - 0.5 * sinV2 * ( abs2( F2 ) + abs2( F3 ) ) );
    I3 = -0.25 * ( abs2( F2 ) - abs2( F3 ) ) * sinV2;
    I4 = real( conj( F1 ) * F2 ) * sinV * 0.5;
    I5 = real( conj( F1 ) * F3 ) * sinV;
    I6 = real( conj( F2 ) * F3 ) * sinV2;
    I7 = imag( conj( F2 ) * F1 ) * sinV;
    I8 = imag( conj( F3 ) * F1 ) * sinV * 0.5;
    I9 = imag( conj( F3 ) * F2 ) * sinV2 * ( -0.5 );

    double coschi = cos( chi );
    double sinchi = sin( chi );
    double sin2chi = 2.0 * sinchi * coschi;
    double cos2chi = 1.0 - 2.0 * sinchi * sinchi;

    double sinL = sqrt( 1. - cosL * cosL );
    double sinL2 = sinL * sinL;
    double sin2L = 2.0 * sinL * cosL;
    double cos2L = 1.0 - 2.0 * sinL2;

    I = I1 + I2 * cos2L + I3 * sinL2 * cos2chi + I4 * sin2L * coschi +
        I5 * sinL * coschi + I6 * cosL + I7 * sinL * sinchi +
        I8 * sin2L * sinchi + I9 * sinL2 * sin2chi;
    return I;
}

void EvtDToKpienu::ResonanceP( const double m, const double q, const double mV,
                               const double mA, const double V_0,
                               const double A1_0, const double A2_0,
                               const double m0, const double width0,
                               const double rBW, double& amplitude,
                               double& delta, EvtComplex& F11, EvtComplex& F21,
                               EvtComplex& F31 ) const
{
    double pKPi = getPStar( mD, m, q );
    double mD2 = mD * mD;
    double m2 = m * m;
    double m02 = m0 * m0;
    double q2 = q * q;
    double mV2 = mV * mV;
    double mA2 = mA * mA;
    double summDm = mD + m;
    double V = V_0 / ( 1.0 - q2 / ( mV2 ) );
    double A1 = A1_0 / ( 1.0 - q2 / ( mA2 ) );
    double A2 = A2_0 / ( 1.0 - q2 / ( mA2 ) );
    double A = summDm * A1;
    double B = 2.0 * mD * pKPi / summDm * V;

    // construct the helicity form factor
    double H0 = 0.5 / ( m * q ) *
                ( ( mD2 - m2 - q2 ) * summDm * A1 -
                  4.0 * ( mD2 * pKPi * pKPi ) / summDm * A2 );
    double Hp = A - B;
    double Hm = A + B;

    // calculate alpha
    double B_Kstar = 2. / 3.;    // B_Kstar = Br(Kstar(892)->k pi)
    double pStar0 = getPStar( m0, mPi, mK );
    double alpha = sqrt( 3. * Pi * B_Kstar / ( pStar0 * width0 ) );

    // construct amplitudes of (non)resonance
    double F = getF1( m, m0, mPi, mK, rBW );
    double width = getWidth1( m, m0, mPi, mK, width0, rBW );

    EvtComplex C( m0 * width0 * F, 0.0 );
    double AA = m02 - m2;
    double BB = -m0 * width;
    EvtComplex amp = C / EvtComplex( AA, BB );
    amplitude = abs( amp );
    delta = atan2( imag( amp ), real( amp ) );

    double alpham2 = alpha * 2.0;
    F11 = amp * alpham2 * q * H0 * root2;
    F21 = amp * alpham2 * q * ( Hp + Hm );
    F31 = amp * alpham2 * q * ( Hp - Hm );
}

void EvtDToKpienu::NRS( const double m, const double q, const double rS,
                        const double rS1, const double a_delta,
                        const double b_delta, const double mA, const double m0,
                        const double width0, double& amplitude, double& delta,
                        EvtComplex& F10 ) const
{
    static const double tmp = ( mK + mPi ) * ( mK + mPi );

    double m2 = m * m;
    double q2 = q * q;
    double mA2 = mA * mA;
    double pKPi = getPStar( mD, m, q );
    double m_K0_1430 = m0;
    double width_K0_1430 = width0;
    double m2_K0_1430 = m_K0_1430 * m_K0_1430;
    double width = getWidth0( m, m_K0_1430, mPi, mK, width_K0_1430 );

    // calculate modul of the amplitude
    double x, Pm;
    if ( m < m_K0_1430 ) {
        x = sqrt( m2 / tmp - 1.0 );
        Pm = 1.0 + rS1 * x;
    } else {
        x = sqrt( m2_K0_1430 / tmp - 1.0 );
        Pm = 1.0 + rS1 * x;
        Pm *= m_K0_1430 * width_K0_1430 /
              sqrt( ( m2_K0_1430 - m2 ) * ( m2_K0_1430 - m2 ) +
                    m2_K0_1430 * width * width );
    }

    // calculate phase of the amplitude
    double pStar = getPStar( m, mPi, mK );
    double delta_bg = atan( 2. * a_delta * pStar /
                            ( 2. + a_delta * b_delta * pStar * pStar ) );
    delta_bg = ( delta_bg > 0 ) ? delta_bg : ( delta_bg + Pi );

    double delta_K0_1430 = atan( m_K0_1430 * width / ( m2_K0_1430 - m2 ) );
    delta_K0_1430 = ( delta_K0_1430 > 0 ) ? delta_K0_1430
                                          : ( delta_K0_1430 + Pi );
    delta = delta_bg + delta_K0_1430;

    EvtComplex ci( cos( delta ), sin( delta ) );
    EvtComplex amp = ci * rS * Pm;
    amplitude = rS * Pm;

    F10 = amp * pKPi * mD / ( 1. - q2 / mA2 );
}

void EvtDToKpienu::ResonanceD( const double m, const double q, const double mV,
                               const double mA, const double TV_0,
                               const double T1_0, const double T2_0,
                               const double m0, const double width0,
                               const double rBW, double& amplitude,
                               double& delta, EvtComplex& F12, EvtComplex& F22,
                               EvtComplex& F32 ) const
{
    double pKPi = getPStar( mD, m, q );
    double mD2 = mD * mD;
    double m2 = m * m;
    double m02 = m0 * m0;
    double q2 = q * q;
    double mV2 = mV * mV;
    double mA2 = mA * mA;
    double summDm = mD + m;
    double TV = TV_0 / ( 1.0 - q2 / ( mV2 ) );
    double T1 = T1_0 / ( 1.0 - q2 / ( mA2 ) );
    double T2 = T2_0 / ( 1.0 - q2 / ( mA2 ) );

    // construct amplitudes of (non)resonance
    double F = getF2( m, m0, mPi, mK, rBW );
    double width = getWidth2( m, m0, mPi, mK, width0, rBW );
    EvtComplex C( m0 * width0 * F, 0.0 );
    double AA = m02 - m2;
    double BB = -m0 * width;

    EvtComplex amp = C / EvtComplex( AA, BB );
    amplitude = abs( amp );
    delta = atan2( imag( amp ), real( amp ) );

    F12 = amp * mD * pKPi / 3. *
          ( ( mD2 - m2 - q2 ) * summDm * T1 - mD2 * pKPi * pKPi / summDm * T2 );
    F22 = amp * root2d3 * mD * m * q * pKPi * summDm * T1;
    F32 = amp * root2d3 * 2. * mD2 * m * q * pKPi * pKPi / summDm * TV;
}

double EvtDToKpienu::getPStar( const double m, const double m1,
                               const double m2 ) const
{
    double s = m * m;
    double s1 = m1 * m1;
    double s2 = m2 * m2;
    double x = s + s1 - s2;
    double t = 0.25 * x * x / s - s1;
    double p;
    if ( t > 0.0 ) {
        p = sqrt( t );
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " Hello, pstar is less than 0.0" << std::endl;
        p = 0.04;
    }
    return p;
}

double EvtDToKpienu::getF1( const double m, const double m0, const double m_c1,
                            const double m_c2, const double rBW ) const
{
    double pStar = getPStar( m, m_c1, m_c2 );
    double pStar0 = getPStar( m0, m_c1, m_c2 );
    double rBW2 = rBW * rBW;
    double pStar2 = pStar * pStar;
    double pStar02 = pStar0 * pStar0;
    double B = 1. / sqrt( 1. + rBW2 * pStar2 );
    double B0 = 1. / sqrt( 1. + rBW2 * pStar02 );
    double F = pStar / pStar0 * B / B0;
    return F;
}

double EvtDToKpienu::getF2( const double m, const double m0, const double m_c1,
                            const double m_c2, const double rBW ) const
{
    double pStar = getPStar( m, m_c1, m_c2 );
    double pStar0 = getPStar( m0, m_c1, m_c2 );
    double rBW2 = rBW * rBW;
    double pStar2 = pStar * pStar;
    double pStar02 = pStar0 * pStar0;
    double B = 1. / sqrt( ( rBW2 * pStar2 - 3. ) * ( rBW2 * pStar2 - 3. ) +
                          9. * rBW2 * pStar2 );
    double B0 = 1. / sqrt( ( rBW2 * pStar02 - 3. ) * ( rBW2 * pStar02 - 3. ) +
                           9. * rBW2 * pStar02 );
    double F = pStar2 / pStar02 * B / B0;
    return F;
}

double EvtDToKpienu::getWidth0( const double m, const double m0,
                                const double m_c1, const double m_c2,
                                const double width0 ) const
{
    double pStar = getPStar( m, m_c1, m_c2 );
    double pStar0 = getPStar( m0, m_c1, m_c2 );
    double width = width0 * pStar / pStar0 * m0 / m;
    return width;
}

double EvtDToKpienu::getWidth1( const double m, const double m0,
                                const double m_c1, const double m_c2,
                                const double width0, const double rBW ) const
{
    double pStar = getPStar( m, m_c1, m_c2 );
    double pStar0 = getPStar( m0, m_c1, m_c2 );
    double F = getF1( m, m0, m_c1, m_c2, rBW );
    double width = width0 * pStar / pStar0 * m0 / m * F * F;
    return width;
}

double EvtDToKpienu::getWidth2( const double m, const double m0,
                                const double m_c1, const double m_c2,
                                const double width0, const double rBW ) const
{
    double pStar = getPStar( m, m_c1, m_c2 );
    double pStar0 = getPStar( m0, m_c1, m_c2 );
    double F = getF2( m, m0, m_c1, m_c2, rBW );
    double width = width0 * pStar / pStar0 * m0 / m * F * F;
    return width;
}

EvtComplex EvtDToKpienu::getCoef( const double rho, const double phi ) const
{
    EvtComplex coef( rho * cos( phi ), rho * sin( phi ) );
    return coef;
}
