/***********************************************************************
 * Software License Agreement (BSD License)
 *
 * Copyright 2019-2020  Dmitry Arkhipkim (arkhipkin@bnl.gov). All rights reserved.
 *   All rights reserved.
 *
 * THE BSD LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/

/** \mainpage kdfinder C++ API documentation
  *  kdfinder is a C++ header-only track finder
  *  mostly optimized for quick track pattern recognition
  *
  *  Author: Dmitry Arkhipkin <arkhipkin@bnl.gov> 2018-2019
  *
  *  kdfinder does not require compiling or installing,
  *  just include kdfinder.hpp and nanoflann.hpp in your code.
  *
  *  See:
  *   - <a href="https://drupal.star.bnl.gov/STAR/comp/db/kdfinder/" >Drupal documentation / TBD </a>
  */

#ifndef  KDFINDER_HPP_
#define  KDFINDER_HPP_

#include "nanoflann.hpp"

#include <algorithm>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <thread>

namespace kdfinder
{

template <typename T>
T rnd_gauss( T mean, T stddev )
{
	T u, v, r, gauss;
	do {
		u = 2 * ((double) rand() / (RAND_MAX)) - 1;
		v = 2 * ((double) rand() / (RAND_MAX)) - 1;
		r = u * u + v * v;
	} while ( r > 1 || r == 0 );
	gauss = u * std::sqrt( -2 * std::log( r ) / r );
	return ( mean + gauss * stddev );
}

template<class T>
class TVector {

	public:
		TVector() : mX1(0), mX2(0), mX3(0) {}
		TVector( T x, T y, T z ) : mX1(x), mX2(y), mX3(z) {}
		TVector( const TVector<T>& v ) : mX1( v.x() ), mX2( v.y() ), mX3( v.z() ) {}
		TVector( const T* a ) : mX1( a[0] ), mX2( a[1] ), mX3( a[2] ) {}
		~TVector() {}

    void setX( T x ) { mX1 = x; }
    void setY( T y ) { mX2 = y; }
    void setZ( T z ) { mX3 = z; }

    void set( T X, T Y, T Z ) { mX1 = X; mX2 = Y; mX3 = Z; }

    void setPhi( T Angle ) {
		double  r = magnitude();
		double th = theta();
		mX1 = r * std::sin( th ) * std::cos( Angle );
		mX2 = r * std::sin( th ) * std::sin( Angle );
	}

    void setTheta( T Angle ) {
		double r  = magnitude();
		double ph = phi();
		mX1 = r*sin(Angle)*cos(ph);
		mX2 = r*sin(Angle)*sin(ph);
		mX3 = r*cos(Angle);
	}

    void setMag( T Mag ) {
		setMagnitude( Mag );
 	}

    void setMagnitude( T r ) {
		double th = theta();
		double ph = phi();
		mX1 = r*sin(th)*cos(ph);
		mX2 = r*sin(th)*sin(ph);
		mX3 = r*cos(th);
	}
    
    const T& x() const { return mX1; }
    const T& y() const { return mX2; }
    const T& z() const { return mX3; }

    const T* xyz() const { return &mX1; }
          T* xyz() { return &mX1; }

	T abs( const TVector<T>& v) { return v.mag(); }

    T theta() const {
		return acos( cosTheta() );
	}

    T cosTheta() const {
		return mX3 / ( mag() + 1e-20 );
	}

    T phi() const {
		return std::atan2( mX2, mX1 );
	}

    T perp() const {
		return std::sqrt( mX1 * mX1 + mX2 * mX2 );
	}

    T perp2() const {
		return mX1 * mX1 + mX2 * mX2;
	}

    T magnitude() const {
		return mag();
	}

    T mag() const {
		return std::sqrt( mX1 * mX1 + mX2 * mX2 + mX3 * mX3 );
	}

    T mag2() const {
		return ( mX1 * mX1 + mX2 * mX2 + mX3 * mX3 );
	}

    T pseudoRapidity() const {
		// change code to more optimal:
		// double m = mag();
		// return 0.5*::log( ( m + z() ) / ( m - z() ) );
		double tmp = std::tan( theta() / 2. ); if ( tmp <= 0. ) { return 1e20; }
		return -std::log( tmp );
	}

	TVector<T>& operator = ( const TVector<T>& v ) {
		mX1 = v.x();  mX2 = v.y();  mX3 = v.z();
		return *this;
	}

    T operator() ( size_t i ) const {
		return (&mX1)[i];
	}

	T operator[] ( size_t i ) const {
		return (&mX1)[i];
	}

    T& operator() ( size_t i ) {
		return (&mX1)[i];
	}

    T& operator[] ( size_t i ) {
		return (&mX1)[i];
	}
    
    T massHypothesis( T mass ) const {
    	return std::sqrt( ( *this ) * ( *this ) + mass * mass );
	}
    
    TVector<T>  unit() const {
		T tmp = mag(); if ( tmp <= 0. ) tmp = 1e-20;
		return *this / tmp;
	}

    TVector<T>  orthogonal() const {
		T X = ( mX1 < 0.0 ) ? -mX1 : mX1;
		T Y = ( mX2 < 0.0 ) ? -mX2 : mX2;
		T Z = ( mX3 < 0.0 ) ? -mX3 : mX3;
		if ( X < Y ) {
			return X < Z ? TVector<T>( 0, mX3, -mX2 ) :  TVector<T>( mX2, -mX1, 0 );
		}
		return  mX2 < mX3 ? TVector<T>( -mX3, 0, mX1 ) :  TVector<T>( mX2, -mX1, 0 );
	}

    void  rotateX( T Angle ) {
		T yPrime = cos( Angle ) * mX2 - sin( Angle ) * mX3;
		T zPrime = sin( Angle ) * mX2 + cos( Angle ) * mX3;
		mX2 = yPrime;
		mX3 = zPrime;
	}

    void  rotateY( T Angle ) {
		T zPrime = cos( Angle ) * mX3 - sin( Angle ) * mX1;
		T xPrime = sin( Angle ) * mX3 + cos( Angle ) * mX1;
		mX1 = xPrime;
		mX3 = zPrime;
	}

    void  rotateZ( T Angle ) {
		T xPrime = cos( Angle ) * mX1 - sin( Angle ) * mX2;
		T yPrime = sin( Angle ) * mX1 + cos( Angle ) * mX2;
		mX1 = xPrime;
		mX2 = yPrime;
	}
    
    TVector<T>  operator - () {
		return TVector<T>( -mX1, -mX2, -mX3 );
	}

    TVector<T>  operator + () {
		return *this;
	}

	TVector<T>& operator *= ( T c ) {
		mX1 *= c; mX2 *= c; mX3 *= c;
		return *this;
	}

    TVector<T>& operator /= ( T c ) {
		mX1 /= c; mX2 /= c; mX3 /= c;
		return *this;
	}

    TVector<T>  pseudoProduct( T X, T Y, T Z ) const {
		return TVector<T>( mX1 * X, mX2 * Y, mX3 * Z );
	}
 
    T angle(const TVector<T>& vec ) const {
		T norm = this->mag2() * vec.mag2(); 
		return norm > 0 ? std::acos( this->dot( vec ) / ( std::sqrt( norm ) ) ) : 0;
	}

    TVector<T> cross( const TVector<T>& v ) const {
		return TVector<T>(
			mX2 * v.z() - mX3 * v.y(),
			mX3 * v.x() - mX1 * v.z(),
			mX1 * v.y() - mX2 * v.x()
		);
	}

    T dot( const TVector<T>& v ) const {
		return mX1 * v.x() + mX2 * v.y() + mX3 * v.z();
	}

    TVector<T> pseudoProduct( const TVector<T>& v ) const {
		return this->pseudoProduct( v.x(), v.y(), v.z() );
	}
    
	bool operator == ( const TVector<T>& v) const {
		return mX1 == v.x() && mX2 == v.y() && mX3 == v.z();
	}

    bool operator != ( const TVector<T>& v ) const {
		return !(*this == v);
	}

	TVector<T>& operator += ( const TVector<T>& v ) {
		mX1 += v.x(); mX2 += v.y(); mX3 += v.z();
		return *this;
	}

    TVector<T>& operator -= ( const TVector<T>& v ) {
		mX1 -= v.x(); mX2 -= v.y(); mX3 -= v.z();
		return *this;
	}

	int valid( T world = 1.e+5 ) const {
		return !bad( world );
	}

	int bad( T world = 1.e+5 ) const {
		for ( size_t i = 0; i < 3; i++ ) {
			if ( !std::isfinite( ( &mX1 )[ i ] ) ) { return 10 + i; }
			if ( std::fabs( ( &mX1 )[ i ]) > world ) { return 20 + i; }
		}
		return 0;
	}

	protected:

    	T mX1, mX2, mX3;

};

template<class T>
inline T abs( const TVector<T>& v ) { return v.mag(); }

template<class T>
inline TVector<T> cross_product( const TVector<T>& v1, const TVector<T>& v2 ) { return v1.cross(v2); }

template<class T>
inline TVector<T>
operator+ ( const TVector<T>& v1, const TVector<T>& v2 ) { return TVector<T>(v1) += v2; }

template<class T>
inline TVector<T>
operator- (const TVector<T>& v1, const TVector<T>& v2) { return TVector<T>(v1) -= v2; }

template<class T>
inline T operator* (const TVector<T>& v1, const TVector<T>& v2 ) { return TVector<T>(v1).dot(v2); }

template<class T>
inline TVector<T> operator* (const TVector<T>& v, double c ) { return TVector<T>(v) *= c; }

template<class T>
inline TVector<T> operator* (T c, const TVector<T>& v ) { return TVector<T>(v) *= c; }

template<class T>
inline TVector<T> operator/ (const TVector<T>& v, double c ) { return TVector<T>(v) /= c; }


template<class T>
class Helix {

public:

    /// curvature, dip angle, phase, origin, h
    Helix( T c, T dip, T phase, const TVector<T>& o, int h = -1 ) { setParameters( c, dip, phase, o, h ); }

	// momentum, origin, magnetic field, charge
	Helix( const TVector<T>& p, const TVector<T>& o, T B, T q ) {
		mH = ( q * B <= 0 ) ? 1 : -1;
		if ( p.y() == 0 && p.x() == 0 ) {
			setPhase( ( M_PI / 4 ) * ( 1 - 2. * mH ) );
		} else {
			setPhase( std::atan2( p.y(), p.x() ) - mH * M_PI / 2 );
			setDipAngle( std::atan2( p.z(), p.perp() ) );
			mOrigin = o;
			setCurvature( std::fabs( ( c_light * nanosecond / meter * q * B / tesla ) / ( abs( p ) / GeV * mCosDipAngle ) / meter ) );
		}
	}

	virtual ~Helix() {}

    T dipAngle()  const { return mDipAngle;  }
    T curvature() const { return mCurvature; }
    T phase()     const { return mPhase; }
	T xcenter()   const { if ( mSingularity ) { return 0; } else { return mOrigin.x() - mCosPhase / mCurvature; } }
	T ycenter()   const { if ( mSingularity ) { return 0; } else { return mOrigin.y() - mSinPhase / mCurvature; } }
	T h()         const { return mH; }

    const TVector<T>& origin() const { return mOrigin; }

	TVector<double> momentum( T B ) const {
		if ( mSingularity ) {
			return ( TVector<T>(0,0,0) );
		}
		T pt = GeV * std::fabs( c_light * nanosecond / meter * B / tesla ) / ( std::fabs( mCurvature ) * meter);
        return ( TVector<T>( pt * std::cos( mPhase + mH * M_PI / 2 ),
			pt * sin( mPhase + mH * M_PI / 2 ),
			pt * tan( mDipAngle ) ) );
	}

	TVector<T> momentumAt(T s, T B) const {
		Helix tmp( *this );
		tmp.moveOrigin( s );
		return tmp.momentum( B );
	}

	int charge( T B ) const {
    	return ( B > 0 ? -mH : mH );
	}

	T geometricSignedDistance( T x, T y) {
		// Geometric signed distance
		T thePath = this->pathLength( x, y );
		TVector<T> DCA2dPosition = this->at( thePath );
		DCA2dPosition.setZ(0);
		TVector<T> position(x,y,0);
		TVector<T> DCAVec = ( DCA2dPosition - position );
		TVector<T> momVec;
		// Deal with straight tracks
		if ( this->mSingularity ) {
			momVec = this->at(1)- this->at(0);
			momVec.setZ(0);
		} else {
			momVec = this->momentumAt( thePath, 1. / tesla );
			momVec.setZ(0);
		}
    
		T cross = DCAVec.x() * momVec.y() - DCAVec.y() * momVec.x();
		T theSign = ( cross >= 0 ) ? 1. : -1.;
		return theSign * DCAVec.perp();
	}

	T curvatureSignedDistance( T x, T y ) {
		if ( this->mSingularity || std::abs( this->mH ) <= 0 ) {
			return ( this->geometricSignedDistance( x, y ) );
		}
		return ( this->geometricSignedDistance( x, y ) ) / ( this->mH );
	}

	T geometricSignedDistance( const TVector<T>& pos ) {
		T sdca2d = this->geometricSignedDistance( pos.x(), pos.y() );
		T theSign = ( sdca2d >= 0 ) ? 1. : -1.;
		return ( this->distance( pos ) ) * theSign;
	}

	T curvatureSignedDistance( const TVector<T>& pos ) {
		T sdca2d = this->curvatureSignedDistance( pos.x(), pos.y() );
		T theSign = ( sdca2d >= 0 ) ? 1. : -1.;
		return ( this->distance( pos ) ) * theSign;
	}

    void setParameters( T c, T dip, T phase, const TVector<T>& o, int h ) {
		mH = ( h >= 0 ) ? 1 : -1;
		mOrigin = o;
		setDipAngle(dip);
		setPhase(phase);
	    setCurvature(c);
		if ( mSingularity && mH == -1 ) { mH = +1; setPhase(mPhase-M_PI); }
	}

    /// coordinates of helix at point s
    T x( T s ) const {
    	if ( mSingularity ) {
			return mOrigin.x() - s * mCosDipAngle * mSinPhase;
		}
		return mOrigin.x() + ( std::cos( mPhase + s * mH * mCurvature * mCosDipAngle ) - mCosPhase ) / mCurvature;
	}

    T y( T s ) const {
		if ( mSingularity ) {
			return mOrigin.y() + s * mCosDipAngle * mCosPhase;
		}
		return mOrigin.y() + ( std::sin( mPhase + s * mH * mCurvature * mCosDipAngle ) - mSinPhase ) / mCurvature;
	}

    T z( T s ) const {
		return mOrigin.z() + s * mSinDipAngle;
	}

    TVector<T>  at( T s ) const {
		return TVector<T>( x(s), y(s), z(s) );
	}

    /// pointing vector of helix at point s
    T cx( T s = 0 ) const {
		if ( mSingularity ) {
			return -mCosDipAngle*mSinPhase;
		}
		return -std::sin( mPhase + s * mH * mCurvature * mCosDipAngle ) * mH * mCosDipAngle;
	}

    T cy( T s = 0 ) const {
		if ( mSingularity ) {
			return mCosDipAngle * mCosPhase;
		}
		return std::cos( mPhase + s * mH * mCurvature * mCosDipAngle ) * mH * mCosDipAngle;
	}

    T cz( T s = 0 ) const {
		return mSinDipAngle;
	}
    
    TVector<T> cat( T s ) const {
		return TVector<T>( cx(s), cy(s), cz(s) );
	}

    /// returns period length of helix
    T period() const {
		if ( mSingularity ) {
            return std::numeric_limits<T>::max();
		}
		return std::fabs( 2 * M_PI / ( mH * mCurvature * mCosDipAngle ) );
	}
    
    /// path length at given r (cylindrical r)
    std::pair<T, T> pathLength( T r ) const {

		std::pair<T,T> value;
    	std::pair<T,T> VALUE(999999999.,999999999.);

		if ( mSingularity ) {
			std::cerr << "   singularity!" << std::endl;
			T t1 = mCosDipAngle * ( mOrigin.x() * mSinPhase - mOrigin.y() * mCosPhase );
			T t12 = mOrigin.y() * mOrigin.y();
			T t13 = mCosPhase * mCosPhase;
			T t15 = r * r;
			T t16 = mOrigin.x() * mOrigin.x();
			T t20 = -mCosDipAngle * mCosDipAngle * ( 2.0 * mOrigin.x() * mSinPhase * mOrigin.y() * mCosPhase +
				t12 - t12 * t13 - t15 + t13 * t16 );
			if ( t20 < 0.) { return VALUE; }
			t20 = std::sqrt( t20 );
			value.first  = ( t1 - t20 ) / ( mCosDipAngle * mCosDipAngle );
			value.second = ( t1 + t20 ) / ( mCosDipAngle * mCosDipAngle );
    	} else {
			//std::cout << "   no singularity!" << std::endl;
			T t1 = mOrigin.y() * mCurvature;
			T t2 = mSinPhase;
			T t3 = mCurvature * mCurvature;
			T t4 = mOrigin.y() * t2;
			T t5 = mCosPhase;
			T t6 = mOrigin.x() * t5;
			T t8 = mOrigin.x() * mOrigin.x();
			T t11 = mOrigin.y() * mOrigin.y();
			T t14 = r * r;
			T t15 = t14 * mCurvature;
			T t17 = t8 * t8;
			T t19 = t11 * t11;
			T t21 = t11 * t3;
			T t23 = t5 * t5;
			T t32 = t14 * t14;
			T t35 = t14 * t3;
			T t38 = 8.0 * t4 * t6 - 4.0 * t1 * t2 * t8 - 4.0 * t11 * mCurvature * t6 +
                 4.0 * t15 * t6 + t17 * t3 + t19 * t3 + 2.0 * t21 * t8 + 4.0 * t8 * t23 -
                 4.0 * t8 * mOrigin.x() * mCurvature * t5 - 4.0 * t11 * t23 -
                 4.0 * t11 * mOrigin.y() * mCurvature * t2 + 4.0 * t11 - 4.0 * t14 +
                 t32 * t3 + 4.0 * t15 * t4 - 2.0 * t35 * t11 - 2.0 * t35 * t8;
			T t40 = ( -t3 * t38 );
			if ( t40 < 0. ) {
				std::cerr << "t40 < 0." << std::endl;
				return VALUE;
			}
			t40 = std::sqrt( t40 );

			T t43 = mOrigin.x() * mCurvature;
			T t45 = 2.0 * t5 - t35 + t21 + 2.0 - 2.0 * t1 * t2 -2.0 * t43 - 2.0 * t43 * t5 + t8 * t3;
			T t46 = mH * mCosDipAngle * mCurvature;

			value.first = ( -mPhase + 2.0 * std::atan( ( -2.0 * t1 + 2.0 * t2 + t40 ) / t45 ) ) / t46;
			value.second = -( mPhase + 2.0 * std::atan( ( 2.0 * t1 - 2.0 * t2 + t40 ) / t45 ) ) / t46;

			T p = period();
			if ( !std::isnan( value.first ) ) {
				if ( std::fabs( value.first - p ) < std::fabs( value.first ) ) { value.first = value.first - p; }
				else if ( std::fabs( value.first + p ) < std::fabs( value.first ) ) { value.first = value.first + p; }
    		} else {
				std::cerr << "value.first = isnan" << std::endl;
			}

			if ( !std::isnan( value.second ) ) {
				if ( std::fabs( value.second - p ) < std::fabs( value.second ) ) { value.second = value.second - p; }
				else if (fabs(value.second+p) < fabs(value.second)) value.second = value.second+p;
    		} else {
				std::cerr << "value.second = isnan" << std::endl;
			}

    	}

    	if ( value.first > value.second ) {
    		std::swap( value.first, value.second );
		}

    	return(value);
	}
    
    /// path length at given r (cylindrical r, cylinder axis at x,y)
    std::pair<T, T> pathLength( T r, T x, T y) {
		//std::cout << "--- pathlength at r, x, y --- " << std::endl;
		T x0 = mOrigin.x();
		T y0 = mOrigin.y();
		mOrigin.setX( x0 - x );
		mOrigin.setY( y0 - y );
		std::pair<T, T> result = this->pathLength( r );
		mOrigin.setX( x0 );
		mOrigin.setY( y0 );
		return result;
	}
    
    /// path length at distance of closest approach to a given point
    T pathLength( const TVector<T>& p, bool scanPeriods = true ) const {
		T s;
		T dx = p.x()-mOrigin.x();
		T dy = p.y()-mOrigin.y();
		T dz = p.z()-mOrigin.z();  
		if ( mSingularity ) {
			s = mCosDipAngle * ( mCosPhase * dy - mSinPhase * dx ) + mSinDipAngle * dz;
		} else {
			const T MaxPrecisionNeeded = 0.0001;
			const int MaxIterations = 100;

			T t34 = mCurvature * mCosDipAngle * mCosDipAngle;
			T t41 = mSinDipAngle * mSinDipAngle;
			T t6, t7, t11, t12, t19;

			s = fudgePathLength( p );

			if ( scanPeriods ) {
				T ds = period();
				int j, jmin = 0;
				T d, dmin = abs(at(s) - p);
				for( j = 1; j < MaxIterations; j++ ) {
					if ( ( d = abs( at( s + j * ds ) - p ) ) < dmin ) {
						dmin = d;
						jmin = j;
					} else {
						break;
					}
            	}
				for( j = -1; -j < MaxIterations; j-- ) {
					if ( ( d = abs( at( s + j * ds ) - p ) ) < dmin ) {
						dmin = d;
						jmin = j;
					} else {
						break;
					}
            	}
				if ( jmin ) { s += jmin * ds; }
        	}

			T sOld = s;
			for ( int i = 0; i < MaxIterations; i++ ) {
				t6  = mPhase + s * mH * mCurvature * mCosDipAngle;
				t7  = std::cos( t6 );
				t11 = dx - ( 1 / mCurvature ) * ( t7 - mCosPhase );
				t12 = std::sin( t6 );
				t19 = dy - ( 1 / mCurvature ) * ( t12 - mSinPhase );
				s -= ( t11 * t12 * mH * mCosDipAngle - t19 * t7 * mH * mCosDipAngle -
					( dz - s * mSinDipAngle ) * mSinDipAngle ) /
					( t12 * t12 * mCosDipAngle * mCosDipAngle + t11 * t7 * t34 +
						t7 * t7 * mCosDipAngle * mCosDipAngle + t19 * t12 * t34 + t41 );
				if ( std::fabs( sOld - s ) < MaxPrecisionNeeded ) { break; }
				sOld = s;
			}
		}
		return s;
	}
    
    /// path length at intersection with plane
    T pathLength( const TVector<T>& r, const TVector<T>& n ) const {

		T s;
        
		if ( mSingularity ) {
			T t = n.z() * mSinDipAngle +
               n.y() * mCosDipAngle * mCosPhase -
               n.x() * mCosDipAngle * mSinPhase;
			if (t == 0) {
				s = NoSolution;
			} else {
				s = ((r - mOrigin)*n)/t;
			}
    	} else {
			const T MaxPrecisionNeeded = 0.0001;
			const int MaxIterations = 20;

			T A = mCurvature * ( ( mOrigin - r ) * n ) - n.x() * mCosPhase - n.y() * mSinPhase;
			T t = mH * mCurvature * mCosDipAngle;
			T u = n.z() * mCurvature * mSinDipAngle;

			T a, f, fp;
			T sOld = s = 0;
			T shiftOld = 0;
			T shift;
			const T angMax = 0.21;
			T deltas = std::fabs( angMax / ( mCurvature * mCosDipAngle ) );

			int i;
			for ( i = 0; i < MaxIterations; i++ ) {
				a  = t * s + mPhase;
				T sina = std::sin( a );
				T cosa = std::cos( a );
				f  = A + n.x() * cosa + n.y() * sina + u*s;
				fp = -n.x() * sina * t + n.y() * cosa * t + u;
				if ( std::fabs( fp ) * deltas <= std::fabs( f ) ) { //too big step
					int sgn = 1;
					if ( fp < 0. ) { sgn = -sgn; }
					if ( f < 0.) { sgn = -sgn; }
					shift = sgn * deltas;
					if ( shift < 0 ) { shift *= 0.9; } // don't get stuck shifting +/-deltas
				} else {
					shift = f/fp;
				}
				s -= shift;
				shiftOld = shift;
				if ( std::fabs( sOld - s ) < MaxPrecisionNeeded ) { break; }
				sOld = s;
			}
			if ( i == MaxIterations ) { return NoSolution; }
    	}
		return s;
	}

    /// path length at distance of closest approach in the xy-plane to a given point
    T pathLength( T x, T y ) const {
		return fudgePathLength( TVector<T>( x, y, 0 ) );
	}

    /// path lengths at dca between two helices 
    std::pair<T, T> pathLengths( const Helix&, T minStepSize = 10 * 0.0001, T minRange = 10 ) const {
		if ( mSingularity != h.mSingularity ) {
        	return std::pair<T,T>( NoSolution, NoSolution );
		}
    
		T s1, s2;
 
		if ( mSingularity ) {
			TVector<T> dv = h.mOrigin - mOrigin;
			TVector<T> a(-mCosDipAngle * mSinPhase, mCosDipAngle * mCosPhase, mSinDipAngle );
			TVector<T> b(-h.mCosDipAngle * h.mSinPhase, h.mCosDipAngle * h.mCosPhase, h.mSinDipAngle );
			T ab = a * b;
			T g  = dv * a;
			T k  = dv * b;
			s2 = ( k - ab * g ) / ( ab * ab - 1. );
			s1 = g + s2 * ab;
			return std::pair<T,T>(s1, s2);
		} else {
			T dx = h.xcenter() - xcenter();
			T dy = h.ycenter() - ycenter();
			T dd = std::sqrt( dx * dx + dy * dy );
			T r1 = 1 / curvature();
			T r2 = 1 / h.curvature();
			T cosAlpha = ( r1 * r1 + dd * dd - r2 * r2 ) / ( 2 * r1 * dd );

			T s;
			T x, y;
			if ( std::fabs( cosAlpha ) < 1 ) {
				T sinAlpha = std::sin( std::acos( cosAlpha ) );
				x = xcenter() + r1 * ( cosAlpha * dx - sinAlpha * dy ) / dd;
				y = ycenter() + r1 * ( sinAlpha * dx + cosAlpha * dy ) / dd;
				s = pathLength( x, y );
				x = xcenter() + r1 * ( cosAlpha * dx + sinAlpha * dy ) / dd;
				y = ycenter() + r1 * ( cosAlpha * dy - sinAlpha * dx ) / dd;
				T a = pathLength(x, y);
				if ( h.distance( at( a ) ) < h.distance( at( s ) ) ) { s = a; }
        	} else {
				int rsign = ( ( r2 - r1 ) > dd ? -1 : 1 );
				x = xcenter() + rsign * r1 * dx / dd;
				y = ycenter() + rsign * r1 * dy / dd;
				s = pathLength(x, y);
			}

			T dmin = h.distance( at( s ) );
			T range = std::max( 2 * dmin, minRange );
			T ds = range / 10;
			T slast = -999999, ss, d;
			s1 = s - range / 2.;
			s2 = s + range / 2.;

			while ( ds > minStepSize ) {
				for ( ss = s1; ss < s2 + ds; ss += ds ) {
					d = h.distance( at( ss ) );
					if ( d < dmin ) {
						dmin = d;
						s = ss;
					}
					slast = ss;
				}
				if ( s == s1 ) {
					d = 0.8 * ( s2 - s1 );
					s1 -= d;
					s2 -= d;
				} else if ( s == slast ) {
					d = 0.8  *( s2 - s1 );
					s1 += d;
					s2 += d;
				} else {
					s1 = s - ds;
					s2 = s + ds;
					ds /= 10;
				}
			}
			return std::pair<T, T>( s, h.pathLength( at( s ) ) );
		}
	}
    
    /// minimal distance between point and helix
    T distance( const TVector<T>& p, bool scanPeriods = true ) const {
		return std::abs( this->at( pathLength( p, scanPeriods ) ) - p );
	}
    
    /// checks for valid parametrization
    bool valid( T world = 1.e+5 ) const {
		return !bad(world);
	}

    int bad( T WorldSize = 1.e+5 ) const {
		int ierr;
    	if ( !std::isfinite( mDipAngle ) )  { return 11; }
		if ( !std::isfinite( mCurvature ) ) { return 12; }
		ierr = mOrigin.bad( WorldSize );
		if ( ierr ) { return 3+ierr*100; }
		if ( std::fabs( mDipAngle ) > 1.58 ) { return 21; }
		T qwe = std::fabs( std::fabs( mDipAngle ) - M_PI / 2 );
		if ( qwe < 1. / WorldSize ) { return 31; }
		if ( std::fabs( mCurvature ) > WorldSize ) { return 22; }
		if ( mCurvature < 0 ) { return 32; }
		if ( std::abs( mH ) != 1 ) { return 24; }
		return 0;
	}
    
    /// move the origin along the helix to s which becomes then s=0
    virtual void moveOrigin( T s ) {
		if ( mSingularity ) {
			mOrigin = at( s );
    	} else {
			TVector<T> newOrigin = at(s);
			T newPhase = std::atan2( newOrigin.y() - ycenter(), newOrigin.x() - xcenter() );
			mOrigin = newOrigin;
			setPhase( newPhase );
		}
	}
    
    static constexpr T NoSolution = 3.e+33;
    
protected:
    Helix() {}
    
    void setCurvature( T val ) {
		if ( val < 0 ) {
			mCurvature = -val;
			mH = -mH;
			setPhase( mPhase + M_PI );
		} else {
			mCurvature = val;
		}                    
		if ( std::fabs( mCurvature ) <= std::numeric_limits<T>::epsilon() ) {
			mSingularity = true; // straight line
		} else {
			mSingularity = false; // curved
		}
	}

    void setPhase( T val ) {
    	mPhase       = val;
    	mCosPhase    = std::cos( mPhase );
    	mSinPhase    = std::sin( mPhase );
    	if ( std::fabs( mPhase ) > M_PI ) {
			mPhase = atan2(mSinPhase, mCosPhase);  // force range [-pi,pi]
		}
	}

    void setDipAngle( T val ) {
		mDipAngle    = val;
		mCosDipAngle = std::cos( mDipAngle );
		mSinDipAngle = std::sin( mDipAngle );
	}

    /// value of S where distance in x-y plane is minimal
	T fudgePathLength( const TVector<T>& p ) const {
		T s;
		T dx = p.x() - mOrigin.x();
		T dy = p.y() - mOrigin.y();
		if ( mSingularity ) {
			s = ( dy * mCosPhase - dx * mSinPhase ) / mCosDipAngle;
		} else {
			s = std::atan2( dy * mCosPhase - dx * mSinPhase, 1 / mCurvature + dx * mCosPhase + dy * mSinPhase ) / ( mH * mCurvature * mCosDipAngle );
		}
		return s;
	}

protected:
    bool mSingularity = false;	// true for straight line case (B=0)
    TVector<T> mOrigin;
    T mDipAngle = 0;
    T mCurvature = 0;
    T mPhase = 0;
    int mH = 0;	// -sign(q*B);

    T mCosDipAngle = 0;
    T mSinDipAngle = 0;
    T mCosPhase = 0;
    T mSinPhase = 0;

	static constexpr T meter = 100.0;
	static constexpr T meter2 = meter * meter;
	static constexpr T second = 1.0;
	static constexpr T nanosecond = 1.e-9;
	static constexpr T GeV = 1.0;
	static constexpr T MeV = 1.e-3;
	static constexpr T volt = 1.e-6 * MeV;
	static constexpr T c_light = 2.99792458e+8 * meter / second; // c_light * meter / sec
	static constexpr T tesla = 1.0;

};



template <class T>
class Data
{

public:
	Data( const std::vector<std::vector<T>>& hits_ ) : meanX(0), meanY(0)
	{
		std::copy( hits_.begin(), hits_.end(), back_inserter( hits ) );
		hits_size = hits.size();
	}

	void means()
	{
		meanX = 0.;
		meanY = 0.;
		for ( size_t i = 0; i < hits_size; i++ ) {
			meanX += hits[i][0];
			meanY += hits[i][1];
		}
		meanX /= hits_size;
		meanY /= hits_size;
	}

	void center(void)
	{
		T sX = 0., sY = 0.;
		for ( size_t i = 0; i < hits_size; i++ ) {
			sX += hits[i][0];
			sY += hits[i][1];
		}
		sX /= hits_size;
		sY /= hits_size;
		for ( size_t i = 0; i < hits_size; i++ ) {
			hits[i][0] -= sX;
			hits[i][1] -= sY;
		}
		meanX = 0.;
		meanY = 0.;
	}

	void scale(void)
	{
		T sXX = 0., sYY = 0., scaling;
		for ( size_t i = 0; i < hits_size; i++ ) {
			sXX += hits[i][0] * hits[i][0];
			sYY += hits[i][1] * hits[i][1];
		}
		scaling = std::sqrt( ( sXX + sYY ) / hits_size / 2.0 );
		for ( size_t i = 0; i < hits_size; i++ ) {
			hits[i][0] /= scaling;
			hits[i][1] /= scaling;
		}
	}

	std::vector<std::vector<T>> hits;
	size_t hits_size;
	T meanX;
	T meanY;
};

template <typename T>
struct Circle {

	Circle() : a(0), b(0), r(0), s(0), g(0), Gx(0), Gy(0), i(0), j(0), chi2(99999.0) {}

	void initFromCircle( Circle<T>* c ) {
		a = c->a;
		b = c->b;
		r = c->r;
		s = c->s;
		g = c->g;
		Gx = c->Gx;
		Gy = c->Gy;
		i = c->i;
		j = c->j;
		chi2 = c->chi2;
	}

	bool is_good()
	{
		return !std::isnan(a) && !isnan(b) && !isnan(r) &&
		       std::isfinite(a) && std::isfinite(b) && std::isfinite(r) &&
		       std::isfinite(j) && j > 0;
	}

	T a;    // x-center
	T b;    // y-center
	T r;    // radius
	T s;    // sigma
	T g;
	T Gx;
	T Gy;
	size_t i; // inner loop iterations
	size_t j; // outer loop iterations
	T chi2;
};

template <typename T>
struct Line {
	Line() : a( std::numeric_limits<T>::infinity() ), b( std::numeric_limits<T>::infinity() ),
			 			siga(0), sigb(0), chi2( std::numeric_limits<T>::infinity() ) {}

	Line( T a_, T b_, T siga_, T sigb_, T chi2_ ) : a( a_ ), b( b_ ), siga( siga_ ), sigb( sigb_ ), chi2( chi2_ ) {}

	bool is_good()
	{
		return !std::isnan(a) && !isnan(b) && !isnan(chi2) &&
		       std::isfinite(a) && std::isfinite(b) && std::isfinite(chi2);
	}
	T a;    // a + b * x
	T b;    //
	T siga; // sigma a
	T sigb; // sigma b
	T chi2;
};

template <typename T>
T DistanceToCircle( Circle<T>* circle, T x, T y )
{
	return std::fabs( std::sqrt( std::pow( x - circle->a, 2 ) + std::pow( y - circle->b, 2 ) ) - circle->r );
}

template <typename T>
T DistanceToCircle( T a, T b, T r, T x, T y )
{
	return std::fabs( std::sqrt( std::pow( x - a, 2 ) + std::pow( y - b, 2 ) ) - r );
}

template <typename T>
T DistanceToLine( T a, T b, T x, T y )
{
	// point-slope line : 1 * y = a + b * x
	// standard line    : b * x - 1 * y + a = 0
	return std::fabs( ( b * x ) + ( -1 * y ) + a ) / std::sqrt( std::pow( b, 2 ) + std::pow( -1, 2 ) );
}


template<class T>
class CircleFit {
            
    public:

        CircleFit() {}

        static Circle<T>* RegularFit( const std::vector<std::vector<T>>& hits, size_t fit = 0 ) {
			switch( fit ) {
				case 2: 
					return CircleFitByTaubin( hits );
					break;
				case 3: 
					return CircleFitByPratt( hits );
					break;
				case 1: 
				default:
					return CircleFitByHyper( hits );
					break;
			}
		}

        static Circle<T>* RobustFit( const std::vector<std::vector<T>>& hits, Circle<T>* circle = 0 /* estimate */, T lambda = 0.01 ) {
			CircleFitByChernovHoussam( hits, circle, lambda );
			return circle;
		}

	private:

		static T Sigma ( Data<T>& data, Circle<T>* circle ) {
			T sum = 0., dx, dy;
			for ( size_t i = 0; i < data.hits_size; i++ ) {
				dx = data.hits[i][0] - circle->a;
				dy = data.hits[i][1] - circle->b;
				sum += std::pow( std::sqrt( dx * dx + dy * dy ) - circle->r, 2 );
			}
			return std::sqrt( sum / data.hits_size );
		}

		static T ChiSqr ( Data<T>& data, Circle<T>* circle ) {
			T sum = 0., dx, dy;
			for ( size_t i = 0; i < data.hits_size; i++ ) {
				dx = data.hits[i][0] - circle->a;
				dy = data.hits[i][1] - circle->b;
				sum += std::pow( std::sqrt( dx * dx + dy * dy ) - circle->r, 2 );
			}
			return ( sum / ( data.hits_size - 3 ) );
		}

		static Circle<T>* CircleFitByHyper( const std::vector<std::vector<T>>& hits ) {
			// From: http://people.cas.uab.edu/~mosya/cl/CPPcircle.html
			size_t i, iter, IterMAX = 99;

			T Xi, Yi, Zi;
			T Mz, Mxy, Mxx, Myy, Mxz, Myz, Mzz, Cov_xy, Var_z;
			T A0, A1, A2, A22;
			T Dy, xnew, x, ynew, y;
			T DET, Xcenter, Ycenter;

			Circle<T>* circle = new Circle<T>();
			circle->a = std::numeric_limits<T>::infinity();
			circle->b = std::numeric_limits<T>::infinity();
			circle->r = std::numeric_limits<T>::infinity();
			circle->s = 0;
			circle->j = 0;

			Data<T> data( hits );
			data.means();

			Mxx = Myy = Mxy = Mxz = Myz = Mzz = 0.;

			for ( i = 0; i < data.hits_size; i++ ) {
				Xi = data.hits[i][0] - data.meanX;   //  centered x-coordinates
				Yi = data.hits[i][1] - data.meanY;   //  centered y-coordinates
				Zi = Xi * Xi + Yi * Yi;
				Mxy += Xi * Yi;
				Mxx += Xi * Xi;
				Myy += Yi * Yi;
				Mxz += Xi * Zi;
				Myz += Yi * Zi;
				Mzz += Zi * Zi;
			}
			Mxx /= data.hits_size;
			Myy /= data.hits_size;
			Mxy /= data.hits_size;
			Mxz /= data.hits_size;
			Myz /= data.hits_size;
			Mzz /= data.hits_size;

			Mz = Mxx + Myy;
			Cov_xy = Mxx * Myy - Mxy * Mxy;
			Var_z = Mzz - Mz * Mz;

			A2 = 4.0 * Cov_xy - 3.0 * Mz * Mz - Mzz;
			A1 = Var_z * Mz + 4.0 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz;
			A0 = Mxz * ( Mxz * Myy - Myz * Mxy ) + Myz * ( Myz * Mxx - Mxz * Mxy ) - Var_z * Cov_xy;
			A22 = A2 + A2;

			for ( x = 0., y = A0, iter = 0; iter < IterMAX; iter++ ) {  // usually, 4-6 iterations are enough
				Dy = A1 + x * ( A22 + 16. * x * x );
				xnew = x - y / Dy;
				if ( ( xnew == x ) || ( !std::isfinite( xnew ) ) ) {
					break;
				}
				ynew = A0 + xnew * ( A1 + xnew * ( A2 + 4.0 * xnew * xnew ) );
				if ( std::abs( ynew ) >= std::abs( y ) ) {
					break;
				}
				x = xnew;
				y = ynew;
			}

			DET = x * x - x * Mz + Cov_xy;
			Xcenter = ( Mxz * ( Myy - x ) - Myz * Mxy ) / DET / 2.0;
			Ycenter = ( Myz * ( Mxx - x ) - Mxz * Mxy ) / DET / 2.0;

			circle->a = Xcenter + data.meanX;
			circle->b = Ycenter + data.meanY;
			circle->r = std::sqrt( Xcenter * Xcenter + Ycenter * Ycenter + Mz - x - x );
			circle->s = Sigma( data, circle );
			circle->j = iter;
			circle->chi2 = ChiSqr( data, circle );

			return circle;
		}

		static Circle<T>* CircleFitByTaubin( const std::vector<std::vector<T>>& hits ) {
			// From: http://people.cas.uab.edu/~mosya/cl/CPPcircle.html
			size_t i, iter, IterMAX = 99;

			T Xi, Yi, Zi;
			T Mz, Mxy, Mxx, Myy, Mxz, Myz, Mzz, Cov_xy, Var_z;
			T A0, A1, A2, A22, A3, A33;
			T Dy, xnew, x, ynew, y;
			T DET, Xcenter, Ycenter;

			Circle<T>* circle = new Circle<T>();

			Data<T> data( hits );
			data.means();

			Mxx = Myy = Mxy = Mxz = Myz = Mzz = 0.;

			for ( i = 0; i < data.hits_size; i++ ) {
				Xi = data.hits[i][0] - data.meanX;   //  centered x-coordinates
				Yi = data.hits[i][0] - data.meanY;   //  centered y-coordinates
				Zi = Xi * Xi + Yi * Yi;

				Mxy += Xi * Yi;
				Mxx += Xi * Xi;
				Myy += Yi * Yi;
				Mxz += Xi * Zi;
				Myz += Yi * Zi;
				Mzz += Zi * Zi;
			}
			Mxx /= data.hits_size;
			Myy /= data.hits_size;
			Mxy /= data.hits_size;
			Mxz /= data.hits_size;
			Myz /= data.hits_size;
			Mzz /= data.hits_size;

			Mz = Mxx + Myy;
			Cov_xy = Mxx * Myy - Mxy * Mxy;
			Var_z = Mzz - Mz * Mz;
			A3 = 4.0 * Mz;
			A2 = -3.0 * Mz * Mz - Mzz;
			A1 = Var_z * Mz + 4.0 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz;
			A0 = Mxz * ( Mxz * Myy - Myz * Mxy ) + Myz * ( Myz * Mxx - Mxz * Mxy ) - Var_z * Cov_xy;
			A22 = A2 + A2;
			A33 = A3 + A3 + A3;

			for ( x = 0., y = A0, iter = 0; iter < IterMAX; iter++ ) {  // usually, 4-6 iterations are enough
				Dy = A1 + x * ( A22 + A33 * x );
				xnew = x - y / Dy;
				if ( ( xnew == x ) || ( !std::isfinite(xnew) ) ) {
					break;
				}
				ynew = A0 + xnew * ( A1 + xnew * ( A2 + xnew * A3 ) );
				if ( std::abs( ynew ) >= std::abs( y ) ) {
					break;
				}
				x = xnew;
				y = ynew;
			}

			DET = x * x - x * Mz + Cov_xy;
			Xcenter = ( Mxz * ( Myy - x ) - Myz * Mxy ) / DET / 2.0;
			Ycenter = ( Myz * ( Mxx - x ) - Mxz * Mxy ) / DET / 2.0;

			circle->a = Xcenter + data.meanX;
			circle->b = Ycenter + data.meanY;
			circle->r = std::sqrt( Xcenter * Xcenter + Ycenter * Ycenter + Mz );
			circle->s = Sigma( data, circle );
			circle->j = iter;
			circle->chi2 = ChiSqr( data, circle );

			return circle;
		}

		static Circle<T>* CircleFitByPratt( const std::vector<std::vector<T>>& hits ) {
			// From: http://people.cas.uab.edu/~mosya/cl/CPPcircle.html
			size_t i, iter, IterMAX = 99;

			T Xi, Yi, Zi;
			T Mz, Mxy, Mxx, Myy, Mxz, Myz, Mzz, Cov_xy, Var_z;
			T A0, A1, A2, A22;
			T Dy, xnew, x, ynew, y;
			T DET, Xcenter, Ycenter;

			Circle<T>* circle = new Circle<T>();

			Data<T> data( hits );
			data.means();

			Mxx = Myy = Mxy = Mxz = Myz = Mzz = 0.;

			for ( i = 0; i < data.hits_size; i++) {
				Xi = data.hits[i][0] - data.meanX;   //  centered x-coordinates
				Yi = data.hits[i][1] - data.meanY;   //  centered y-coordinates
				Zi = Xi * Xi + Yi * Yi;

				Mxy += Xi * Yi;
				Mxx += Xi * Xi;
				Myy += Yi * Yi;
				Mxz += Xi * Zi;
				Myz += Yi * Zi;
				Mzz += Zi * Zi;
			}
			Mxx /= data.hits_size;
			Myy /= data.hits_size;
			Mxy /= data.hits_size;
			Mxz /= data.hits_size;
			Myz /= data.hits_size;
			Mzz /= data.hits_size;

			Mz = Mxx + Myy;
			Cov_xy = Mxx * Myy - Mxy * Mxy;
			Var_z = Mzz - Mz * Mz;

			A2 = 4.0 * Cov_xy - 3.0 * Mz * Mz - Mzz;
			A1 = Var_z * Mz + 4.0 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz;
			A0 = Mxz * ( Mxz * Myy - Myz * Mxy ) + Myz * ( Myz * Mxx - Mxz * Mxy ) - Var_z * Cov_xy;
			A22 = A2 + A2;

			for ( x = 0., y = A0, iter = 0; iter < IterMAX; iter++ ) {  // usually, 4-6 iterations are enough
				Dy = A1 + x * ( A22 + 16. * x * x );
				xnew = x - y / Dy;
				if ( ( xnew == x ) || ( !std::isfinite( xnew ) ) ) {
					break;
				}
				ynew = A0 + xnew * ( A1 + xnew * ( A2 + 4.0 * xnew * xnew ) );
				if ( std::abs( ynew ) >= std::abs( y ) ) {
					break;
				}
				x = xnew;
				y = ynew;
			}

			DET = x * x - x * Mz + Cov_xy;
			Xcenter = ( Mxz * ( Myy - x ) - Myz * Mxy ) / DET / 2.0;
			Ycenter = ( Myz * ( Mxx - x ) - Mxz * Mxy ) / DET / 2.0;

			circle->a = Xcenter + data.meanX;
			circle->b = Ycenter + data.meanY;
			circle->r = std::sqrt( Xcenter * Xcenter + Ycenter * Ycenter + Mz + x + x );
			circle->s = Sigma( data, circle );
			circle->j = iter;
			circle->chi2 = ChiSqr( data, circle );

			return circle;
		}

		static T pythag( T a, T b ) {
			T absa = std::abs( a ), absb = std::abs( b );
			if ( absa > absb ) {
				return absa * std::sqrt( 1.0 + std::pow( absb / absa, 2 ) );
			}
			return ( absb == 0.0 ? 0.0 : absb * std::sqrt( 1.0 + std::pow( absa / absb, 2 ) ) );
		}

		static T OptimalRadius( Data<T>& data, Circle<T>& circle ) {
			T Mr = 0., dx, dy;
			for ( size_t i = 0; i < data.hits_size; i++ ) {
				dx = data.hits[i][0] - circle.a;
				dy = data.hits[i][1] - circle.b;
				Mr += std::sqrt( dx * dx + dy * dy );
			}
			return Mr / data.hits_size;
		}

		static void eigen2x2( T a, T b, T c, T& d1, T& d2, T& Vx, T& Vy ) {
			T disc, f;
			disc = pythag( a - b, 2.0 * c );
			d1 = ( a + b > 0. ) ? ( a + b + disc ) / 2.0 : ( a + b - disc ) / 2.0;
			d2 = ( a * b - c * c ) / d1;
			if ( std::abs( a - d1)  > std::abs( b - d1 ) ) {
				if ( ( f = pythag( c, d1 - a ) ) == 0. ) {
					Vx = 1.0;
					Vy = 0.;
					return;
				} else {
					Vx = c / f;
					Vy = ( d1 - a ) / f;
				}
			} else {
				if ( ( f = pythag( c, d1 - b ) ) == 0. ) {
					Vx = 1.0;
					Vy = 0.;
					return;
				} else {
					Vx = ( d1 - b ) / f;
					Vy = c / f;
				}
			}
			return;
		}

		static T SigmaWithLargeCircleOption( Data<T>& data, Circle<T>& circle ) {
			int i, n = data.hits_size;
			T sum = 0., dx, dy, r;
			std::vector<T> D(n);
			T LargeCircle = 10.0, a0, b0, del, s, c, x, y, z, p, t, g, W, Z;
			if ( std::abs( circle.a ) < LargeCircle && std::abs( circle.b ) < LargeCircle ) {
				for ( i = 0; i < n; i++ ) {
					dx = data.hits[i][0] - circle.a;
					dy = data.hits[i][1] - circle.b;
					D[i] = std::sqrt( dx * dx + dy * dy );
					sum += D[i];
				}
				r = sum / n;
				for ( sum = 0., i = 0; i < n; i++ ) {
					sum += std::pow( D[i] - r, 2 );
				}
				return sum / n;
			} else {
				a0 = circle.a - data.meanX;
				b0 = circle.b - data.meanY;
				del = 1.0 / std::sqrt( a0 * a0 + b0 * b0 );
				s = b0 * del;
				c = a0 * del;
				for ( W = Z = 0., i = 0; i < n; i++) {
					x = data.hits[i][0] - data.meanX;
					y = data.hits[i][1] - data.meanY;
					z = x * x + y * y;
					p = x * c + y * s;
					t = del * z - 2.0 * p;
					g = t / ( 1.0 + std::sqrt( 1.0 + del * t ) );
					W += ( z + p * g ) / ( 2.0 + del * g );
					Z += z;
				}
				W /= n;
				Z /= n;
				return Z - W * ( 2.0 + del * del * W );
			}
		}

		static void GradientHessian( Data<T>& data, Circle<T>& circle, T& F1, T& F2, T& A11, T& A22, T& A12 ) {

			int i, n = data.hits_size;
			T LargeCircle = 10.0, dx, dy, r, u, v, Mr, Mu, Mv, Muu, Mvv, Muv, Muur, Mvvr, Muvr;
			T a0, b0, del, dd, s, c, x, y, a, b, z, p, t, w, g, g1, gg1, gg2;
			T X,Y,R,U,V,Tr,W,AA,BB,AB,AG,BG,GG,UUR,VVR,UVR;

			if ( std::abs( circle.a ) < LargeCircle && std::abs( circle.b ) < LargeCircle) {
				for ( Mr = Mu = Mv = Muu = Mvv = Muv = Muur = Mvvr = Muvr = 0., i = 0; i < n; i++ ) {
					dx = data.hits[i][0] - circle.a;
					dy = data.hits[i][1] - circle.b;
					r = std::sqrt( dx * dx + dy * dy );
					u = dx / r;
					v = dy / r;
					Mr += r;
					Mu += u;
					Mv += v;
					Muu += u * u;
					Mvv += v * v;
					Muv += u * v;
					Muur += u * u / r;
					Mvvr += v * v / r;
					Muvr += u * v / r;
				}
				Mr /= n;
				Mu /= n;
				Mv /= n;
				Muu /= n;
				Mvv /= n;
				Muv /= n;
				Muur /= n;
				Mvvr /= n;
				Muvr /= n;

				F1 = circle.a + Mu * Mr - data.meanX;
				F2 = circle.b + Mv * Mr - data.meanY;

				A11 = 1.0 - Mu * Mu - Mr * Mvvr;
				A22 = 1.0 - Mv * Mv - Mr * Muur;
				A12 = -Mu * Mv + Mr * Muvr;
			} else {
				a0 = circle.a - data.meanX;
				b0 = circle.b - data.meanY;
				del = 1.0 / std::sqrt( a0 * a0 + b0 * b0 );
				dd = del * del;
				s = b0 * del;
				c = a0 * del;
				for ( X = Y = R = Tr = W = AA = BB = AB = AG = BG = GG = 0., i = 0; i < n; i++ ) {
					x = data.hits[i][0] - data.meanX;
					y = data.hits[i][1] - data.meanY;
					z = x * x + y * y;
					p = x * c + y * s;
					t = 2.0 * p - del * z;
					w = std::sqrt( 1.0 - del * t );
					g = -t / ( 1.0 + w );
					g1 = 1.0 / ( 1.0 + del * g );
					gg1 = g * g1;
					gg2 = g / ( 2.0 + del * g );
					a = ( x + g * c ) / w;
					b = ( y + g * s ) / w;
					X += x * gg1;
					Y += y * gg1;
					R += z + t * gg2;
					Tr += t * gg1;
					W += t * gg1 * gg2;
					AA += a * a * g1;
					BB += b * b * g1;
					AB += a * b * g1;
					AG += a * gg1;
					BG += b * gg1;
					GG += g * gg1;
				}
				X /= n;
				Y /= n;
				R /= n;
				Tr /= n;
				W /= n;
				AA /= n;
				BB /= n;
				AB /= n;
				AG /= n;
				BG /= n;
				GG /= n;

				U = ( Tr - del * W ) * c / 2.0 - X + R * c / 2.0;
				V = ( Tr - del * W ) * s /2.0 - Y + R * s / 2.0;

				F1 = del * ( ( dd * R * U - del * W * c + Tr * c ) / 2.0 - X );
				F2 = del * ( ( dd * R * V - del * W * s + Tr * s ) / 2.0 - Y );

				UUR = ( ( GG - R / 2.0 ) * c + 2.0 * ( AG - U ) ) * c + AA;
				VVR = ( ( GG - R / 2.0 ) * s + 2.0 * ( BG - V ) ) * s + BB;
				UVR = ( ( GG - R / 2.0 ) * c + ( AG - U ) ) * s + ( BG - V ) * c + AB;

				A11 = dd * ( U * ( 2.0 * c - dd * U ) - R * s * s / 2.0 - VVR * ( 1.0 + dd * R / 2.0 ) );
				A22 = dd * ( V * ( 2.0 * s - dd * V ) - R * c * c / 2.0 - UUR * ( 1.0 + dd * R / 2.0 ) );
				A12 = dd * ( U * s + V * c + R * s * c / 2.0 - dd * U * V + UVR * ( 1.0 + dd * R / 2.0 ) );
			}
		}

		static int CircleFitByChernovHoussam( const std::vector<std::vector<T>>& hits, Circle<T>* circleIni = 0, T LambdaIni = 0.01 ) {

			if ( circleIni == 0 ) {
				circleIni = RegularFit( hits );
			}

			Data<T> data( hits );
			int i, n = data.hits_size, iter, inner, IterMAX = 200, check_line = 1, code;

			T lambda;
			T F1,F2,A11,A22,A12,dX,dY,Mxx,Myy,Mxy,Mxxy,dx,dy;
			T d1,d2,dmin=1.0,Vx,Vy,dL1,dL2,VLx = 0,VLy = 0,aL = 0,bL = 0,R = 0,G1,G2,sBest,gBest,AB,progress;

			T ParLimit2 = 100.;
			T epsilon = 1.e+7 * std::numeric_limits<T>::epsilon();
			T factor1 = 32., factor2 = 32.;
			T ccc = 0.4;
			T factorUp = 10., factorDown = 0.1;

			Circle<T> Old, New;

			data.means();

			// New = circleIni;
			New.initFromCircle( circleIni );

			New.s = SigmaWithLargeCircleOption( data, New );
			GradientHessian( data, New, F1, F2, A11, A22, A12 );
			New.Gx = F1;
			New.Gy = F2;
			New.g = std::sqrt( F1 * F1 + F2 * F2 );

			lambda = LambdaIni;
			iter = inner = code = 0;
			sBest = gBest = progress = std::numeric_limits<T>::max();;

			do { // NextIteration
				bool break_outer = false;

				if ( iter > 0) {
					progress = ( std::abs( New.a - Old.a ) + std::abs( New.b - Old.b ) ) / 
						( std::pow( Old.a, 2 ) + std::pow( Old.b, 2 ) + 1.0 );
				}

				Old = New;
				if ( ++iter > IterMAX ) {
					break;
				}

				eigen2x2( A11, A22, A12, d1, d2, Vx, Vy );
				dmin = ( d1 < d2 ) ? d1 : d2;

				AB = std::sqrt( std::pow( Old.a, 2 ) + std::pow( Old.b, 2 ) ) + 1.0;

				if ( ( Old.g < factor1 * std::numeric_limits<T>::epsilon() ) && ( progress < epsilon ) ) {
					break;
				}

				if ( ( Old.s >= sBest ) && ( Old.g >= gBest ) ) {
					break;
				}

				if ( sBest > Old.s ) {
					sBest = Old.s;
				}
				if ( gBest > Old.g ) {
					gBest = Old.g;
				}

				G1 = Vx * F1 + Vy * F2;
				G2 = Vx * F2 - Vy * F1;

				do { // try_again

					if ( lambda < std::abs( G1 ) / AB / ccc - d1 ) {
						lambda = std::abs( G1 ) / AB / ccc - d1;
					}
					if ( lambda < std::abs( G2 ) / AB / ccc - d2 ) {
						lambda = std::abs( G2 ) / AB / ccc - d2;
					}

					dX = Old.Gx * ( Vx * Vx / ( d1 + lambda ) + Vy * Vy / ( d2 + lambda ) ) + Old.Gy * Vx * Vy * ( 1.0 / ( d1 + lambda ) - 1.0 / ( d2 + lambda ) );
					dY = Old.Gx * Vx * Vy * ( 1.0 / ( d1 + lambda ) - 1.0 / ( d2 + lambda ) ) + Old.Gy * ( Vx * Vx / ( d2 + lambda ) + Vy * Vy / ( d1 + lambda ) );

					New.a = Old.a - dX;
					New.b = Old.b - dY;

					if ( ( New.a == Old.a ) && ( New.b == Old.b ) ) {
						break;
					}

					if ( std::abs( New.a ) > ParLimit2 || std::abs( New.b ) > ParLimit2 ) {
						if ( check_line ) {
							for ( Mxx = Myy = Mxy = 0., i = 0; i < n; i++) {
								dx = data.hits[i][0] - data.meanX;
								dy = data.hits[i][1] - data.meanY;
								Mxx += dx * dx;
								Myy += dy * dy;
								Mxy += dx * dy;
							}

							eigen2x2( Mxx, Myy, Mxy, dL1, dL2, VLx, VLy );

							for ( Mxxy = 0., i = 0; i < n; i++ ) {
								dx = ( data.hits[i][0] - data.meanX ) * VLx + ( data.hits[i][1] - data.meanY ) * VLy;
								dy = ( data.hits[i][1] - data.meanY ) * VLx - ( data.hits[i][0] - data.meanX ) * VLy;
								Mxxy += dx * dx * dy;
							}

							R = ( Mxxy > 0. ) ? ParLimit2 : -ParLimit2;
							aL = -VLy * R;
							bL =  VLx * R;
							check_line = 0;
						}

						if ( ( New.a * VLy - New.b * VLx ) * R > 0. ) {
							New.a = aL;
							New.b = bL;
							New.s = SigmaWithLargeCircleOption( data, New );
							GradientHessian( data, New, F1, F2, A11, A22, A12 );
							New.Gx = F1;
							New.Gy = F2;
							New.g = std::sqrt( F1 * F1 + F2 * F2 );
							lambda = LambdaIni;
							sBest = gBest = std::numeric_limits<T>::max();;
							break;
						}
					}

					New.s = SigmaWithLargeCircleOption( data, New );
					GradientHessian( data, New, F1, F2, A11, A22, A12 );
					New.Gx = F1;
					New.Gy = F2;
					New.g = std::sqrt( F1 * F1 + F2 * F2 );

					if ( New.s < sBest * ( 1.0 + factor2 * std::numeric_limits<T>::epsilon() ) ) {
						lambda *= factorDown;
						break;
					} else {
						if ( ++inner > IterMAX ) {
							break_outer = true;
							break;
						}
						lambda = factorUp*lambda;
						continue;
					}

				} while( 1 ); // try_again loop

				if ( break_outer ) {
					break;    // break out of two loops
				}

			} while( 1 ); // NextIteration loop

			Old.r = OptimalRadius( data, Old );
			Old.i = iter;
			Old.j = inner;

			circleIni->initFromCircle( &Old );
			circleIni->chi2 = ChiSqr( data, circleIni );

			if ( iter  > IterMAX ) {
				code = 1;
			}
			if ( inner > IterMAX ) {
				code = 2;
			}
			if ( ( dmin <= 0. ) && ( code == 0 ) ) {
				code = 3;
			}

			return code;

		} /* CircleFitByChernovHoussam */

};


template<class T>
class LinearFit {
            
    public:
                
        LinearFit() {}
        
        static Line<T>* RegularFit( const std::vector<std::vector<T>>& hits ) {
            // From: Numerical Recipes in C++: The Art of Scientific Computing 3rd Ed.
            size_t n = hits.size();

            T ss, sx = 0., sy = 0., st2 = 0., t, sxoss, sigdat = 0.,
                  a = 0., b = 0., siga = 0., sigb = 0., chi2 = 0.;
            for ( size_t i = 0; i < n; i++ ) {
                sx += hits[i][0];
                sy += hits[i][1];
            }
            ss = n;
            sxoss = sx / ss;
            for ( size_t i = 0; i < n; i++ ) {
                t = hits[i][0] - sxoss;
                st2 += t * t;
                b += t * hits[i][1];
            }
            b /= st2;
            a = ( sy - sx * b ) / ss;
            siga = std::sqrt( ( 1.0 + sx * sx / ( ss * st2 ) ) / ss );
            sigb = std::sqrt( 1.0 / st2 );
            for ( size_t i = 0; i < n; i++ ) {
                chi2 += std::pow( hits[i][1] - a - b * hits[i][0], 2 );
            }
            if ( n > 2 ) {
                sigdat = std::sqrt( chi2 / ( n - 2 ) );
            }
            siga *= sigdat;
            sigb *= sigdat;

			chi2 = chi2 / ( n - 2 );

            Line<T>* line = new Line<T>( a, b, siga, sigb, chi2 );

            return line;
        }
        
        static Line<T>* RobustFit( const std::vector<std::vector<T>>& hits, Line<T>* line = 0 /* estimate */ ) {
            // From: Numerical Recipes in C++: The Art of Scientific Computing 3rd Ed.
            size_t n = hits.size();
            T abdev = 0;
            T b1, b2, f, f1, f2;
            T a, b, sigb, chi2;

            if ( line == 0 ) {
                line = RegularFit( hits );
            }
            
            a = line->a;
            b = line->b;
            sigb = line->sigb;
            chi2 = line->chi2;

            // robust fit
            b1 = b;
            f1 = rofunc( a, b1, abdev, hits );
            if ( sigb > 0.0 ) {
                b2 = b + SIGN( 3.0 * sigb, f1 );
                f2 = rofunc( a, b2, abdev, hits );
                if ( b2 == b1 ) {
                    abdev /= n;
                    return 0;
                }
                while ( f1 * f2 > 0.0 ) {
                    b = b2 + 1.6 * ( b2 - b1 );
                    b1 = b2;
                    f1 = f2;
                    b2 = b;
                    f2 = rofunc( a, b2, abdev, hits );
                }
                sigb = 0.01 * sigb;
                while ( std::abs( b2 - b1 ) > sigb ) {
                    b = b1 + 0.5 * ( b2 - b1 );
                    if ( b == b1 || b == b2 ) break;
                    f = rofunc( a, b, abdev, hits );
                    if ( f * f1 >= 0.0 ) {
                        f1 = f;
                        b1 = b;
                    } else {
                        f2 = f;
                        b2 = b;
                    }
                }
            }
            abdev /= n;
            line->b = b;
            line->sigb = sigb * 100.0; // restore original sigb or not?
            line->chi2 = chi2;

            return line;
        }
        
    private:

        static T rofunc( T& a, const T b, T& abdev, const std::vector<std::vector<T>>& hits ) {
            size_t n = hits.size();
            const T EPS = std::numeric_limits<T>::epsilon();
			size_t j;
            T d, sum = 0.0;
            std::vector<T> arr( n );
            for ( j = 0; j < n; j++ ) { arr[j] = hits[j][1] - b * hits[j][0]; }
            if ( ( n & 1 ) == 1) {
                a = select( ( n - 1 ) >> 1, arr );
            } else {
                j = n >> 1;
                a = 0.5 * ( select( j - 1, arr ) + select( j, arr ) );
            }
            abdev = 0.0;
            for ( j = 0; j < n; j++ ) {
                d = hits[j][1] - ( b * hits[j][0] + a );
                abdev += std::abs(d);
                if ( hits[j][1] != 0.0 ) { d /= std::abs( hits[j][1] ); }
                if ( std::abs(d) > EPS ) { sum += (d >= 0.0 ? hits[j][0] : -hits[j][0] ); }
            }
            return sum;
        }
        static T select( const int k, std::vector<T>& arr ) {
            int i, ir, j, l, mid, n = arr.size();
            T a;
            l = 0;
            ir = n-1;
            for (;;) {
                if (ir <= l + 1) {
                    if (ir == l + 1 && arr[ir] < arr[l]) {
                        SWAP( arr[l], arr[ir] );
                    }
                    return arr[k];
                } else {
                    mid = ( l + ir ) >> 1;
                    SWAP( arr[mid], arr[l+1] );
                    if ( arr[l] > arr[ir] ) {
                        SWAP( arr[l], arr[ir] );
                    }
                    if ( arr[l+1] > arr[ir] ) {
                        SWAP( arr[l+1], arr[ir] );
                    }
                    if ( arr[l] > arr[l+1] ) {
                        SWAP( arr[l], arr[l+1] );
                    }
                    i = l + 1;
                    j = ir;
                    a = arr[ l + 1 ];
                    for (;;) {
                        do { i++; } while ( arr[i] < a );
                        do { j--; } while ( arr[j] > a );
                        if ( j < i ) { break; }
                        SWAP( arr[i], arr[j] );
                    }
                    arr[ l + 1 ] = arr[j];
                    arr[ j ] = a;
                    if ( j >= k ) { ir = j - 1; }
                    if ( j <= k ) { l = i; }
                }
            }
        }
        
        static inline T SQR(const T a) { return a * a; }
        static inline T SIGN( const T& a, const T& b ) { return b >= 0 ? ( a >= 0 ? a : -a ) : ( a >= 0 ? -a : a ); }
        static inline float SIGN( const float& a, const double& b ) { return b >= 0 ? ( a >= 0 ? a : -a ) : ( a >= 0 ? -a : a ); }
        static inline float SIGN( const double& a, const float& b ) { return (float)( b >= 0 ? ( a >= 0 ? a : -a ) : ( a >= 0 ? -a : a) ); }
        static inline void SWAP( T& a, T& b ) { T dum = a; a = b; b = dum; }

};


template <typename T>
struct KDPointCloud {
	std::vector<std::vector<T> >  pts;
	inline size_t kdtree_get_point_count() const
	{
		return pts.size();
	}
	inline T kdtree_distance(const T *p1, const size_t idx_p2,size_t /*size*/) const
	{
		const T d0=p1[0]-pts[idx_p2][0];
		const T d1=p1[1]-pts[idx_p2][1];
		const T d2=p1[2]-pts[idx_p2][2];
		return d0*d0+d1*d1+d2*d2;
	}
	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim==0) return pts[idx][0];
		else if (dim==1) return pts[idx][1];
		else return pts[idx][2];
	}
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /*bb*/) const
	{
		return false;
	}
};

template <typename T>
struct KDTriplet {
	T dist1;
	T dist2;
	T angle;
	size_t n1;
	size_t n2;
};

template <class T>
class TrackCandidate
{

public:

	TrackCandidate( std::vector<std::vector<T>>& hits, T B ): mCircle(0), mLine(0), mB(B), mMinR(0), mMaxR(0) {
		// copy hits to local container
		std::copy( hits.begin(), hits.end(), back_inserter( mHits ) );
		// check if radius of the last hit is closer to 0,0
		// if    so - reverse hits
		size_t n = hits.size();
		if ( std::sqrt( mHits[0][0] * mHits[0][0] + mHits[0][1]*hits[0][1] ) > std::sqrt( mHits[ n - 1][0] * mHits[ n - 1][0] + mHits[ n - 1][1] * mHits[n - 1][1] ) ) {
			std::reverse( mHits.begin(), mHits.end() );
		}
		calcMinMaxR();
	}

	~TrackCandidate() {
		mHits.clear();
		if ( mCircle ) {
			delete mCircle;
			mCircle = 0;
		}
		if ( mLine ) {
			delete mLine;
			mLine = 0;
		}
	}

	bool isFitted() const {
		return mCircle && mLine && mCircle->is_good() && mLine->is_good();
	}

	void refit() {
		radiusFit();
		if ( mCircle->is_good() ) {
			szFit();
		}
	}

	size_t nhits() const {
		return mHits.size();
	}

	std::vector<T>& getFirstHit() {
		return mHits[0];
	}

	std::vector<T>& getLastHit() {
		return mHits.back();
	}

	std::vector<std::vector<T>>& getHits() {
		return mHits;
	}

	void deleteHits() {
		mHits.clear();
		std::vector<std::vector<T>>().swap(mHits);
		if ( mCircle ) {
			delete mCircle;
			mCircle = 0;
		}
		if ( mLine ) {
			delete mLine;
			mLine = 0;
		}
	}

	T sign() const {
		T s = getS( 3 );
		return ( s * mB ) < 0 ? -1. : 1;
	}

	T momentum() const {
		return std::sqrt( std::pow( Pt(), 2 ) + std::pow( Pl(), 2 ) );
	}

	T Pt() const {
		return 0.3 * std::fabs( mB ) * mCircle->r * 0.01;
	}

	T Pl() const {
		T Pl = Pt() * mLine->b;
		if ( ( mHits[1][2] - mHits[0][2] ) > 0 ) {
			if ( Pl < 0 ) {
				Pl *= -1;
			}
		} else {
			if ( Pl > 0 ) {
				Pl *= -1;
			}
		}
		return Pl;
	}

	// momentum for i-th hit
	std::vector<T> getMomForHit( size_t i ) {
		T pt = Pt();
		T pl = Pl();
		std::vector<T>& myHit = mHits[i];
		T phi = calcAlpha( i );
		T ptX = -pt * std::sin( phi );
		T ptY = +pt * std::cos( phi );
		if ( i > 0 ) {
			std::vector<T>& myHitBefore = mHits[ i - 1 ];
			T difX = myHit[0] - myHitBefore[0];
			T difY = myHit[1] - myHitBefore[1];
			if ( ( ptX * difX + ptY * difY ) < 0 ) {
				ptX *= -1;
				ptY *= -1;
			}
		} else if ( !mHits.empty() ) {
			std::vector<T>& myHitAfter = mHits[ 1 ];
			T difX = myHitAfter[0] - myHit[0];
			T difY = myHitAfter[1] - myHit[1];
			if ( ( ptX * difX + ptY * difY ) < 0 ) {
				ptX *= -1;
				ptY *= -1;
			}
		}
		std::vector<T> result(3);
		result[0] = ptX;
		result[1] = ptY;
		result[2] = pl;
		return result;
	}

	std::vector<T> getPosForHit( size_t i ) {
		TVector<T> v0( mCircle->r, 0, 0 );
		v0.rotateZ( calcAlpha(i) );
		v0.setZ( calcZPosByS( getS( i ) ) );
		std::vector<T> result(3);
		result[0] = v0.x() + mCircle->a;
		result[1] = v0.y() + mCircle->b;
		result[2] = v0.z();
		return result;
	}

	std::vector<T>& getHit( size_t i ) {
		return mHits[i];
	}

	T approxLength() {
		// FIXME: need proper solution here
		T length = 0;
		if ( isFitted() ) {
			std::vector<T> p = getMomForHit(0);
			std::vector<T> o = getPosForHit(0);
			std::vector<T> oL = getPosForHit( mHits.size() - 1 );
			Helix<T> helix( TVector<T>( p[0], p[1], p[2] ), TVector<T>( o[0], o[1], o[2] ), mB, sign() );
			length = helix.pathLength( TVector<T>( oL[0], oL[1], oL[2] ), false );
		} else {
			for ( size_t i = 1, ilen = mHits.size(); i < ilen; i++ ) {
				length += std::sqrt( std::pow( mHits[i][0] - mHits[i-1][0], 2 ) +
			                     std::pow( mHits[i][1] - mHits[i-1][1], 2 ) +
			                     std::pow( mHits[i][2] - mHits[i-1][2], 2 ) );
			}
		}
		return length;
	}

	T getB() const {
		return mB;
	}

	T minR() const {
		return mMinR;
	}

	T maxR() const {
		return mMaxR;
	}

	T radius() const {
		return mCircle->r;
	}

	T radius_err() const {
		return mCircle->s;
	}

	T tanl() const {
		return mLine->b;
	}

	T tanl_err() const {
		return mLine->sigb;
	}

	T dip() const {
		return std::cos( std::atan( mLine->b ) );
	}

	T distanceToCircle( T x, T y ) {
		return DistanceToCircle<T>( mCircle, x, y );
	}

	T distanceToCircle( const std::vector<T>& hit ) {
		return DistanceToCircle<T>( mCircle, hit[0], hit[1] );
	}

	void mergeCandidate( TrackCandidate<T>* candidate ) {
		// default sorting minR -> maxR for both tracks
		if ( minR() > candidate->maxR() ) {
			// add hits to front
			std::vector<std::vector<T>>& hits = candidate->getHits(); // FIXME: reserve extra capacity for hits
			std::reverse( hits.begin(), hits.end() );
			std::reverse( mHits.begin(), mHits.end() );
			std::copy( hits.begin(), hits.end(), back_inserter( mHits ) );
			std::reverse( mHits.begin(), mHits.end() );
		} else {
			// add hits to back
			std::vector<std::vector<T>>& hits = candidate->getHits(); // FIXME: reserve extra capacity for hits
			std::copy( hits.begin(), hits.end(), back_inserter( mHits ) );
		}
		candidate->deleteHits();
		refit();
		calcMinMaxR();
	}

	void print() {
		std::cout << "--- track params ---" << std::endl;
		if ( mCircle && mCircle->is_good() ) {
			std::cout << "   x0: " << mCircle->a << ", y0: " << mCircle->b << ", r: " << mCircle->r << ", sig: " << mCircle->s << std::endl;
		} else {
			std::cout << "   xy fit either not done or failed" << std::endl;
		}
		if ( mLine && mLine->is_good() ) {
			std::cout << "   z0: " << mLine->a << ", tanl: " << mLine->b << ", sig_z0: " << mLine->siga << ", sig_tanl: " << mLine->sigb << ", chi2: " << mLine->chi2 << std::endl;
		} else {
			std::cout << "   sz fit either not done or failed" << std::endl;
		}
		if ( mCircle->is_good() && mLine->is_good() ) {
			std::cout << "   Length: " << approxLength() << ", Pt: " << Pt() << ", Pl: " << Pl() << ", mom: " << momentum() << std::endl;
		}
	}

	void print_xy() {
		std::cout << "--- xy track params ---" << std::endl;
		std::cout << "    x0: " << mCircle->a << ", y0: " << mCircle->b << ", r: " << mCircle->r << ", sig: " << mCircle->s << std::endl;
		for( size_t i = 0; i < mHits.size(); i++ ) {
			std::cout << i << "-th hit, x: " << mHits[i][0] << ", y: " << mHits[i][1] << std::endl;
		}	
	}

	void print_sz() {
		std::cout << "--- sz track params ---" << std::endl;
		for( size_t i = 0; i < mHits.size(); i++ ) {
			std::cout << i << "-th hit, s: " << getS( i ) << ", z: " << mHits[i][2] << std::endl;
		}	
	}

	Circle<T>* getCircleFit() { return mCircle; }

private:

	void calcMinMaxR() {
		if ( mHits.empty() ) { return; }
		// min-max radius, used in track merging
		T radius = std::sqrt( mHits[0][0]*mHits[0][0] + mHits[0][1]*mHits[0][1] );
		mMinR = mMaxR = radius;
		for ( size_t i = 1, ilen = mHits.size(); i < ilen; i++ ) {
			radius = std::sqrt( mHits[i][0]*mHits[i][0] + mHits[i][1]*mHits[i][1] );
			if ( radius < mMinR ) {
				mMinR = radius;
			}
			if ( radius > mMaxR ) {
				mMaxR = radius;
			}
		}
	}

	void radiusFit() {
		if ( mCircle ) {
			delete mCircle;
		}
		// approximate fit first (fast!)
		mCircle = CircleFit<T>::RegularFit( mHits );

		// check for outliers, refit with outliers weighted down
		if ( mCircle->chi2 > 5.0 ) {
			CircleFit<T>::RobustFit( mHits, mCircle );
		}
	}

	void szFit() {
		// linear fit in sz
		if ( mLine ) {
			delete mLine;
		}
		if ( !mCircle || !mCircle->is_good() ) {
			return;    // can't fit sz without xy fit
		}
		std::vector<std::vector<T>> szhits;
		szhits.reserve( mHits.size() );
		for ( size_t i = 0, ilen = mHits.size(); i < ilen; i++ ) {
			std::vector<T> szhit(2);
			szhit[0] = getS( i );   // S-angle for i-th hit
			szhit[1] = mHits[i][2]; // z-coordinate for i-th hit
			szhits.emplace_back( szhit );
		}

		// approximate fit first (fast!)
		mLine = LinearFit<T>::RegularFit( szhits );

		// check chi2 for possible outliers, refit with outliers weighted down
		if ( mLine->chi2 > 5.0 ) {
			LinearFit<T>::RobustFit( szhits, mLine );
		}

	}

	T getS( size_t i ) const {
		T cx0 = mCircle->a, cy0 = mCircle->b;
		T dphi = phi_mpi_pi( std::atan2( mHits[0][1] - cy0, mHits[0][0] - cx0 ) - std::atan2( mHits[i][1] - cy0, mHits[i][0] - cx0 ) );
		return mCircle->r * dphi;
	}

	T calcAlpha( size_t i ) {
		return ( M_PI + std::atan2( -( mHits[i][1] - mCircle->b ), -( mHits[i][0] - mCircle->a ) ) );
	}

	T calcZPosByS( T s ) {
		return ( s * mLine->b + mLine->a );
	}

	T phi_mpi_pi( T x ) const {
		while ( x >= M_PI ) { x -= M_PI * 2.0; }
		while ( x < -M_PI ) { x += M_PI * 2.0; }
		return x;
	}

	/* members of the class */

	Circle<T>* mCircle; // xy fit results
	Line<T>* mLine;    // sz fit results
	std::vector<std::vector<T>> mHits;
	T mB;
	T mMinR;
	T mMaxR;
};

template <typename T>
bool pointsort ( std::vector<T> a, std::vector<T> b )
{
	return ( ( a[0]*a[0] + a[1]*a[1] ) > ( b[0]*b[0]+b[1]*b[1] ) );
}

template <typename T>
bool tripletsort ( KDTriplet<T> a, KDTriplet<T> b )
{
	return ( ( a.dist1 + a.dist2 ) < ( b.dist1 + b.dist2 ) );
}

template <typename T>
T angle_between_vectors( T x1, T y1, T z1, T x2, T y2, T z2 )
{
	T ad = std::sqrt( x1*x1 + y1*y1 + z1*z1 ), bd = std::sqrt( x2*x2 + y2*y2 + z2*z2 );
	T c0 = ( x1 / ad + x2 / bd ), c1 = ( y1 / ad + y2 / bd ), c2 = ( z1 / ad + z2 / bd );
	T d0 = ( x1 / ad - x2 / bd ), d1 = ( y1 / ad - y2 / bd ), d2 = ( z1 / ad - z2 / bd );
	return (T)( 2.0 * std::atan( std::sqrt( d0*d0 + d1*d1 + d2*d2 ) / std::sqrt( c0*c0 + c1*c1 + c2*c2 ) ) );
}

template <typename T> 
void make_triplets( size_t triplet_begin, size_t triplet_end, std::vector<std::vector<double>>& hits,
				nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<T, KDPointCloud<T> >, KDPointCloud<T>,3>& index,
				std::vector<std::vector<KDTriplet<T>>>& triplets,
				T search_radius = 10, T search_angle = M_PI / 8 ) {

	search_radius *= search_radius; // KD-tree is using search_radius^2
	nanoflann::SearchParams params;
	params.sorted = true;
	std::vector<std::pair<size_t,T> > ret_matches;
	double query_pt[3] = { 0., 0., 0.};
	size_t nMatches = 0;
	T dist1, dist2, angle;
	T current_search_radius, current_search_angle;

	for ( size_t h = triplet_begin; h < triplet_end; h++ ) {

		std::vector<T>& hit = hits[ h ];
		current_search_radius = search_radius;
		current_search_angle = search_angle;

		query_pt[0] = hits[h][0];
		query_pt[1] = hits[h][1];
		query_pt[2] = hits[h][2];
		nMatches = index.radiusSearch( &query_pt[0], current_search_radius, ret_matches, params );

		if ( nMatches < 5 ) {
			// if too few matches, increase search radius, retry search
			current_search_radius += 0.5 * search_radius;
			current_search_angle += M_PI / 32;
			nMatches = index.radiusSearch( &query_pt[0], current_search_radius, ret_matches, params );
		}

		if ( nMatches < 5 ) {
			// if STILL too few matches, increase search radius AGAIN, retry search
			current_search_radius += 0.5 * search_radius;
			current_search_angle += M_PI / 32;
			nMatches = index.radiusSearch( &query_pt[0], current_search_radius, ret_matches, params );
		}

		if ( nMatches < 3 ) {
			continue;    // well, we tried, but can't make any triplets
		}

		for ( size_t i = 1; i < nMatches; i++ ) { // skip original hit

			if ( !triplets[h].empty() && i > 6 ) { break; } // optimization: reduces time to process large events by 1/3

			std::vector<T>& hit1 = hits[ ret_matches[i].first ];
			dist1 = ret_matches[i].second; // hit distance to the point
			for ( size_t j = i + 1; j < nMatches; j++ ) {
				if ( !triplets[h].empty() && i > 6 ) { break; } // optimization: reduces time to process large events by 1/3

				std::vector<T>& hit2 = hits[ ret_matches[j].first ];
				dist2 = ret_matches[j].second;
				angle = M_PI - angle_between_vectors( hit1[0] - hit[0], hit1[1] - hit[1], hit1[2] - hit[2],
				                                      hit2[0] - hit[0], hit2[1] - hit[1], hit2[2] - hit[2] );
				if ( angle > current_search_angle ) {
					continue;    // discard curvy triplets
				}

				KDTriplet<T> t;
				t.dist1 = dist1;
				t.dist2 = dist2;
				t.angle = angle;
				t.n1 = ret_matches[i].first;
				t.n2 = ret_matches[j].first;
				triplets[h].push_back( t );
			}
		}

		// sort triplets by smallest distance
		if ( !triplets.empty() ) {
			std::sort( triplets[h].begin(), triplets[h].end(), tripletsort<T> );
		}
	}

}

template <typename T> 
std::vector<std::vector<std::vector<T> > > find_tracks( std::vector<std::vector<T>>& points, std::vector<std::vector<T>>& unused_hits,
        T search_radius = 10, T search_angle = M_PI / 8, size_t min_track_size = 10, 
		size_t nthreads = 1,
		bool stats = false )
{

	auto start0 = std::chrono::system_clock::now(); // total time

	const size_t TOTALHITS = points.size();

	// map of used hits
	std::vector<bool> used_hits( TOTALHITS, false );

	// data points storage container
	KDPointCloud<T> cloud;

	// generate / import points:
	auto start1 = std::chrono::system_clock::now(); // data import time
	cloud.pts.resize( TOTALHITS );
	for ( size_t i = 0, ilen = TOTALHITS; i < ilen; i++ ) {
		cloud.pts[i] = points[i];
	}
	auto end1 = std::chrono::system_clock::now(); // data import time

	// alias points as hits
	std::vector<std::vector<double>>& hits = cloud.pts;

	// declare kdtree index
	nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<T, KDPointCloud<T> >,
	          KDPointCloud<T>,3> index(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));

	// build kdtree index
	auto start3 = std::chrono::system_clock::now(); // kd-index build time
	index.buildIndex();
	auto end3 = std::chrono::system_clock::now(); // kd-index build time


	// make triplets
	auto start4 = std::chrono::system_clock::now(); // make triplets time

	std::vector<std::vector<KDTriplet<T>>> triplets( TOTALHITS );
	if ( nthreads > 1 ) {
		std::vector<std::thread> t(nthreads);
		size_t dnthreads = TOTALHITS / nthreads, triplet_begin = 0, triplet_end = 0;
		for ( size_t i = 0; i < nthreads; i++ ) {
			triplet_begin = i * dnthreads;
			triplet_end = ( i + 1 ) * dnthreads;
			if ( triplet_end > TOTALHITS ) { triplet_end = TOTALHITS; }
			if ( i == ( nthreads - 1 ) ) { triplet_end = TOTALHITS; }
			t[i] = std::thread( make_triplets<T>, triplet_begin, triplet_end, std::ref( hits ),
						std::ref( index ), std::ref( triplets ), search_radius, search_angle );
		}
    	for ( size_t i = 0; i < nthreads; i++ ) {
	        t[i].join();
    	}
	} else {
		make_triplets<T>( 0, TOTALHITS, hits, index, triplets, search_radius, search_angle );
	}

	auto end4 = std::chrono::system_clock::now(); // make triplets time

	// count triplets
	size_t triplet_ctr = 0;
	if ( stats ) {
		for ( size_t i = 0, ilen = triplets.size(); i < ilen; i++ ) {
			triplet_ctr += triplets[i].size();
		}
	}

	// reconstruct track candidates
	auto start5 = std::chrono::system_clock::now(); // construct track candidates

	std::vector<std::vector<size_t> > reco_tracks;
	size_t n0, n1;
	bool found;
	for ( size_t i = 0; i < TOTALHITS; i++ ) {

		if ( used_hits[i] != false ) {
			continue;
		}
		if ( triplets[i].size() == 0 ) {
			continue;
		}

		KDTriplet<T>& triplet = triplets[i][0];
		std::vector<size_t> track;
		track.emplace_back( i );
		used_hits[i] = true;

		if ( used_hits[ triplet.n1 ] == false ) {
			n0 = i;
			n1 = triplet.n1;
			do {
				found = false;
				track.emplace( track.begin(), n1 );
				used_hits[ n1 ] = true;
				if ( triplets[ n1 ].size() == 0 ) { // reached hit with no triplets
					break; // break out of do-while
				}
				// search for a good triplet
				for ( size_t t = 0; t < triplets[ n1 ].size(); t++ ) {
					if ( triplets[ n1 ][ t ].n1 == n0 && used_hits[ triplets[ n1 ][ t ].n2 ] == false ) {
						n0 = n1;
						n1 = triplets[ n1 ][ t ].n2;
						found = true;
						break; // found candidate via n2, break out of for-loop
					} else if ( triplets[ n1 ][ t ].n2 == n0 && used_hits[ triplets[ n1 ][ t ].n1 ] == false ) {
						n0 = n1;
						n1 = triplets[ n1 ][ t ].n1;
						found = true;
						break; // found candidate via n1, break out of for-loop
					}
				}
			} while( found == true );
		}

		if ( used_hits[ triplet.n2 ] == false ) {
			n0 = i;
			n1 = triplet.n2;
			do {
				found = false;
				track.emplace_back( n1 );
				used_hits[ n1 ] = true;
				if ( triplets[ n1 ].size() == 0 ) { // reached hit with no triplets
					break; // break out of do-while
				}
				// search for a good triplet
				for ( size_t t = 0; t < triplets[ n1 ].size(); t++ ) {
					if ( triplets[ n1 ][ t ].n1 == n0 && used_hits[ triplets[ n1 ][ t ].n2 ] == false ) {
						n0 = n1;
						n1 = triplets[ n1 ][ t ].n2;
						found = true;
						break; // found candidate via n2, break out of for-loop
					} else if ( triplets[ n1 ][ t ].n2 == n0 && used_hits[ triplets[ n1 ][ t ].n1 ] == false ) {
						n0 = n1;
						n1 = triplets[ n1 ][ t ].n1;
						found = true;
						break; // found candidate via n1, break out of for-loop
					}
				}
			} while( found == true );
		}

		if ( ( TOTALHITS < 1000 && track.size() >= 5 ) || track.size() >= min_track_size ) {
			// record track, keep hits marked as used
			reco_tracks.emplace_back( track );
		} else {
			// release unused hits
			for ( size_t i = 0, ilen = track.size(); i < ilen; i++ ) {
				used_hits[ track[i] ] = false;
			}
		}

	}
	auto end5 = std::chrono::system_clock::now();

	// convert hit ids to xyz

	auto start6 = std::chrono::system_clock::now();
	std::vector<std::vector<std::vector<T> > > final_tracks;
	final_tracks.reserve(reco_tracks.size());
	for( size_t i = 0, ilen = reco_tracks.size(); i < ilen; i++ ) {
		std::vector<std::vector<T> > final_track;
		for ( size_t j = 0, jlen = reco_tracks[i].size(); j < jlen; j++ ) {
			final_track.emplace_back( hits[ reco_tracks[i][j] ] );
		}
		final_tracks.emplace_back( final_track );
	}
	auto end6 = std::chrono::system_clock::now();


	size_t un_hits = 0, us_hits = 0;
	for ( size_t i = 0, ilen = used_hits.size(); i < ilen; i++ ) {
		if ( used_hits[i] ) {
			us_hits += 1;
		} else {
			un_hits += 1;
			unused_hits.emplace_back( hits[i] );
		}
	}

	auto end0 = std::chrono::system_clock::now(); // total time

	if ( stats ) {

		double us_pct = roundf( ( (double)us_hits / (double)used_hits.size() * 100.0 )  * 1000 ) / 1000;
		double un_pct = roundf( ( (double)un_hits / (double)used_hits.size() * 100.0 )  * 1000 ) / 1000;

		std::cout << "---------- KDfinder stats ----------\n";
		std::chrono::duration<double> elapsed_seconds1 = end1 - start1;
		std::cout << "     data import:    " << elapsed_seconds1.count() << "s\n";
		std::chrono::duration<double> elapsed_seconds3 = end3 - start3;
		std::cout << "     kd index build: " << elapsed_seconds3.count() << "s\n";
		std::chrono::duration<double> elapsed_seconds4 = end4 - start4;
		std::cout << "     triplet build:  " << elapsed_seconds4.count() << "s\n";
		std::chrono::duration<double> elapsed_seconds5 = end5 - start5;
		std::cout << "     track reco:     " << elapsed_seconds5.count() << "s\n";
		std::chrono::duration<double> elapsed_seconds6 = end6 - start6;
		std::cout << "     hit converter:  " << elapsed_seconds6.count() << "s\n";
		std::chrono::duration<double> elapsed_seconds0 = end0 - start0;
		std::cout << "     = total time:   " << elapsed_seconds0.count() << "s\n";
		std::cout << " --- objects --- \n";
		std::cout << "     triplets created: " << triplet_ctr << std::endl;
		std::cout << "     tracks   created: " << final_tracks.size() << std::endl;
		std::cout << " --- hits --- \n";
		std::cout << "     hits   used: " << us_hits << "\t\t or " << us_pct << "%" << "\n";
		std::cout << "     hits unused: " << un_hits << "\t\t or " << un_pct << "%" << "\n";
	}

	return final_tracks;
}

template <typename T>
std::vector<std::vector<std::vector<T> > > find_tracks_iterative( std::vector<std::vector<T>>& hits, std::vector<std::vector<T>>& unused_hits,
        T search_radius1 = 10, T search_angle1 = M_PI / 8, size_t min_track_size1 = 10,
        T search_radius2 = 12, T search_angle2 = M_PI / 8, size_t min_track_size2 = 6,
		size_t nthreads = 1,
        bool stats = false )
{
	// allows to make two iterations => gets better results in reasonable time

	std::vector<std::vector<std::vector<T> > > tracks = kdfinder::find_tracks<T>( hits, unused_hits,
	        search_radius1, search_angle1, min_track_size1, nthreads, stats );
	if ( unused_hits.size() > 10 ) {
		hits.clear();
		hits.insert( hits.begin(), std::make_move_iterator(unused_hits.begin()), std::make_move_iterator(unused_hits.end()) );
		unused_hits.erase( unused_hits.begin(), unused_hits.end() );
		std::vector<std::vector<std::vector<T> > > extra_tracks = kdfinder::find_tracks<T>( hits, unused_hits,
		        search_radius2, search_angle2, min_track_size2, nthreads, stats );
		if ( !extra_tracks.empty() ) {
			tracks.insert( tracks.end(), std::make_move_iterator(extra_tracks.begin()), std::make_move_iterator(extra_tracks.end()) );
			extra_tracks.erase(extra_tracks.begin(), extra_tracks.end());
		}
	}

	return tracks;
}

template <typename T>
std::vector<TrackCandidate<T>*> get_track_candidates( std::vector<std::vector<std::vector<T> > >& tracks, T B, bool stats = false )
{
	auto begin1 = std::chrono::system_clock::now();
	
	std::vector<TrackCandidate<T>*> candidates;
	for ( size_t i = 0, ilen = tracks.size(); i < ilen; i++ ) {
		TrackCandidate<T>* trk = new TrackCandidate<T>( tracks[i], B );
		trk->refit();
		if ( trk->isFitted() ) {
			candidates.emplace_back( trk );
		}
	}

	auto end1 = std::chrono::system_clock::now();

	if ( stats ) {
		std::cout << "---------- KDcandidates stats ----------\n";
		std::chrono::duration<double> elapsed_seconds1 = end1 - begin1;
		std::cout << "   conversion time : " << elapsed_seconds1.count() << "s\n";
		std::cout << "   input tracks    : " << tracks.size() << "\n";
		std::cout << "   candidates      : " << candidates.size() << "\n";
	}

	return candidates;
}

template <typename T>
bool candidatesort ( const TrackCandidate<T>* a, const TrackCandidate<T>* b )
{
	return ( a->nhits() > b->nhits() );
}

template <typename T>
bool candidatesortradius ( const TrackCandidate<T>* a, const TrackCandidate<T>* b )
{
	return ( a->minR() > b->minR() );
}

template <typename T>
bool ismergedcandidate( const TrackCandidate<T>* o )
{
	return ( o->nhits() == 0 );
}

template <typename T>
std::vector<TrackCandidate<T>*> merge_track_candidates(
		std::vector<TrackCandidate<T>*>& candidates,
        T c_tanl = 5.0 /* n sigma */, T c_xy = 5.0 /* n sigma */, T c_dist = 60.0 /* 3D distance between hits */,
		bool stats = false )
{
	auto start = std::chrono::system_clock::now();
	size_t tracks_before_merge = candidates.size();
	if ( candidates.size() < 2 ) {
		return candidates;    // need 2+ tracks to check for merging possibility
	}

	std::sort( candidates.begin(), candidates.end(), candidatesortradius<T> );
	
	for ( size_t i = 0, ilen = candidates.size(); i < ilen; i++ ) {

		if ( !candidates[i]->nhits() ) { continue; }
		if ( !candidates[i]->isFitted() ) { continue; }
		if ( ( candidates[i]->maxR() - candidates[i]->minR() ) < 10.0 ) { continue; } 

		for ( size_t j = i + 1, jlen = candidates.size(); j < jlen; j++ ) {
			if ( !candidates[i]->nhits() ) { break; }
			if ( !candidates[j]->nhits() ) { continue; }
			if ( !candidates[j]->isFitted() ) { continue; }
			if ( candidates[j]->maxR() > candidates[i]->minR() ) { continue; }
			if ( ( candidates[j]->maxR() - candidates[j]->minR() ) < 10.0 ) { continue; } 

			T distTanlErr = candidates[i]->nhits() > candidates[j]->nhits() ? candidates[i]->tanl_err() : candidates[j]->tanl_err();
        	if ( std::fabs( candidates[i]->tanl() - candidates[j]->tanl() ) > ( c_tanl * distTanlErr ) ) { continue; }

			T distToCircle = candidates[i]->nhits() > candidates[j]->nhits() ? 
				std::fabs( candidates[i]->distanceToCircle( candidates[j]->getLastHit() ) ) : std::fabs( candidates[j]->distanceToCircle( candidates[i]->getFirstHit() ) );
			T circleErr = candidates[i]->nhits() > candidates[j]->nhits() ? candidates[i]->radius_err() : candidates[j]->radius_err();
        	if ( std::fabs( distToCircle ) > ( c_xy * circleErr ) ) { continue; }

			std::vector<T>& hit1 = candidates[i]->getFirstHit();
			std::vector<T>& hit2 = candidates[j]->getLastHit();
			T dist = std::sqrt( std::pow(hit1[0] - hit2[0],2) + std::pow(hit1[1] - hit2[1],2) + std::pow(hit1[2] - hit2[2],2) );
			if ( dist > c_dist ) { continue; }

			candidates[i]->mergeCandidate( candidates[j] );
			// don't break out of the loop yet, there might be more candidates coming
		}
	}

	candidates.erase( std::remove_if( candidates.begin(), candidates.end(), ismergedcandidate<T> ), candidates.end() );

	auto stop = std::chrono::system_clock::now();
	if ( stats ) {
		size_t tracks_after_merge = candidates.size();
		std::cout << "---------- KDmerger stats ----------\n";
		std::chrono::duration<double> elapsed_seconds = stop - start;
		std::cout << "   tracks before: " << tracks_before_merge << "\n";
		std::cout << "   tracks merged: " << ( tracks_before_merge - tracks_after_merge ) << " = " << roundf( (double)( tracks_before_merge - tracks_after_merge ) / (double)tracks_before_merge * 100.0 * 100 ) / 100.0 << "%\n";
		std::cout << "   tracks after : " << tracks_after_merge << "\n";
		std::cout << "   time: " << elapsed_seconds.count() << "s\n";
	}

	return candidates;
}

template <typename T>
bool elementsort ( const std::tuple<T,Helix<T>*,size_t> a, const std::tuple<T,Helix<T>*,size_t> b )
{
	return ( std::get<0>( a ) < std::get<0>( b ) );
}

template <typename T>
bool vertexsort ( const std::pair< std::vector<T>, std::vector<size_t> > a, const std::pair< std::vector<T>, std::vector<size_t> > b )
{
	return ( ( std::get<1>( a ) ).size() > ( std::get<1>( b ) ).size() );
}


template <typename T>
std::vector< std::pair< std::vector<T>, std::vector<size_t> > >
	find_vertex_seeds( std::vector<TrackCandidate<T>*>& candidates, T x0 = 0, T y0 = 0,
	 			T c_z_dist = 0.5 /* cm */, T c_xy_dist = 2.0 /* cm */, T c_min_tracks = 3, bool stats = false ) {

	auto start = std::chrono::system_clock::now();

	std::vector< std::pair< std::vector<T>, std::vector<size_t> > > vertices;

	std::vector< std::tuple<T, Helix<T>*, size_t > > elements;
	std::vector< std::tuple<T, Helix<T>*, size_t > > sequence;
	Helix<T>* helix;
	T dca_z, dca_xy;

	size_t stat_candidates = candidates.size(),
			stat_elements = 0,
			stat_good_sequences = 0,
			stat_bad_sequences = 0;

	for ( size_t i = 0, ilen = candidates.size(); i < ilen; i++ ) {
		// make helix out of a candidate, get DCA of ( 0, 0 ), push to elements
        std::vector<T> o = candidates[i]->getPosForHit(0);
        std::vector<T> p = candidates[i]->getMomForHit(0);
		helix = new Helix<T>( TVector<T>( p[0], p[1], p[2] ), TVector<T>( o[0], o[1], o[2] ), candidates[i]->getB(), candidates[i]->sign() );
		T s = helix->pathLength( x0, y0 );
		TVector<T> pos = helix->at( s );
		dca_z = pos.z();
		dca_xy = pos.perp();
		if ( dca_xy > c_xy_dist ) { continue; }
		elements.push_back( std::make_tuple( dca_z, helix, i ) );
	}
	stat_elements = elements.size();

	// sort elements by dca_z
	std::sort( elements.begin(), elements.end(), elementsort<T> );

	// scan z-direction, find sequences
	for ( size_t i = 0, ilen = elements.size(); i < ilen; i++ ) {

		if ( sequence.empty() ) {
			// push helix and candidate to sequence
			sequence.push_back( elements[i] );

		} else {
			// get last element of sequence, check z distance to the previous element in sequence
			T dz = std::fabs( std::get<0>( elements[i] ) - std::get<0>( sequence.back() ) );

			if ( dz > c_z_dist ) {
				// stop sequence

				if ( sequence.size() >= c_min_tracks ) {
					// calculate average x, y, z of vertex, push mean(x,y,z) and candidates to vertices array
					std::vector<T> mean(6,0);
					std::vector<size_t> track_ids;
					track_ids.reserve( sequence.size() );

					// mean:
					for ( size_t k = 0, klen = sequence.size(); k < klen; k++ ) {
						Helix<T>* helix = std::get<1>( sequence[k] );
						T s = helix->pathLength( x0, y0 );
						TVector<T> pos = helix->at( s );
						mean[0] += pos.x();
						mean[1] += pos.y();
						mean[2] += pos.z();
						track_ids.push_back( std::get<2>( sequence[k] ) );
					}
					mean[0] /= sequence.size();
					mean[1] /= sequence.size();
					mean[2] /= sequence.size();

					// deviation:
					for ( size_t k = 0, klen = sequence.size(); k < klen; k++ ) {
						Helix<T>* helix = std::get<1>( sequence[k] );
						T s = helix->pathLength( 0, 0 );
						TVector<T> pos = helix->at( s );
						mean[3] += std::fabs( mean[0] - pos.x() );
						mean[4] += std::fabs( mean[1] - pos.y() );
						mean[5] += std::fabs( mean[2] - pos.z() );
					}
					mean[3] /= sequence.size();
					mean[4] /= sequence.size();
					mean[5] /= sequence.size();

					vertices.push_back( std::make_pair( mean, track_ids ) );
					stat_good_sequences += 1;
				} else {
					stat_bad_sequences += 1;
				}

				// cleanup => delete all sequence ptrs, clear sequence
				for ( auto it = sequence.begin() ; it != sequence.end(); ++it ) { delete std::get<1>( *it ); }
				sequence.clear();

			} else {
				// start sequence => add element to sequence
				sequence.push_back( elements[i] );

			}
		}
	}

	// sort vertices by number of tracks per vertex
	std::sort( vertices.begin(), vertices.end(), vertexsort<T> );

	auto stop = std::chrono::system_clock::now();

	if ( stats ) {
		std::cout << "---------- KDvertexer stats ----------\n";
		std::cout << "   input candidates: " << stat_candidates << "\n";
		std::cout << "   proj. candidates: " << stat_elements << "\n";
		std::cout << "   good pre-vertices: " << stat_good_sequences << "\n";
		std::cout << "   weak pre-vertices: " << stat_bad_sequences << "\n";
		std::chrono::duration<double> elapsed_seconds = stop - start;
		std::cout << "   time: " << elapsed_seconds.count() << "s\n";
	}

	return vertices;
}

template <typename T>
size_t get_track_color( T pt )
{
	T r = 0, g = 0, b = 0, p = std::fabs( pt ), maxp = 4.5, colval = std::min( 1.0, p / maxp ), colvaltimes4 = colval * 4.0;
	if ( colval < 0.25 ) {
		b = g = colvaltimes4;
		b += 1.0 - colvaltimes4;
	} else if ( colval < 0.5 ) {
		b = g = 1.0 - ( colvaltimes4 - 1.0 );
		g += colvaltimes4 - 1.0;
	} else if ( colval < 0.75 ) {
		g = r = colvaltimes4 - 2.0;
		g += 1.0 - ( colvaltimes4 - 2.0 );
	} else {
		g = r = 1.0 - ( colvaltimes4 - 3.0 );
		r += colvaltimes4 - 3.0;
	}
	if ( rand() % 1000 < 10 ) {
		r = 1.0;
	}
	return ( (int( r * 255) & 0xff) << 16 ) + ( ( int( g * 255 ) & 0xff) << 8 ) + (int( b * 255 ) & 0xff );
}


template <typename T>
std::string export_candidates_json( std::vector<TrackCandidate<T>*>& candidates, std::vector<std::vector<T> >& hits,
                                    size_t run = 0, size_t event_number = 0, size_t evt_time = 0, double B = -0.5 )
{

	std::stringstream ofs;

	ofs << "{ \n"
	    << " \"EVENT\": {  \n"
	    << "  \"runid\": " << run << ", \n"
	    << "  \"evtid\": "<< event_number << ", \n"
	    << "  \"time\": "<< evt_time <<", \n"
	    << "  \"B\": "<< B <<" \n"
	    << " }, \n";

	ofs << "\"META\": { \n"
	    << " \"HITS\": { \n"
	    << "    \"TPC\": { \n"
	    << "       \"type\": \"3D\", \n"
	    << "       \"options\": { \n"
	    << "           \"size\": 2, \n"
	    << "           \"color\": 16777215 \n"
	    << "       }\n"
	    << "   }\n"
	    << " },\n"
	    << " \"TRACKS\": { \n"
	    << "    \"TPC\": { \n"
	    << "      \"r_min\": 0,\n"
	    << "      \"r_max\": 2000,\n"
	    << "      \"size\": 2, \n"
	    << "      \"thickness\": 2 \n"
	    << "    }\n"
	    << "  }\n"
	    << "},\n"
	    << " \"TRACKS\": {\n"
	    << "   \"TPC\": [\n";

	for ( int i = 0, ilen = candidates.size(); i < ilen; i++ ) {
		std::vector<T> hit = candidates[i]->getPosForHit(0);
		std::vector<T> mom = candidates[i]->getMomForHit(0);
		size_t nhits = candidates[i]->nhits();
		T pt = candidates[i]->Pt();
		T sign = candidates[i]->sign();
		T length = candidates[i]->approxLength();
		size_t color = get_track_color<T>( pt );
		ofs << "{ \"color\": "<< color <<", \"pt\": "<< pt <<", \"xyz\":[" << hit[0] << "," << hit[1] << "," << hit[2] << "], \"pxyz\":[" << mom[0] << "," << mom[1] << "," << mom[2] << "],\"l\":" << length << ",\"nh\":" << nhits << ",\"q\":" << sign << "}";
		if ( i != ( ilen - 1 ) ) {
			ofs << ",";
		} else {
			ofs << "\n";
		}
	}

	ofs << "    ]\n"
	    << "  },\n"
	    << " \"HITS\": {\n"
	    << "   \"TPC\": [\n";

	for ( int i = 0, ilen = hits.size(); i < ilen; i++ ) {
		ofs << "[ " << hits[i][0] << "," << hits[i][1] << "," << hits[i][2] << " ]";
		if ( i != ( ilen - 1 ) ) {
			ofs << ",";
		} else {
			ofs << "\n";
		}
	}

	ofs << "   ]\n"
	    << " }\n"
	    << "}\n";

	return ofs.str();
}


template <typename T>
std::string export_json( std::vector<std::vector<std::vector<T> > >& tracks, std::vector<std::vector<T> >& hits,
                         size_t run = 0, size_t event_number = 0, size_t evt_time = 0, double B = -0.5 )
{

	std::stringstream ofs;

	ofs << "{ \n"
	    << " \"EVENT\": {  \n"
	    << "  \"runid\": " << run << ", \n"
	    << "  \"evtid\": "<< event_number << ", \n"
	    << "  \"time\": "<< evt_time <<", \n"
	    << "  \"B\": "<< B <<" \n"
	    << " }, \n";

	ofs << "\"META\": { \n"
	    << " \"HITS\": { \n"
	    << "    \"TPC\": { \n"
	    << "       \"type\": \"3D\", \n"
	    << "       \"options\": { \n"
	    << "           \"size\": 2, \n"
	    << "           \"color\": 16777215 \n"
	    << "       }\n"
	    << "   }\n"
	    << " },\n"
	    << " \"TRACKS\": { \n"
	    << "    \"TPC\": { \n"
	    << "      \"size\": 2, \n"
	    << "      \"thickness\": 2, \n"
	    << "      \"color\": 255 \n"
	    << "    }\n"
	    << "  }\n"
	    << "},\n"
	    << " \"TRACKS\": {\n"
	    << "   \"TPC\": [\n";

	for ( int i = 0, ilen = tracks.size(); i < ilen; i++ ) {
		ofs << "{ \"nh\": "<< tracks[i].size() <<", \"pts\": [ ";
		for ( int j = 0, jlen = tracks[i].size(); j < jlen; j++ ) {
			ofs << "[ " << tracks[i][j][0] << "," << tracks[i][j][1] << "," << tracks[i][j][2] << " ]";
			if ( j != ( jlen - 1 ) ) {
				ofs << ",";
			}
		}
		ofs << " ] }\n";
		if ( i != ( ilen - 1 ) ) {
			ofs << ",";
		}
	}

	ofs << "    ]\n"
	    << "  },\n"
	    << " \"HITS\": {\n"
	    << "   \"TPC\": [\n";

	for ( int i = 0, ilen = hits.size(); i < ilen; i++ ) {
		ofs << "[ " << hits[i][0] << "," << hits[i][1] << "," << hits[i][2] << " ]";
		if ( i != ( ilen - 1 ) ) {
			ofs << ",";
		} else {
			ofs << "\n";
		}
	}

	ofs << "   ]\n"
	    << " }\n"
	    << "}\n";

	return ofs.str();
}

template <typename T>
std::string export_candidates_json_old( std::vector<TrackCandidate<T>*>& candidates, std::vector<ulong>& trigger_ids,
                                        size_t run = 0, size_t event_number = 0, size_t evt_time = 0, T B = -0.5, bool round = false )
{

	std::stringstream ofs;
	ofs << "{"
	    << "\"R\": " << run << ", "
	    << "\"Evt\": " << event_number << ", "
	    << "\"B\": " << B << ","
	    << "\"tm\": " << evt_time << ", "
	    << "\"trig\": [";

	for ( size_t i = 0; i < trigger_ids.size(); i++ ) {
		if ( i != 0 ) {
			ofs << ",";
		}
		ofs << trigger_ids[i];
	}

	ofs << "],\n"
	    << "\"tracks\": [\n";

	for ( int i = 0, ilen = candidates.size(); i < ilen; i++ ) {
		std::vector<T> hit = candidates[i]->getPosForHit(0);
		std::vector<T> mom = candidates[i]->getMomForHit(0);
		size_t nhits = candidates[i]->nhits();
		T pt = candidates[i]->Pt();
		T sign = candidates[i]->sign();
		T length = candidates[i]->approxLength();
		size_t color = get_track_color<T>( pt );

		if ( round ) {
			hit[0] = std::floor( hit[0] * 10 ) / 10;
			hit[1] = std::floor( hit[0] * 10 ) / 10;
			hit[2] = std::floor( hit[0] * 10 ) / 10;
			mom[0] = std::floor( mom[0] * 1000) / 1000;
			mom[1] = std::floor( mom[1] * 1000) / 1000;
			mom[2] = std::floor( mom[2] * 1000) / 1000;
			length = std::floor( length * 100 ) / 100;
		}

		ofs << "{ \"color\": "<< color <<", \"pt\": "<< pt <<", \"xyz\":[" << hit[0] << "," << hit[1] << "," << hit[2] << "], \"pxyz\":[" << mom[0] << "," << mom[1] << "," << mom[2] << "],\"l\":" << length << ",\"nh\":" << nhits << ",\"q\":" << sign << "}";
		if ( i != ( ilen - 1 ) ) {
			ofs << ",";
		} else {
			ofs << "\n";
		}
	}

	ofs << "] }";
	return ofs.str();
}

} /* namespace kdfinder */

#endif /* KDFINDER_H_ */

