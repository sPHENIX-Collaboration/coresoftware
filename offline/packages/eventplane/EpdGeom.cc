/*************************************
 *
 * \author Mike Lisa
 * \date 4 Jan 2018
 *
 * \description:
 *  class defining geometrical aspects of an EPD Tile
 *  (position of center, RandomPointOnTile(), etc.)
 *
 * The user may pass the PP/TT/SN _or_ the uniqueID to
 *   most functions.  No option to pass StEpdHit object,
 *   because we want to avoid StObject-dependence.
 *
 *
 * Adapted for use in sPHENIX
 *   by Brennan Schaefer
 *   11 July 2022
 *
 * The common convension for coordinate systems is used here,
 *   where 'x' points out from the center of the collider ring.
 *   A simple mapping from STAR to sPHENIX is a rotation by 90
 *   degrees clockwise viewed from the sky looking down at RHIC
 *
 *   The angle PHI is zero in the x-direction
 *
 * The STAR-EPD-East plane is mapped to sPHENIX-EPD-South
 * The STAR-EPD-West plane is mapped to sPHENIX-EPD-North
 *
 * Abbreviations
 *   PP = position;
 *   TT = tilenumber;
 *   SN = southnorth;
 *   SS = 30 deg, super sector of two sectors
 *
 *************************************/

#include "EpdGeom.h"
#include <TMath.h>
#include <phool/PHRandomSeed.h>

EpdGeom::EpdGeom()
  : mPP(0)
  , mTT(0)
  , mSN(0)
  , mPhiCenter()
  , mRmin()
  , mRmax()
  , mRave()
{
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  m_Seed = PHRandomSeed();  // fixed seed is handled in this function
  gsl_rng_set(m_RandomGenerator, m_Seed);

  InitializeGeometry();
}

EpdGeom::~EpdGeom()
{
  gsl_rng_free(m_RandomGenerator);
}

void EpdGeom::InitializeGeometry()
{
  //  First, save the phi values of the centers of all tiles.
  //  I am aware that this is a it wasteful, since all the even (or odd) numbered
  //  tiles of a given position all have the same phi.  But you'll go nuts otherwise.
  //  Better to waste a little memory than to get confused.
  //
  double DeltaPhiSS = 30.0 * M_PI / 180.0;  // 30 degree supersectors
  short SN = 0;                               // South
  for (int PP = 1; PP < 13; PP++)
  {
    double phiSS = M_PI / 2.0 - (PP - 0.5) * DeltaPhiSS;
    if (phiSS < 0.0)
    {
      phiSS += 2.0 * M_PI;
    }
    mPhiCenter[PP - 1][0][SN] = phiSS;
    for (int TT = 2; TT < 32; TT += 2)
    {  // EVENS
      mPhiCenter[PP - 1][TT - 1][SN] = phiSS - DeltaPhiSS / 4.0;
    }
    for (int TT = 3; TT < 32; TT += 2)
    {  // ODDS
      mPhiCenter[PP - 1][TT - 1][SN] = phiSS + DeltaPhiSS / 4.0;
    }
  }
  SN = 1;  // North
  for (int PP = 1; PP < 13; PP++)
  {
    double phiSS = M_PI / 2.0 + (PP - 0.5) * DeltaPhiSS;
    if (phiSS > 2.0 * M_PI)
    {
      phiSS -= 2.0 * M_PI;
    }
    mPhiCenter[PP - 1][0][SN] = phiSS;
    for (int TT = 2; TT < 32; TT += 2)
    {  // EVENS
      mPhiCenter[PP - 1][TT - 1][SN] = phiSS + DeltaPhiSS / 4.0;
    }
    for (int TT = 3; TT < 32; TT += 2)
    {  // ODDS
      mPhiCenter[PP - 1][TT - 1][SN] = phiSS - DeltaPhiSS / 4.0;
    }
  }

  //  Now the inner, outer, and average radius of a _ROW_
  const double RowHeight[16] = {4.4, 4.4, 4.4, 5.53, 5.53, 5.53, 5.53, 5.53, 5.53, 5.53, 5.53, 5.53, 5.53, 5.53, 5.53, 5.53};
  double Rminimum = 4.6;  // the distance from beamline to the inner edge of tile 1
  mRmin[0] = Rminimum;    // row 1 (tiles 1)

  for (int irow = 1; irow < 16; irow++)
  {
    mRmin[irow] = mRmin[irow - 1] + RowHeight[irow - 1];
    mRmax[irow - 1] = mRmin[irow];
  }
  mRmax[15] = mRmin[15] + RowHeight[15];
  for (int irow = 0; irow < 16; irow++)
  {
    mRave[irow] = 0.5 * (mRmin[irow] + mRmax[irow]);
  }
}

void EpdGeom::SetPpTtSn(short uniqueID)
{
  mPP = abs(uniqueID / 100);
  mTT = abs(uniqueID % 100);
  mSN = (uniqueID > 0) ? +1 : -1;
}

double EpdGeom::GetZwheel()
{
  const double z_EPD = 316.0;  // Distance (cm) of EPD from center of TPC, in the z-direction
  return z_EPD * mSN;
}

//------------------------------------------------
short EpdGeom::Row(short uniqueID)
{
  SetPpTtSn(uniqueID);
  return this->Row();
}

short EpdGeom::Row(short PP, short TT, short SN)
{
  mPP = PP;
  mTT = TT;
  mSN = SN;
  return this->Row();
}

short EpdGeom::Row()
{
  return mTT / 2 + 1;
}

//------------------------------------------------
TVector3 EpdGeom::TileCenter(short uniqueID)
{
  SetPpTtSn(uniqueID);
  return this->TileCenter();
}

TVector3 EpdGeom::TileCenter(short PP, short TT, short SN)
{
  mPP = PP;
  mTT = TT;
  mSN = SN;
  return this->TileCenter();
}

TVector3 EpdGeom::TileCenter()
{
  //  double Rmin,Rmax;
  //  GetRminRmax(&Rmin,&Rmax);
  double ZZ = this->GetZwheel();
  TVector3 cent(mRave[this->Row() - 1], 0.0, ZZ);
  //  cent.SetXYZ(0.5*(Rmin+Rmax),0.0,ZZ);
  //  cent.RotateZ(this->GetPhiCenter());
  int sn = (mSN > 0) ? 1 : 0;
  cent.RotateZ(mPhiCenter[mPP - 1][mTT - 1][sn]);
  return cent;
}

//-----------------------------------------------------
TVector3 EpdGeom::RandomPointOnTile(short uniqueID)
{
  SetPpTtSn(uniqueID);
  return this->RandomPointOnTile();
}

TVector3 EpdGeom::RandomPointOnTile(short PP, short TT, short SN)
{
  mPP = PP;
  mTT = TT;
  mSN = SN;
  return this->RandomPointOnTile();
}

TVector3 EpdGeom::RandomPointOnTile()
{
  double GapWidth = 0.08;  // one half of the glue gap width
  double Aparam = 2.0 * tan(7.5 * M_PI / 180.0);
  double Bparam = -2.0 * GapWidth;

  double ZZ = this->GetZwheel();
  short  RR = this->Row();
  double Rmin = mRmin[RR - 1];
  double Rmax = mRmax[RR - 1];

  double Xmin = Rmin + GapWidth;
  double Xmax = Rmax - GapWidth;

  if (1 == RR)
  {
    Xmin -= 2.0 * GapWidth;
  }  // no glue on the "inside" of tile 1
  if (16 == RR)
  {
    Xmax += GapWidth;
  }  // no glue on "outside" of TT30,31

  // the reason for this next command is that Tile 01 is a pain in the neck.
  // I didn't figure out an easy way to get a random point inside the pentagon,
  // so I make the outer radius a little too big.  Then I get a point, and if
  // it doesn't fit in the tile, I try again until it's finally there.

  if (1 == RR)
  {
    Xmax += GapWidth;
  }

  double A = Aparam;
  if (1 == RR)
  {
    A = A * 2.0;
  }

  double gamma = 0.5 * A * pow(Xmin, 2) + Bparam * Xmin;
  double alpha = 0.5 * A * pow(Xmax, 2) + Bparam * Xmax - gamma;

  double q = gsl_rng_uniform_pos(RandomGenerator());
  double XX = (sqrt(pow(Bparam, 2) + 2.0 * A * (alpha * q + gamma)) - Bparam) / A;
  q = gsl_rng_uniform_pos(RandomGenerator());
  double DeltaY = A * XX + Bparam;
  double YY = (q - 0.5) * DeltaY;

  TVector3 Point(XX, YY, ZZ);
  //  Point.RotateZ(this->GetPhiCenter());
  int sn = (mSN > 0) ? 1 : 0;
  Point.RotateZ(mPhiCenter[mPP - 1][mTT - 1][sn]);

  // if this is Tile 01, there's the possibility that the point does
  // not fit into the tile after all, so check and if it doesn't
  // then try again.

  if (1 == RR)
  {
    if (!(this->IsInTile(Point.X(), Point.Y())))
    {
      return this->RandomPointOnTile();
    }
  }

  return Point;
}

//----------------------------------------------------
void EpdGeom::GetCorners(short uniqueID, int* nCorners, double* x, double* y)
{
  SetPpTtSn(uniqueID);
  GetCorners(nCorners, x, y);
}

void EpdGeom::GetCorners(short position, short tilenumber, short southnorth, int* nCorners, double* x, double* y)
{
  mPP = position;
  mTT = tilenumber;
  mSN = southnorth;
  GetCorners(nCorners, x, y);
}

void EpdGeom::GetCorners(int* nCorners, double* xc, double* yc)
{
  double x[5];
  double y[5];
  // we provide the user five corners.  For tiles 2-31, the fifth "corner" is junk.
  // only tile 1 is a pentagon
  double OpeningAngle = 7.5 * M_PI / 180.0;
  double GapWidth = 0.08;  // gap between tiles / 2
  short RR = this->Row();
  double Rmin = mRmin[RR - 1];
  double Rmax = mRmax[RR - 1];
  //  GetRminRmax(&Rmin,&Rmax);
  if (1 == RR)
  {
    *nCorners = 5;
    double xtmp[3], ytmp[3];
    xtmp[0] = Rmin;
    ytmp[0] = +Rmin * tan(OpeningAngle);
    xtmp[1] = Rmax;
    ytmp[1] = +Rmax * tan(OpeningAngle);
    xtmp[2] = Rmax;
    ytmp[2] = -Rmax * tan(OpeningAngle);
    for (int ic = 0; ic < 3; ic++)
    {
      x[ic] = xtmp[ic] * cos(OpeningAngle) - ytmp[ic] * sin(OpeningAngle);
      y[ic] = +xtmp[ic] * sin(OpeningAngle) + ytmp[ic] * cos(OpeningAngle);
    }
    y[0] -= GapWidth;
    y[1] -= GapWidth;
    x[1] -= GapWidth;
    x[2] -= GapWidth;
    x[3] = x[1];
    y[3] = -y[1];
    x[4] = x[0];
    y[4] = -y[0];
  }
  else
  {
    *nCorners = 4;
    x[0] = Rmin + GapWidth;
    y[0] = +Rmin * tan(OpeningAngle) - GapWidth;
    x[1] = Rmax - GapWidth;
    y[1] = +Rmax * tan(OpeningAngle) - GapWidth;
    x[2] = Rmax - GapWidth;
    y[2] = -Rmax * tan(OpeningAngle) + GapWidth;
    x[3] = Rmin + GapWidth;
    y[3] = -Rmin * tan(OpeningAngle) + GapWidth;
    x[4] = -999;
    y[4] = -999;  // unused for TT!=1

    if (16 == RR)
    {  // there is no glue "outside" TT30,31
      x[1] += GapWidth;
      x[2] += GapWidth;
    }
  }
  //  double phi = this->GetPhiCenter();
  int sn = (mSN > 0) ? 1 : 0;
  double phi = mPhiCenter[mPP - 1][mTT - 1][sn];
  for (int icorn = 0; icorn < (*nCorners); icorn++)
  {
    xc[icorn] = +x[icorn] * cos(phi) - y[icorn] * sin(phi);
    yc[icorn] = +x[icorn] * sin(phi) + y[icorn] * cos(phi);
  }
}

//---------------------------------------------------------------------
bool EpdGeom::IsInTile(short uniqueID, double x, double y)
{
  SetPpTtSn(uniqueID);
  return this->IsInTile(x, y);
}

bool EpdGeom::IsInTile(short position, short tilenumber, short southnorth, double x, double y)
{
  mPP = position;
  mTT = tilenumber;
  mSN = southnorth;
  return this->IsInTile(x, y);
}

bool EpdGeom::IsInTile(double x, double y)
{
  double PolygonX[6];
  double PolygonY[6];
  int numberOfCorners;
  this->GetCorners(&numberOfCorners, PolygonX, PolygonY);
  PolygonX[numberOfCorners] = PolygonX[0];
  PolygonY[numberOfCorners] = PolygonY[0];
  return TMath::IsInside(x, y, numberOfCorners + 1, PolygonX, PolygonY);
}
