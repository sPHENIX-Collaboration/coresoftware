#include "PHG4TpcCentralMembrane.h"

#include <phparameter/PHParameterInterface.h>  // for PHParameterInterface

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for get_volume_id
#include <g4main/PHG4Hitv1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TVector3.h>
#include <TString.h>

#include <iostream>  // for operator<<, endl, basi...

class PHCompositeNode;

namespace
{
  // unique detector id for all direct lasers
  static const int detId = PHG4HitDefs::get_volume_id("PHG4TpcCentralMembrane");
  
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
  
  template<class T> inline T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

}  // namespace

//_____________________________________________________________
// all distances in mm, all angles in rad
// class that generates stripes and dummy hit coordinates
// stripes have width of one mm, length of one pad width, and are centered in middle of sector gaps
PHG4TpcCentralMembrane::PHG4TpcCentralMembrane(const std::string& name)
  : SubsysReco(name)
  , PHParameterInterface(name)
{
  InitializeParameters();

  // set to 1.0 mm for all else
  for (int j = 0; j < nRadii; j++)
  {
    for (int i = 0; i < nStripes_R1; i++)
    {
      str_width_R1_e[i][j] = 1.0 * mm;
      str_width_R1[i][j] = 1.0 * mm;
    }
    for (int i = 0; i < nStripes_R2; i++)
    {
      str_width_R2[i][j] = 1.0 * mm;
    }
    for (int i = 0; i < nStripes_R3; i++)
    {
      str_width_R3[i][j] = 1.0 * mm;
    }
  }

  CalculateVertices(nStripes_R1, nPads_R1, R1_e, spacing_R1_e, x1a_R1_e, y1a_R1_e, x1b_R1_e, y1b_R1_e, x2a_R1_e, y2a_R1_e, x2b_R1_e, y2b_R1_e, x3a_R1_e, y3a_R1_e, x3b_R1_e, y3b_R1_e, padfrac_R1, str_width_R1_e, widthmod_R1_e, nGoodStripes_R1_e, keepUntil_R1_e, nStripesIn_R1_e, nStripesBefore_R1_e);
  CalculateVertices(nStripes_R1, nPads_R1, R1, spacing_R1, x1a_R1, y1a_R1, x1b_R1, y1b_R1, x2a_R1, y2a_R1, x2b_R1, y2b_R1, x3a_R1, y3a_R1, x3b_R1, y3b_R1, padfrac_R1, str_width_R1, widthmod_R1, nGoodStripes_R1, keepUntil_R1, nStripesIn_R1, nStripesBefore_R1);
  CalculateVertices(nStripes_R2, nPads_R2, R2, spacing_R2, x1a_R2, y1a_R2, x1b_R2, y1b_R2, x2a_R2, y2a_R2, x2b_R2, y2b_R2, x3a_R2, y3a_R2, x3b_R2, y3b_R2, padfrac_R2, str_width_R2, widthmod_R2, nGoodStripes_R2, keepUntil_R2, nStripesIn_R2, nStripesBefore_R2);
  CalculateVertices(nStripes_R3, nPads_R3, R3, spacing_R3, x1a_R3, y1a_R3, x1b_R3, y1b_R3, x2a_R3, y2a_R3, x2b_R3, y2b_R3, x3a_R3, y3a_R3, x3b_R3, y3b_R3, padfrac_R3, str_width_R3, widthmod_R3, nGoodStripes_R3, keepUntil_R3, nStripesIn_R3, nStripesBefore_R3);
}

//______________________________________________________
PHG4TpcCentralMembrane::~PHG4TpcCentralMembrane()
{
  for (auto&& hit : PHG4Hits)
  {
    delete hit;
  }
}

//_____________________________________________________________
int PHG4TpcCentralMembrane::InitRun(PHCompositeNode* /* topNode */)
{
  // setup parameters
  UpdateParametersWithMacro();
  electrons_per_stripe = get_int_param("electrons_per_stripe");
  electrons_per_gev = get_double_param("electrons_per_gev");

  std::cout << "PHG4TpcCentralMembrane::InitRun - electrons_per_stripe: " << electrons_per_stripe << std::endl;
  std::cout << "PHG4TpcCentralMembrane::InitRun - electrons_per_gev " << electrons_per_gev << std::endl;
  
  // reset g4hits
  for (auto&& hit : PHG4Hits)
  {
    delete hit;
  }

  PHG4Hits.clear();

  /*
   * utility function to
   * - duplicate generated G4Hit to cover both sides of the central membrane
   * - adjust hit time and z,
   * - insert in container
   */
  auto adjust_hits = [&](PHG4Hit* source) {
    // adjust time to account for central membrane delay
    source->set_t(0, m_centralMembraneDelay);
    source->set_t(1, m_centralMembraneDelay);

    // assign to positive side
    source->set_z(0, 1.);
    source->set_z(1, 1.);
    PHG4Hits.push_back(source);

    // clone
    // assign to negative side and insert in list
    auto copy = new PHG4Hitv1(source);
    copy->set_z(0, -1.);
    copy->set_z(1, -1.);
    PHG4Hits.push_back(copy);
  };

  // loop over petalID
  for (int i = 0; i < 18; i++)
  {
    // loop over radiusID
    for (int j = 0; j < 8; j++)
    {
      // loop over stripeID
      for (int k = 0; k < nGoodStripes_R1_e[j]; k++)
      {
        adjust_hits(GetPHG4HitFromStripe(i, 0, j, k, electrons_per_stripe));
      }

      // loop over stripeID
      for (int k = 0; k < nGoodStripes_R1[j]; k++)
      {
        adjust_hits(GetPHG4HitFromStripe(i, 1, j, k, electrons_per_stripe));
      }

      // loop over stripeID
      for (int k = 0; k < nGoodStripes_R2[j]; k++)
      {
        adjust_hits(GetPHG4HitFromStripe(i, 2, j, k, electrons_per_stripe));
      }

      // loop over stripeID
      for (int k = 0; k < nGoodStripes_R3[j]; k++)
      {
        adjust_hits(GetPHG4HitFromStripe(i, 3, j, k, electrons_per_stripe));
      }
    }
  }

  m_eventNum = 0;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
int PHG4TpcCentralMembrane::process_event(PHCompositeNode* topNode)
{
  
  if(m_eventNum % m_eventModulo != 0)
  {
    if(Verbosity()) std::cout << "Event " << m_eventNum << " will not generate CM hits" << std::endl;
    m_eventNum++;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // load g4hit container
  auto g4hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hitcontainer)
  {
    std::cout << PHWHERE << "Could not locate g4 hit node " << hitnodename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // copy all hits from G4hits vector into container
  for (const auto& hit : PHG4Hits)
  {
    auto copy = new PHG4Hitv1(hit);
    g4hitcontainer->AddHit(detId, copy);
  }

  m_eventNum++;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
void PHG4TpcCentralMembrane::SetDefaultParameters()
{
  // same gas parameters as in PHG4TpcElectronDrift::SetDefaultParameters

  // Data on gasses @20 C and 760 Torr from the following source:
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // diffusion and drift velocity for 400kV for NeCF4 50/50 from calculations:
  // http://skipper.physics.sunysb.edu/~prakhar/tpc/HTML_Gases/split.html
  static constexpr double Ne_dEdx = 1.56;    // keV/cm
  static constexpr double CF4_dEdx = 7.00;   // keV/cm
  static constexpr double Ne_NTotal = 43;    // Number/cm
  static constexpr double CF4_NTotal = 100;  // Number/cm
  static constexpr double Tpc_NTot = 0.5 * Ne_NTotal + 0.5 * CF4_NTotal;
  static constexpr double Tpc_dEdx = 0.5 * Ne_dEdx + 0.5 * CF4_dEdx;
  static constexpr double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;

  // number of electrons per deposited GeV in TPC gas
  set_default_double_param("electrons_per_gev", Tpc_ElectronsPerKeV * 1000000.);

  /// mean number of electrons per stripe
  set_default_int_param("electrons_per_stripe", 100);
}

//_____________________________________________________________
void PHG4TpcCentralMembrane::CalculateVertices(
    int nStripes, int nPads,
    const std::array<double, nRadii>& R,
    std::array<double, nRadii>& spacing,
    double x1a[][nRadii], double y1a[][nRadii],
    double x1b[][nRadii], double y1b[][nRadii],
    double x2a[][nRadii], double y2a[][nRadii],
    double x2b[][nRadii], double y2b[][nRadii],
    double x3a[][nRadii], double y3a[][nRadii],
    double x3b[][nRadii], double y3b[][nRadii],
    double padfrac,
    double str_width[][nRadii],
    const std::array<double, nRadii>& widthmod,
    std::array<int, nRadii>& nGoodStripes,
    const std::array<int, nRadii>& keepUntil,
    std::array<int, nRadii>& nStripesIn,
    std::array<int, nRadii>& nStripesBefore)
{
  const double phi_module = M_PI / 6.0;  // angle span of a module
  const int pr_mult = 3;                 // multiples of intrinsic resolution of pads
  const int dw_mult = 8;                 // multiples of diffusion width
  const double diffwidth = 0.6 * mm;     // diffusion width
  const double adjust = 0.015;           //arbitrary angle to center the pattern in a petal

  double theta = 0.0;
  //center coords
  double cx[nStripes][nRadii];
  double cy[nStripes][nRadii];
  //corner coords
  /* double tempX1a[nStripes][nRadii], tempY1a[nStripes][nRadii];
  double tempX1b[nStripes][nRadii], tempY1b[nStripes][nRadii];
  double tempX2a[nStripes][nRadii], tempY2a[nStripes][nRadii];
  double tempX2b[nStripes][nRadii], tempY2b[nStripes][nRadii];
  double rotatedX1a[nStripes][nRadii], rotatedY1a[nStripes][nRadii];
  double rotatedX1b[nStripes][nRadii], rotatedY1b[nStripes][nRadii];
  double rotatedX2a[nStripes][nRadii], rotatedY2a[nStripes][nRadii];
  double rotatedX2b[nStripes][nRadii], rotatedY2b[nStripes][nRadii]; */

  //calculate spacing first:
  for (int i = 0; i < nRadii; i++)
  {
    spacing[i] = 2.0 * ((dw_mult * diffwidth / R[i]) + (pr_mult * phi_module / nPads));
  }

  //vertex calculation
  for (int j = 0; j < nRadii; j++)
  {
    int i_out = 0;
    for (int i = keepThisAndAfter[j]; i < keepUntil[j]; i++)
    {
      if (j % 2 == 0)
      {
        theta = i * spacing[j] + (spacing[j] / 2) - adjust;
        cx[i_out][j] = R[j] * cos(theta);
        cy[i_out][j] = R[j] * sin(theta);
      }
      else
      {
        theta = (i + 1) * spacing[j] - adjust;
        cx[i_out][j] = R[j] * cos(theta);
        cy[i_out][j] = R[j] * sin(theta);
      }

      TVector3 corner[4];
      corner[0].SetXYZ(-padfrac + arc_r, -(widthmod[j] * str_width[i][j]) / 2, 0);  //"1a" = length of the pad, but not including the arc piece
      corner[1].SetXYZ(padfrac - arc_r, -(widthmod[j] * str_width[i][j]) / 2, 0);   //"1b" = length of the pad, but not including the arc piece
      corner[2].SetXYZ(-padfrac + arc_r, (widthmod[j] * str_width[i][j]) / 2, 0);   //"2a" = length of the pad, but not including the arc piece
      corner[3].SetXYZ(padfrac - arc_r, (widthmod[j] * str_width[i][j]) / 2, 0);    //"2b" = length of the pad, but not including the arc piece

      TVector3 rotatedcorner[4];
      for (int n = 0; n < 4; n++)
      {
        rotatedcorner[n] = corner[n];
        rotatedcorner[n].RotateZ(theta);
      }

      x1a[i_out][j] = rotatedcorner[0].X() + cx[i_out][j];
      x1b[i_out][j] = rotatedcorner[1].X() + cx[i_out][j];
      x2a[i_out][j] = rotatedcorner[2].X() + cx[i_out][j];
      x2b[i_out][j] = rotatedcorner[3].X() + cx[i_out][j];

      y1a[i_out][j] = rotatedcorner[0].Y() + cy[i_out][j];
      y1b[i_out][j] = rotatedcorner[1].Y() + cy[i_out][j];
      y2a[i_out][j] = rotatedcorner[2].Y() + cy[i_out][j];
      y2b[i_out][j] = rotatedcorner[3].Y() + cy[i_out][j];

      /* x1a[i_out][j] = cx[i_out][j] - padfrac + arc_r;
      y1a[i_out][j] = cy[i_out][j] - str_width/2;
      x1b[i_out][j] = cx[i_out][j] + padfrac - arc_r;
      y1b[i_out][j] = cy[i_out][j] - str_width/2;
      x2a[i_out][j] = cx[i_out][j] - padfrac + arc_r;
      y2a[i_out][j] = cy[i_out][j] + str_width/2;
      x2b[i_out][j] = cx[i_out][j] + padfrac - arc_r;
      y2b[i_out][j] = cy[i_out][j] + str_width/2;
      
      tempX1a[i_out][j] = x1a[i_out][j] - cx[i_out][j];
      tempY1a[i_out][j] = y1a[i_out][j] - cy[i_out][j];
      tempX1b[i_out][j] = x1b[i_out][j] - cx[i_out][j];
      tempY1b[i_out][j] = y1b[i_out][j] - cy[i_out][j];
      tempX2a[i_out][j] = x2a[i_out][j] - cx[i_out][j];
      tempY2a[i_out][j] = y2a[i_out][j] - cy[i_out][j];
      tempX2b[i_out][j] = x2b[i_out][j] - cx[i_out][j];
      tempY2b[i_out][j] = y2b[i_out][j] - cy[i_out][j];

      rotatedX1a[i_out][j] = tempX1a[i_out][j]*cos(theta) - tempY1a[i_out][j]*sin(theta);
      rotatedY1a[i_out][j] = tempX1a[i_out][j]*sin(theta) + tempY1a[i_out][j]*cos(theta);
      rotatedX1b[i_out][j] = tempX1b[i_out][j]*cos(theta) - tempY1b[i_out][j]*sin(theta);
      rotatedY1b[i_out][j] = tempX1b[i_out][j]*sin(theta) + tempY1b[i_out][j]*cos(theta);
      rotatedX2a[i_out][j] = tempX2a[i_out][j]*cos(theta) - tempY2a[i_out][j]*sin(theta);
      rotatedY2a[i_out][j] = tempX2a[i_out][j]*sin(theta) + tempY2a[i_out][j]*cos(theta);
      rotatedX2b[i_out][j] = tempX2b[i_out][j]*cos(theta) - tempY2b[i_out][j]*sin(theta);
      rotatedY2b[i_out][j] = tempX2b[i_out][j]*sin(theta) + tempY2b[i_out][j]*cos(theta);*/

      /* x1a[i_out][j] = rotatedX1a[i_out][j] + cx[i_out][j];
      y1a[i_out][j] = rotatedY1a[i_out][j] + cy[i_out][j];
      x1b[i_out][j] = rotatedX1b[i_out][j] + cx[i_out][j];
      y1b[i_out][j] = rotatedY1b[i_out][j] + cy[i_out][j];
      x2a[i_out][j] = rotatedX2a[i_out][j] + cx[i_out][j];
      y2a[i_out][j] = rotatedY2a[i_out][j] + cy[i_out][j];
      x2b[i_out][j] = rotatedX2b[i_out][j] + cx[i_out][j];
      y2b[i_out][j] = rotatedY2b[i_out][j] + cy[i_out][j]; */

      x3a[i_out][j] = (x1a[i_out][j] + x2a[i_out][j]) / 2.0;
      y3a[i_out][j] = (y1a[i_out][j] + y2a[i_out][j]) / 2.0;
      x3b[i_out][j] = (x1b[i_out][j] + x2b[i_out][j]) / 2.0;
      y3b[i_out][j] = (y1b[i_out][j] + y2b[i_out][j]) / 2.0;

      i_out++;

      nStripesBefore_R1_e[0] = 0;

      nStripesIn[j] = keepUntil[j] - keepThisAndAfter[j];
      if (j == 0)
      {
        nStripesBefore[j] = 0;
      }
      else
      {
        nStripesBefore[j] = nStripesIn[j - 1] + nStripesBefore[j - 1];
      }
      nStripesBefore_R1_e[0] = 0;
    }
    nGoodStripes[j] = i_out;
  }
}

int PHG4TpcCentralMembrane::SearchModule(int /*nStripes*/,
                                         const double x1a[][nRadii], const double x1b[][nRadii],
                                         const double x2a[][nRadii], const double x2b[][nRadii],
                                         const double y1a[][nRadii], const double y1b[][nRadii],
                                         const double y2a[][nRadii], const double y2b[][nRadii],
                                         const double x3a[][nRadii], const double y3a[][nRadii],
                                         const double x3b[][nRadii], const double y3b[][nRadii],
                                         double x, double y, const std::array<int, nRadii>& nGoodStripes) const
{
  int c = 0;

  for (int j = 0; j < nRadii; j++)
  {
    for (int i = 0; i < nGoodStripes[j]; i++)
    {
      if (((y1a[i][j] > y) != (y2a[i][j] > y) && (x < (x2a[i][j] - x1a[i][j]) * (y - y1a[i][j]) / (y2a[i][j] - y1a[i][j]) + x1a[i][j])))
        c = !c;
      if (((y1b[i][j] > y) != (y1a[i][j] > y) && (x < (x1a[i][j] - x1b[i][j]) * (y - y1b[i][j]) / (y1a[i][j] - y1b[i][j]) + x1b[i][j])))
        c = !c;
      if (((y2b[i][j] > y) != (y1b[i][j] > y) && (x < (x1b[i][j] - x2b[i][j]) * (y - y2b[i][j]) / (y1b[i][j] - y2b[i][j]) + x2b[i][j])))
        c = !c;
      if (((y2a[i][j] > y) != (y2b[i][j] > y) && (x < (x2b[i][j] - x2a[i][j]) * (y - y2a[i][j]) / (y2b[i][j] - y2a[i][j]) + x2a[i][j])))
        c = !c;

      //check inside arcs
      if (c == 0)
      {
        if (((x - x3a[i][j]) * (x - x3a[i][j]) + (y - y3a[i][j]) * (y - y3a[i][j])) <= arc_r * arc_r)
        {
          c = !c;
        }
        else if (((x - x3b[i][j]) * (x - x3b[i][j]) + (y - y3b[i][j]) * (y - y3b[i][j])) <= arc_r * arc_r)
        {
          c = !c;
        }
      }
    }
  }
  return c;
}

int PHG4TpcCentralMembrane::getSearchResult(double xcheck, double ycheck) const
{
  const double phi_petal = M_PI / 9.0;  // angle span of one petal
  const double end_R1_e = 312.0 * mm;   // arbitrary radius between R1_e and R1
  const double end_R1 = 408.0 * mm;     // arbitrary radius between R1 and R2
  const double end_R2 = 580.0 * mm;     // arbitrary radius between R2 and R3

  double r, phi, phimod, xmod, ymod;

  r = sqrt(xcheck * xcheck + ycheck * ycheck);
  phi = atan(ycheck / xcheck);
  if ((xcheck < 0.0) && (ycheck > 0.0))
  {
    phi = phi + M_PI;
  }
  else if ((xcheck > 0.0) && (ycheck < 0.0))
  {
    phi = phi + 2.0 * M_PI;
  }

  phimod = fmod(phi, phi_petal);
  xmod = r * cos(phimod);
  ymod = r * sin(phimod);

  int result = 0;

  if (r <= end_R1_e)
  {
    result = SearchModule(nStripes_R1, x1a_R1_e, x1b_R1_e, x2a_R1_e, x2b_R1_e, y1a_R1_e, y1b_R1_e, y2a_R1_e, y2b_R1_e, x3a_R1_e, y3a_R1_e, x3b_R1_e, y3b_R1_e, xmod, ymod, nGoodStripes_R1_e);
  }
  else if ((r > end_R1_e) && (r <= end_R1))
  {
    result = SearchModule(nStripes_R1, x1a_R1, x1b_R1, x2a_R1, x2b_R1, y1a_R1, y1b_R1, y2a_R1, y2b_R1, x3a_R1, y3a_R1, x3b_R1, y3b_R1, xmod, ymod, nGoodStripes_R1);
  }
  else if ((r > end_R1) && (r <= end_R2))
  {
    result = SearchModule(nStripes_R2, x1a_R2, x1b_R2, x2a_R2, x2b_R2, y1a_R2, y1b_R2, y2a_R2, y2b_R2, x3a_R2, y3a_R2, x3b_R2, y3b_R2, xmod, ymod, nGoodStripes_R2);
  }
  else if ((r > end_R2) && (r <= end_CM))
  {
    result = SearchModule(nStripes_R3, x1a_R3, x1b_R3, x2a_R3, x2b_R3, y1a_R3, y1b_R3, y2a_R3, y2b_R3, x3a_R3, y3a_R3, x3b_R3, y3b_R3, xmod, ymod, nGoodStripes_R3);
  }

  return result;
}

PHG4Hit* PHG4TpcCentralMembrane::GetPHG4HitFromStripe(int petalID, int moduleID, int radiusID, int stripeID, int nElectrons) const
{                                       //this function generates a PHG4 hit using coordinates from a stripe
  const double phi_petal = M_PI / 9.0;  // angle span of one petal
  PHG4Hit* hit;
  TVector3 dummyPos0, dummyPos1;

  //could put in some sanity checks here but probably not necessary since this is only really used within the class
  //petalID ranges 0-17, module ID 0-3, stripeID varies - nGoodStripes for each module
  //radiusID ranges 0-7

  //from phg4tpcsteppingaction.cc
  hit = new PHG4Hitv1();
  hit->set_layer(-1);  // dummy number
  //here we set the entrance values in cm
  if (moduleID == 0)
  {
    hit->set_x(0, x3a_R1_e[stripeID][radiusID] / cm);
    hit->set_y(0, y3a_R1_e[stripeID][radiusID] / cm);
  }
  else if (moduleID == 1)
  {
    hit->set_x(0, x3a_R1[stripeID][radiusID] / cm);
    hit->set_y(0, y3a_R1[stripeID][radiusID] / cm);
  }
  else if (moduleID == 2)
  {
    hit->set_x(0, x3a_R2[stripeID][radiusID] / cm);
    hit->set_y(0, y3a_R2[stripeID][radiusID] / cm);
  }
  else if (moduleID == 3)
  {
    hit->set_x(0, x3a_R3[stripeID][radiusID] / cm);
    hit->set_y(0, y3a_R3[stripeID][radiusID] / cm);
  }
  hit->set_z(0, 0.0 / cm);

  // check if you need to rotate coords to another petal
  if (petalID > 0)
  {
    dummyPos0.SetXYZ(hit->get_x(0), hit->get_y(0), hit->get_z(0));
    dummyPos0.RotateZ(petalID * phi_petal);
    hit->set_x(0, dummyPos0.X());
    hit->set_y(0, dummyPos0.Y());
  }

  // TODO: use actual stripe direction for the momentum
  hit->set_px(1, 500.0);
  hit->set_py(1, 500.0);
  hit->set_pz(1, 500.0);

  // time in ns
  hit->set_t(0, 0);

  //set and save the track ID
  hit->set_trkid(-1);  // dummy number

  // here we just update the exit values, it will be overwritten
  // for every step until we leave the volume or the particle
  // ceases to exist
  if (moduleID == 0)
  {
    hit->set_x(1, x3b_R1_e[stripeID][radiusID] / cm);
    hit->set_y(1, y3b_R1_e[stripeID][radiusID] / cm);
  }
  else if (moduleID == 1)
  {
    hit->set_x(1, x3b_R1[stripeID][radiusID] / cm);
    hit->set_y(1, y3b_R1[stripeID][radiusID] / cm);
  }
  else if (moduleID == 2)
  {
    hit->set_x(1, x3b_R2[stripeID][radiusID] / cm);
    hit->set_y(1, y3b_R2[stripeID][radiusID] / cm);
  }
  else if (moduleID == 3)
  {
    hit->set_x(1, x3b_R3[stripeID][radiusID] / cm);
    hit->set_y(1, y3b_R3[stripeID][radiusID] / cm);
  }
  hit->set_z(1, 0.0 / cm);

  // check if you need to rotate coords to another petal
  if (petalID > 0)
  {
    dummyPos1.SetXYZ(hit->get_x(1), hit->get_y(1), hit->get_z(1));
    dummyPos1.RotateZ(petalID * phi_petal);
    hit->set_x(1, dummyPos1.X());
    hit->set_y(1, dummyPos1.Y());
  }

  // TODO: use actual stripe direction for the momentum
  hit->set_px(1, 500.0);
  hit->set_py(1, 500.0);
  hit->set_pz(1, 500.0);

  hit->set_t(1, 0);  // dummy number, nanosecond

  // calculate deposited energy corresponding to number of electrons per stripe
  const double edep = nElectrons / electrons_per_gev;
  hit->set_edep(edep);
  hit->set_eion(edep);

  /*
  if (hit->get_edep()){ //print out hits
    double rin = sqrt(hit->get_x(0) * hit->get_x(0) + hit->get_y(0) * hit->get_y(0));
    double rout = sqrt(hit->get_x(1) * hit->get_x(1) + hit->get_y(1) * hit->get_y(1));
    std::cout << "Added Tpc g4hit with rin, rout = " << rin << "  " << rout
	 << " g4hitid " << hit->get_hit_id() << std::endl;
    std::cout << " xin " << hit->get_x(0)
	 << " yin " << hit->get_y(0)
	 << " zin " << hit->get_z(0)
	 << " rin " << rin
	 << std::endl;
    std::cout << " xout " << hit->get_x(1)
	 << " yout " << hit->get_y(1)
	 << " zout " << hit->get_z(1)
	 << " rout " << rout
	 << std::endl;
    std::cout << " xav " << (hit->get_x(1) + hit->get_x(0)) / 2.0
	 << " yav " << (hit->get_y(1) + hit->get_y(0)) / 2.0
	 << " zav " << (hit->get_z(1) + hit->get_z(0)) / 2.0
	 << " rav " << (rout + rin) / 2.0
	 << std::endl;
  }
  */

  return hit;
}

int PHG4TpcCentralMembrane::getStripeID(double xcheck, double ycheck) const
{
  //check if point came from stripe then see which stripe it is
  //213 stripes in a petal, 18 petals, ntotstripes = 3834
  int result;
  int fullID = -1;
  //const double adjust = 0.015; //arbitrary angle to center the pattern in a petal
  const double phi_petal = M_PI / 9.0;  // angle span of one petal

  // check if in a stripe
  result = getSearchResult(xcheck, ycheck);

  // find which stripe
  if (result == 1)
  {
    //std::cout << "on a stripe" << std::endl;
    //convert coords to radius n angle
    double r = sqrt(xcheck * xcheck + ycheck * ycheck);
    double phi = atan(ycheck / xcheck);
    if ((xcheck < 0.0) && (ycheck > 0.0))
    {
      phi = phi + M_PI;
    }
    else if ((xcheck > 0.0) && (ycheck < 0.0))
    {
      phi = phi + 2.0 * M_PI;
    }
    //get angle within first petal
    double phimod = fmod(phi, phi_petal);
    double xmod = r * cos(phimod);
    double ymod = r * sin(phimod);

    int petalID = phi / phi_petal;

    int phiID = 0;
    for (int j = 0; j < nRadii; j++)
    {
      if (((R1_e[j] - padfrac_R1) < r) && (r < (R1_e[j] + padfrac_R1)))
      {  // check if radius is in stripe
        int rID = j;
        std::cout << "rID: " << rID << std::endl;
        //'angle' is to the center of a stripe
        for (int i = 0; i < nGoodStripes_R1_e[j]; i++)
        {
          //if (j % 2 == 0){
          //theta = i*spacing[j];
          //angle = theta + (spacing[j]/2) - adjust;
          // look at distance from center line of stripe
          // if distance from x,y to center line < str_width
          // calculate slope n then do dist

          double m = (y3b_R1_e[i][j] - y3a_R1_e[i][j]) / (x3b_R1_e[i][j] - x3a_R1_e[i][j]);
          /*std::cout << "y2: " << y3b_R1_e[i][j] << std::endl;
    std::cout << "y1: " << y3a_R1_e[i][j] << std::endl;
    std::cout << "x2: " << x3b_R1_e[i][j] << std::endl;
    std::cout << "x1: " << x3a_R1_e[i][j] << std::endl;
    std::cout << "xc: " << xcheck << std::endl;
    std::cout << "yc: " << ycheck << std::endl;
	  std::cout << "m: " << m << std::endl;  */
          //std::cout << fabs((-m)*xcheck + ycheck) << std::endl;
          double dist = fabs((-m) * xmod + ymod) / sqrt(1 + m * m);
          //std::cout << "dist:" << dist << std::endl;
          if (dist < ((widthmod_R1_e[j] * str_width_R1_e[i][j]) / 2.0))
          {
            phiID = i;
            //std::cout << "phiID: " << phiID << std::endl;
          }
        }

        std::cout << "nStripesBefore: " << nStripesBefore_R1_e[j] << std::endl;
        fullID = petalID * nStripesPerPetal + nStripesBefore_R1_e[j] + phiID;
        //std::cout << "fullID: " << fullID << std::endl;
      }
      else if (((R1[j] - padfrac_R1) < r) && (r < (R1[j] + padfrac_R1)))
      {
        //std::cout << "R1" << std::endl;
        for (int i = 0; i < nGoodStripes_R1[j]; i++)
        {
          // look at distance from center line of stripe
          double m = (y3b_R1[i][j] - y3a_R1[i][j]) / (x3b_R1[i][j] - x3a_R1[i][j]);
          double dist = fabs(m * xmod - ymod) / sqrt(1 + m * m);
          if (dist < ((widthmod_R1[j] * str_width_R1[i][j]) / 2.0))
          {
            phiID = i;
          }
        }

        fullID = petalID * nStripesPerPetal + nStripesBefore_R1[j] + phiID;
      }
      else if (((R2[j] - padfrac_R2) < r) && (r < (R2[j] + padfrac_R2)))
      {
        //std::cout << "R2" << std::endl;
        for (int i = 0; i < nGoodStripes_R2[j]; i++)
        {
          // look at distance from center line of stripe
          double m = (y3b_R2[i][j] - y3a_R2[i][j]) / (x3b_R2[i][j] - x3a_R2[i][j]);
          double dist = fabs(m * xmod - ymod) / sqrt(1 + m * m);
          if (dist < ((widthmod_R2[j] * str_width_R2[i][j]) / 2.0))
          {
            phiID = i;
          }
        }

        fullID = petalID * nStripesPerPetal + nStripesBefore_R2[j] + phiID;
      }
      else if (((R3[j] - padfrac_R3) < r) && (r < (R3[j] + padfrac_R3)))
      {
        //std::cout << "R3" << std::endl;
        for (int i = 0; i < nGoodStripes_R3[j]; i++)
        {
          // look at distance from center line of stripe
          double m = (y3b_R3[i][j] - y3a_R3[i][j]) / (x3b_R3[i][j] - x3a_R3[i][j]);
          double dist = fabs(m * xmod - ymod) / sqrt(1 + m * m);
          if (dist < ((widthmod_R3[j] * str_width_R3[i][j]) / 2.0))
          {
            phiID = i;
          }
        }

        fullID = petalID * nStripesPerPetal + nStripesBefore_R3[j] + phiID;
      }
    }
  }
  else
  {
    fullID = -1;
  }

  return fullID;
}
