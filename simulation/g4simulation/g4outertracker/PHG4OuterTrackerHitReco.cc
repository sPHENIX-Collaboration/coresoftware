// this is the new trackbase version 

#include "PHG4OuterTrackerHitReco.h"

#include <outertracker/CylinderGeomOuterTracker.h>
#include <outertracker/OuterTrackerDefs.h>
#include <outertracker/OuterTrackerHit.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>                          // for TrkrHit
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <g4detectors/PHG4CylinderGeom.h>               // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <phparameter/PHParameterContainerInterface.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                         // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                               // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                             // for PHObject
#include <phool/PHTimeServer.h>                         // for PHTimeServer
#include <phool/PHTimer.h>                              // for PHTimer
#include <phool/getClass.h>
#include <phool/phool.h>                                // for PHWHERE

#include <TVector3.h>                                   // for TVector3, ope...

#include <cmath>
#include <cstdlib>
#include <cstring>                                     // for memset
#include <iostream>
#include <memory>                                       // for allocator_tra...
#include <vector>                                       // for vector

using namespace std;

PHG4OuterTrackerHitReco::PHG4OuterTrackerHitReco(const string &name)
  : SubsysReco(name)
  , PHParameterContainerInterface(name)
  , detector(name)
  , _timer(PHTimeServer::get()->insert_new(name))
  , chkenergyconservation(0)
{
  SetDefaultParameters();  // sets default timing window
  //memset(nbins, 0, sizeof(nbins));

  if (Verbosity() > 0)
    cout << "Creating PHG4OuterTrackerHitReco for name = " << name << endl;
}

int PHG4OuterTrackerHitReco::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }
  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
  {
    cout << "Could not locate g4 hit node " << hitnodename << endl;
    exit(1);
  }

  geonodename = "CYLINDERGEOM_" + detector;
  PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  if (!geo)
  {
    cout << "Could not locate geometry node " << geonodename << endl;
    exit(1);
  }
  if (Verbosity() > 0)
  {
    geo->identify();
  }

  TrkrHitSetContainer *hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if(!hitsetcontainer)
    {
      PHNodeIterator dstiter(dstNode);
      PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
      if (!DetNode)
	{
	  DetNode = new PHCompositeNode("TRKR");
	  dstNode->addNode(DetNode);
	}

      hitsetcontainer = new TrkrHitSetContainer();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(hitsetcontainer, "TRKR_HITSET", "PHObject");
      DetNode->addNode(newNode);
    }

  TrkrHitTruthAssoc *hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");
  if(!hittruthassoc)
    {
      PHNodeIterator dstiter(dstNode);
      PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
      if (!DetNode)
	{
	  DetNode = new PHCompositeNode("TRKR");
	  dstNode->addNode(DetNode);
	}

      hittruthassoc = new TrkrHitTruthAssoc();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(hittruthassoc, "TRKR_HITTRUTHASSOC", "PHObject");
      DetNode->addNode(newNode);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4OuterTrackerHitReco::process_event(PHCompositeNode *topNode)
{
  //cout << PHWHERE << "Entering process_event for PHG4OuterTrackerHitReco" << endl;

  _timer.get()->restart();
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
  {
    cout << "Could not locate g4 hit node " << hitnodename << endl;
    exit(1);
  }

  geonodename = "CYLINDERGEOM_" + detector;
  PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  if (!geo)
  {
    cout << "Could not locate geometry node " << geonodename << endl;
    exit(1);
  }

  // Get the TrkrHitSetContainer node
  TrkrHitSetContainer *trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if(!trkrhitsetcontainer)
    {
      cout << "Could not locate TRKR_HITSET node, quit! " << endl;
      exit(1);
    }

  // Get the TrkrHitTruthAssoc node
  TrkrHitTruthAssoc *hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if(!hittruthassoc)
    {
      cout << "Could not locate TRKR_HITTRUTHASSOC node, quit! " << endl;
      exit(1);
    }

  // loop over all of the layers in the g4hit container
  PHG4HitContainer::LayerIter layer;
  pair<PHG4HitContainer::LayerIter, PHG4HitContainer::LayerIter> layer_begin_end = g4hit->getLayers();
  for (layer = layer_begin_end.first; layer != layer_begin_end.second; ++layer)
  {
    //cout << "---------- PHG4OuterTrackerHitReco:  Looping over layers " << endl;

    // loop over the hits in this layer
    PHG4HitContainer::ConstIterator hiter;
    PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);

    // we need the geometry object for this layer
    CylinderGeomOuterTracker *layergeom = dynamic_cast<CylinderGeomOuterTracker *>(geo->GetLayerGeom(*layer));
    if (!layergeom)
      exit(1);

    if (Verbosity() > 2)
      layergeom->identify();

    // Get some layer parameters for later use
    double xpixw = layergeom->get_average_radius() * layergeom->get_pixel_phi_spacing();
    double xpixw_half = xpixw / 2.0;
    double zpixw = layergeom->get_pixel_z_spacing();
    double zpixw_half = zpixw / 2.0;
    int maxNX = layergeom->get_num_pixels_phi();
    int maxNZ = layergeom->get_num_pixels_z();

    // Now loop over all g4 hits for this layer
    for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      if (Verbosity() > 4)
	{
	  cout << "From PHG4OuterTrackerHitReco: Call hit print method: " << endl;
	  hiter->second->print();
	}
      // checking ADC timing integration window cut
      if (Verbosity() > 1)
        cout << " layer " << *layer << " t0 " << hiter->second->get_t(0) << " t1 " << hiter->second->get_t(1)
             << " tmin " << tmin << " tmax " << tmax
             << endl;
      if (hiter->second->get_t(0) > tmax) continue;
      if (hiter->second->get_t(1) < tmin) continue;

      // get_property_int(const PROPERTY prop_id) const {return INT_MIN;}

      TVector3 local_in(hiter->second->get_local_x(0), hiter->second->get_local_y(0), hiter->second->get_local_z(0));
      TVector3 local_out(hiter->second->get_local_x(1), hiter->second->get_local_y(1), hiter->second->get_local_z(1));
      TVector3 midpoint((local_in.X() + local_out.X()) / 2.0, (local_in.Y() + local_out.Y()) / 2.0, (local_in.Z() + local_out.Z()) / 2.0);

      if (Verbosity() > 4)
      {
        cout << "  world entry point position: " << hiter->second->get_x(0) << " " << hiter->second->get_y(0) << " " << hiter->second->get_z(0) << endl;
        cout << "  world exit point position: " << hiter->second->get_x(1) << " " << hiter->second->get_y(1) << " " << hiter->second->get_z(1) << endl;
        cout << "  G4 local coords of entry point " << hiter->second->get_local_x(0) << " " << hiter->second->get_local_y(0) << " " << hiter->second->get_local_z(0) << endl;
        cout << "  G4 local coords of exit point " << hiter->second->get_local_x(1) << " " << hiter->second->get_local_y(1) << " " << hiter->second->get_local_z(1) << endl;
        cout << endl;
      }

      // Get the pixel number of the entry location
      int pixel_phi_index_in, pixel_z_index_in;
      layergeom->find_pixel_index_values(local_in[0], local_in[1], local_in[2], pixel_phi_index_in, pixel_z_index_in);
      // Get the pixel number of the exit location
      int pixel_phi_index_out, pixel_z_index_out;
      layergeom->find_pixel_index_values(local_out[0], local_out[1], local_out[2], pixel_phi_index_out, pixel_z_index_out);

      if (Verbosity() > 0)
	{ 
	  cout << "entry pixel phi index " << pixel_phi_index_in << " exit pixel phi_index " << pixel_phi_index_out << endl;
	  cout << "exit pixel phi index " << pixel_phi_index_out << " exit pixel phi_index " << pixel_phi_index_out << endl;
	}

      //vector<int> vpixel;
      vector<int> vxbin;
      vector<int> vzbin;
      vector<double> vlen;
      vector<pair<double, double> > venergy;
      //double trklen = 0.0;

      //===================================================
      // OK, now we have found which sensor the hit is in, extracted the hit
      // position in local coordinates,  and found the pixel numbers of the
      // entry point and exit point

      //====================================================
      // Beginning of charge sharing implementation
      //    Find tracklet line inside sensor
      //    Divide tracklet line into n segments (vary n until answer stabilizes)
      //    Find centroid of each segment
      //    Diffuse charge at each centroid
      //    Apportion charge between neighboring pixels
      //    Add the pixel energy contributions from different track segments together
      //====================================================

      // Because this detector is a cylinder centered on zero, the best coordinate system is cylindrical (r,phi,z)
      // In local coords, the pixels are flat squares in the space of (phi*radius, z) and their thickness is in the radial direction
      // So we want to express our segment location in (r, phi, z) coords

      double R_in = sqrt(pow(local_in.X(), 2)+pow(local_in.Y(), 2)) ;
      double R_out = sqrt(pow(local_out.X(), 2)+pow(local_out.Y(), 2)) ;
      double Phi_in = atan2(local_in.Y(), local_in.X());
      double Phi_out = atan2(local_out.Y(), local_out.X());
      double Z_in = local_in.Z();
      double Z_out = local_out.Z();

      double path_dR, path_dPhi, path_dZ;
      path_dR = R_out - R_in;
      path_dPhi = Phi_out - Phi_in;
      path_dZ = Z_out - Z_in;

      // See figure 7.3 of the thesis by  Lucasz Maczewski (arXiv:10053.3710) for diffusion simulations in a MAPS epitaxial layer
      // The diffusion widths below were inspired by those plots, corresponding to where the probability drops off to 1/3 of the peak value
      // However note that we make the simplifying assumption that the probability distribution is flat within this diffusion width,
      // while in the simulation it is not
      double diffusion_width_max = 25.0e-04;  // maximum diffusion radius 35 microns, in cm
      double diffusion_width_min = 8.0e-04;   // minimum diffusion radius 12 microns, in cm

      double ydrift_max = path_dR;
      int nsegments = 4;

      // we want to make a list of all pixels possibly affected by this hit
      // we take the entry and exit locations in local coordinates, and build
      // a rectangular array of pixels that encompasses both, with "nadd" pixels added all around

      int xbin_in = pixel_phi_index_in;
      int zbin_in = pixel_z_index_in;
      int xbin_out = pixel_phi_index_out;
      int zbin_out = pixel_z_index_out;

      int xbin_max, xbin_min;
      int nadd = 2;
      if (xbin_in > xbin_out)
      {
        xbin_max = xbin_in + nadd;
        xbin_min = xbin_out - nadd;
      }
      else
      {
        xbin_max = xbin_out + nadd;
        xbin_min = xbin_in - nadd;
      }

      int zbin_max, zbin_min;
      if (zbin_in > zbin_out)
      {
        zbin_max = zbin_in + nadd;
        zbin_min = zbin_out - nadd;
      }
      else
      {
        zbin_max = zbin_out + nadd;
        zbin_min = zbin_in - nadd;
      }

      // need to check that values of xbin and zbin are within the valid range
      if (xbin_min < 0) xbin_min = 0;
      if (zbin_min < 0) zbin_min = 0;
      if (xbin_max > maxNX) xbin_max = maxNX;
      if (zbin_max > maxNZ) zbin_max = maxNZ;

      if (Verbosity() > 1)
      {
        cout << " xbin_in " << xbin_in << " xbin_out " << xbin_out << " xbin_min " << xbin_min << " xbin_max " << xbin_max << endl;
        cout << " zbin_in " << zbin_in << " zbin_out " << zbin_out << " zbin_min " << zbin_min << " zbin_max " << zbin_max << endl;
      }

      // skip this hit if it involves an unreasonable  number of pixels
      // this skips it if either the xbin or ybin range traversed is greater than 8 (for 8 adding two pixels at each end makes the range 12)
      if (xbin_max - xbin_min > 12 || zbin_max - zbin_min > 12)
        continue;

      // this hit is skipped earlier if this dimensioning would be exceeded
      double pixenergy[13][13] = {};  // init to 0
      double pixeion[13][13] = {};    // init to 0

      // Loop over track segments and diffuse charge at each segment location, collect energy in pixels
      for (int i = 0; i < nsegments; i++)
      {
        // Find the tracklet segment location
        // If there are n segments of equal length, we want 2*n intervals
        // The 1st segment is centered at interval 1, the 2nd at interval 3, the nth at interval 2n -1
        double interval = 2 * (double) i + 1;
        double frac = interval / (double) (2 * nsegments);

	double segR = R_in + path_dR * frac;
	double segPhiR = segR * (Phi_in + path_dPhi * frac);
	double segZ = Z_in + path_dZ * frac;

        //  Find the distance to the back of the sensor from the segment location
        // That projection changes only the value of y
        double ydrift = R_out - segR;

        // Caculate the charge diffusion over this drift distance
        // increases from diffusion width_min to diffusion_width_max
        double ydiffusion_radius = diffusion_width_min + (ydrift / ydrift_max) * (diffusion_width_max - diffusion_width_min);

        if (Verbosity() > 10)
          cout << " segment " << i
               << " interval " << interval
               << " frac " << frac
               << " local_in.X " << local_in.X()
               << " local_in.Z " << local_in.Z()
               << " local_in.Y " << local_in.Y()
               << " path_dR " << path_dR
               << " path_dPhi " << path_dPhi
               << " path_dZ " << path_dZ
               << " segR " << segR
               << " segZ " << segZ
               << " segPhiR " << segPhiR
               << " ydrift " << ydrift
               << " ydrift_max " << ydrift_max
               << " ydiffusion_radius " << ydiffusion_radius
               << endl;

        // Now find the area of overlap of the diffusion circle with each pixel and apportion the energy
        for (int ix = xbin_min; ix <= xbin_max; ix++)
        {
          for (int iz = zbin_min; iz <= zbin_max; iz++)
          {

	    double phicenter, zcenter, pix_radius;
	    layergeom->find_pixel_center(ix, iz, phicenter, zcenter);
	    pix_radius = layergeom->get_average_radius();
	    //cout << " trying pixel phi center " << phicenter << " pixel z center " << zcenter << endl; 
            // note that (x1,z1) is the top left corner, (x2,z2) is the bottom right corner of the pixel - circle_rectangle_intersection expects this ordering
            double x1 = pix_radius * phicenter - xpixw_half;  // lower limit in phi direction
            double z1 = zcenter + zpixw_half;  //z upper limit
            double x2 = pix_radius * phicenter + xpixw_half;  // upper limit in phi direction
            double z2 = zcenter - zpixw_half;  // z lower limit
	    //cout << " xpixw_half " << xpixw_half << " zpixw_half " << zpixw_half << " x1, z1, x2, z2 " << x1 << ", "<< z1 << ", " << x2 << ", " << z2 << endl;
            // here segvec.X and segvec.Z are the center of the circle, and diffusion_radius is the circle radius
            // circle_rectangle_intersection returns the overlap area of the circle and the pixel. It is very fast if there is no overlap.
            double pixarea_frac = circle_rectangle_intersection(x1, z1, x2, z2, segPhiR, segZ, ydiffusion_radius) / (M_PI * pow(ydiffusion_radius, 2));
            // assume that the energy is deposited uniformly along the tracklet length, so that this segment gets the fraction 1/nsegments of the energy
            double this_pix_energy = pixarea_frac * hiter->second->get_edep() / (float) nsegments;
            pixenergy[ix - xbin_min][iz - zbin_min] += this_pix_energy;
            if (hiter->second->has_property(PHG4Hit::prop_eion))
            {
              pixeion[ix - xbin_min][iz - zbin_min] += pixarea_frac * hiter->second->get_eion() / (float) nsegments;
            }
            if (this_pix_energy > 0 && Verbosity() > 5)
            {
              cout << "     xbin " << ix << " zbin " << iz
                   << " pixel_area fraction of circle " << pixarea_frac << " accumulated pixel energy " << pixenergy[ix - xbin_min][iz - zbin_min]
                   << endl;
            }
          }
        }
      }  // end loop over segments

      // now we have the energy deposited in each pixel, summed over all tracklet segments. We make a vector of all pixels with non-zero energy deposited
      for (int ix = xbin_min; ix <= xbin_max; ix++)
      {
        for (int iz = zbin_min; iz <= zbin_max; iz++)
        {
          if (pixenergy[ix - xbin_min][iz - zbin_min] > 0.0)
          {
            //int pixnum = layergeom->get_pixel_number_from_xbin_zbin(ix, iz);
            //vpixel.push_back(pixnum);
            vxbin.push_back(ix);
            vzbin.push_back(iz);
            pair<double, double> tmppair = make_pair(pixenergy[ix - xbin_min][iz - zbin_min], pixeion[ix - xbin_min][iz - zbin_min]);
            venergy.push_back(tmppair);
            if (Verbosity() > 1)
              cout << " Added pixel with  xbin " << ix << " zbin " << iz << " to vectors with energy " << pixenergy[ix - xbin_min][iz - zbin_min] << endl;
          }
        }
      }

      //===================================
      // End of charge sharing implementation
      //===================================

      // loop over all fired cells for this g4hit and add them to the TrkrHitSet
      if(Verbosity() > 5) cout << " number of hits in OuterTracker layer " << *layer << " is " << vxbin.size() << endl;
      for (unsigned int i1 = 0; i1 < vxbin.size(); i1++)  // loop over all fired cells
      {
	// This is the new storage object version
	//====================================

	// We need to create the TrkrHitSet if not already made - each TrkrHitSet should correspond to a chip for the OuterTracker	
	TrkrDefs::hitsetkey hitsetkey = OuterTrackerDefs::genHitSetKey(*layer);
	// Use existing hitset or add new one if needed
	TrkrHitSetContainer::Iterator hitsetit = trkrhitsetcontainer->findOrAddHitSet(hitsetkey);

	// generate the key for this hit
	TrkrDefs::hitkey hitkey = OuterTrackerDefs::genHitKey(vzbin[i1], vxbin[i1]);
	// See if this hit already exists
	TrkrHit *hit = nullptr;
	hit = hitsetit->second->getHit(hitkey);
	if(!hit)
	  {
	    // Otherwise, create a new one
	    hit = new OuterTrackerHit();
	    hitsetit->second->addHitSpecificKey(hitkey, hit);
	  }

	// Either way, add the energy to it
	hit->addEnergy(venergy[i1].first);
	
	// now we update the TrkrHitTruthAssoc map - the map contains <hitsetkey, std::pair <hitkey, g4hitkey> >
	// There is only one TrkrHit per pixel, but there may be multiple g4hits
	// How do we know how much energy from PHG4Hit went into TrkrHit? We don't, have to sort it out in evaluator to save memory

	// How do we check if this association already exists?
	if(Verbosity() > 5) cout << "PHG4OuterTrackerHitReco: adding association entry for hitkey " << hitkey << " and g4hitkey " << hiter->first << endl; 
	hittruthassoc->addAssoc(hitsetkey, hitkey, hiter->first);
	
      } // end loop over hit cells
    }  // end loop over g4hits for this layer

  } // end loop over layers

  // print the list of entries in the association table
  if(Verbosity() > 2)
  {
    cout << "From PHG4OuterTrackerHitReco: " << endl;
    hittruthassoc->identify();
  }

  if (chkenergyconservation)
  {
    CheckEnergy(topNode);
  }
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4OuterTrackerHitReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4OuterTrackerHitReco::CheckEnergy(PHCompositeNode *topNode)
{
  return 0;
}

bool PHG4OuterTrackerHitReco::lines_intersect(
    double ax,
    double ay,
    double bx,
    double by,
    double cx,
    double cy,
    double dx,
    double dy,
    double *rx,  // intersection point (output)
    double *ry)
{
  // Find if a line segment limited by points A and B
  // intersects line segment limited by points C and D.
  // First check if an infinite line defined by A and B intersects
  // segment (C,D). If h is from 0 to 1 line and line segment intersect
  // Then check in intersection point is between C and D

  double ex = bx - ax;  // E=B-A
  double ey = by - ay;
  double fx = dx - cx;  // F=D-C
  double fy = dy - cy;
  double px = -ey;  // P
  double py = ex;

  double bottom = fx * px + fy * py;  // F*P
  double gx = ax - cx;                // A-C
  double gy = ay - cy;
  double top = gx * px + gy * py;  // G*P

  double h = 99999.;
  if (bottom != 0.)
  {
    h = top / bottom;
  }

  //intersection point R = C + F*h
  if (h > 0. && h < 1.)
  {
    *rx = cx + fx * h;
    *ry = cy + fy * h;
    //cout << "      line/segment intersection coordinates: " << *rx << " " << *ry << endl;
    if ((*rx > ax && *rx > bx) || (*rx < ax && *rx < bx) || (*ry < ay && *ry < by) || (*ry > ay && *ry > by))
    {
      //cout << "       NO segment/segment intersection!" << endl;
      return false;
    }
    else
    {
      //cout << "       segment/segment intersection!" << endl;
      return true;
    }
  }

  return false;
}

bool PHG4OuterTrackerHitReco::line_and_rectangle_intersect(
    double ax,
    double ay,
    double bx,
    double by,
    double cx,
    double cy,
    double dx,
    double dy,
    double *rr  // length of the line segment inside the rectangle (output)
)
{
  // find if a line isegment limited by points (A,B)
  // intersects with a rectangle defined by two
  // corner points (C,D) two other points are E and F
  //   E--------D
  //   |        |
  //   |        |
  //   C--------F

  if (cx > dx || cy > dy)
  {
    cerr << "ERROR: Bad rectangle definition!" << endl;
    return false;
  }

  double ex = cx;
  double ey = dy;
  double fx = dx;
  double fy = cy;
  double rx = 99999.;
  double ry = 99999.;

  vector<double> vx;
  vector<double> vy;

  bool i1 = lines_intersect(ax, ay, bx, by, cx, cy, fx, fy, &rx, &ry);
  if (i1)
  {
    vx.push_back(rx);
    vy.push_back(ry);
  }
  bool i2 = lines_intersect(ax, ay, bx, by, fx, fy, dx, dy, &rx, &ry);
  if (i2)
  {
    vx.push_back(rx);
    vy.push_back(ry);
  }
  bool i3 = lines_intersect(ax, ay, bx, by, ex, ey, dx, dy, &rx, &ry);
  if (i3)
  {
    vx.push_back(rx);
    vy.push_back(ry);
  }
  bool i4 = lines_intersect(ax, ay, bx, by, cx, cy, ex, ey, &rx, &ry);
  if (i4)
  {
    vx.push_back(rx);
    vy.push_back(ry);
  }

  //cout << "Rectangle intersections: " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;
  //cout << "Number of intersections = " << vx.size() << endl;

  *rr = 0.;
  if (vx.size() == 2)
  {
    *rr = sqrt((vx[0] - vx[1]) * (vx[0] - vx[1]) + (vy[0] - vy[1]) * (vy[0] - vy[1]));
    //  cout << "Length of intersection = " << *rr << endl;
  }
  if (vx.size() == 1)
  {
    // find which point (A or B) is within the rectangle
    if (ax > cx && ay > cy && ax < dx && ay < dy)  // point A is inside the rectangle
    {
      //cout << "Point A is inside the rectangle." << endl;
      *rr = sqrt((vx[0] - ax) * (vx[0] - ax) + (vy[0] - ay) * (vy[0] - ay));
    }
    if (bx > cx && by > cy && bx < dx && by < dy)  // point B is inside the rectangle
    {
      //cout << "Point B is inside the rectangle." << endl;
      *rr = sqrt((vx[0] - bx) * (vx[0] - bx) + (vy[0] - by) * (vy[0] - by));
    }
  }

  if (i1 || i2 || i3 || i4)
  {
    return true;
  }
  return false;
}

double PHG4OuterTrackerHitReco::circle_rectangle_intersection(double x1, double y1, double x2, double y2, double mx, double my, double r)
{
  // Find the area of overlap of a circle and rectangle
  // Calls sA, which uses an analytic formula to determine the integral of the circle between limits set by the corners of the rectangle

  // move the rectangle to the frame where the circle is at (0,0)
  x1 -= mx;
  x2 -= mx;
  y1 -= my;
  y2 -= my;

  if (Verbosity() > 100)
  {
    cout << " mx " << mx << " my " << my << " r " << r << " x1 " << x1 << " x2 " << x2 << " y1 " << y1 << " y2 " << y2 << endl;
    cout << " sA21 " << sA(r, x2, y1)
         << " sA11 " << sA(r, x1, y1)
         << " sA22 " << sA(r, x2, y2)
         << " sA12 " << sA(r, x1, y2)
         << endl;
  }

  return sA(r, x2, y1) - sA(r, x1, y1) - sA(r, x2, y2) + sA(r, x1, y2);
}

double PHG4OuterTrackerHitReco::sA(double r, double x, double y)
{
  // Uses analytic formula for the integral of a circle between limits set by the corner of a rectangle
  // It is called repeatedly to find the overlap area between the circle and rectangle
  // I found this code implementing the integral on a web forum called "ars technica",
  // https://arstechnica.com/civis/viewtopic.php?t=306492
  // posted by "memp"

  double a;

  if (x < 0)
  {
    return -sA(r, -x, y);
  }

  if (y < 0)
  {
    return -sA(r, x, -y);
  }

  if (x > r)
  {
    x = r;
  }

  if (y > r)
  {
    y = r;
  }

  if (x * x + y * y > r * r)
  {
    a = r * r * asin(x / r) + x * sqrt(r * r - x * x) + r * r * asin(y / r) + y * sqrt(r * r - y * y) - r * r * M_PI_2;

    a *= 0.5;
  }
  else
  {
    a = x * y;
  }

  return a;
}

void PHG4OuterTrackerHitReco::set_timing_window(const int detid, const double tmin_in, const double tmax_in)
{
  // now replace it with the new value
  tmin = tmin_in;
  tmax = tmax_in;

  cout << "PHG4OuterTrackerHitReco: Set OuterTracker timing window parameters from macro for layer = " << detid << " to tmin = " << tmin << " tmax = " << tmax << endl;

  return;
}

void PHG4OuterTrackerHitReco::SetDefaultParameters()
{
  cout << "PHG4OuterTrackerHitReco: Setting OuterTracker timing window defaults to tmin = -50 and  tmax = 50 ns" << endl;
  tmin = -50.0;
  tmax = 50.0;

  return;
}
