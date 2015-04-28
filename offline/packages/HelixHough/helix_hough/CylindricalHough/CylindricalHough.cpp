#include "CylindricalHough.h"
#include "ZHough_Cylindrical.h"


using namespace std;


bool CylindricalHough::intersect_circles(bool hel, double startx, double starty, double rad_det, double rad_trk, double cx, double cy, double& x, double& y)
{
  double cx_det = -vertex_x;
  double cy_det = -vertex_y;
  
  double d2 = ((cx-cx_det)*(cx-cx_det) + (cy-cy_det)*(cy-cy_det));
  double d = sqrt(d2);
  if(d > (rad_det + rad_trk))
  {
    return false;
  }
  if(d < fabs(rad_det - rad_trk))
  {
    return false;
  }
  
  double r2 = rad_trk*rad_trk;
  
  double d_inv = 1./d;
  double R2 = rad_det*rad_det;
  double a = 0.5*(R2 - r2 + d2)*d_inv;
  double h = a*d_inv;
  double P2x = cx_det + (cx-cx_det)*h;
  double P2y = cy_det + (cy-cy_det)*h;;
  h = sqrt(R2 - a*a);
  
  double ux = -(cy-cy_det)*d_inv;
  double uy = (cx-cx_det)*d_inv;
  double P3x1 = P2x + ux*h;
  double P3y1 = P2y + uy*h;
  ux = -ux;
  uy = -uy;
  double P3x2 = P2x + ux*h;
  double P3y2 = P2y + uy*h;
  
  double d1_2 = (startx - P3x1)*(startx - P3x1) + (starty - P3y1)*(starty - P3y1);
  double d2_2 = (startx - P3x2)*(startx - P3x2) + (starty - P3y2)*(starty - P3y2);
  
  if(d1_2 < d2_2)
  {
    x = P3x1;
    y = P3y1;
  }
  else
  {
    x = P3x2;
    y = P3y2;
  }
  
  return true;
}


CylindricalHough::CylindricalHough(vector<float>& detrad, unsigned int inv_radius_nbin, unsigned int center_angle_nbin, unsigned int dca_origin_nbin, CircleResolution& min_resolution, CircleResolution& max_resolution, CircleRange& range, unsigned int z0_nbin, unsigned int theta_nbin, ZResolution& minzres, ZResolution& maxzres, ZRange& zrange, double sxy, double sz) : CircleHough(inv_radius_nbin, center_angle_nbin, dca_origin_nbin, min_resolution, max_resolution, range, z0_nbin, theta_nbin, minzres, maxzres, zrange, sxy, sz, false), vertex_x(0.), vertex_y(0.), vertex_z(0.), phicut(0.1)
{
  for(unsigned int i=0;i<detrad.size();++i)
  {
    detector_radii.push_back(detrad[i]);
  }
  init_ZHough(z0_nbin, theta_nbin, minzres, maxzres, zrange);
}


CylindricalHough::~CylindricalHough()
{
  
}


void CylindricalHough::init_ZHough(int z0_nbin, unsigned int theta_nbin, ZResolution& minzres, ZResolution& maxzres, ZRange& zrange)
{
  if(zhough!=NULL){delete zhough;}
  zhough = new ZHough_Cylindrical(z0_nbin, theta_nbin, minzres, maxzres, zrange, sigma_xy, sigma_z);
  ((ZHough_Cylindrical*)(zhough))->setNLayers(detector_radii.size());
  vector<double> lxy, lz;
  for(unsigned int i=0;i<detector_radii.size();++i)
  {
    lxy.push_back(0.005/3.);
    lz.push_back(0.0425/3.);
  }
  ((ZHough_Cylindrical*)(zhough))->setLayerResolution(lxy, lz);
  ((ZHough_Cylindrical*)(zhough))->setVertexResolution(0.002, 0.01);
  
  zhough->setCircleHough(this);
}


void CylindricalHough::setLayerResolution(vector<double>& lxy, vector<double>& lz)
{
  ((ZHough_Cylindrical*)(zhough))->setLayerResolution(lxy, lz);
}


void CylindricalHough::setVertexResolution(double vxy, double vz)
{
  ((ZHough_Cylindrical*)(zhough))->setVertexResolution(vxy, vz);
}


void CylindricalHough::customFindHelicesInit(vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits, unsigned int min_zhits, unsigned int max_zhits, double chi2_cut, float xydiffcut, vector<SimpleTrack3D>& tracks, unsigned int maxtracks)
{
  AngleIndexList templist;
  angle_list.clear();
  angle_list.assign(detector_radii.size(), templist);
  for(unsigned int i=0;i<hits.size();i++)
  {
    AngleIndexPair temppair(atan2(hits[i].get_y(), hits[i].get_x()),hits[i].index);
    angle_list[hits[i].layer].addPair( temppair );
  }
}


//params:
//0 <--> xydiffcut
//1 <--> chi2_cut
//2 <--> min_hits
//3 <--> max_kappa_cut
void CylindricalHough::addHits(unsigned int zlevel, vector<SimpleTrack3D>& temptracks, vector<SimpleTrack3D>& tracks, vector<float>& params, int tracks_per_hit, float z_cut)
{
  vector<double> chi2_hit;
//   float xydiffcut = params[0];
  float chi2_cut = params[1];
  unsigned int min_hits = (unsigned int)(params[2]);
  float max_kappa_cut = params[3];
  
  for(unsigned int trk=0;trk<temptracks.size();++trk)
  {
    if(temptracks[trk].hits.size()<3 && using_vertex==false){continue;}
    else if(temptracks[trk].hits.size()<2 && using_vertex==true){continue;}
    zhough->fitTrack(temptracks[trk], chi2_hit);
    SimpleTrack3D temptrack = temptracks[trk];
    
    float cx =0.0; // center of rotation, x
    float cy =0.0; // center of rotation, y
    float r  =0.0; // radius of rotation
    
    if(fabs(temptracks[trk].kappa) < 1.0e-8){temptracks[trk].kappa = 1.0e-8;}
    // radius or something very straight
    r = 1./temptracks[trk].kappa;
    
    // center of rotation
    cx = (temptracks[trk].d+r)*cos(temptracks[trk].phi);
    cy = (temptracks[trk].d+r)*sin(temptracks[trk].phi);
    
    // counter-clockwise is represented by hel=true
    double v1x = temptracks[trk].hits[0].get_x();
    double v1y = temptracks[trk].hits[0].get_y();
    double v1z = temptracks[trk].hits[0].get_z();
    double v2x = temptracks[trk].hits.back().get_x();
    double v2y = temptracks[trk].hits.back().get_y();
    double diffx = v2x-v1x;
    double diffy = v2y-v1y;
    double cross = diffx*cy - diffy*cx;
    bool helicity=true;
    if(cross<0.){helicity=false;}
    
    
    vector<bool> layer_used;
    layer_used.assign(angle_list.size(), false);
    for(unsigned int h=0;h<temptracks[trk].hits.size();++h)
    {
      layer_used[temptracks[trk].hits[h].layer]=true;
    }
    vector<int> unused_layers;
    for(unsigned int ll=0;ll<layer_used.size();++ll)
    {
      if(layer_used[ll]==false){unused_layers.push_back(ll);}
    }
    
    for(unsigned int ll=0;ll<unused_layers.size();ll++)
    {
      if(((unused_layers.size() - ll) + temptracks[trk].hits.size()) < min_hits){break;}
      
      int layer = unused_layers[ll];
      int closest_layer = temptracks[trk].hits[0].layer;
      double startx = v1x;
      double starty = v1y;
      double startz = v1z;
      for(int cl=1;cl<(int)(temptracks[trk].hits.size());cl++)
      {
        if( fabs(layer - temptracks[trk].hits[cl].layer) < fabs(layer - closest_layer) )
        {
          closest_layer = temptracks[trk].hits[cl].layer;
          startx = temptracks[trk].hits[cl].get_x();
          starty = temptracks[trk].hits[cl].get_y();
          startz = temptracks[trk].hits[cl].get_z();
        }
      }
      
      double x_intersect=0.;
      double y_intersect=0.;
      double phi_error = phicut;
      
      double xy_d_cut = 3.*((ZHough_Cylindrical*)(zhough))->getLayerResolution_xy(layer);
      phi_error += xy_d_cut/(2.*M_PI*detector_radii[layer]);
//       if(phi_error < 0.02){phi_error = 0.02;}
      
      bool try_helicity = helicity;
//       if(layer < temptracks[trk].hits[0].layer){try_helicity = !helicity;}
      
      bool intersected = intersect_circles(try_helicity, startx, starty, detector_radii[layer], r, cx, cy, x_intersect, y_intersect);
      if( ( intersected==false ) || ( x_intersect != x_intersect ) || ( y_intersect != y_intersect ) )
      {
        continue;
      }
      double intersect_phi = atan2(y_intersect, x_intersect);
      if(intersect_phi < 0.){intersect_phi += 2.*M_PI;}
      double k = temptracks[trk].kappa;
      double D = sqrt((startx-x_intersect)*(startx-x_intersect) + (starty-y_intersect)*(starty-y_intersect));
      double dzdl = temptracks[trk].dzdl;
      double s=0.;
      double z_intersect=0.;
      if(0.5*k*D > 0.01)
      {
        double v = 0.5*k*D;
        if(v >= 0.999999){v = 0.999999;}
        s = 2.*asin(v)/k;
      }
      else
      {
        double temp1 = k*D*0.5;temp1*=temp1;
        double temp2 = D*0.5;
        s += 2.*temp2;
        temp2*=temp1;
        s += temp2/3.;
        temp2*=temp1;
        s += (3./20.)*temp2;
        temp2*=temp1;
        s += (5./56.)*temp2;
      }
      double dz = sqrt(s*s*dzdl*dzdl/(1. - dzdl*dzdl));
      if(layer < closest_layer){dz=-dz;}
      if(dzdl>0.){z_intersect = startz + dz;}
      else{z_intersect = startz - dz;}
      if(z_intersect != z_intersect)
      {
        continue;
      }
      
      vector<AngleIndexPair*> hit_candidates;
      angle_list[layer].getRangeList(intersect_phi, phi_error, hit_candidates);
      
      for(unsigned int hit_cand=0;hit_cand<hit_candidates.size();hit_cand++)
      {
//         //has hit been used in too many tracks?
//         if(used_vec[hit_candidates[hit_cand]->index] >= tracks_per_hit){continue;}
        //is this hit close enough in z?
        
        double z_error = z_cut + 3.*((ZHough_Cylindrical*)(zhough))->getLayerResolution_z(layer);
        if(fabs( ((*(hits_vec[0]))[hit_candidates[hit_cand]->index]).get_z() - z_intersect ) >= z_error )
        {
          continue;
        }
        
        temptracks[trk].hits.push_back((*(hits_vec[0]))[hit_candidates[hit_cand]->index]);
        //fit the track with the hit candidate added on
        double chi2 = zhough->fitTrack(temptracks[trk], chi2_hit);
        //if this hit doesn't belong, then get rid of it
        if( (chi2 != chi2) || (chi2 > chi2_cut) || (temptracks[trk].kappa > max_kappa_cut))
        {
          temptracks[trk] = temptrack;
        }
        else
        {
          temptrack = temptracks[trk];
          break;
        }
      }
    }
    if(temptracks[trk].hits.size() >= min_hits)
    {
      // loop over all tracks
      for(unsigned int itrack=0; itrack<tracks.size();itrack++)
      {
        // sort all hit indexes within the track
        sort (tracks[itrack].hits.begin(), tracks[itrack].hits.end(),SimpleHit3D_LessThan);
      }
      
      // now sort the tracks themselves
      sort (tracks.begin(), tracks.end(), SimpleTrack3D_LessThan);
      
      sort (temptracks[trk].hits.begin(), temptracks[trk].hits.end(),SimpleHit3D_LessThan);
      
      
      if(binary_search(tracks.begin(), tracks.end(), temptracks[trk], SimpleTrack3D_LessThan) == false)
      {
        // add to the output list
        tracks.push_back(temptracks[trk]);
        // add the hits to the global usage list
        for(unsigned int h=0;h<temptracks[trk].hits.size();h++)
        {
          // doesn't this double count the already used hits??? -MPM
          // no, actually this is the first time the hits get marked -Theo K
          // so the usage in ZHough during the track construction is not passed up??? -MPM
          used_vec[temptracks[trk].hits[h].index]++;
        }
      }
    }
    
  }
  
//   //---------------------------
//   // Eliminate duplicate tracks
//   //---------------------------
//   
//   // loop over all tracks
//   for(unsigned int itrack=0; itrack<tracks.size();itrack++)
//   {
//     // sort all hit indexes within the track
//     sort (tracks[itrack].hits.begin(), tracks[itrack].hits.end(),SimpleHit3D_LessThan);
//   }
//   
//   // now sort the tracks themselves
//   sort (tracks.begin(), tracks.end(), SimpleTrack3D_LessThan);
//   
//   // remove the duplicate tracks which will be neighbors after the sorting above
//   vector<SimpleTrack3D>::iterator it = unique (tracks.begin(), tracks.end(), SimpleTrack3D_Equality);
//   
//   // resize the vector to eliminate the duplicates
//   tracks.resize( it - tracks.begin() );
  
}


