#include "ZHough_Cylindrical.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include "CircleHough.h"
#include <Eigen/LU>
#include <Eigen/Core>


using namespace std;
using namespace Eigen;


ZHough_Cylindrical::ZHough_Cylindrical(unsigned int z0_nbin, unsigned int theta_nbin, ZResolution& min_resolution, ZResolution& max_resolution, ZRange& range, double sxy, double sz) : ZHough(z0_nbin, theta_nbin, min_resolution, max_resolution, range, sxy, sz), nlayers(6)
{
  
}


ZHough_Cylindrical::~ZHough_Cylindrical()
{
  
}


void ZHough_Cylindrical::findTracks(unsigned int zoomlevel, float xydiffcut, unsigned int max_hits, unsigned int tracks_per_hit, double chi2_cut, float max_kappa_cut, vector<SimpleTrack3D>& tracks, vector<SimpleHit3D>& hits, float min_k, float max_k, float min_phi, float max_phi, float min_d, float max_d, float min_z0, float max_z0, float min_theta, float max_theta, vector<float>& params)
{
  if(_verbosity > 0)
  {
    cout<<"findTracks"<<endl;
    cout<<min_k<<" "<<max_k<<" "<<min_phi<<" "<<max_phi<<" "<<min_d<<" "<<max_d<<" "<<min_z0<<" "<<max_z0<<" "<<min_theta<<" "<<max_theta<<endl;
  }
  
  vector<SimpleHit3D> comb_hits = (*(hits_vec[zoomlevel+1]));
  if(using_vertex==false)
  {
    findTracksCombo_noVertex(comb_hits, zoomlevel, xydiffcut, max_hits, tracks_per_hit, chi2_cut, max_kappa_cut, tracks, hits, min_k, max_k, min_phi, max_phi, min_d, max_d, min_z0, max_z0, min_theta, max_theta, params);
  }
  else
  {
    findTracksCombo_withVertex(comb_hits, zoomlevel, xydiffcut, max_hits, tracks_per_hit, chi2_cut, max_kappa_cut, tracks, hits, min_k, max_k, min_phi, max_phi, min_d, max_d, min_z0, max_z0, min_theta, max_theta, params);
  }
  
}


void ZHough_Cylindrical::findTracksCombo_noVertex(vector<SimpleHit3D>& comb_hits,  unsigned int zoomlevel, float xydiffcut, unsigned int max_hits, unsigned int tracks_per_hit, double chi2_cut, float max_kappa_cut, vector<SimpleTrack3D>& tracks, vector<SimpleHit3D>& hits, float min_k, float max_k, float min_phi, float max_phi, float min_d, float max_d, float min_z0, float max_z0, float min_theta, float max_theta, vector<float>& params)
{
  unsigned int req_hits = (unsigned int)(params[2]);
  
  vector<unsigned int> index_holder;
  index_holder.assign(comb_hits.size(), 0);
  for(unsigned int i=0;i<comb_hits.size();++i)
  {
    index_holder[i]=comb_hits[i].index;
    comb_hits[i].index = i;
  }
  
  double phidiff;
  vector<double> chi2_hits;
  double chi2=0.;
  vector<bool> use_here;
  use_here.assign(comb_hits.size(), true);
  SimpleTrack3D temp_track;
  temp_track.hits.resize(3, comb_hits[0]);
  
  //loop over combinations of 3 hits, each from a different layer
  for(unsigned int i1=0;i1<comb_hits.size();++i1)
  {
    if(use_here[i1]==false){continue;}
    if(used_vec[index_holder[i1]] >= tracks_per_hit){continue;}
    temp_track.hits[0] = comb_hits[i1];
    for(unsigned int i2=(i1+1);i2<comb_hits.size();++i2)
    {
      if(use_here[i2]==false){continue;}
      if(used_vec[index_holder[i1]] >= tracks_per_hit){continue;}
      if(used_vec[index_holder[i2]] >= tracks_per_hit){continue;}
      if(comb_hits[i2].layer == comb_hits[i1].layer){continue;}
      phidiff = temp_track.hits[0].get_phi() - comb_hits[i2].get_phi();
      if(phidiff < -M_PI){phidiff += 2.*M_PI;}
      if(phidiff > M_PI){phidiff -= 2.*M_PI;}
      if(fabs(phidiff) > M_PI/3.){continue;}
      temp_track.hits[1] = comb_hits[i2];
      float xdiff2 = (temp_track.hits[0].get_x() - temp_track.hits[1].get_x());
      xdiff2*=xdiff2;
      float ydiff2 = (temp_track.hits[0].get_y() - temp_track.hits[1].get_y());
      ydiff2*=ydiff2;
      float xydiff2 = xdiff2 + ydiff2;
      float zdiff2 = (temp_track.hits[0].get_z() - temp_track.hits[1].get_z());
      zdiff2*=zdiff2;
      for(unsigned int i3=(i2+1);i3<comb_hits.size();++i3)
      {
        if(use_here[i3]==false){continue;}
        if(used_vec[index_holder[i1]] >= tracks_per_hit){continue;}
        if(used_vec[index_holder[i2]] >= tracks_per_hit){continue;}
        if(used_vec[index_holder[i3]] >= tracks_per_hit){continue;}
        if((comb_hits[i3].layer == comb_hits[i2].layer) || (comb_hits[i3].layer == comb_hits[i1].layer)){continue;}
        phidiff = temp_track.hits[1].get_phi() - comb_hits[i3].get_phi();
        if(phidiff < -M_PI){phidiff += 2.*M_PI;}
        if(phidiff > M_PI){phidiff -= 2.*M_PI;}
        if(fabs(phidiff) > M_PI/3.){continue;}
        temp_track.hits[2] = comb_hits[i3];
        float xdiff2_2 = (temp_track.hits[1].get_x() - temp_track.hits[2].get_x());
        xdiff2_2*=xdiff2_2;
        float ydiff2_2 = (temp_track.hits[1].get_y() - temp_track.hits[2].get_y());
        ydiff2_2*=ydiff2_2;
        float xydiff2_2 = xdiff2_2 + ydiff2_2;
        float zdiff2_2 = (temp_track.hits[1].get_z() - temp_track.hits[2].get_z());
        zdiff2_2*=zdiff2_2;
        
        if(fabs(xydiff2*zdiff2_2/(zdiff2*xydiff2_2) - 1.) > 0.2){continue;}
        
        temp_track.hits.resize(3, comb_hits[0]);//make sure there are only 3 hits in this track
        chi2 = fitTrack(temp_track, chi2_hits);
        if((chi2 < chi2_cut)&&(temp_track.kappa < max_kappa_cut))
        {
          //we found a track candidate.  see if there are more hits in comb_hits to add to this track
          vector<SimpleTrack3D> ttracks;
          SimpleTrack3D ttrack = temp_track;
          for(unsigned int j=0;j<temp_track.hits.size();++j)
          {
            ttrack.hits[j].index = index_holder[ttrack.hits[j].index];
            ttrack.hits[j].index = index_mapping[ttrack.hits[j].index];
          }
          ttracks.push_back(ttrack);
          unsigned int start_size = tracks.size();
          chough->addHits((unsigned int)(params[4]), ttracks, tracks, params, tracks_per_hit, params[5]);
          if(tracks.size() > start_size)
          {
            if(tracks.back().hits.size() >= req_hits)
            {
              unsigned int count = 0;
              for(unsigned int j=0;j<temp_track.hits.size();++j)
              {
                use_here[temp_track.hits[j].index] = false;
              }
              for(unsigned int j=0;j<temp_track.hits.size();++j)
              {
                used_vec[index_holder[temp_track.hits[j].index]]++;
                count++;
              }
              for(unsigned int j=temp_track.hits.size();j<tracks.back().hits.size();++j)
              {
                for(unsigned int i=0;i<(*(hits_vec[zoomlevel+1])).size();++i)
                {
                  if(tracks.back().hits[j].index == index_mapping[(*(hits_vec[zoomlevel+1]))[i].index])
                  {
                    used_vec[(*(hits_vec[zoomlevel+1]))[i].index]++;
                    count++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  //restore indexes
  for(unsigned int i=0;i<comb_hits.size();++i)
  {
    comb_hits[i].index = index_holder[i];
  }
}


void ZHough_Cylindrical::findTracksCombo_withVertex(vector<SimpleHit3D>& comb_hits,  unsigned int zoomlevel, float xydiffcut, unsigned int max_hits, unsigned int tracks_per_hit, double chi2_cut, float max_kappa_cut, vector<SimpleTrack3D>& tracks, vector<SimpleHit3D>& hits, float min_k, float max_k, float min_phi, float max_phi, float min_d, float max_d, float min_z0, float max_z0, float min_theta, float max_theta, vector<float>& params)
{
  unsigned int req_hits = (unsigned int)(params[2]);
  
  vector<unsigned int> index_holder;
  index_holder.assign(comb_hits.size(), 0);
  for(unsigned int i=0;i<comb_hits.size();++i)
  {
    index_holder[i]=comb_hits[i].index;
    comb_hits[i].index = i;
  }
  
  double phidiff;
  vector<double> chi2_hits;
  double chi2=0.;
  vector<bool> use_here;
  use_here.assign(comb_hits.size(), true);
  SimpleTrack3D temp_track;
  temp_track.hits.resize(2, comb_hits[0]);
  
  //loop over combinations of 2 hits, each from a different layer
  for(unsigned int i1=0;i1<comb_hits.size();++i1)
  {
    if(use_here[i1]==false){continue;}
    if(used_vec[index_holder[i1]] >= tracks_per_hit){continue;}
    temp_track.hits[0] = comb_hits[i1];
    
    float xdiff2_v = (temp_track.hits[0].get_x());
    xdiff2_v*=xdiff2_v;
    float ydiff2_v = (temp_track.hits[0].get_y());
    ydiff2_v*=ydiff2_v;
    float xydiff2_v = xdiff2_v + ydiff2_v;
    float zdiff2_v = (temp_track.hits[0].get_z());
    zdiff2_v*=zdiff2_v;
    
    for(unsigned int i2=(i1+1);i2<comb_hits.size();++i2)
    {
      if(use_here[i1]==false){continue;}
      if(use_here[i2]==false){continue;}
      if(comb_hits[i2].layer == comb_hits[i1].layer){continue;}
      if(used_vec[index_holder[i1]] >= tracks_per_hit){continue;}
      if(used_vec[index_holder[i2]] >= tracks_per_hit){continue;}
      phidiff = temp_track.hits[0].get_phi() - comb_hits[i2].get_phi();
      if(phidiff < -M_PI){phidiff += 2.*M_PI;}
      if(phidiff > M_PI){phidiff -= 2.*M_PI;}
      if(fabs(phidiff) > M_PI/3.){continue;}
      temp_track.hits[1] = comb_hits[i2];
      float xdiff2 = (temp_track.hits[0].get_x() - temp_track.hits[1].get_x());
      xdiff2*=xdiff2;
      float ydiff2 = (temp_track.hits[0].get_y() - temp_track.hits[1].get_y());
      ydiff2*=ydiff2;
      float xydiff2 = xdiff2 + ydiff2;
      float zdiff2 = (temp_track.hits[0].get_z() - temp_track.hits[1].get_z());
      zdiff2*=zdiff2;
      if(fabs(xydiff2*zdiff2_v/(zdiff2*xydiff2_v) - 1.) > 0.2){continue;}
      
      temp_track.hits.resize(2, comb_hits[0]);//make sure there are only 2 hits in this track
      chi2 = fitTrack(temp_track, chi2_hits);
      if((chi2 < chi2_cut)&&(temp_track.kappa < max_kappa_cut))
      {
        //we found a track candidate.  see if there are more hits in comb_hits to add to this track
        vector<SimpleTrack3D> ttracks;
        SimpleTrack3D ttrack = temp_track;
        for(unsigned int j=0;j<temp_track.hits.size();++j)
        {
          ttrack.hits[j].index = index_holder[ttrack.hits[j].index];
          ttrack.hits[j].index = index_mapping[ttrack.hits[j].index];
        }
        ttracks.push_back(ttrack);
        unsigned int start_size = tracks.size();
        chough->addHits((unsigned int)(params[4]), ttracks, tracks, params, tracks_per_hit, params[5]);
        if(tracks.size() > start_size)
        {
          if(tracks.back().hits.size() >= req_hits)
          {
            unsigned int count = 0;
            for(unsigned int j=0;j<temp_track.hits.size();++j)
            {
              use_here[temp_track.hits[j].index] = false;
            }
            for(unsigned int j=0;j<temp_track.hits.size();++j)
            {
              used_vec[index_holder[temp_track.hits[j].index]]++;
              count++;
            }
            for(unsigned int j=temp_track.hits.size();j<tracks.back().hits.size();++j)
            {
              for(unsigned int i=0;i<(*(hits_vec[zoomlevel+1])).size();++i)
              {
                if(tracks.back().hits[j].index == index_mapping[(*(hits_vec[zoomlevel+1]))[i].index])
                {
                  used_vec[(*(hits_vec[zoomlevel+1]))[i].index]++;
                  count++;
                }
              }
            }
          }
        }
      }
    }
  }
}



double ZHough_Cylindrical::fitTrack(SimpleTrack3D& track, vector<double>& chi2_hit)
{
  if(using_vertex == true)
  {
    track.hits.push_back(SimpleHit3D(0.,0.,0., 0, 0));
  }
  
  chi2_hit.clear();
  chi2_hit.resize(track.hits.size(), 0.);
  
  MatrixXd y = MatrixXd::Zero(track.hits.size(), 1);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    y(i, 0) = ( pow(track.hits[i].get_x(),2) + pow(track.hits[i].get_y(),2) );
    if((using_vertex==true ) && (i == (track.hits.size() - 1))){y(i, 0) /= vertex_sigma_xy;}
    else{y(i, 0) /= layer_xy_resolution[track.hits[i].layer];}
  }
  
  MatrixXd X = MatrixXd::Zero(track.hits.size(), 3);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    X(i, 0) = track.hits[i].get_x();
    X(i, 1) = track.hits[i].get_y();
    X(i, 2) = -1.;
    if((using_vertex==true ) && (i == (track.hits.size() - 1)))
    {
      X(i, 0) /= vertex_sigma_xy;
      X(i, 1) /= vertex_sigma_xy;
      X(i, 2) /= vertex_sigma_xy;
    }
    else
    {
      X(i, 0) /= layer_xy_resolution[track.hits[i].layer];
      X(i, 1) /= layer_xy_resolution[track.hits[i].layer];
      X(i, 2) /= layer_xy_resolution[track.hits[i].layer];
    }
  }
  
  MatrixXd Xt = X.transpose();
  
  MatrixXd prod = Xt*X;
  MatrixXd inv = prod.fullPivLu().inverse();
  
  MatrixXd beta = inv*Xt*y;
  
  double cx = beta(0,0)*0.5;
  double cy = beta(1,0)*0.5;
  double r = sqrt(cx*cx + cy*cy - beta(2,0));
  
  double phi = atan2(cy, cx);
  double d = sqrt(cx*cx + cy*cy) - r;
  double k = 1./r;
  
  MatrixXd diff = y - (X*beta);
  MatrixXd chi2 = (diff.transpose())*diff;
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    chi2_hit[i] += diff(i,0)*diff(i,0);
  }
  
  double dx = d*cos(phi);
  double dy = d*sin(phi);
  
  MatrixXd y2 = MatrixXd::Zero(track.hits.size(), 1);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    y2(i,0) = track.hits[i].get_z();
    if((using_vertex==true ) && (i == (track.hits.size() - 1))){y2(i, 0) /= vertex_sigma_z;}
    else{y2(i, 0) /= layer_z_resolution[track.hits[i].layer];}
  }
  
  MatrixXd X2 = MatrixXd::Zero(track.hits.size(), 2);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    double D = sqrt( pow(dx - track.hits[i].get_x(), 2) + pow(dy - track.hits[i].get_y(),2));
    double s = 0.0;
    
    if(0.5*k*D > 0.1)
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
    
    X2(i,0) = s;  
    X2(i,1) = 1.0;
    
    if((using_vertex==true ) && (i == (track.hits.size() - 1)))
    {
      X2(i, 0) /= vertex_sigma_z;
      X2(i, 1) /= vertex_sigma_z;
    }
    else
    {
      X2(i, 0) /= layer_z_resolution[track.hits[i].layer];
      X2(i, 1) /= layer_z_resolution[track.hits[i].layer];
    }
  }
  
  MatrixXd Xt2 = X2.transpose();
  MatrixXd prod2 = Xt2*X2;
  MatrixXd inv2 = prod2.fullPivLu().inverse();
  MatrixXd beta2 = inv2*Xt2*y2;
  
  MatrixXd diff2 = y2 - (X2*beta2);
  MatrixXd chi2_z = (diff2.transpose())*diff2;
  for(unsigned int i=0;i<track.hits.size();i++) 
  {
    chi2_hit[i] += diff2(i,0)*diff2(i,0);
  }
  
  double z0 = beta2(1,0);
  double dzdl = beta2(0,0)/sqrt(1. + beta2(0,0)*beta2(0,0));
  
  track.phi = phi;
  track.d = d;
  track.kappa = k;
  track.dzdl = dzdl;
  track.z0 = z0;
//   cout<<"fit params = "<<k<<" "<<phi<<" "<<d<<" "<<" "<<z0<<" "<<dzdl<<endl;
  
  if(track.kappa!=0.)
  {
    r=1./track.kappa;
  }
  else
  {
    r=1.0e10;
  }
  
  cx = (track.d+r)*cos(track.phi);
  cy = (track.d+r)*sin(track.phi);
  
  double chi2_tot = 0.;
  for(unsigned int h=0;h<track.hits.size();h++)
  {
    double dx1 = track.hits[h].get_x() - cx;
    double dy1 = track.hits[h].get_y() - cy;
    
    double dx2 = track.hits[h].get_x() + cx;
    double dy2 = track.hits[h].get_y() + cy;
    
    double xydiff1 = sqrt(dx1*dx1 + dy1*dy1) - r;
    double xydiff2 = sqrt(dx2*dx2 + dy2*dy2) - r;
    double xydiff = xydiff2;
    if(fabs(xydiff1) < fabs(xydiff2)){ xydiff = xydiff1; }
    
    double ls_xy = layer_xy_resolution[track.hits[h].layer];
    double ls_z = layer_z_resolution[track.hits[h].layer];
    if((using_vertex == true) && (h == (track.hits.size() - 1)))
    {
      ls_xy = vertex_sigma_xy;
      ls_z = vertex_sigma_z;
    }
    
    chi2_hit[h] = 0.;
    chi2_hit[h] += xydiff*xydiff/(ls_xy*ls_xy);
    chi2_hit[h] += diff2(h,0)*diff2(h,0);
//     chi2_hit[h] += diff2(h,0)*diff2(h,0)/(ls_z*ls_z);
    
    chi2_tot += chi2_hit[h];
  }
  
  unsigned int deg_of_freedom = 2*track.hits.size() - 5;
  
  if(using_vertex == true)
  {
    track.hits.pop_back();
    chi2_hit.pop_back();
  }
  
//   for(unsigned int i=0;i<track.hits.size();++i)
//   {
//     cout<<"chi^2 ( "<<track.hits[i].layer<<" ) = "<<chi2_hit[i]<<endl;
//   }
//   cout<<endl;
  
  return (chi2_tot)/((double)(deg_of_freedom));
}

