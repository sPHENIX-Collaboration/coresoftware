#include "ThreeHitSeedGrower.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <Eigen/LU>
#include <Eigen/Core>


using namespace std;
using namespace Eigen;


ThreeHitSeedGrower::ThreeHitSeedGrower(vector<float>& detrad, unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0, HelixResolution& min_resolution, HelixResolution& max_resolution, HelixRange& range) : HelixHough(n_phi, n_d, n_k, n_dzdl, n_z0, min_resolution, max_resolution, range), using_vertex(false), vertex_sigma_xy(0.002), vertex_sigma_z(0.005), chi2_cut(3.)
{
	unsigned int n_layers = detrad.size();
	_max_hits = n_layers + 2;
  for(unsigned int i=0;i<n_layers;++i)
  {
    detector_radii.push_back(detrad[i]);
  }
}


void ThreeHitSeedGrower::findTracks(vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks)
{cout<<"findTracks::entered with "<<hits.size()<<" hits"<<endl;
  vector<double> chi2_hits;
	double chi2;
  SimpleTrack3D temp_track;
  temp_track.hits.resize(3, hits[0]);
  for(unsigned int i1=0;i1<hits.size();++i1)
  {
    temp_track.hits[0] = hits[i1];
    for(unsigned int i2=(i1+1);i2<hits.size();++i2)
			{//TO DO: ADD PHICUT
      if( (hits[i2].layer == hits[i1].layer)){continue;}
      temp_track.hits[1] = hits[i2];
      for(unsigned int i3=(i2+1);i3<hits.size();++i3)
      {
        if((hits[i3].layer == hits[i2].layer) || (hits[i3].layer == hits[i1].layer)){continue;}
        temp_track.hits[2] = hits[i3];
				
				chi2 = fitTrack(temp_track, chi2_hits);
				if(chi2 > chi2_cut){continue;}
				vector<unsigned int> nhit_layer;
				nhit_layer.assign(detector_radii.size(),0);
				for(unsigned int ihit = 0; ihit<3;ihit++)
					{
						nhit_layer[temp_track.hits[ihit].layer]+=1;
					}
				if(GrowTrack(temp_track, nhit_layer, hits, tracks, 0) == false)
					{//refit if track wasn't grown.
						vector<unsigned int> tempcomb;
						tempcomb.assign(3,0);
						tempcomb[0] = temp_track.hits[0].index;
						tempcomb[1] = temp_track.hits[1].index;
						tempcomb[2] = temp_track.hits[2].index;
						sort(tempcomb.begin(), tempcomb.end());
						set<vector<unsigned int> >::iterator it = combos.find(tempcomb);
						if(it != combos.end()){continue;}
						combos.insert(tempcomb);
						chi2 = fitTrack(temp_track, chi2_hits);
						tracks.push_back(temp_track);
						cout<<"findTrack::added a 3 hit track"<<endl;
					}
				
      }
    }
  }
	cout<<"leaving findTrack"<<endl;
}
//this is my attempt at implementation of 3 hit seeding adding hits up to 8 total
//plan for this function:
//will figure out which layers need new hits the most from nhit_layer
//incrementally adds a hit from hits and fits it
//tests to see if the chi squared is still in bound 
//recursively calls grow track on newly accepted track and adds it to tracks.
//
bool ThreeHitSeedGrower::GrowTrack(SimpleTrack3D seed_track, vector<unsigned int>& nhit_layer, vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks, unsigned int c_hit){
	unsigned int init_size = seed_track.hits.size();
	vector<double> chi2_hits;
	double chi2;
	bool hit_added;
	while(c_hit < hits.size())
		{
			hit_added = addOneHit(seed_track, nhit_layer, c_hit, hits);
			c_hit+=1;
			if(hit_added == false)
				continue;
			unsigned int final_size = seed_track.hits.size();
			if(final_size < _max_hits)
				{
					if(GrowTrack(seed_track, nhit_layer, hits, tracks, c_hit) == false)
						//attempt to continue growing track. if fails, reset fitting parameters
						{
							chi2 = fitTrack(seed_track,chi2_hits);
						}
				}
			vector<unsigned int> tempcomb;
			tempcomb.assign(seed_track.hits.size(),0);
			for(unsigned int ihit = 0; ihit < seed_track.hits.size();ihit++)
				tempcomb[ihit] = seed_track.hits[ihit].index;
			sort(tempcomb.begin(), tempcomb.end());
			set<vector<unsigned int> >::iterator it = combos.find(tempcomb);
			if(it != combos.end()){
				seed_track.hits.pop_back();
				continue;
			}        
			combos.insert(tempcomb);
			tracks.push_back(seed_track);		
			cout<<"GrowTrack::added a track of size: "<<seed_track.hits.size()<<endl;
			return true;
		}
	return false;
}

bool ThreeHitSeedGrower::addOneHit(SimpleTrack3D & seed_track, vector<unsigned int> & nhit_layer, unsigned int c_hit, vector<SimpleHit3D>& hits) 
{//cout<<"addOneHit::newly entered"<<endl;
	for(unsigned int ihit =0; ihit<seed_track.hits.size();ihit++)
		{
			//	cout<<"addOneHit::hit: "<<ihit<<" index = "<<seed_track.hits[ihit].index<<endl;
			if(seed_track.hits[ihit].index == hits[c_hit].index)
			return false;
		}
	//	cout<<"addOneHit::hit: "<<c_hit<<" passed with index = "<<hits[c_hit].index<<endl;
  //this conditional logic specifically applies to the VTX, maybe modularize this later.                                   
	if(nhit_layer[hits[c_hit].layer] == 0 || (nhit_layer[hits[c_hit].layer] == 1  && hits[c_hit].layer >1)){
		vector<double> chi2_hits;
		double chi2;
		seed_track.hits.push_back(hits[c_hit]);
		chi2 = fitTrack(seed_track,chi2_hits);
		if(chi2 <= chi2_cut)
			{
				return true;
				nhit_layer[hits[c_hit].layer]+=1;;
			}
		else
			seed_track.hits.pop_back();
	}
	return false;
}

void ThreeHitSeedGrower::finalize(vector<SimpleTrack3D>& input, vector<SimpleTrack3D>& output)
{
  unsigned int nt = input.size();
  for(unsigned int i=0;i<nt;++i)
  {
    output.push_back(input[i]);
  }
}


double ThreeHitSeedGrower::fitTrack(SimpleTrack3D & track, vector<double>& chi2_hit)
{
  if(using_vertex == true)
  {
    track.hits.push_back(SimpleHit3D(0.,0., 0.,0., 0.,0., 0, 0));
  }
  
  chi2_hit.clear();
  chi2_hit.resize(track.hits.size(), 0.);
  
  MatrixXf y = MatrixXf::Zero(track.hits.size(), 1);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    y(i, 0) = ( pow(track.hits[i].x,2) + pow(track.hits[i].y,2) );
    if((using_vertex==true ) && (i == (track.hits.size() - 1))){y(i, 0) /= vertex_sigma_xy;}
    else{y(i, 0) /= layer_xy_resolution[track.hits[i].layer];}
  }
  
  MatrixXf X = MatrixXf::Zero(track.hits.size(), 3);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    X(i, 0) = track.hits[i].x;
    X(i, 1) = track.hits[i].y;
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
  
  MatrixXf Xt = X.transpose();
  
  MatrixXf prod = Xt*X;
  MatrixXf inv = prod.fullPivLu().inverse();
  
  MatrixXf beta = inv*Xt*y;
  
  float cx = beta(0,0)*0.5;
  float cy = beta(1,0)*0.5;
  float r = sqrt(cx*cx + cy*cy - beta(2,0));
  
  float phi = atan2(cy, cx);
  float d = sqrt(cx*cx + cy*cy) - r;
  float k = 1./r;
  
  MatrixXf diff = y - (X*beta);
  MatrixXf chi2 = (diff.transpose())*diff;
  
  float dx = d*cos(phi);
  float dy = d*sin(phi);
  
  MatrixXf y2 = MatrixXf::Zero(track.hits.size(), 1);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    y2(i,0) = track.hits[i].z;
    if((using_vertex==true ) && (i == (track.hits.size() - 1))){y2(i, 0) /= vertex_sigma_z;}
    else{y2(i, 0) /= layer_z_resolution[track.hits[i].layer];}
  }
  
  MatrixXf X2 = MatrixXf::Zero(track.hits.size(), 2);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    float D = sqrt( pow(dx - track.hits[i].x, 2) + pow(dy - track.hits[i].y,2));
    float s = 0.0;
    
    if(0.5*k*D > 0.1)
    {
      float v = 0.5*k*D;
      if(v >= 0.999999){v = 0.999999;}
      s = 2.*asin(v)/k;
    }
    else
    {
      float temp1 = k*D*0.5;temp1*=temp1;
      float temp2 = D*0.5;
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
  
  MatrixXf Xt2 = X2.transpose();
  MatrixXf prod2 = Xt2*X2;
  MatrixXf inv2 = prod2.fullPivLu().inverse();
  MatrixXf beta2 = inv2*Xt2*y2;
  
  MatrixXf diff2 = y2 - (X2*beta2);
  MatrixXf chi2_z = (diff2.transpose())*diff2;
  
  float z0 = beta2(1,0);
  float dzdl = beta2(0,0)/sqrt(1. + beta2(0,0)*beta2(0,0));
  
  track.phi = phi;
  track.d = d;
  track.kappa = k;
  track.dzdl = dzdl;
  track.z0 = z0;
  
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
  
  float chi2_tot = 0.;
  for(unsigned int h=0;h<track.hits.size();h++)
  {
    float dx1 = track.hits[h].x - cx;
    float dy1 = track.hits[h].y - cy;
    
    float dx2 = track.hits[h].x + cx;
    float dy2 = track.hits[h].y + cy;
    
    float xydiff1 = sqrt(dx1*dx1 + dy1*dy1) - r;
    float xydiff2 = sqrt(dx2*dx2 + dy2*dy2) - r;
    float xydiff = xydiff2;
    if(fabs(xydiff1) < fabs(xydiff2)){ xydiff = xydiff1; }
    
    float ls_xy = layer_xy_resolution[track.hits[h].layer];
    if((using_vertex == true) && (h == (track.hits.size() - 1)))
    {
      ls_xy = vertex_sigma_xy;
    }
    
    chi2_hit[h] = 0.;
    chi2_hit[h] += xydiff*xydiff/(ls_xy*ls_xy);
    chi2_hit[h] += diff2(h,0)*diff2(h,0);
    
    chi2_tot += chi2_hit[h];
  }
  
  unsigned int deg_of_freedom = 2*track.hits.size() - 5;
  
  if(using_vertex == true)
  {
    track.hits.pop_back();
    chi2_hit.pop_back();
  }
    
  return (chi2_tot)/((double)(deg_of_freedom));
}


