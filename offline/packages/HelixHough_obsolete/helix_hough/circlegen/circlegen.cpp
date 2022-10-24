// a. dion - author of code
// j. nagle - added more documentation (07/29/2011)
//
// routine generates a ROOT output with M hits in cylindrical detector layers
// (with radii and number of layers defined below in "vector<double> radii")
// from N (as user input) tracks following a helical pattern from a uniform solenoidal 
// magnetic field.
//
// now with autogen and make, puts executable in install/bin
//
// recent addition for smearing of position resolution of hits (as test)
// see boolean "resolution_smear" below
//
// j.nagle - added ability to include noise hits
//

#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include <sstream>
#include <iostream>
#include <vector>

using namespace std;

bool intersect_circles(bool hel, double startx, double starty, double rad_det, double rad_trk, double cx, double cy, double& x, double& y)
{
  double d2 = (cx*cx + cy*cy);
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
  double P2x = cx*h;
  double P2y = cy*h;
  h = sqrt(R2 - a*a);
  
  double ux = -cy*d_inv;
  double uy = cx*d_inv;
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


void print_usage()
{
  cout<<"usage: circle_gen <# of events> <max # of tracks per event> <output file name>"<<endl;
}


int main(int argc, char** argv)
{
  stringstream ss;
  unsigned int nmaxtracks=0;
  unsigned int nevents=0;

  bool all_same_track_number = true;

  // inclusion of noise hits by these fake tracks
  bool faketracks  = true;
  int  nfaketracks = 250;

  bool resolution_smear = true;
  double smear_x =  (70.0e-4/sqrt(12.))/sqrt(2.0);  //  70 microns as example
  double smear_y =  (70.0e-4/sqrt(12.))/sqrt(2.0);  //  70 microns as example
  double smear_z = 425.0e-4/sqrt(12.);  // 425 microns as example
  
  

  // how many layers to include....  
  //===============================
//    sPHENIX style with 6 layers
  int nlayers = 6;
  vector<double> radii;
  radii.assign(nlayers,0.);
  radii[0]=2.5;
  radii[1]=5.0;
  radii[2]=10.0;
  radii[3]=14.0;
  radii[4]=40.0;
  radii[5]=60.0;
  
  vector<double> smear_x_layer;smear_x_layer.assign(nlayers,0);
  vector<double> smear_y_layer;smear_y_layer.assign(nlayers,0);
  vector<double> smear_z_layer;smear_z_layer.assign(nlayers,0);
  smear_x_layer[0] = (50.0e-4/sqrt(12.))/sqrt(2.0);
  smear_y_layer[0] = (50.0e-4/sqrt(12.))/sqrt(2.0);
  smear_z_layer[0] = (425.0e-4/sqrt(12.));
  smear_x_layer[1] = (50.0e-4/sqrt(12.))/sqrt(2.0);
  smear_y_layer[1] = (50.0e-4/sqrt(12.))/sqrt(2.0);
  smear_z_layer[1] = (425.0e-4/sqrt(12.));
  smear_x_layer[2] = (80.0e-4/sqrt(12.))/sqrt(2.0);
  smear_y_layer[2] = (80.0e-4/sqrt(12.))/sqrt(2.0);
  smear_z_layer[2] = (1000.0e-4/sqrt(12.));
  smear_x_layer[3] = (80.0e-4/sqrt(12.))/sqrt(2.0);
  smear_y_layer[3] = (80.0e-4/sqrt(12.))/sqrt(2.0);
  smear_z_layer[3] = (1000.0e-4/sqrt(12.));
  smear_x_layer[4] = (100.0e-4/sqrt(12.))/sqrt(2.0);
  smear_y_layer[4] = (100.0e-4/sqrt(12.))/sqrt(2.0);
  smear_z_layer[4] = (10000.0e-4/sqrt(12.));
  smear_x_layer[5] = (100.0e-4/sqrt(12.))/sqrt(2.0);
  smear_y_layer[5] = (100.0e-4/sqrt(12.))/sqrt(2.0);
  smear_z_layer[5] = (10000.0e-4/sqrt(12.));
  

  // current SVX
//   int nlayers = 4;
//   vector<double> radii;
//   radii.assign(nlayers,0.);
//   radii[0]=2.5;
//   radii[1]=5.0;
//   radii[2]=10.0;
//   radii[3]=14.0;
  cout << "Circle Gen Code : running with layers = " << nlayers << " and resolution smearing = " << resolution_smear << endl;

  double kappa_min=0.001;
  double kappa_max=0.03;
  double d_min=-0.2;
  double d_max=0.2;
  double phi_min=0.;
  double phi_max=6.28318530717958623;
  double dzdl_min = -0.5;
  double dzdl_max = 0.5;
  double z0_min=-5.;
  double z0_max=5.;
    
  ss.clear();ss.str("");
  ss<<argv[1];
  ss>>nevents;
  if(nevents==0){print_usage();return -1;}
  ss.clear();ss.str("");
  ss<<argv[2];
  ss>>nmaxtracks;
  if(nmaxtracks==0){print_usage();return -1;}
  ss.clear();ss.str("");
  ss<<argv[3];
  if(ss.str() == ""){print_usage();return -1;}
  string filename = ss.str();
  
  unsigned int track=0;
  unsigned int nhits=0;
  double x_hits[12];
  double y_hits[12];
  double z_hits[12];
  int layer[12];
  double trk_kappa=0.;
  double trk_d=0.;
  double trk_phi=0.;
  double trk_dzdl=0.;
  double trk_z0=0.;
  unsigned char charge=0;
  TFile* ofile = new TFile(filename.c_str(), "recreate");
  TTree* etree = new TTree("events", "event list");
  TTree* otree = new TTree("tracks", "tracks in a single event");
  otree->Branch("track", &track, "track/i");
  otree->Branch("nhits", &nhits, "nhits/i");
  otree->Branch("x_hits", x_hits, "x_hits[nhits]/D");
  otree->Branch("y_hits", y_hits, "y_hits[nhits]/D");
  otree->Branch("z_hits", z_hits, "z_hits[nhits]/D");
  otree->Branch("layer", layer, "layer[nhits]/i");
  otree->Branch("trk_kappa", &trk_kappa, "trk_kappa/D");
  otree->Branch("trk_d", &trk_d, "trk_d/D");
  otree->Branch("trk_phi", &trk_phi, "trk_phi/D");
  otree->Branch("trk_dzdl", &trk_dzdl, "trk_dzdl/D");
  otree->Branch("trk_z0", &trk_z0, "trk_z0/D");
  otree->Branch("charge", &charge, "charge/b");
  etree->Branch("tracklist", &otree);
  
  double secondary_proportion=0.1;
  
  
  TRandom3 rand;
  double kappa=0.;
  double d=0.;
  double phi=0.;
  double dzdl=0.;
  double z0=0.;
  
  for(unsigned int ev=0;ev<nevents;ev++)
  {
    otree->Reset();
    unsigned int ntracks = (unsigned int)(rand.Uniform(1., nmaxtracks));
    if (all_same_track_number) ntracks = (unsigned int) nmaxtracks;

    for(unsigned int i=0;i<ntracks;i++)
    {
      kappa = rand.Uniform(kappa_min, kappa_max);
      d=0.;
      phi = rand.Uniform(phi_min, phi_max);
      dzdl = rand.Uniform(dzdl_min, dzdl_max);
      z0 = 0.;
      bool hel=false;
      double temp = rand.Uniform(-1., 1.);
      if(temp > 0.){hel = true;}
      
      double ux=0.;
      double uy=0.;
      double tanphi=0.;
      if((fabs(phi) < 7.85398163397448279e-01) || (fabs(phi - 3.1415926535897931) < 7.85398163397448279e-01))
      {
        tanphi = tan(phi);
        ux = sqrt(1./(1. + tanphi*tanphi));
        if(phi > 1.57079632679489656 && phi < 4.71238898038468967){ux=-ux;}
        uy = ux*tanphi;
      }
      else
      {
        tanphi = tan(1.57079632679489656 - phi);//cotangent
        uy = sqrt(1./(tanphi*tanphi + 1.));
        if(phi > 3.1415926535897931){uy=-uy;}
        ux = uy*tanphi;
      }
      double dx = ux*d;
      double dy = uy*d;
      double dispx = ux/kappa;
      double dispy = uy/kappa;
      double cx = dx + dispx;
      double cy = dy + dispy;
      
      double x=0.;
      double y=0.;
      double z=0.;
      
      double startx=dx;
      double starty=dy;
      double startz=z0;
      nhits=0;
      for(unsigned int l=0;l<radii.size();l++)
      {
        if(intersect_circles(hel, startx, starty, radii[l], 1./kappa, cx, cy, x, y) == true)
        {
          layer[nhits]=l;
          // could smear these by some resolution factor here...
	  x_hits[nhits]=x;
          y_hits[nhits]=y;
	  if (resolution_smear) {
	    x_hits[nhits] += rand.Gaus(0.0,smear_x_layer[l]); // in cm units
	    y_hits[nhits] += rand.Gaus(0.0,smear_y_layer[l]); // in cm units
	  }
          
          double k = kappa;
          double D = sqrt((startx-x)*(startx-x) + (starty-y)*(starty-y));
          
          double s=0.;
          if(0.5*k*D > 0.1)
          {
            double v = 0.5*k*D;
            if(v >= 0.999999){v=0.999999;}
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
          if(dzdl>0.){z = startz + dz;}
          else{z = startz - dz;}
          
          z_hits[nhits]=z;
	  if (resolution_smear) {
	    z_hits[nhits] += rand.Gaus(0.0,smear_z_layer[l]); // in cm units
	  }
          
          
          if((x_hits[nhits] == x_hits[nhits]) && (y_hits[nhits] == y_hits[nhits]) && (z_hits[nhits] == z_hits[nhits]))
          {
            nhits++;
            startx = x;
            starty = y;
            startz = z;
          }
          else
          {
            break;
          }
        }
      }

      if(nhits>0)
      {
        double v1x = x_hits[0];
        double v1y = y_hits[0];
        double v2x = x_hits[1];
        double v2y = y_hits[1];
        
        double diffx = v2x-v1x;
        double diffy = v2y-v1y;
        
        double cross = diffx*cy - diffy*cx;
        
        int cross_sign = 1;
        if(cross < 0.){cross_sign=-1;}
//         cout<<hel<<" "<<cross_sign<<endl;
        
        
        trk_kappa=kappa;
        trk_d=d;
        trk_phi=phi;
        trk_dzdl=dzdl;
        trk_z0=z0;
        charge=0;
	if(hel==true){charge=1;}
        track=i;
//         for(int nnh=0;nnh<nhits;++nnh)
//         {
//           cout<<x_hits[nnh]<<" "<<y_hits[nnh]<<" "<<z_hits[nnh]<<endl;
//         }
        otree->Fill();
      }

    } // end loop over tracks
    

    // how about noise tracks (i.e. made up of random hits from nlayers different tracks (so same spatial distribution)
    // j.nagle - add switch for this features....
    if (faketracks) {

    for(unsigned int j=0;j<nfaketracks;j++) {

      

      for(unsigned int i=0;i<nlayers;i++)
	{
          nhits=0;
	  kappa = rand.Uniform(kappa_min, kappa_max);
	  d=0.;
	  phi = rand.Uniform(phi_min, phi_max);
	  dzdl = rand.Uniform(dzdl_min, dzdl_max);
	  z0 = 0.;
	  bool hel=false;
	  double temp = rand.Uniform(-1., 1.);
	  if(temp > 0.){hel = true;}
	  
	  double ux=0.;
	  double uy=0.;
	  double tanphi=0.;
	  if((fabs(phi) < 7.85398163397448279e-01) || (fabs(phi - 3.1415926535897931) < 7.85398163397448279e-01))
	    {
	      tanphi = tan(phi);
	      ux = sqrt(1./(1. + tanphi*tanphi));
	      if(phi > 1.57079632679489656 && phi < 4.71238898038468967){ux=-ux;}
	      uy = ux*tanphi;
	    }
	  else
	    {
	      tanphi = tan(1.57079632679489656 - phi);//cotangent
	      uy = sqrt(1./(tanphi*tanphi + 1.));
	      if(phi > 3.1415926535897931){uy=-uy;}
	      ux = uy*tanphi;
	    }
	  double dx = ux*d;
	  double dy = uy*d;
	  double dispx = ux/kappa;
	  double dispy = uy/kappa;
	  double cx = dx + dispx;
	  double cy = dy + dispy;
	  
	  double x=0.;
	  double y=0.;
	  double z=0.;
	  
	  double startx=dx;
	  double starty=dy;
	  double startz=z0;

	  // only check one layer for each of these tracks	  
	  for(unsigned int l=i;l<i+1;l++)
	    {
	      if(intersect_circles(hel, startx, starty, radii[l], 1./kappa, cx, cy, x, y) == true)
		{
		  layer[nhits]=l;
		  // could smear these by some resolution factor here...
		  x_hits[nhits]=x;
		  y_hits[nhits]=y;
		  if (resolution_smear) {
		    x_hits[nhits] += rand.Gaus(0.0,smear_x); // in cm units
		    y_hits[nhits] += rand.Gaus(0.0,smear_y); // in cm units
		  }
		  
		  double k = kappa;
		  double D = sqrt((startx-x)*(startx-x) + (starty-y)*(starty-y));
		  
		  double s=0.;
		  if(0.5*k*D > 0.1)
		    {
		      s = 2.*asin(0.5*k*D)/k;
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
		  if(dzdl>0.){z = startz + dz;}
		  else{z = startz - dz;}
		  
		  z_hits[nhits]=z;
		  if (resolution_smear) {
		    z_hits[nhits] += rand.Gaus(0.0,smear_z); // in cm units
		  }
		  
		  nhits++;
		  
		  startx = x;
		  starty = y;
		  startz = z;
		}
	    }
	    
	    if(nhits>0)
            {
              trk_kappa=kappa;
              trk_d=d;
              trk_phi=phi;
              trk_dzdl=dzdl;
              trk_z0=z0;
              charge=0;
              //if(hel==true){charge=1;}
              charge=1; // fake track anyway so no real charge
              track=otree->GetEntries();
              otree->Fill();
            }
	    
	    

	} // end loop over tracks leaving hits in each layer

      

    } // end loop over fake tracks
    
    } // end if faketracks

    /*
    ntracks = (unsigned int)((secondary_proportion)*((double)ntracks));

    //now the secondaries
    for(unsigned int i=0;i<ntracks;i++)
    {
      kappa = rand.Uniform(kappa_min, kappa_max);
      d = rand.Uniform(d_min, d_max);
      phi = rand.Uniform(phi_min, phi_max);
      dzdl = rand.Uniform(dzdl_min, dzdl_max);
      z0 = rand.Uniform(z0_min, z0_max);
      bool hel=false;
      double temp = rand.Uniform(-1., 1.);
      if(temp > 0.){hel = true;}
      
      double ux=0.;
      double uy=0.;
      double tanphi=0.;
      if((fabs(phi) < 7.85398163397448279e-01) || (fabs(phi - 3.1415926535897931) < 7.85398163397448279e-01))
      {
        tanphi = tan(phi);
        ux = sqrt(1./(1. + tanphi*tanphi));
        if(phi > 1.57079632679489656 && phi < 4.71238898038468967){ux=-ux;}
        uy = ux*tanphi;
      }
      else
      {
        tanphi = tan(1.57079632679489656 - phi);//cotangent
        uy = sqrt(1./(tanphi*tanphi + 1.));
        if(phi > 3.1415926535897931){uy=-uy;}
        ux = uy*tanphi;
      }
      double dx = ux*d;
      double dy = uy*d;
      double dispx = ux/kappa;
      double dispy = uy/kappa;
      double cx = dx + dispx;
      double cy = dy + dispy;
      
      double x=0.;
      double y=0.;
      double z=0.;
      
      double startx=dx;
      double starty=dy;
      double startz=z0;
      nhits=0;
      for(unsigned int l=0;l<radii.size();l++)
      {
        if(intersect_circles(hel, startx, starty, radii[l], 1./kappa, cx, cy, x, y) == true)
        {
          layer[nhits]=l;
          x_hits[nhits]=x;
          y_hits[nhits]=y;
          
          double k = kappa;
          double D = sqrt((startx-x)*(startx-x) + (starty-y)*(starty-y));
          
          double s=0.;
          if(0.5*k*D > 0.1)
          {
            s = 2.*asin(0.5*k*D)/k;
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
          if(dzdl>0.){z = startz + dz;}
          else{z = startz - dz;}
          
          z_hits[nhits]=z;
          
          nhits++;
          
          startx = x;
          starty = y;
          startz = z;
        }
      }
      if(nhits>0)
      {
        trk_kappa=kappa;
        trk_d=d;
        trk_phi=phi;
        trk_dzdl=dzdl;
        trk_z0=z0;
        charge=0;
        if(hel==true){charge=1;}
        track=i;
        otree->Fill();
      }
    }
    */
    
    
    etree->Fill();
  }
  
  
  etree->Write();
  ofile->Close();
  ofile->Delete();
  
  
  
  return 0;
}


