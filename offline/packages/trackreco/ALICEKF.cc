#include "ALICEKF.h"

#include "GPUTPCTrackLinearisation.h"
#include "GPUTPCTrackParam.h"

#include <trackbase/TrkrCluster.h>
#include <Geant4/G4SystemOfUnits.hh>

#include "TFile.h"
#include "TNtuple.h"

//#define _DEBUG_

#if defined(_DEBUG_)
#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void)0
#endif

#define LogError(exp) if(Verbosity()>0) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) if(Verbosity()>0) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

using namespace std;
typedef vector<TrkrDefs::cluskey> keylist;

bool ALICEKF::checknan(double val, std::string name, int num)
{
  if(std::isnan(val))
  {
    if(Verbosity()>0) std::cout << "WARNING: " << name << " is NaN for seed " << num << ". Aborting this seed.\n";
  }
  return std::isnan(val);
}

double ALICEKF::get_Bz(double x, double y, double z)
{
  if(_use_const_field) return 1.4;
  double p[4] = {x*cm,y*cm,z*cm,0.*cm};
  double bfield[3];
  _B->GetFieldValue(p,bfield);
  return bfield[2]/tesla;
}

vector<SvtxTrack_v1> ALICEKF::ALICEKalmanFilter(vector<keylist> trackSeedKeyLists,bool use_nhits_limit)
{
//  TFile* f = new TFile("/sphenix/u/mjpeters/macros_hybrid/detectors/sPHENIX/pull.root", "RECREATE");
//  TNtuple* ntp = new TNtuple("pull","pull","cx:cy:cz:xerr:yerr:zerr:tx:ty:tz:layer:xsize:ysize:phisize:phierr:zsize");
  vector<SvtxTrack_v1> seeds_vector;
  int nseeds = 0;
  if(Verbosity()>0) std::cout << "min clusters per track: " << _min_clusters_per_track << "\n";
  for(vector<keylist>::iterator trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    if(trackKeyChain->size()<2) continue;
    if(use_nhits_limit && trackKeyChain->size() < _min_clusters_per_track) continue;
    // get starting cluster from key
    TrkrCluster* startCluster = _cluster_map->findCluster(trackKeyChain->at(0));
    // Transform sPHENIX coordinates into ALICE-compatible coordinates
    double x0 = startCluster->getPosition(0);
    double y0 = startCluster->getPosition(1);
    double z0 = startCluster->getPosition(2);
    LogDebug("Initial (x,y,z): (" << x0 << "," << y0 << "," << z0 << ")" << endl);
    // ALICE x coordinate = distance from beampipe
    double alice_x0 = sqrt(x0*x0+y0*y0);
    double alice_y0 = 0;
    double alice_z0 = z0;
    // Initialize track and linearisation
    GPUTPCTrackParam trackSeed;
    trackSeed.InitParam();
    trackSeed.SetX(alice_x0);
    trackSeed.SetY(alice_y0);
    trackSeed.SetZ(alice_z0);
    double x = x0;
    double y = y0;
    #if defined(_DEBUG_)
    double z = z0;
    double alice_x = sqrt(x0*x0+y0*y0);
    #endif
    double trackCartesian_x = 0.;
    double trackCartesian_y = 0.;
    double trackCartesian_z = 0.;
    // Pre-set momentum-based parameters to improve numerical stability
    TrkrCluster* SecondCluster = _cluster_map->findCluster(trackKeyChain->at(1));
    double second_x = SecondCluster->getPosition(0);
    double second_y = SecondCluster->getPosition(1);
    double second_z = SecondCluster->getPosition(2);
    double first_phi = atan2(y0,x0);
    double second_alice_x = second_x*cos(first_phi)+second_y*sin(first_phi);
    double delta_alice_x = second_alice_x - alice_x0;
    //double second_alice_y = (second_x/cos(first_phi)-second_y/sin(first_phi))/(sin(first_phi)/cos(first_phi)+cos(first_phi)/sin(first_phi));
    double second_alice_y = -second_x*sin(first_phi)+second_y*cos(first_phi);
    double init_SinPhi = second_alice_y / sqrt(delta_alice_x*delta_alice_x + second_alice_y*second_alice_y);
    double delta_z = second_z - z0;
    double init_DzDs = delta_z / sqrt(delta_alice_x*delta_alice_x + second_alice_y*second_alice_y);
    trackSeed.SetSinPhi(init_SinPhi);
    LogDebug("Set initial SinPhi to " << init_SinPhi << endl);
    trackSeed.SetDzDs(init_DzDs);
    LogDebug("Set initial DzDs to " << init_DzDs << endl);
    
    // get initial pt estimate
    std::vector<std::pair<double,double>> pts;
    for(int c=0;c<trackKeyChain->size();++c)
    {
      TrkrCluster* cl = _cluster_map->findCluster(trackKeyChain->at(c));
      pts.push_back(std::make_pair(cl->getX(),cl->getY()));
    }
    double R;
    double x_center;
    double y_center;
    CircleFitByTaubin(pts,R,x_center,y_center);
    if(Verbosity()>1) cout << "circle fit parameters: R=" << R << ", X0=" << x_center << ", Y0=" << y_center << endl;
    double init_QPt = 1./(0.3*R/100.*get_Bz(x0,y0,z0));
    // determine charge
    double phi_first = atan2(y0,x0);
    if(Verbosity()>1) cout << "phi_first: " << phi_first << endl;
    double phi_second = atan2(second_y,second_x);
    if(Verbosity()>1) cout << "phi_second: " << phi_second << endl;
    double dphi = phi_second - phi_first;
    if(Verbosity()>1) cout << "dphi: " << dphi << endl;
    if(dphi>M_PI) dphi = 2*M_PI - dphi;
    if(Verbosity()>1) cout << "corrected dphi: " << dphi << endl;
    if(dphi<0) init_QPt *= -1;
    if(Verbosity()>0) cout << "initial QPt: " << init_QPt << endl;
    trackSeed.SetQPt(init_QPt);

    GPUTPCTrackLinearisation trackLine(trackSeed);

    LogDebug(endl << endl << "------------------------" << endl << "seed size: " << trackKeyChain->size() << endl << endl << endl);
    int cluster_ctr = 1;
    bool aborted = false;
    // starting at second cluster, perform track propagation
    std::vector<double> cx;
    std::vector<double> cy;
    std::vector<double> cz;
    std::vector<double> tx;
    std::vector<double> ty;
    std::vector<double> tz;
    std::vector<double> xerr;
    std::vector<double> yerr;
    std::vector<double> zerr;
    std::vector<double> layer;
    std::vector<double> xsize;
    std::vector<double> ysize;
    std::vector<double> phisize;
    std::vector<double> phierr;
    std::vector<double> zsize;
    for(keylist::iterator clusterkey = next(trackKeyChain->begin()); clusterkey != trackKeyChain->end(); ++clusterkey)
    {
      LogDebug("cluster " << cluster_ctr << " -> " << cluster_ctr + 1 << endl);
      LogDebug("this cluster (x,y,z) = (" << x << "," << y << "," << z << ")" << endl);
      LogDebug("layer " << (int)TrkrDefs::getLayer(*clusterkey) << endl);
      // get cluster from key
      TrkrCluster* nextCluster = _cluster_map->findCluster(*clusterkey);
      // find ALICE x-coordinate
      double nextCluster_x = nextCluster->getPosition(0);
      double nextCluster_xerr = sqrt(nextCluster->getError(0,0));
      double nextCluster_y = nextCluster->getPosition(1);
      double nextCluster_yerr = sqrt(nextCluster->getError(1,1));
      double nextCluster_z = nextCluster->getPosition(2);
      double nextCluster_zerr = sqrt(nextCluster->getError(2,2));
      // rotate track coordinates to match orientation of next cluster
      double newPhi = atan2(nextCluster_y,nextCluster_x);
      LogDebug("new phi = " << newPhi << endl);
      double oldPhi = atan2(y,x);
      LogDebug("old phi = " << oldPhi << endl);
      double alpha = newPhi - oldPhi;
      LogDebug("alpha = " << alpha << endl);
      if(!trackSeed.Rotate(alpha,trackLine,_max_sin_phi))
      {
        LogWarning("Rotate failed! Aborting for this seed...\n");
        aborted = true;
        break;
      }
      double nextAlice_x = nextCluster_x*cos(newPhi)+nextCluster_y*sin(newPhi);
      LogDebug("track coordinates (ALICE) after rotation: (" << trackSeed.GetX() << "," << trackSeed.GetY() << "," << trackSeed.GetZ() << ")" << endl);
      LogDebug("Transporting from " << alice_x << " to " << nextAlice_x << "...");
      GPUTPCTrackParam::GPUTPCTrackFitParam fp;
      trackSeed.CalculateFitParameters(fp);
  //    for(int i=1;i<=10;i++)
  //    {
        double track_x = trackSeed.GetX()*cos(newPhi)-trackSeed.GetY()*sin(newPhi);
        double track_y = trackSeed.GetX()*sin(newPhi)+trackSeed.GetY()*cos(newPhi);
        double track_z = trackSeed.GetZ();
        if(!trackSeed.TransportToX(nextAlice_x,_Bzconst*get_Bz(track_x,track_y,track_z),_max_sin_phi)) // remember: trackLine was here
        {
          LogWarning("Transport failed! Aborting for this seed...\n");
          aborted = true;
          break;
        }
  //    }
      // convert ALICE coordinates to sPHENIX cartesian coordinates, for debugging
      double predicted_alice_x = trackSeed.GetX();
      LogDebug("new track ALICE x = " << trackSeed.GetX() << endl);
      double predicted_alice_y = trackSeed.GetY();
      LogDebug("new track ALICE y = " << trackSeed.GetY() << endl);
      double predicted_z = trackSeed.GetZ();
      LogDebug("new track z = " << trackSeed.GetZ() << endl);
      double cos_phi = x/sqrt(x*x+y*y);
      LogDebug("cos_phi = " << cos_phi << endl);
      double sin_phi = y/sqrt(x*x+y*y);
      LogDebug("sin phi = " << sin_phi << endl);
      trackCartesian_x = predicted_alice_x*cos_phi-predicted_alice_y*sin_phi;
      trackCartesian_y = predicted_alice_x*sin_phi+predicted_alice_y*cos_phi;
      trackCartesian_z = predicted_z;
      LogDebug("Track transported to (x,y,z) = (" << trackCartesian_x << "," << trackCartesian_y << "," << trackCartesian_z << ")" << endl);
      LogDebug("Next cluster is at (x,y,z) = (" << nextCluster_x << "," << nextCluster_y << "," << nextCluster_z << ")" << endl);
      LogDebug("Cluster errors: (" << nextCluster_xerr << ", " << nextCluster_yerr << ", " << nextCluster_zerr << ")" << endl);
      LogDebug("track coordinates (ALICE) after rotation: (" << trackSeed.GetX() << "," << trackSeed.GetY() << "," << trackSeed.GetZ() << ")" << endl);
      //double nextCluster_alice_y = (nextCluster_x/cos(newPhi) - nextCluster_y/sin(newPhi))/(tan(newPhi)+1./tan(newPhi));
      //double nextCluster_alice_y = 0.;
      double nextCluster_alice_y = -nextCluster_x*sin(newPhi)+nextCluster_y*cos(newPhi);
      LogDebug("next cluster ALICE y = " << nextCluster_alice_y << endl);
      double y2_error = nextCluster->getError(0,0)*sin(newPhi)*sin(newPhi)+2*nextCluster->getError(0,1)*cos(newPhi)*sin(newPhi)+nextCluster->getError(1,1)*cos(newPhi)*cos(newPhi);
      double z2_error = nextCluster_zerr*nextCluster_zerr;
      LogDebug("track ALICE SinPhi = " << trackSeed.GetSinPhi() << endl);
      // Apply Kalman filter
      if(!trackSeed.Filter(nextCluster_alice_y,nextCluster_z,y2_error,z2_error,_max_sin_phi))
      {
	LogError("Kalman filter failed for seed " << nseeds << "! Aborting for this seed..." << endl);
        aborted = true;
        break;
      }
      #if defined(_DEBUG_)
      double track_pt = 1./trackSeed.GetQPt();
      double track_pY = track_pt*trackSeed.GetSinPhi();
      double track_pX = sqrt(track_pt*track_pt-track_pY*track_pY);
      double track_px = track_pX*cos(newPhi)-track_pY*sin(newPhi);
      double track_py = track_pX*sin(newPhi)+track_pY*cos(newPhi);
      double track_pz = track_pt*trackSeed.GetDzDs();
      double track_pterr = sqrt(trackSeed.GetErr2QPt())/(trackSeed.GetQPt()*trackSeed.GetQPt());
      #endif
      LogDebug("track pt = " << track_pt << " +- " << track_pterr << endl);
      LogDebug("track ALICE p = (" << track_pX << ", " << track_pY << ", " << track_pz << ")" << endl);
      LogDebug("track p = (" << track_px << ", " << track_py << ", " << track_pz << ")" << endl);
      x = nextCluster_x;
      y = nextCluster_y;
      #if defined(_DEBUG_)
      z = nextCluster_z;
      alice_x = nextAlice_x;
      #endif
      ++cluster_ctr;
      //if(cluster_ctr>10)
      {
        cx.push_back(nextCluster_x);
        cy.push_back(nextCluster_y);
        cz.push_back(nextCluster_z);
        tx.push_back(trackCartesian_x);
        ty.push_back(trackCartesian_y);
        tz.push_back(trackCartesian_z);
        xerr.push_back(nextCluster_xerr);
        yerr.push_back(nextCluster_yerr);
        zerr.push_back(nextCluster_zerr);
        layer.push_back(TrkrDefs::getLayer(*clusterkey));
        xsize.push_back(sqrt(nextCluster->getSize(0,0)));
        ysize.push_back(sqrt(nextCluster->getSize(1,1)));
        phisize.push_back(nextCluster->getPhiSize());
        phierr.push_back(nextCluster->getPhiError());
        zsize.push_back(nextCluster->getZSize());
      }
    }
    if(aborted) continue;
    if(Verbosity()>0) cout << "finished track\n";
/*
    // transport to beamline
//    float old_phi = atan2(y,x);
    float trackX = trackSeed.GetX();
    for(int i=99;i>=0;i--)
    {
      if(!trackSeed.TransportToX(i/100.*trackX,trackLine,_Bz,_max_sin_phi))
      {
        LogWarning("Transport failed! Aborting for this seed...\n");
        aborted = true;
        break;
      }
//      float new_phi = atan2(trackSeed.GetX()*sin(old_phi)+trackSeed.GetY()*cos(old_phi),trackSeed.GetX()*cos(old_phi)-trackSeed.GetY()*sin(old_phi));
//      if(!trackSeed.Rotate(new_phi-old_phi,trackLine,_max_sin_phi))
//      {
//        LogWarning("Rotate failed! Aborting for this seed...\n");
//        aborted = true;
//        break;
//      }
//      old_phi = new_phi;
    }
    if(aborted) continue;
    cout << "transported to beamline\n";
    // find nearest vertex
    double beamline_X = trackSeed.GetX();
    double beamline_Y = trackSeed.GetY();
*/
    double track_phi = atan2(y,x);
/*
    double beamline_x = beamline_X*cos(track_phi)-beamline_Y*sin(track_phi);
    double beamline_y = beamline_X*sin(track_phi)+beamline_Y*cos(track_phi);
    double beamline_z = trackSeed.GetZ();
    double min_dist = 1e9;
    int best_vtx = -1;
    for(int i=0;i<_vertex_x.size();++i)
    {
      double delta_x = beamline_x-_vertex_x[i];
      double delta_y = beamline_y-_vertex_y[i];
      double delta_z = beamline_z-_vertex_z[i];
      double dist = sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
      if(dist<min_dist)
      {
        min_dist = dist;
        best_vtx = i;
      }
    }
    cout << "best vtx:\n";
    cout << "("<<_vertex_x[best_vtx]<<","<<_vertex_y[best_vtx]<<","<<_vertex_z[best_vtx]<<")\n";
    // Fit to vertex point
    double vertex_phi = atan2(_vertex_y[best_vtx],_vertex_x[best_vtx]);
    cout << "vertex_phi: " << vertex_phi << "\n";
    cout << "track_phi: " << track_phi << "\n";
//    double alpha = vertex_phi - track_phi;
    // Here's where we need to be careful about the vertex position.
    // Most clusters are at roughly the same spatial phi, with only a little rotation required between them.
    // This is no longer guaranteed for the vertex - its phi could be anywhere,
    // including on the opposite side of the origin.
    // If it ends up on the opposite side, then we need to transport to *negative* radius in order to get close to it.
    // We will simplify this condition to abs(alpha)>pi/2, which assumes that (innermost TPC cluster R) >> (vertex R).

    bool crosses_origin = false;
    if(alpha<-M_PI/4)
    {
      while(alpha<-M_PI/4) alpha += M_PI/2;
      crosses_origin = true;
    }
    if(alpha>M_PI/4)
    {
      while(alpha>M_PI/4) alpha -= M_PI/2;
      crosses_origin = true;
    }
    if(crosses_origin) cout << "bad\n";
    cout << "alpha: " << alpha << "\n";

    if(!trackSeed.Rotate(alpha,trackLine,_max_sin_phi))
    {
      LogWarning("Rotate failed! Aborting for this seed...\n");
      aborted = true;
      continue;
    }
    LogDebug("ALICE coordinates after rotation: (" << trackSeed.GetX() << ", " << trackSeed.GetY() << ", " << trackSeed.GetZ() << ")\n");
    cout << "rotated to vertex\n";

    double vertex_X = sqrt(_vertex_x[best_vtx]*_vertex_x[best_vtx]+_vertex_y[best_vtx]*_vertex_y[best_vtx]);
    if(crosses_origin) vertex_X = -vertex_X;
    if(!trackSeed.TransportToX(vertex_X,trackLine,_Bz,_max_sin_phi))
    {
      LogWarning("Transport failed! Aborting for this seed...\n");
      aborted = true;
      continue;
    }
    LogDebug("Track transported to (x,y,z) = (" << trackSeed.GetX()*cos(vertex_phi)-trackSeed.GetY()*sin(vertex_phi) << "," << trackSeed.GetX()*sin(vertex_phi)+trackSeed.GetY()*cos(vertex_phi) << "," << trackSeed.GetZ() << ")" << endl);
    LogDebug("Next cluster is at (x,y,z) = (" << _vertex_x[best_vtx] << "," << _vertex_y[best_vtx] << "," << _vertex_z[best_vtx] << ")" << endl);

    double vertex_Y = -_vertex_x[best_vtx]*sin(vertex_phi)+_vertex_y[best_vtx]*cos(vertex_phi);
    cout << "vertex Y: " << vertex_Y << "\n";
    cout << "transported to vertex\n";
    double vertex_Yerr = -_vertex_xerr[best_vtx]*sin(vertex_phi)+_vertex_yerr[best_vtx]*cos(vertex_phi);
    cout << "vertex Y err: " << vertex_Yerr << "\n";

    if(!trackSeed.Filter(vertex_Y,_vertex_z[best_vtx],vertex_Yerr*vertex_Yerr,_vertex_zerr[best_vtx]*_vertex_zerr[best_vtx],_max_sin_phi))
    {
      cout << "filter failed\n";
      if (Verbosity() >= 1)
        LogError("Kalman filter failed for seed " << nseeds << "! Aborting for this seed..." << endl);
      aborted = true;
      continue;
    }
*/
    double track_pt = 1./trackSeed.GetQPt();
    #if defined(_DEBUG_)
    double track_pY = track_pt*trackSeed.GetSinPhi();
    double track_pX = sqrt(track_pt*track_pt-track_pY*track_pY);
    double track_px = track_pX*cos(track_phi)-track_pY*sin(track_phi);
    double track_py = track_pX*sin(track_phi)+track_pY*cos(track_phi);
    double track_pz = track_pt*trackSeed.GetDzDs();
    #endif
    double track_pterr = sqrt(trackSeed.GetErr2QPt())/(trackSeed.GetQPt()*trackSeed.GetQPt());
    // If Kalman filter doesn't do its job (happens often with short seeds), use the circle-fit estimate as the central value
    if(track_pterr>fabs(track_pt)) track_pt = 1./init_QPt;
    LogDebug("track pt = " << track_pt << " +- " << track_pterr << endl);
    LogDebug("track ALICE p = (" << track_pX << ", " << track_pY << ", " << track_pz << ")" << endl);
    LogDebug("track p = (" << track_px << ", " << track_py << ", " << track_pz << ")" << endl);

/*    
    if(cluster_ctr!=1 && !trackSeed.CheckNumericalQuality())
    {
      cout << "ERROR: Track seed failed numerical quality check before conversion to sPHENIX coordinates! Skipping this one.\n";
      aborted = true;
      continue;
    } 
*/    
    //    pt:z:dz:phi:dphi:c:dc
    // Fill NT with track parameters
    // double StartEta = -log(tan(atan(z0/sqrt(x0*x0+y0*y0))));
    if(aborted) continue;
//    double track_pt = fabs( 1./(trackSeed.GetQPt()));
    if(checknan(track_pt,"pT",nseeds)) continue;
//    double track_pterr = sqrt(trackSeed.GetErr2QPt())/(trackSeed.GetQPt()*trackSeed.GetQPt());
    if(checknan(track_pterr,"pT err",nseeds)) continue;
    LogDebug("Track pterr = " << track_pterr << endl);
    double track_x = trackSeed.GetX()*cos(track_phi)-trackSeed.GetY()*sin(track_phi);
    double track_y = trackSeed.GetX()*sin(track_phi)+trackSeed.GetY()*cos(track_phi);
    double track_z = trackSeed.GetZ();
    if(checknan(track_z,"z",nseeds)) continue;
    double track_zerr = sqrt(trackSeed.GetErr2Z());
    if(checknan(track_zerr,"zerr",nseeds)) continue;
    double last_cluster_phierr = _cluster_map->findCluster(trackKeyChain->back())->getPhiError();
    // phi error assuming error in track radial coordinate is zero
    double track_phierr = sqrt(pow(last_cluster_phierr,2)+(pow(trackSeed.GetX(),2)*trackSeed.GetErr2Y()) / 
      pow(pow(trackSeed.GetX(),2)+pow(trackSeed.GetY(),2),2));
    if(checknan(track_phierr,"phierr",nseeds)) continue;
    LogDebug("Track phi = " << atan2(track_py,track_px) << endl);
    LogDebug("Track phierr = " << track_phierr << endl);
    double track_curvature = trackSeed.GetKappa(_Bzconst*get_Bz(track_x,track_y,track_z));
    if(checknan(track_curvature,"curvature",nseeds)) continue;
    double track_curverr = sqrt(trackSeed.GetErr2QPt())*_Bzconst*get_Bz(track_x,track_y,track_z);
    if(checknan(track_curverr,"curvature error",nseeds)) continue;
    SvtxTrack_v1 track;
    track.set_id(nseeds);
//    track.set_vertex_id(_vertex_ids[best_vtx]);
    for (unsigned int j = 0; j < trackKeyChain->size(); ++j)
    {
      track.insert_cluster_key(trackKeyChain->at(j));
    }
    track.set_chisq(trackSeed.GetChi2());
    track.set_ndf(trackSeed.GetNDF());
    int track_charge = 0;
    if(trackSeed.GetQPt()<0) track_charge = -1 * _fieldDir;
    else track_charge = 1 * _fieldDir;
    track.set_charge(track_charge);
    double s = sin(track_phi);
    double c = cos(track_phi);
    double p = trackSeed.GetSinPhi();
    // TrkrCluster *cl = _cluster_map->findCluster(trackKeyChain->at(0));
    track.set_x(trackSeed.GetX()*c-trackSeed.GetY()*s);//_vertex_x[best_vtx]);  //track.set_x(cl->getX());
    track.set_y(trackSeed.GetX()*s+trackSeed.GetY()*c);//_vertex_y[best_vtx]);  //track.set_y(cl->getY());
    track.set_z(trackSeed.GetZ());//_vertex_z[best_vtx]);  //track.set_z(cl->getZ());
    if(Verbosity()>0) cout << "x " << track.get_x() << "\n";
    if(Verbosity()>0) cout << "y " << track.get_y() << "\n";
    if(Verbosity()>0) cout << "z " << track.get_z() << "\n";
    if(checknan(p,"ALICE sinPhi",nseeds)) continue;
    double d = trackSeed.GetDzDs();
    if(checknan(d,"ALICE dz/ds",nseeds)) continue;
    double pY = track_pt*p;
    double pX = sqrt(track_pt*track_pt-pY*pY);
    track.set_px(pX*c-pY*s);
    track.set_py(pX*s+pY*c);
    track.set_pz(track_pt * trackSeed.GetDzDs()); 
    const double* cov = trackSeed.GetCov();
    bool cov_nan = false;
    for(int i=0;i<15;i++)
    {
      if(checknan(cov[i],"covariance element "+std::to_string(i),nseeds)) cov_nan = true;
    }
    if(cov_nan) continue;
    // make this into an actual Eigen matrix
    Eigen::Matrix<double,5,5> ecov;
    ecov(0,0)=cov[0];
    ecov(0,1)=cov[1];
    ecov(0,2)=cov[2];
    ecov(0,3)=cov[3];
    ecov(0,4)=cov[4];
    ecov(1,1)=cov[5];
    ecov(1,2)=cov[6];
    ecov(1,3)=cov[7];
    ecov(1,4)=cov[8];
    ecov(2,2)=cov[9];
    ecov(2,3)=cov[10];
    ecov(2,4)=cov[11];
    ecov(3,3)=cov[12];
    ecov(3,4)=cov[13];
    ecov(4,4)=cov[14];
    // symmetrize
    ecov(1,0)=ecov(0,1);
    ecov(2,0)=ecov(0,2);
    ecov(3,0)=ecov(0,3);
    ecov(4,0)=ecov(0,4);
    ecov(2,1)=ecov(1,2);
    ecov(3,1)=ecov(1,3);
    ecov(4,1)=ecov(1,4);
    ecov(3,2)=ecov(2,3);
    ecov(4,2)=ecov(2,4);
    ecov(4,3)=ecov(3,4);
    // make rotation matrix based on the following:
    // x = X*cos(track_phi) - Y*sin(track_phi)
    // y = X*sin(track_phi) + Y*cos(track_phi)
    // z = Z
    // pY = pt*sinphi
    // pX = sqrt(pt**2 - pY**2)
    // px = pX*cos(track_phi) - pY*sin(track_phi)
    // py = pX*sin(track_phi) + pY*cos(track_phi)
    // pz = pt*(dz/ds)
    Eigen::Matrix<double,6,5> J;
    J(0,0) = -s; // dx/dY
    J(0,1) = 0.; // dx/dZ
    J(0,2) = 0.; // dx/d(sinphi)
    J(0,3) = 0.; // dx/d(dz/ds)
    J(0,4) = 0.; // dx/d(Q/pt)

    J(1,0) = c;  // dy/dY
    J(1,1) = 0.; // dy/dZ
    J(1,2) = 0.; // dy/d(sinphi)
    J(1,3) = 0.; // dy/d(dz/ds)
    J(1,4) = 0.; // dy/d(Q/pt)

    J(2,0) = 0.; // dz/dY
    J(2,1) = 1.; // dz/dZ
    J(2,2) = 0.; // dz/d(sinphi)
    J(2,3) = 0.; // dz/d(dz/ds)
    J(2,4) = 0.; // dz/d(Q/pt)

    J(3,0) = 0.; // dpx/dY
    J(3,1) = 0.; // dpx/dZ
    J(3,2) = -track_pt*(p*c/sqrt(1-p*p)+s); // dpx/d(sinphi)
    J(3,3) = 0.; // dpx/d(dz/ds)
    J(3,4) = track_pt*track_pt*track_charge*(p*s-c*sqrt(1-p*p)); // dpx/d(Q/pt)

    J(4,0) = 0.; // dpy/dY
    J(4,1) = 0.; // dpy/dZ
    J(4,2) = track_pt*(c-p*s/sqrt(1-p*p)); // dpy/d(sinphi)
    J(4,3) = 0.; // dpy/d(dz/ds)
    J(4,4) = -track_pt*track_pt*track_charge*(p*c+s*sqrt(1-p*p)); // dpy/d(Q/pt)

    J(5,0) = 0.; // dpz/dY
    J(5,1) = 0.; // dpz/dZ
    J(5,2) = 0.; // dpz/d(sinphi)
    J(5,3) = track_pt; // dpz/d(dz/ds)
    J(5,4) = -track_pt*track_pt*track_charge*d; // dpz/d(Q/pt)
    bool cov_rot_nan = false;
    for(int i=0;i<6;i++)
    {
      for(int j=0;j<5;j++)
      {
        if(checknan(J(i,j),"covariance rotator element ("+std::to_string(i)+","+std::to_string(j)+")",nseeds))
        {
          cov_rot_nan = true;
          continue;
        }
      }
    }
    if(cov_rot_nan) continue;

    // the heavy lifting happens here
    Eigen::Matrix<double,6,6> scov = J*ecov*J.transpose();
    
    // fill SvtxTrack covariance matrix with results
    for(int i=0;i<6;i++)
    {
      for(int j=0;j<6;j++)
      {
        track.set_error(i, j, scov(i,j));
      }
    }
/*
    // Proceed with the absolutely hellish coordinate transformation of the covariance matrix.
    // Derived from:
    // 1) Taking the Jacobian of the conversion from (Y,Z,SinPhi,DzDs,Q/Pt) to (x,y,z,px,py,pz)
    // 2) Computing (Jacobian)*(ALICE covariance matrix)*(transpose of Jacobian)
    track.set_error(0, 0, cov[0]*s*s);
    track.set_error(0, 1, -cov[0]*c*s);
    track.set_error(0, 2, -cov[1]*s);
    track.set_error(0, 3, cov[2]*s*s/q-cov[4]*s*(-c/(q*q)+p*s/(q*q)));
    track.set_error(0, 4, -cov[2]*c*s/q-cov[4]*s*(-c*p/(q*q)-s/(q*q)));
    track.set_error(0, 5, cov[4]*d*s/(q*q)-cov[3]*s/q);
    track.set_error(1, 1, cov[0]*c*c);
    track.set_error(1, 2, cov[1]*c);
    track.set_error(1, 3, -cov[2]*c*s/q+cov[4]*c*(-c/(q*q)+p*s/(q*q)));
    track.set_error(1, 4, cov[2]*c*c/q+cov[4]*c*(-c*p/(q*q)-s/(q*q)));
    track.set_error(1, 5, cov[4]*d*c/(q*q)+cov[3]*c/q);
    track.set_error(2, 2, cov[5]);
    track.set_error(2, 3, -cov[6]*s/q+cov[8]*(-c/(q*q)+p*s/(q*q)));
    track.set_error(2, 4, cov[6]*c/q+cov[8]*(-c*p/(q*q)-s/(q*q)));
    track.set_error(2, 5, -cov[8]*d/(q*q)+cov[7]/q);
    track.set_error(3, 3, cov[9]*s*s/(q*q)-cov[11]*(-c/(q*q*q)+p*s/(q*q*q)) + (-c/(q*q)+p*s/(q*q))*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(3, 4, -cov[9]*c*s/(q*q)+cov[11]*(-c/(q*q*q)+p*s/(q*q*q)) + (-c*p/(q*q)-s/(q*q))*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(3, 5, -cov[10]*s/(q*q)+cov[13]/q*(-c/(q*q)+p*s/(q*q))-d/(q*q)*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(4, 4, c/q*(c/q*cov[9]+cov[11]*(-c*p/(q*q)-s/(q*q)))+(-c*p/(q*q)-s/(q*q))*(c/q*cov[11]+cov[14]*(-c*p/(q*q)-s/(q*q))));
    track.set_error(4, 5, cov[10]*c/(q*q)+cov[13]/q*(-c*p/(q*q)-s/(q*q))-d/(q*q)*(c/q*cov[11]+cov[14]*(-c*p/(q*q)-s/(q*q))));
    track.set_error(5, 5, -d/(q*q)*(-d*cov[14]/(q*q)+cov[13]/q)-d*cov[13]/(q*q*q)+cov[12]/(q*q));
    // symmetrize covariance
    track.set_error(1, 0, track.get_error(0, 1));
    track.set_error(2, 0, track.get_error(0, 2));
    track.set_error(3, 0, track.get_error(0, 3));
    track.set_error(4, 0, track.get_error(0, 4));
    track.set_error(5, 0, track.get_error(0, 5));
    track.set_error(2, 1, track.get_error(1, 2));
    track.set_error(3, 1, track.get_error(1, 3));
    track.set_error(4, 1, track.get_error(1, 4));
    track.set_error(5, 1, track.get_error(1, 5));
    track.set_error(3, 2, track.get_error(2, 3));
    track.set_error(4, 2, track.get_error(2, 4));
    track.set_error(5, 2, track.get_error(2, 5));
    track.set_error(4, 3, track.get_error(3, 4));
    track.set_error(5, 3, track.get_error(3, 5));
    track.set_error(5, 4, track.get_error(4, 5));
*/

    if(!covIsPosDef(track))
    {
      repairCovariance(track);
    }
/*
    for(int w=0;w<cx.size();w++)
    {
      ntp->Fill(cx[w],cy[w],cz[w],xerr[w],yerr[w],zerr[w],tx[w],ty[w],tz[w],layer[w],xsize[w],ysize[w],phisize[w],phierr[w],zsize[w]);
    }
    cx.clear();
    cy.clear();
    cz.clear();
    tx.clear();
    ty.clear();
    tz.clear();
    xerr.clear();
    yerr.clear();
    zerr.clear();
    layer.clear();
    xsize.clear();
    ysize.clear();
    phisize.clear();
    phierr.clear();
    zsize.clear();
*/
    seeds_vector.push_back(track);
    ++nseeds;
  }
//  f->cd();
//  ntp->Write();
//  f->Close();
  if(Verbosity()>0) cout << "number of seeds: " << nseeds << "\n";
  return seeds_vector;
}

Eigen::Matrix<double,6,6> ALICEKF::getEigenCov(SvtxTrack_v1 &track)
{
  Eigen::Matrix<double,6,6> cov;
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<6;j++)
    {
      cov(i,j) = track.get_error(i,j);
    }
  }
  return cov;
}

bool ALICEKF::covIsPosDef(SvtxTrack_v1 &track)
{
  // put covariance matrix into Eigen container
  Eigen::Matrix<double,6,6> cov = getEigenCov(track);
  // attempt Cholesky decomposition
  Eigen::LLT<Eigen::Matrix<double,6,6>> chDec(cov);
  // if Cholesky decomposition does not exist, matrix is not positive definite
  return (chDec.info() != Eigen::NumericalIssue);
}

void ALICEKF::repairCovariance(SvtxTrack_v1 &track)
{
  // find closest positive definite matrix
  Eigen::Matrix<double,6,6> cov = getEigenCov(track);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6>> solver(cov);
  Eigen::Matrix<double,6,1> D = solver.eigenvalues();
  Eigen::Matrix<double,6,6> Q = solver.eigenvectors();
  Eigen::Matrix<double,6,1> Dp = D.cwiseMax(1e-15);
  Eigen::Matrix<double,6,6> Z = Q*Dp.asDiagonal()*Q.transpose();
  // updates covariance matrix
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<6;j++)
    {
      track.set_error(i,j,Z(i,j));
    }
  }
}
void ALICEKF::CircleFitByTaubin (std::vector<std::pair<double,double>> points, double &R, double &X0, double &Y0)
/*  
      Circle fit to a given set of data points (in 2D)
      This is an algebraic fit, due to Taubin, based on the journal article
      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
                  Space Curves Defined By Implicit Equations, With 
                  Applications To Edge And Range Image Segmentation",
                  IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
*/
{
  int iter,IterMAX=99;
  
  double Mz,Mxy,Mxx,Myy,Mxz,Myz,Mzz,Cov_xy,Var_z;
  double A0,A1,A2,A22,A3,A33;
  double x,y;
  double DET,Xcenter,Ycenter;
  
  // Compute x- and y- sample means   
  double meanX = 0;
  double meanY = 0;
  double weight = 0;
  for(unsigned int i = 0; i < points.size(); ++i)
    {
      meanX += points[i].first;
      meanY += points[i].second;
      weight++;
    }
  meanX /= weight;
  meanY /= weight;

  //     computing moments 
  
  Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;
  
  for (unsigned int i=0; i<points.size(); i++)
    {
      double Xi = points[i].first - meanX;   //  centered x-coordinates
      double Yi = points[i].second - meanY;   //  centered y-coordinates
      double Zi = Xi*Xi + Yi*Yi;
      
      Mxy += Xi*Yi;
      Mxx += Xi*Xi;
      Myy += Yi*Yi;
      Mxz += Xi*Zi;
      Myz += Yi*Zi;
      Mzz += Zi*Zi;
    }
  Mxx /= weight;
  Myy /= weight;
  Mxy /= weight;
  Mxz /= weight;
  Myz /= weight;
  Mzz /= weight;
  
  //  computing coefficients of the characteristic polynomial
  
  Mz = Mxx + Myy;
  Cov_xy = Mxx*Myy - Mxy*Mxy;
  Var_z = Mzz - Mz*Mz;
  A3 = 4*Mz;
  A2 = -3*Mz*Mz - Mzz;
  A1 = Var_z*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
  A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
  A22 = A2 + A2;
  A33 = A3 + A3 + A3;
  
  //    finding the root of the characteristic polynomial
  //    using Newton's method starting at x=0  
  //    (it is guaranteed to converge to the right root)
  
  for (x=0.,y=A0,iter=0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
    {
      double Dy = A1 + x*(A22 + A33*x);
      double xnew = x - y/Dy;
      if ((xnew == x)||(!isfinite(xnew))) break;
      double ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
      if (fabs(ynew)>=fabs(y))  break;
      x = xnew;  y = ynew;
    }
  
  //  computing parameters of the fitting circle
  
  DET = x*x - x*Mz + Cov_xy;
  Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
  Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;
  
  //  assembling the output
  
  X0 = Xcenter + meanX;
  Y0 = Ycenter + meanY;
  R = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz);
}
