#ifndef PHG4MAPSCELLRECO_H
#define PHG4MAPSCELLRECO_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <string>
#include <map>
#include <vector>

class PHCompositeNode;
class PHG4CylinderCell;

class PHG4MapsCellReco : public SubsysReco
{
 public:

  PHG4MapsCellReco(const std::string &name);

  virtual ~PHG4MapsCellReco(){}
  
  //! module initialization
  int InitRun(PHCompositeNode *topNode);
  
    //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);
  
  void Detector(const std::string &d) {detector = d;}
  void checkenergy(const int i=1) {chkenergyconservation = i;}

 protected:

  bool lines_intersect(
		       double ax,
		       double ay,
		       double bx,
		       double by,
		       double cx,
		       double cy,
		       double dx,
		       double dy,
		       double* rx, // intersection point (output)
		       double* ry
		       );
  
  bool  line_and_rectangle_intersect(
				     double ax,
				     double ay,
				     double bx,
				     double by,
				     double cx,
				     double cy,
				     double dx,
				     double dy,
				     double* rr // length of the line segment inside the rectangle (output)
				     );

  double circle_square_intersection(double x1, double y1,  double x2,  double y2,  double mx, double my,  double r);
  double sA(double r, double x, double y) ;
  
  //void set_size(const int i, const double sizeA, const int sizeB, const int what);
  int CheckEnergy(PHCompositeNode *topNode);
  static std::pair<double, double> get_etaphi(const double x, const double y, const double z);
  static double get_eta(const double radius, const double z);
  std::map<int, int>  binning;
  std::map<int, std::pair <double,int> > cell_size; // cell size in eta/nslats
  std::map<int, std::pair <double,double> > zmin_max; // zmin/zmax for each layer for faster lookup
  std::map<int, double> etastep;
  std::string detector;
  std::string hitnodename;
  std::string cellnodename;
  std::string geonodename;
  std::string seggeonodename;
  PHTimeServer::timer _timer;
  int nbins[2];
  int chkenergyconservation;

  std::map<unsigned long long, PHG4CylinderCell*> celllist;  // This map holds the hit cells
};

#endif
