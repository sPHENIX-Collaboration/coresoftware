#ifndef __CYLINDRICALHOUGH__
#define __CYLINDRICALHOUGH__

#include "CircleHough.h"
#include <cmath>
#include <iostream>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


class AngleIndexPair
{
  public:
    AngleIndexPair(float ang, unsigned int idx) : angle(ang), index(idx)
    {
      float twopi = 2.*M_PI;
      int a = (int)(angle/twopi);
      angle -= a*twopi;
      while(angle < 0.){angle += twopi;}
      while(angle >= twopi){angle -= twopi;}
    }
    ~AngleIndexPair(){}
    
    bool operator<(const AngleIndexPair& other) const
    {
      return angle<other.angle;
    }
    
    static float absDiff(float angle1, float angle2)
    {
      float diff = ( angle1 - angle2 );
      while(diff > M_PI){diff -= 2.*M_PI;}
      while(diff < -M_PI){diff += 2.*M_PI;}
      diff = fabs(diff);
      return diff;
    }
    
    float angle;
    unsigned int index;
};


class AngleIndexList
{
  public:
    AngleIndexList() : sorted(false) {}
    ~AngleIndexList(){}
    
    void addPair(AngleIndexPair& angind)
    {
      sorted=false;
      vec.push_back(angind);
    }
    
    
    void getRangeListSimple(float angle, float error, std::vector<AngleIndexPair*>& result)
    {
      result.clear();
      
      for(unsigned int i=0;i<vec.size();i++)
      {
        if(AngleIndexPair::absDiff(angle, vec[i].angle) <= error)
        {
          result.push_back(&(vec[i]));
        }
      }
    }
    
    void getRangeList(float angle, float error, std::vector<AngleIndexPair*>& result)
    {
      float twopi = 2.*M_PI;
      int a = (int)(angle/twopi);
      angle -= a*twopi;
      while(angle < 0.){angle += twopi;}
      while(angle >= twopi){angle -= twopi;}
      
      if(vec.size() <= 4){return getRangeListSimple(angle, error, result);}
      
      result.clear();
      
      unsigned int closest = findClosest(angle);
      //first, traverse upward
      unsigned int current = closest;
      unsigned int lowest = 0;
      unsigned int highest = vec.size()-1;
      while(true)
      {
        if(AngleIndexPair::absDiff(angle, vec[current].angle) <= error)
        {
          result.push_back(&(vec[current]));
          current = (current+1)%(vec.size());
          if(current==closest){break;}
        }
        else
        {
          break;
        }
      }
      
      if(closest==0){return;}
      
      //now, traverse downward
      if(current <= closest)
      {
        lowest=current;
      }
      else
      {
        highest=current;
      }
      current = closest-1;
      while(true)
      {
        if(AngleIndexPair::absDiff(angle, vec[current].angle) <= error)
        {
          result.push_back(&(vec[current]));
          if( (current==lowest) || (current==highest) ){break;}
          current = ((current + vec.size()) - 1)%(vec.size());
        }
        else
        {
          break;
        }
      }
    }
    
  private:
    unsigned int findClosestSimple(float angle, unsigned int lower, unsigned int upper)
    {
      unsigned int closest = lower;
      float diff = AngleIndexPair::absDiff(vec[closest].angle, angle);
      for(unsigned int i=(lower+1);i<=upper;i++)
      {
        float tempdiff = AngleIndexPair::absDiff(vec[i].angle, angle);
        if( tempdiff < diff )
        {
          closest = i;
          diff = tempdiff;
        }
      }
      
      return closest;
    }
    
    
    
    unsigned int findClosest(float angle)
    {
      if(vec.size() <= 4){return findClosestSimple(angle, 0, vec.size()-1);}
      
      if(sorted==false)
      {
        std::sort(vec.begin(), vec.end());
        sorted=true;
      }
      
      unsigned int lower = 0;
      unsigned int upper = vec.size() - 1;
      unsigned int middle = vec.size()/2;
      while(true)
      {
        if((upper - lower) <= 4){return findClosestSimple(angle, lower, upper);}
        
        if(angle <= vec[middle].angle)
        {
          upper = middle;
          middle = (lower + upper)/2;
        }
        else
        {
          lower = middle;
          middle = (lower + upper)/2;
        }
      }
    }
    
    std::vector<AngleIndexPair> vec;
    bool sorted;
};


class CylindricalHough : public CircleHough
{
  public:
    CylindricalHough(std::vector<float>& detrad, unsigned int inv_radius_nbin, unsigned int center_angle_nbin, unsigned int dca_origin_nbin, CircleResolution& min_resolution, CircleResolution& max_resolution, CircleRange& range, unsigned int z0_nbin, unsigned int theta_nbin, ZResolution& minzres, ZResolution& maxzres, ZRange& zrange, double sxy=70.e-4, double sz=500.e-4);
    ~CylindricalHough();
    
    void customFindHelicesInit(std::vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits, unsigned int min_zhits, unsigned int max_zhits, double chi2_cut, float xydiffcut, std::vector<SimpleTrack3D>& tracks, unsigned int maxtracks=0);
    
    void addHits(unsigned int zlevel, std::vector<SimpleTrack3D>& temptracks, std::vector<SimpleTrack3D>& tracks, std::vector<float>& params, int tracks_per_hit, float z_cut);
    
    void init_ZHough(int z0_nbin, unsigned int theta_nbin, ZResolution& minzres, ZResolution& maxzres, ZRange& zrange);
    
    void setVertex(double vx, double vy, double vz)
    {
      vertex_x = vx;
      vertex_y = vy;
      vertex_z = vz;
    }
    
    
    void setPhiCut(double pc){phicut=pc;}
    
    bool intersect_circles(bool hel, double startx, double starty, double rad_det, double rad_trk, double cx, double cy, double& x, double& y);
    
    void setLayerResolution(std::vector<double>& lxy, std::vector<double>& lz);
    void setVertexResolution(double vxy, double vz);
    
  private:
    std::vector<AngleIndexList> angle_list;
    std::vector<float> detector_radii;
    double vertex_x, vertex_y, vertex_z;
    double phicut;
    
};




#endif
