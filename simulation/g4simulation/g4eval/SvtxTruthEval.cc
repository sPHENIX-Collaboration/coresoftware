#include "SvtxTruthEval.h"

#include "BaseTruthEval.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <g4main/PHG4Particle.h>

#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrDefs.h>

#include <tpc/TpcDefs.h>
#include <intt/InttDefs.h>
#include <mvtx/MvtxDefs.h>
#include <micromegas/MicromegasDefs.h>

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>               // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <mvtx/CylinderGeom_Mvtx.h>

#include <intt/CylinderGeomIntt.h>


#include <phool/getClass.h>
#include <phool/phool.h>                                // for PHWHERE

#include <TVector3.h>

#include <cassert>
#include <cmath>                                       // for sqrt, NAN, fabs
#include <cstdlib>                                     // for abs
#include <iostream>
#include <map>
#include <set>
#include <utility>

using namespace std;

SvtxTruthEval::SvtxTruthEval(PHCompositeNode* topNode)
  : _basetrutheval(topNode)
  , _truthinfo(nullptr)
  , _strict(false)
  , _verbosity(0)
  , _errors(0)
  , _do_cache(true)
  , _cache_all_truth_hits()
  , _cache_all_truth_hits_g4particle()
  , _cache_all_truth_clusters_g4particle()
  , _cache_get_innermost_truth_hit()
  , _cache_get_outermost_truth_hit()
  , _cache_get_primary_particle_g4hit()
{
  get_node_pointers(topNode);
  iclus = 0;
}

SvtxTruthEval::~SvtxTruthEval()
{
  if (_verbosity > 0)
  {
    if ((_errors > 0) || (_verbosity > 1))
    {
      cout << "SvtxTruthEval::~SvtxTruthEval() - Error Count: " << _errors << endl;
    }
  }
}

void SvtxTruthEval::next_event(PHCompositeNode* topNode)
{
  _cache_all_truth_hits.clear();
  _cache_all_truth_hits_g4particle.clear();
  _cache_all_truth_clusters_g4particle.clear();
  _cache_get_innermost_truth_hit.clear();
  _cache_get_outermost_truth_hit.clear();
  _cache_get_primary_particle_g4hit.clear();

  _basetrutheval.next_event(topNode);

  get_node_pointers(topNode);
}

/// \todo this copy may be too expensive to call a lot...
std::set<PHG4Hit*> SvtxTruthEval::all_truth_hits()
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    if (!_cache_all_truth_hits.empty())
    {
      return _cache_all_truth_hits;
    }
  }

  // since the SVTX can be composed of two different trackers this is a
  // handy function to spill out all the g4hits from both "detectors"

  std::set<PHG4Hit*> truth_hits;

  // loop over all the g4hits in the cylinder layers
  if (_g4hits_svtx)
  {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_svtx->getHits().first;
         g4iter != _g4hits_svtx->getHits().second;
         ++g4iter)
    {
      PHG4Hit* g4hit = g4iter->second;
      truth_hits.insert(g4hit);
    }
  }

  // loop over all the g4hits in the ladder layers
  if (_g4hits_tracker)
  {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_tracker->getHits().first;
         g4iter != _g4hits_tracker->getHits().second;
         ++g4iter)
    {
      PHG4Hit* g4hit = g4iter->second;
      truth_hits.insert(g4hit);
    }
  }

  // loop over all the g4hits in the maps layers
  if (_g4hits_maps)
  {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_maps->getHits().first;
         g4iter != _g4hits_maps->getHits().second;
         ++g4iter)
    {
      PHG4Hit* g4hit = g4iter->second;
      truth_hits.insert(g4hit);
    }
  }

 // loop over all the g4hits in the maicromegas layers
  if (_g4hits_mms)
  {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_mms->getHits().first;
         g4iter != _g4hits_mms->getHits().second;
         ++g4iter)
    {
      PHG4Hit* g4hit = g4iter->second;
      truth_hits.insert(g4hit);
    }
  }

  if (_do_cache) _cache_all_truth_hits = truth_hits;

  return truth_hits;
}

std::set<PHG4Hit*> SvtxTruthEval::all_truth_hits(PHG4Particle* particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits_g4particle.find(particle);
    if (iter != _cache_all_truth_hits_g4particle.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all the g4hits in the cylinder layers
  if (_g4hits_svtx)
  {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_svtx->getHits().first;
         g4iter != _g4hits_svtx->getHits().second;
         ++g4iter)
    {
      PHG4Hit* g4hit = g4iter->second;
      if (!is_g4hit_from_particle(g4hit, particle)) continue;
      truth_hits.insert(g4hit);
    }
  }

  // loop over all the g4hits in the ladder layers
  if (_g4hits_tracker)
  {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_tracker->getHits().first;
         g4iter != _g4hits_tracker->getHits().second;
         ++g4iter)
    {
      PHG4Hit* g4hit = g4iter->second;
      if (!is_g4hit_from_particle(g4hit, particle)) continue;
      truth_hits.insert(g4hit);
    }
  }

  // loop over all the g4hits in the maps ladder layers
  if (_g4hits_maps)
  {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_maps->getHits().first;
         g4iter != _g4hits_maps->getHits().second;
         ++g4iter)
    {
      PHG4Hit* g4hit = g4iter->second;
      if (!is_g4hit_from_particle(g4hit, particle)) continue;
      truth_hits.insert(g4hit);
    }
  }

  // loop over all the g4hits in the micromegas layers
  if (_g4hits_mms)
  {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_mms->getHits().first;
         g4iter != _g4hits_mms->getHits().second;
         ++g4iter)
    {
      PHG4Hit* g4hit = g4iter->second;
      if (!is_g4hit_from_particle(g4hit, particle)) continue;
      truth_hits.insert(g4hit);
    }
  }

  if (_do_cache) _cache_all_truth_hits_g4particle.insert(make_pair(particle, truth_hits));

  return truth_hits;
}

std::map<unsigned int, std::shared_ptr<TrkrCluster> > SvtxTruthEval::all_truth_clusters(PHG4Particle* particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::map<unsigned int, std::shared_ptr<TrkrCluster> >();
  }

  if (_strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++_errors;
    return std::map<unsigned int, std::shared_ptr<TrkrCluster> >();
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, std::map<unsigned int, std::shared_ptr<TrkrCluster> > >::iterator iter =
        _cache_all_truth_clusters_g4particle.find(particle);
    if (iter != _cache_all_truth_clusters_g4particle.end())
    {
      return iter->second;
    }
  }

  if(_verbosity > 0)
    cout << PHWHERE << " Truth clustering for particle " << particle->get_track_id() << endl;;

  // get all g4hits for this particle
  std::set<PHG4Hit*> g4hits = all_truth_hits(particle);

  float ng4hits = g4hits.size();
  if(ng4hits == 0)
    return std::map<unsigned int, std::shared_ptr<TrkrCluster> >();

  // container for storing truth clusters
  //std::map<unsigned int, TrkrCluster*> truth_clusters;
  std::map<unsigned int, std::shared_ptr<TrkrCluster>> truth_clusters;

  // convert truth hits for this particle to truth clusters in each layer
  // loop over layers
  unsigned int layer;
  for(layer = 0; layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc +_nlayers_mms; ++layer)
    {
      float gx = NAN;
      float gy = NAN;
      float gz = NAN;
      float gt = NAN;
      float gedep = NAN;

      std::vector<PHG4Hit*> contributing_hits;
      std::vector<double> contributing_hits_energy;
      std::vector<std::vector<double>> contributing_hits_entry;
      std::vector<std::vector<double>> contributing_hits_exit;
      LayerClusterG4Hits(g4hits, contributing_hits, contributing_hits_energy, contributing_hits_entry, contributing_hits_exit, layer, gx, gy, gz, gt, gedep);
      if(!(gedep > 0)) continue;

      // we have the cluster in this layer from this truth particle
      // add the cluster to a TrkrCluster object
      TrkrDefs::cluskey ckey;
      if(layer >= _nlayers_maps + _nlayers_intt && layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)  // in TPC
	{
	  unsigned int side = 0;
	  if(gz > 0) side = 1;
	  // need dummy sector here
	  unsigned int sector = 0;
	  ckey = TpcDefs::genClusKey(layer, sector, side, iclus);
	}
      else if(layer < _nlayers_maps)  // in MVTX
	{
	  unsigned int stave = 0;
	  unsigned int chip = 0;
	  ckey = MvtxDefs::genClusKey(layer, stave, chip, iclus);
	}
      else if(layer >= _nlayers_maps && layer < _nlayers_maps  + _nlayers_intt)  // in INTT
	{
	  // dummy ladder and phi ID
	  unsigned int ladderzid = 0;
	  unsigned int ladderphiid = 0;
	  ckey = InttDefs::genClusKey(layer, ladderzid, ladderphiid,iclus);
	}
      else if(layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc)    // in MICROMEGAS
	{
	  unsigned int tile = 0;
	  MicromegasDefs::SegmentationType segtype;
	  segtype  =  MicromegasDefs::SegmentationType::SEGMENTATION_PHI;
	  TrkrDefs::hitsetkey hkey = MicromegasDefs::genHitSetKey(layer, segtype, tile);
  	  ckey = MicromegasDefs::genClusterKey(hkey, iclus);
	}
      else
	{
	  std::cout << PHWHERE << "Bad layer number: " << layer << std::endl;
	  continue;
	}

      std::shared_ptr<TrkrClusterv1> clus(new TrkrClusterv1());
      clus->setClusKey(ckey);
      iclus++;

      // estimate cluster ADC value
      unsigned int adc_value = getAdcValue(gedep);
      //std::cout << " gedep " << gedep << " adc_value " << adc_value << std::endl;
      clus->setAdc(adc_value);
      clus->setPosition(0, gx);
      clus->setPosition(1, gy);
      clus->setPosition(2, gz);
      clus->setGlobal();

      // record which g4hits contribute to this truth cluster
      for(unsigned int i=0; i< contributing_hits.size(); ++i)
	{
	  _truth_cluster_truth_hit_map.insert(make_pair(ckey,contributing_hits[i]));
	}

      // Estimate the size of the truth cluster
      float g4phisize = NAN;
      float g4zsize = NAN;
      G4ClusterSize(layer, contributing_hits_entry, contributing_hits_exit, g4phisize, g4zsize);

      for(int i1=0;i1<3;++i1)
	for(int i2=0;i2<3;++i2)
	{
	  clus->setSize(i1, i2, 0.0);
	  clus->setError(i1, i2, 0.0);
	}
      clus->setError(0,0,gedep);  // stores truth energy
      clus->setSize(1, 1, g4phisize);
      clus->setSize(2, 2, g4zsize);
      clus->setError(1, 1, g4phisize/sqrt(12));
      clus->setError(2, 2, g4zsize/sqrt(12.0));

      truth_clusters.insert(make_pair(layer, clus));

    }  // end loop over layers for this particle

  if (_do_cache) _cache_all_truth_clusters_g4particle.insert(make_pair(particle, truth_clusters));

  return truth_clusters;
}

void SvtxTruthEval::LayerClusterG4Hits(std::set<PHG4Hit*> truth_hits, std::vector<PHG4Hit*> &contributing_hits, std::vector<double> &contributing_hits_energy, std::vector<std::vector<double>> &contributing_hits_entry, std::vector<std::vector<double>> &contributing_hits_exit, float layer, float &x, float &y, float &z,  float &t, float &e)
{
  bool use_geo = true;

  // Given a set of g4hits, cluster them within a given layer of the TPC

  float gx = 0.0;
  float gy = 0.0;
  float gz = 0.0;
  float gr = 0.0;
  float gt = 0.0;
  float gwt = 0.0;

  if (layer >= _nlayers_maps + _nlayers_intt && layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)   // in TPC
    {
      //cout << "layer = " << layer << " _nlayers_maps " << _nlayers_maps << " _nlayers_intt " << _nlayers_intt << endl;

      // This calculates the truth cluster position for the TPC from all of the contributing g4hits from a g4particle, typically 2-4 for the TPC
      // Complicated, since only the part of the energy that is collected within a layer contributes to the position
      //===============================================================================

      PHG4CylinderCellGeom* GeoLayer = _tpc_geom_container->GetLayerCellGeom(layer);
      // get layer boundaries here for later use
      // radii of layer boundaries
      float rbin = GeoLayer->get_radius() - GeoLayer->get_thickness() / 2.0;
      float rbout = GeoLayer->get_radius() + GeoLayer->get_thickness() / 2.0;

      if(_verbosity > 0)
	cout << " TruthEval::LayerCluster hits for layer  " << layer << " with rbin " << rbin << " rbout " << rbout << endl;

      // we do not assume that the truth hits know what layer they are in
      for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
	   iter != truth_hits.end();
	   ++iter)
	{

	  PHG4Hit* this_g4hit = *iter;
	  float rbegin = sqrt(this_g4hit->get_x(0) * this_g4hit->get_x(0) + this_g4hit->get_y(0) * this_g4hit->get_y(0));
	  float rend = sqrt(this_g4hit->get_x(1) * this_g4hit->get_x(1) + this_g4hit->get_y(1) * this_g4hit->get_y(1));
	  //cout << " Eval: g4hit " << this_g4hit->get_hit_id() <<  " layer " << layer << " rbegin " << rbegin << " rend " << rend << endl;

	  // make sure the entry point is at lower radius
	  float xl[2];
	  float yl[2];
	  float zl[2];

	  if (rbegin < rend)
	    {
	      xl[0] = this_g4hit->get_x(0);
	      yl[0] = this_g4hit->get_y(0);
	      zl[0] = this_g4hit->get_z(0);
	      xl[1] = this_g4hit->get_x(1);
	      yl[1] = this_g4hit->get_y(1);
	      zl[1] = this_g4hit->get_z(1);
	    }
	  else
	    {
	      xl[0] = this_g4hit->get_x(1);
	      yl[0] = this_g4hit->get_y(1);
	      zl[0] = this_g4hit->get_z(1);
	      xl[1] = this_g4hit->get_x(0);
	      yl[1] = this_g4hit->get_y(0);
	      zl[1] = this_g4hit->get_z(0);
	      swap(rbegin, rend);
	      //cout << "swapped in and out " << endl;
	    }

	  // check that the g4hit is not completely outside the cluster layer. Just skip this g4hit if it is
	  if (rbegin < rbin && rend < rbin)
	    continue;
	  if (rbegin > rbout && rend > rbout)
	    continue;


	  if(_verbosity > 0)
	    {
	      cout << "     keep g4hit with rbegin " << rbegin << " rend " << rend
		   << "         xbegin " <<  xl[0] << " xend " << xl[1]
		   << " ybegin " << yl[0] << " yend " << yl[1]
		   << " zbegin " << zl[0] << " zend " << zl[1]
		   << endl;
	    }


	  float xin = xl[0];
	  float yin = yl[0];
	  float zin = zl[0];
	  float xout = xl[1];
	  float yout = yl[1];
	  float zout = zl[1];

	  float t = NAN;

	  if (rbegin < rbin)
	    {
	      // line segment begins before boundary, find where it crosses
	      t = line_circle_intersection(xl, yl, zl, rbin);
	      if (t > 0)
		{
		  xin = xl[0] + t * (xl[1] - xl[0]);
		  yin = yl[0] + t * (yl[1] - yl[0]);
		  zin = zl[0] + t * (zl[1] - zl[0]);
		}
	    }

	  if (rend > rbout)
	    {
	      // line segment ends after boundary, find where it crosses
	      t = line_circle_intersection(xl, yl, zl, rbout);
	      if (t > 0)
		{
		  xout = xl[0] + t * (xl[1] - xl[0]);
		  yout = yl[0] + t * (yl[1] - yl[0]);
		  zout = zl[0] + t * (zl[1] - zl[0]);
		}
	    }

	  double rin = sqrt(xin*xin + yin*yin);
	  double rout = sqrt(xout*xout + yout*yout);

	  // we want only the fraction of edep inside the layer
	  double efrac =  this_g4hit->get_edep() * (rout - rin) / (rend - rbegin);
	  gx += (xin + xout) * 0.5 * efrac;
	  gy += (yin + yout) * 0.5 * efrac;
	  gz += (zin + zout) * 0.5 * efrac;
	  gt += this_g4hit->get_avg_t() * efrac;
	  gr += (rin + rout) * 0.5 * efrac;
	  gwt += efrac;

	  if(_verbosity > 0)
	    {
	      cout << "      rin  " << rin << " rout " << rout
		   << " xin " << xin << " xout " << xout << " yin " << yin << " yout " << yout << " zin " << zin << " zout " << zout
		   << " edep " << this_g4hit->get_edep()
		   << " this_edep " <<  efrac << endl;
	      cout << "              xavge " << (xin+xout) * 0.5 << " yavge " << (yin+yout) * 0.5 << " zavge " << (zin+zout) * 0.5 << " ravge " << (rin+rout) * 0.5
		   << endl;
	    }

	  // Capture entry and exit points
	  std::vector<double> entry_loc;
	  entry_loc.push_back(xin);
	  entry_loc.push_back(yin);
	  entry_loc.push_back(zin);
	  std::vector<double> exit_loc;
	  exit_loc.push_back(xout);
	  exit_loc.push_back(yout);
	  exit_loc.push_back(zout);

	  // this_g4hit is inside the layer, add it to the vectors
	  contributing_hits.push_back(this_g4hit);
	  contributing_hits_energy.push_back( this_g4hit->get_edep() * (zout - zin) / (zl[1] - zl[0]) );
	  contributing_hits_entry.push_back(entry_loc);
	  contributing_hits_exit.push_back(exit_loc);

	}  // loop over contributing hits

      if(gwt == 0)
	{
	  e = gwt;
	  return;  // will be discarded
	}

      gx /= gwt;
      gy /= gwt;
      gz /= gwt;
      gr = (rbin + rbout) * 0.5;
      gt /= gwt;

      if(_verbosity > 0)
	{
	  cout << " weighted means:   gx " << gx << " gy " << gy << " gz " << gz << " gr " << gr << " e " << gwt << endl;
	}

      if(use_geo)
	{
	  // The energy weighted values above have significant scatter due to fluctuations in the energy deposit from Geant
	  // Calculate the geometric mean positions instead
	  float rentry = 999.0;
	  float xentry = 999.0;
	  float yentry = 999.0;
	  float zentry = 999.0;
	  float rexit = - 999.0;
	  float xexit = -999.0;
	  float yexit = -999.0;
	  float zexit = -999.0;

	  for(unsigned int ientry = 0; ientry < contributing_hits_entry.size(); ++ientry)
	    {
	      float tmpx = contributing_hits_entry[ientry][0];
	      float tmpy = contributing_hits_entry[ientry][1];
	      float tmpr = sqrt(tmpx*tmpx + tmpy*tmpy);

	      if(tmpr < rentry)
		{
		  rentry =  tmpr;
		  xentry = contributing_hits_entry[ientry][0];
		  yentry = contributing_hits_entry[ientry][1];
		  zentry = contributing_hits_entry[ientry][2];
		}

	      tmpx = contributing_hits_exit[ientry][0];
	      tmpy = contributing_hits_exit[ientry][1];
	      tmpr = sqrt(tmpx*tmpx + tmpy*tmpy);

	      if(tmpr > rexit)
		{
		  rexit =  tmpr;
		  xexit = contributing_hits_exit[ientry][0];
		  yexit = contributing_hits_exit[ientry][1];
		  zexit = contributing_hits_exit[ientry][2];
		}
	    }

	  float geo_r = (rentry+rexit)*0.5;
	  float geo_x = (xentry+xexit)*0.5;
	  float geo_y = (yentry+yexit)*0.5;
	  float geo_z = (zentry+zexit)*0.5;

	  if(rexit > 0)
	    {
	      gx = geo_x;
	      gy = geo_y;
	      gz = geo_z;
	      gr = geo_r;
	    }


	  if(_verbosity > 0)
	    {
	      cout << "      rentry  " << rentry << " rexit " << rexit
		   << " xentry " << xentry << " xexit " << xexit << " yentry " << yentry << " yexit " << yexit << " zentry " << zentry << " zexit " << zexit << endl;

	      cout  << " geometric means: geo_x " << geo_x << " geo_y " << geo_y << " geo_z " << geo_z  << " geo r " << geo_r <<  " e " << gwt << endl << endl;
	    }
	}

    }  // if TPC
  else
    {
      // not TPC, one g4hit per cluster
      for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
	   iter != truth_hits.end();
	   ++iter)
	{

	  PHG4Hit* this_g4hit = *iter;

	  if(this_g4hit->get_layer() != (unsigned int) layer) continue;

	  gx = this_g4hit->get_avg_x();
	  gy = this_g4hit->get_avg_y();
	  gz = this_g4hit->get_avg_z();
	  gt = this_g4hit->get_avg_t();
	  gwt += this_g4hit->get_edep();

	  // Capture entry and exit points
	  std::vector<double> entry_loc;
	  entry_loc.push_back(this_g4hit->get_x(0));
	  entry_loc.push_back(this_g4hit->get_y(0));
	  entry_loc.push_back(this_g4hit->get_z(0));
	  std::vector<double> exit_loc;
	  exit_loc.push_back(this_g4hit->get_x(1));
	  exit_loc.push_back(this_g4hit->get_y(1));
	  exit_loc.push_back(this_g4hit->get_z(1));

	  // this_g4hit is inside the layer, add it to the vectors
	  contributing_hits.push_back(this_g4hit);
	  contributing_hits_energy.push_back( this_g4hit->get_edep() );
	  contributing_hits_entry.push_back(entry_loc);
	  contributing_hits_exit.push_back(exit_loc);
	}
    }  // not TPC

  x = gx;
  y = gy;
  z = gz;
  t = gt;
  e = gwt;

  return;
}

void SvtxTruthEval::G4ClusterSize(unsigned int layer, std::vector<std::vector<double>> contributing_hits_entry,std::vector<std::vector<double>> contributing_hits_exit, float &g4phisize, float &g4zsize)
{

  // sort the contributing g4hits in radius
  double inner_radius = 100.;
  double inner_x = NAN;
  double inner_y = NAN;
  double inner_z = NAN;;

  double outer_radius = 0.;
  double outer_x = NAN;
  double outer_y = NAN;
  double outer_z = NAN;

  for(unsigned int ihit=0;ihit<contributing_hits_entry.size(); ++ihit)
    {
      double rad1 = sqrt(pow(contributing_hits_entry[ihit][0], 2) + pow(contributing_hits_entry[ihit][1], 2));
      if(rad1 < inner_radius)
	{
	  inner_radius = rad1;
	  inner_x = contributing_hits_entry[ihit][0];
	  inner_y = contributing_hits_entry[ihit][1];
	  inner_z = contributing_hits_entry[ihit][2];
	}

      double rad2 = sqrt(pow(contributing_hits_exit[ihit][0], 2) + pow(contributing_hits_exit[ihit][1], 2));
      if(rad2 > outer_radius)
	{
	  outer_radius = rad2;
	  outer_x = contributing_hits_exit[ihit][0];
	  outer_y = contributing_hits_exit[ihit][1];
	  outer_z = contributing_hits_exit[ihit][2];
	}
    }

  double inner_phi =  atan2(inner_y, inner_x);
  double outer_phi =  atan2(outer_y, outer_x);
  double avge_z = (outer_z + inner_z) / 2.0;

  // Now fold these with the expected diffusion and shaping widths
  // assume spread is +/- equals this many sigmas times diffusion and shaping when extending the size
  double sigmas = 2.0;

  double radius = (inner_radius + outer_radius)/2.;
  if(radius > 28 && radius < 80)  // TPC
    {
      PHG4CylinderCellGeom*layergeom = _tpc_geom_container->GetLayerCellGeom(layer);

      double tpc_length = 211.0;  // cm
      double drift_velocity = 8.0 / 1000.0;  // cm/ns

      // Phi size
      //======
      double diffusion_trans =  0.006;  // cm/SQRT(cm)
      double phidiffusion = diffusion_trans * sqrt(tpc_length / 2. - fabs(avge_z));

      double added_smear_trans = 0.085; // cm
      double gem_spread = 0.04;  // 400 microns

      if(outer_phi < inner_phi) swap(outer_phi, inner_phi);

      // convert diffusion from cm to radians
      double g4max_phi =  outer_phi + sigmas * sqrt(  pow(phidiffusion, 2) + pow(added_smear_trans, 2) + pow(gem_spread, 2) ) / radius;
      double g4min_phi =  inner_phi - sigmas * sqrt(  pow(phidiffusion, 2) + pow(added_smear_trans, 2) + pow(gem_spread, 2) ) / radius;

      // find the bins containing these max and min z edges
      unsigned int phibinmin = layergeom->get_phibin(g4min_phi);
      unsigned int phibinmax = layergeom->get_phibin(g4max_phi);
      unsigned int phibinwidth = phibinmax - phibinmin + 1;
      g4phisize = (double) phibinwidth * layergeom->get_phistep() * layergeom->get_radius();

      // Z size
      //=====
      double g4max_z = 0;
      double g4min_z = 0;

      outer_z = fabs(outer_z);
      inner_z = fabs(inner_z);

      double diffusion_long = 0.015;  // cm/SQRT(cm)
      double zdiffusion = diffusion_long * sqrt(tpc_length / 2. - fabs(avge_z)) ;
      double zshaping_lead = 32.0 * drift_velocity;  // ns * cm/ns = cm
      double zshaping_tail = 48.0 * drift_velocity;
      double added_smear_long = 0.105;  // cm

      // largest z reaches gems first, make that the outer z
      if(outer_z < inner_z) swap(outer_z, inner_z);
      g4max_z = outer_z  + sigmas*sqrt(pow(zdiffusion,2) + pow(added_smear_long,2) + pow(zshaping_lead, 2));
      g4min_z = inner_z  -  sigmas*sqrt(pow(zdiffusion,2) + pow(added_smear_long,2) + pow(zshaping_tail, 2));

      // find the bins containing these max and min z edges
      unsigned int binmin = layergeom->get_zbin(g4min_z);
      unsigned int binmax = layergeom->get_zbin(g4max_z);
      if(binmax < binmin) swap(binmax, binmin);
      unsigned int binwidth = binmax - binmin + 1;

      // multiply total number of bins that include the edges by the bin size
      g4zsize = (double) binwidth * layergeom->get_zstep();
    }
  else if(radius > 5 && radius < 20)  // INTT
    {
      // All we have is the position and layer number

      CylinderGeomIntt *layergeom = dynamic_cast<CylinderGeomIntt *>(_intt_geom_container->GetLayerGeom(layer));

      // inner location
      double world_inner[3] = {inner_x, inner_y, inner_z};
      TVector3 world_inner_vec = {inner_x, inner_y, inner_z};

      int segment_z_bin, segment_phi_bin;
      layergeom->find_indices_from_world_location(segment_z_bin, segment_phi_bin, world_inner);

      TVector3 local_inner_vec =  layergeom->get_local_from_world_coords(segment_z_bin, segment_phi_bin, world_inner_vec);
      double yin = local_inner_vec[1];
      double zin = local_inner_vec[2];
      int strip_y_index, strip_z_index;
      layergeom->find_strip_index_values(segment_z_bin, yin, zin, strip_y_index, strip_z_index);

	// outer location
      double world_outer[3] = {outer_x, outer_y, outer_z};
      TVector3 world_outer_vec = {outer_x, outer_y, outer_z};

      layergeom->find_indices_from_world_location(segment_z_bin, segment_phi_bin, world_outer);

      TVector3 local_outer_vec =  layergeom->get_local_from_world_coords(segment_z_bin, segment_phi_bin, world_outer_vec);
      double yout = local_outer_vec[1];
      double zout = local_outer_vec[2];
      int strip_y_index_out, strip_z_index_out;
      layergeom->find_strip_index_values(segment_z_bin, yout, zout, strip_y_index_out, strip_z_index_out);

      int strips = abs(strip_y_index_out - strip_y_index) + 1;
      int cols = abs(strip_z_index_out - strip_z_index) + 1;


      double strip_width = (double) strips * layergeom->get_strip_y_spacing(); // cm
      double strip_length = (double) cols * layergeom->get_strip_z_spacing(); // cm

      g4phisize = strip_width;
      g4zsize = strip_length;

      /*
      if(Verbosity() > 1)
	cout << " INTT: layer " << layer << " strips " << strips << " strip pitch " <<  layergeom->get_strip_y_spacing() << " g4phisize "<< g4phisize
	     << " columns " << cols << " strip_z_spacing " <<  layergeom->get_strip_z_spacing() << " g4zsize " << g4zsize << endl;
      */
    }
  else if(radius > 80)  // MICROMEGAS
    {
      // made up for now
      g4phisize = 300e-04;
      g4zsize = 300e-04;
    }
  else  // MVTX
    {
      unsigned int stave, stave_outer;
      unsigned int chip, chip_outer;
      int row, row_outer;
      int column, column_outer;

      // add diffusion to entry and exit locations
      double max_diffusion_radius = 25.0e-4;  // 25 microns
      double min_diffusion_radius = 8.0e-4;  // 8 microns

      CylinderGeom_Mvtx *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(_mvtx_geom_container->GetLayerGeom(layer));

      TVector3 world_inner = {inner_x, inner_y, inner_z};
      std::vector<double> world_inner_vec = { world_inner[0], world_inner[1], world_inner[2] };
      layergeom->get_sensor_indices_from_world_coords(world_inner_vec, stave, chip);
      TVector3 local_inner = layergeom->get_local_from_world_coords(stave, chip, world_inner);

      TVector3 world_outer = {outer_x, outer_y, outer_z};
      std::vector<double> world_outer_vec = { world_outer[0], world_outer[1], world_outer[2] };
      layergeom->get_sensor_indices_from_world_coords(world_outer_vec, stave_outer, chip_outer);
      TVector3 local_outer = layergeom->get_local_from_world_coords(stave_outer, chip_outer, world_outer);

      double diff =  max_diffusion_radius * 0.6;  // factor of 0.6 gives decent agreement with low occupancy reco clusters
      if(local_outer[0] < local_inner[0])
	diff = -diff;
      local_outer[0] += diff;
      local_inner[0] -= diff;

      double diff_outer = min_diffusion_radius * 0.6;
      if(local_outer[2] < local_inner[2])
	diff_outer = -diff_outer;
      local_outer[2] += diff_outer;
      local_inner[2] -= diff_outer;

      layergeom->get_pixel_from_local_coords(local_inner, row, column);
      layergeom->get_pixel_from_local_coords(local_outer, row_outer, column_outer);

      if(row_outer < row) swap(row_outer, row);
      unsigned int rows = row_outer - row + 1;
      g4phisize = (double) rows * layergeom->get_pixel_x();

      if(column_outer < column) swap(column_outer, column);
      unsigned int columns = column_outer - column + 1;
      g4zsize = (double) columns * layergeom->get_pixel_z();

      /*
      if(Verbosity() > 1)
	cout << " MVTX: layer " << layer << " rows " << rows << " pixel x " <<  layergeom->get_pixel_x() << " g4phisize "<< g4phisize
	     << " columns " << columns << " pixel_z " <<  layergeom->get_pixel_z() << " g4zsize " << g4zsize << endl;
      */

    }

}

std::set<PHG4Hit*> SvtxTruthEval::get_truth_hits_from_truth_cluster(TrkrDefs::cluskey ckey)
{
  std::set<PHG4Hit *> g4hit_set;


  std::pair<std::multimap<TrkrDefs::cluskey, PHG4Hit*>::iterator,
	    std::multimap<TrkrDefs::cluskey,PHG4Hit*>::iterator>  ret =   _truth_cluster_truth_hit_map.equal_range(ckey);

  for(std::multimap<TrkrDefs::cluskey, PHG4Hit*>::iterator jter = ret.first; jter != ret.second; ++jter)
    {
      g4hit_set.insert(jter->second);
    }

  return g4hit_set;
}

PHG4Hit* SvtxTruthEval::get_innermost_truth_hit(PHG4Particle* particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++_errors;
    return nullptr;
  }

  PHG4Hit* innermost_hit = nullptr;
  float innermost_radius = FLT_MAX;

  std::set<PHG4Hit*> truth_hits = all_truth_hits(particle);
  for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
       iter != truth_hits.end();
       ++iter)
  {
    PHG4Hit* candidate = *iter;
    float x = candidate->get_x(0);  // use entry points
    float y = candidate->get_y(0);  // use entry points
    float r = sqrt(x * x + y * y);
    if (r < innermost_radius)
    {
      innermost_radius = r;
      innermost_hit = candidate;
    }
  }

  return innermost_hit;
}

PHG4Hit* SvtxTruthEval::get_outermost_truth_hit(PHG4Particle* particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++_errors;
    return nullptr;
  }

  PHG4Hit* outermost_hit = nullptr;
  float outermost_radius = FLT_MAX * -1.0;

  if (_do_cache)
  {
    std::map<PHG4Particle*, PHG4Hit*>::iterator iter =
        _cache_get_outermost_truth_hit.find(particle);
    if (iter != _cache_get_outermost_truth_hit.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits = all_truth_hits(particle);
  for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
       iter != truth_hits.end();
       ++iter)
  {
    PHG4Hit* candidate = *iter;
    float x = candidate->get_x(1);  // use exit points
    float y = candidate->get_y(1);  // use exit points
    float r = sqrt(x * x + y * y);
    if (r > outermost_radius)
    {
      outermost_radius = r;
      outermost_hit = candidate;
    }
  }
  if (_do_cache) _cache_get_outermost_truth_hit.insert(make_pair(particle, outermost_hit));

  return outermost_hit;
}

PHG4Particle* SvtxTruthEval::get_particle(PHG4Hit* g4hit)
{
  return _basetrutheval.get_particle(g4hit);
}

int SvtxTruthEval::get_embed(PHG4Particle* particle)
{
  return _basetrutheval.get_embed(particle);
}

PHG4VtxPoint* SvtxTruthEval::get_vertex(PHG4Particle* particle)
{
  return _basetrutheval.get_vertex(particle);
}

bool SvtxTruthEval::is_primary(PHG4Particle* particle)
{
  return _basetrutheval.is_primary(particle);
}

PHG4Particle* SvtxTruthEval::get_primary_particle(PHG4Hit* g4hit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(g4hit);
  }
  else if (!g4hit)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<PHG4Hit*, PHG4Particle*>::iterator iter =
        _cache_get_primary_particle_g4hit.find(g4hit);
    if (iter != _cache_get_primary_particle_g4hit.end())
    {
      return iter->second;
    }
  }

  PHG4Particle* primary = _basetrutheval.get_primary_particle(g4hit);

  if (_do_cache) _cache_get_primary_particle_g4hit.insert(make_pair(g4hit, primary));

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
  }

  return primary;
}

PHG4Particle* SvtxTruthEval::get_primary_particle(PHG4Particle* particle)
{
  return _basetrutheval.get_primary_particle(particle);
}

PHG4Particle* SvtxTruthEval::get_particle(const int trackid)
{
  return _basetrutheval.get_particle(trackid);
}

bool SvtxTruthEval::is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle)
{
  return _basetrutheval.is_g4hit_from_particle(g4hit, particle);
}

bool SvtxTruthEval::are_same_particle(PHG4Particle* p1, PHG4Particle* p2)
{
  return _basetrutheval.are_same_particle(p1, p2);
}

bool SvtxTruthEval::are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2)
{
  return _basetrutheval.are_same_vertex(vtx1, vtx2);
}

void SvtxTruthEval::get_node_pointers(PHCompositeNode* topNode)
{
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  _g4hits_mms = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");
  _g4hits_svtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  _g4hits_tracker = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  _g4hits_maps = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");

  _mms_geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  _tpc_geom_container = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  _intt_geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  _mvtx_geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  return;
}

bool SvtxTruthEval::has_node_pointers()
{
  if (_strict)
    assert(_truthinfo);
  else if (!_truthinfo)
    return false;

  if (_strict)
    assert(_g4hits_mms || _g4hits_svtx || _g4hits_tracker || _g4hits_maps);
  else if (!_g4hits_mms && !_g4hits_svtx && !_g4hits_tracker && !_g4hits_maps)
    return false;

  return true;
}

float SvtxTruthEval::line_circle_intersection(float x[], float y[], float z[], float radius)
{
  // parameterize the line in terms of t (distance along the line segment, from 0-1) as
  // x = x0 + t * (x1-x0); y=y0 + t * (y1-y0); z = z0 + t * (z1-z0)
  // parameterize the cylinder (centered at x,y = 0,0) as  x^2 + y^2 = radius^2,   then
  // (x0 + t*(x1-z0))^2 + (y0+t*(y1-y0))^2 = radius^2
  // (x0^2 + y0^2 - radius^2) + (2x0*(x1-x0) + 2y0*(y1-y0))*t +  ((x1-x0)^2 + (y1-y0)^2)*t^2 = 0 = C + B*t + A*t^2
  // quadratic with:  A = (x1-x0)^2+(y1-y0)^2 ;  B = 2x0*(x1-x0) + 2y0*(y1-y0);  C = x0^2 + y0^2 - radius^2
  // solution: t = (-B +/- sqrt(B^2 - 4*A*C)) / (2*A)

  float A = (x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]);
  float B = 2.0 * x[0] * (x[1] - x[0]) + 2.0 * y[0] * (y[1] - y[0]);
  float C = x[0] * x[0] + y[0] * y[0] - radius * radius;
  float tup = (-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
  float tdn = (-B - sqrt(B * B - 4.0 * A * C)) / (2.0 * A);

  // The limits are 0 and 1, but we allow a little for floating point precision
  float t;
  if (tdn >= -0.0e-4 && tdn <= 1.0004)
    t = tdn;
  else if (tup >= -0.0e-4 && tup <= 1.0004)
    t = tup;
  else
  {
    cout << PHWHERE << "   **** Oops! No valid solution for tup or tdn, tdn = " << tdn << " tup = " << tup << endl;
    cout << "   radius " << radius << " rbegin " << sqrt(x[0] * x[0] + y[0] * y[0]) << " rend " << sqrt(x[1] * x[1] + y[1] * y[1]) << endl;
    cout << "   x0 " << x[0] << " x1 " << x[1] << endl;
    cout << "   y0 " << y[0] << " y1 " << y[1] << endl;
    cout << "   z0 " << z[0] << " z1 " << z[1] << endl;
    cout << "   A " << A << " B " << B << " C " << C << endl;

    t = -1;
  }

  return t;
}

unsigned int SvtxTruthEval::getAdcValue(double gedep)
{
  // see TPC digitizer for algorithm

  // drift electrons per GeV of energy deposited in the TPC
  double Ne_dEdx = 1.56;   // keV/cm
  double CF4_dEdx = 7.00;  // keV/cm
  double Ne_NTotal = 43;    // Number/cm
  double CF4_NTotal = 100;  // Number/cm
  double Tpc_NTot = 0.5*Ne_NTotal + 0.5*CF4_NTotal;
  double Tpc_dEdx = 0.5*Ne_dEdx + 0.5*CF4_dEdx;
  double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  double electrons_per_gev = Tpc_ElectronsPerKeV * 1e6;

  double gem_amplification = 1400; // GEM output electrons per drifted electron
  double input_electrons = gedep * electrons_per_gev * gem_amplification;

  // convert electrons after GEM to ADC output
  double ChargeToPeakVolts = 20;
  double ADCSignalConversionGain = ChargeToPeakVolts * 1.60e-04 * 2.4;  // 20 (or 30) mV/fC * fC/electron * scaleup factor
  double adc_input_voltage = input_electrons * ADCSignalConversionGain;  // mV, see comments above
  unsigned int adc_output = (unsigned int) (adc_input_voltage * 1024.0 / 2200.0);  // input voltage x 1024 channels over 2200 mV max range
  if (adc_output > 1023) adc_output = 1023;

  return adc_output;
}
