/*!
 *  \file       PHTpcEventExporter.cc
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */
#include "PHTpcEventExporter.h"

#include "Track.h"  // for Track

#include "externals/kdfinder.hpp"  // for get_track_color, TrackCa...

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <phool/PHLog.h>

#include <log4cpp/CategoryStream.hh>  // for CategoryStream

#include <TVector3.h>  // for TVector3

#include <cstddef>  // for size_t
#include <map>      // for _Rb_tree_const_iterator
#include <ostream>  // for operator<<, basic_ostream
#include <utility>  // for pair

PHTpcEventExporter::PHTpcEventExporter()
{
}

void PHTpcEventExporter::exportEvent(TrkrClusterContainer* cluster_map, TrkrHitSetContainer* hitsets,std::vector<kdfinder::TrackCandidate<double>*> candidates,
                                     double B, const std::string& filename)
{
  // export hits + seeds to json for checks
  LOG_DEBUG("tracking.PHTpcEventExporter.exportEvent") << "exporting event, hits + seeds => file " << filename;

  std::stringstream ofs;

  ofs << "{ \n"
      << " \"EVENT\": {  \n"
      << "  \"runid\": " << 0 << ", \n"
      << "  \"evtid\": " << 0 << ", \n"
      << "  \"time\": " << 0 << ", \n"
      << "  \"B\": " << B << " \n"
      << " }, \n";

  ofs << "\"META\": { \n"
      << " \"HITS\": { \n"
      << "    \"TPC\": { \n"
      << "       \"type\": \"3D\", \n"
      << "       \"options\": { \n"
      << "           \"size\": 2, \n"
      << "           \"color\": 16777215 \n"
      << "       }\n"
      << "   }\n"
      << " },\n"
      << " \"TRACKS\": { \n"
      << "    \"TPC\": { \n"
      << "      \"r_min\": 0,\n"
      << "      \"r_max\": 800,\n"
      << "      \"size\": 2, \n"
      << "      \"thickness\": 2 \n"
      << "    }\n"
      << "  }\n"
      << "},\n"
      << " \"TRACKS\": {\n"
      << "   \"TPC\": [\n";

  for (int i = 0, ilen = candidates.size(); i < ilen; i++)
  {
    std::vector<double> hit = candidates[i]->getPosForHit(0);
    std::vector<double> mom = candidates[i]->getMomForHit(0);
    size_t nhits = candidates[i]->nhits();
    double pt = candidates[i]->Pt();
    double sign = candidates[i]->sign();
    double length = candidates[i]->approxLength();
    size_t color = kdfinder::get_track_color<double>(pt);
    ofs << "{ \"color\": " << color << ", \"pt\": " << pt << ", \"xyz\":[" << hit[0] << "," << hit[1] << "," << hit[2]
        << "], \"pxyz\":[" << mom[0] << "," << mom[1] << "," << mom[2]
        << "],\"l\":" << length << ",\"nh\":" << nhits << ",\"q\":" << sign << "}";
    if (i != (ilen - 1))
    {
      ofs << ",";
    }
    else
    {
      ofs << "\n";
    }
  }

  ofs << "    ]\n"
      << "  },\n"
      << " \"HITS\": {\n"
      << "   \"TPC\": [\n";
  
  auto hitsetrange = hitsets->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (auto hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr){
    std::string separator = "";
    auto range = cluster_map->getClusters(hitsetitr->first);
    for( auto it = range.first; it != range.second; ++it )
      {
	TrkrCluster* cluster = it->second;
	if ((std::pow((double) cluster->getPosition(0), 2) +
	     std::pow((double) cluster->getPosition(1), 2) +
	     std::pow((double) cluster->getPosition(2), 2)) > (25.0 * 25.0))
	  {
	    ofs << "[ " << cluster->getPosition(0) << "," << cluster->getPosition(1) << "," << cluster->getPosition(2) << " ]";
	    separator = ",";
	  }
	else
	  {
	    separator = "";
	  }
	++it;
	if (it != range.second)
	  {
	    ofs << separator;
	  }
	else
	  {
	    ofs << "\n";
	    break;
	  }
      }
  }
  
  ofs << "  ]\n"
      << " }\n"
      << "}\n";
  
  std::ofstream ofile;
  ofile.open(filename);
  ofile << ofs.str();
  ofile.close();
}

void PHTpcEventExporter::exportEvent(TrkrClusterContainer* cluster_map, TrkrHitSetContainer* hitsets, std::vector<PHGenFit2::Track*> gtracks,
                                     double B, const std::string& filename)
{
  // export hits + reco-d tracks to json for checks
  LOG_DEBUG("tracking.PHTpcEventExporter.exportEvent") << "exporting event, hits + GenFit tracks => file " << filename;

  std::stringstream ofs;

  ofs << "{ \n"
      << " \"EVENT\": {  \n"
      << "  \"runid\": " << 0 << ", \n"
      << "  \"evtid\": " << 0 << ", \n"
      << "  \"time\": " << 0 << ", \n"
      << "  \"B\": " << B << " \n"
      << " }, \n";

  ofs << "\"META\": { \n"
      << " \"HITS\": { \n"
      << "    \"TPC\": { \n"
      << "       \"type\": \"3D\", \n"
      << "       \"options\": { \n"
      << "           \"size\": 2, \n"
      << "           \"color\": 16777215 \n"
      << "       }\n"
      << "   }\n"
      << " },\n"
      << " \"TRACKS\": { \n"
      << "    \"TPC\": { \n"
      << "      \"r_min\": 0,\n"
      << "      \"r_max\": 800,\n"
      << "      \"size\": 2, \n"
      << "      \"thickness\": 2 \n"
      << "    }\n"
      << "  }\n"
      << "},\n"
      << " \"TRACKS\": {\n"
      << "   \"TPC\": [\n";

  TVector3 pos, mom;
  double charge, length, pt;
  int nhits;
  size_t color;
  for (int i = 0, ilen = gtracks.size(); i < ilen; i++)
  {
    bool rc = true;
    try
    {
      rc = gtracks[i]->get_track_info(pos, mom, charge, nhits, length);
    }
    catch (...)
    {
      // GenFit tracklength exception?
      continue;
    }
    if (!rc)
    {
      LOG_WARN("tracking.PHTpcEventExporter.exportEvent") << "got track with broken track info, id: " << i;
      continue;
    }
    pt = mom.Perp();
    color = kdfinder::get_track_color<double>(pt);
    ofs << "{ \"color\": " << color << ", \"pt\": " << pt << ", \"xyz\":[" << pos.X() << "," << pos.Y() << "," << pos.Z()
        << "], \"pxyz\":[" << mom.X() << "," << mom.Y() << "," << mom.Z() << "],\"l\":"
        << length << ",\"nh\":" << nhits << ",\"q\":" << (charge < 0 ? -1 : 1) << "}";
    if (i != (ilen - 1))
    {
      ofs << ",";
    }
    else
    {
      ofs << "\n";
    }
  }

  ofs << "    ]\n"
      << "  },\n"
      << " \"HITS\": {\n"
      << "   \"TPC\": [\n";
  auto hitsetrange = hitsets->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (auto hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr){
    std::string separator = "";
    auto range = cluster_map->getClusters(hitsetitr->first);
    for( auto it = range.first; it != range.second; ++it )
      {
	TrkrCluster* cluster = it->second;
	if ((std::pow((double) cluster->getPosition(0), 2) +
	     std::pow((double) cluster->getPosition(1), 2) +
	     std::pow((double) cluster->getPosition(2), 2)) > (25.0 * 25.0))
	  {
	    ofs << "[ " << cluster->getPosition(0) << "," << cluster->getPosition(1) << "," << cluster->getPosition(2) << " ]";
	    separator = ",";
	  }
	else
	  {
	    separator = "";
	  }
	++it;
	if (it != range.second)
	  {
	    ofs << separator;
	  }
	else
	  {
	    ofs << "\n";
	    break;
	  }
      }
  }

  ofs << "  ]\n"
      << " }\n"
      << "}\n";

  std::ofstream ofile;
  ofile.open(filename);
  ofile << ofs.str();
  ofile.close();
}

void PHTpcEventExporter::exportEvent(std::vector<PHGenFit2::Track*> gtracks,
                                     double B, const std::string& filename)
{
  // export final reco-d tracks to json purely for live display purposes (no hits = much smaller file size)
  LOG_DEBUG("tracking.PHTpcEventExporter.exportEvent") << "exporting event, just GenFit tracks => file " << filename;

  std::stringstream ofs;

  ofs << "{ \n"
      << " \"EVENT\": {  \n"
      << "  \"runid\": " << 0 << ", \n"
      << "  \"evtid\": " << 0 << ", \n"
      << "  \"time\": " << 0 << ", \n"
      << "  \"B\": " << B << " \n"
      << " }, \n";

  ofs << "\"META\": { \n"
      << " \"HITS\": { \n"
      << "    \"TPC\": { \n"
      << "       \"type\": \"3D\", \n"
      << "       \"options\": { \n"
      << "           \"size\": 2, \n"
      << "           \"color\": 16777215 \n"
      << "       }\n"
      << "   }\n"
      << " },\n"
      << " \"TRACKS\": { \n"
      << "    \"TPC\": { \n"
      << "      \"r_min\": 0,\n"
      << "      \"r_max\": 800,\n"
      << "      \"size\": 2, \n"
      << "      \"thickness\": 2 \n"
      << "    }\n"
      << "  }\n"
      << "},\n"
      << " \"TRACKS\": {\n"
      << "   \"TPC\": [\n";

  TVector3 pos, mom;
  double charge, length, pt;
  int nhits;
  size_t color;
  for (int i = 0, ilen = gtracks.size(); i < ilen; i++)
  {
    gtracks[i]->get_track_info(pos, mom, charge, nhits, length);
    pt = mom.Perp();
    color = kdfinder::get_track_color<double>(pt);
    ofs << "{ \"color\": " << color << ", \"pt\": " << pt << ", \"xyz\":[" << pos.X() << "," << pos.Y() << "," << pos.Z()
        << "], \"pxyz\":[" << mom.X() << "," << mom.Y() << "," << mom.Z() << "],\"l\":"
        << length << ",\"nh\":" << nhits << ",\"q\":" << (charge < 0 ? -1 : 1) << "}";
    if (i != (ilen - 1))
    {
      ofs << ",";
    }
    else
    {
      ofs << "\n";
    }
  }

  ofs << "    ]\n"
      << "  }\n"
      << "}\n";

  std::ofstream ofile;
  ofile.open(filename);
  ofile << ofs.str();
  ofile.close();
}
