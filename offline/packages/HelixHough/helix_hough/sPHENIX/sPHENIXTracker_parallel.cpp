#include "sPHENIXTracker.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <sys/time.h>


using namespace std;
using namespace Eigen;
using namespace SeamStress;


void sPHENIXTracker::initSplitting(vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits)
{
  initEvent(hits, min_hits);
  (*(hits_vec[0])) = hits;
  zoomranges.clear();
  for(unsigned int z=0;z<=max_zoom;z++)
  {
    zoomranges.push_back(top_range);
  }
}


void sPHENIXTracker::findHelicesParallelOneHelicity(vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits, vector<SimpleTrack3D>& tracks)
{
  unsigned int hits_per_thread = (hits.size() + 2*nthreads)/nthreads;
  unsigned int pos=0;
  while(pos < hits.size())
  {
    for(unsigned int i=0;i<nthreads;++i)
    {
      if(pos>=hits.size()){break;}
      for(unsigned int j=0;j<hits_per_thread;++j)
      {
        if(pos>=hits.size()){break;}
        split_input_hits[i].push_back(hits[pos]);
        pos+=1;
      }
    }
  }
  for(unsigned int i=0;i<nthreads;++i)
  {
    thread_trackers[i]->setTopRange(top_range);
    thread_trackers[i]->initSplitting(split_input_hits[i], thread_min_hits, thread_max_hits);
  }
  pins->sewStraight(&sPHENIXTracker::splitHitsParallelThread, nthreads);
  thread_ranges.clear();
  thread_hits.clear();
  
  unsigned int nbins = split_output_hits[0]->size();
  for(unsigned int b=0;b<nbins;++b)
  {
    thread_ranges.push_back( (*(split_ranges[0]))[b] );
    thread_hits.push_back(vector<SimpleHit3D>());
    for(unsigned int i=0;i<nthreads;++i)
    {
      for(unsigned int j=0;j<(*(split_output_hits[i]))[b].size();++j)
      {
        thread_hits.back().push_back( (*(split_output_hits[i]))[b][j] );
      }
    }
  }
  
  pins->sewStraight(&sPHENIXTracker::findHelicesParallelThread, nthreads);
}


void sPHENIXTracker::findHelicesParallel(vector<SimpleHit3D>& hits, unsigned int min_hits, unsigned int max_hits, vector<SimpleTrack3D>& tracks)
{
  thread_min_hits = min_hits;
  thread_max_hits = max_hits;
  
  for(unsigned int i=0;i<nthreads;++i)
  {
    thread_tracks[i].clear();
    thread_trackers[i]->clear();
    if(cluster_start_bin!=0){thread_trackers[i]->setClusterStartBin(cluster_start_bin-1);}
    else{thread_trackers[i]->setClusterStartBin(0);}
  }
  
  initSplitting(hits, min_hits, max_hits);
  
  if(separate_by_helicity==true)
  {
    for(unsigned int i=0;i<nthreads;++i)
    {
      thread_trackers[i]->setSeparateByHelicity(true);
      thread_trackers[i]->setOnlyOneHelicity(true);
      thread_trackers[i]->setHelicity(true);
      split_output_hits[i]->clear();
      split_input_hits[i].clear();
    }
    findHelicesParallelOneHelicity(hits, min_hits, max_hits, tracks);
    
    for(unsigned int i=0;i<nthreads;++i)
    {
      thread_trackers[i]->setSeparateByHelicity(true);
      thread_trackers[i]->setOnlyOneHelicity(true);
      thread_trackers[i]->setHelicity(false);
      split_output_hits[i]->clear();
      split_input_hits[i].clear();
    }
    findHelicesParallelOneHelicity(hits, min_hits, max_hits, tracks);
  }
  else
  {
    for(unsigned int i=0;i<nthreads;++i)
    {
      thread_trackers[i]->setSeparateByHelicity(false);
      thread_trackers[i]->setOnlyOneHelicity(false);
      split_output_hits[i]->clear();
      split_input_hits[i].clear();
    }
    
    findHelicesParallelOneHelicity(hits, min_hits, max_hits, tracks);
  }
  
  vector<SimpleTrack3D> temp_tracks;
  for(unsigned int i=0;i<nthreads;++i)
  {
    vector<HelixKalmanState>* states = &(thread_trackers[i]->getKalmanStates());
    for(unsigned int j=0;j<thread_tracks[i].size();++j)
    {
      track_states.push_back((*states)[j]);
      temp_tracks.push_back(thread_tracks[i][j]);
    }
  }
  finalize(temp_tracks, tracks);
}


void sPHENIXTracker::splitHitsParallelThread(void* arg)
{
  unsigned long int w = (*((unsigned long int *)arg));
  thread_trackers[w]->splitIntoBins(thread_min_hits, thread_max_hits, *(split_ranges[w]), *(split_output_hits[w]), 0);
}


void sPHENIXTracker::findHelicesParallelThread(void* arg)
{
  unsigned long int w = (*((unsigned long int *)arg));
  
  for(unsigned int i=w;i<thread_ranges.size();i+=nthreads)
  {
    if(thread_hits[i].size() == 0){continue;}
    thread_trackers[w]->setTopRange(thread_ranges[i]);
    thread_trackers[w]->findHelices(thread_hits[i], thread_min_hits, thread_max_hits, thread_tracks[w]);
  }
}
