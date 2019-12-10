#include "HelixHough.h"

#include "HelixRange.h"
#include "SimpleHit3D.h"

#include <memory>
#include <vector>

using namespace std;

void HelixHough::splitIntoBins(unsigned int min_hits, unsigned int max_hits,
                               vector<HelixRange>& ranges,
                               vector<vector<SimpleHit3D> >& split_hits,
                               unsigned int zoomlevel) {
  index_mapping.clear();
  index_mapping.resize(hits_vec[zoomlevel]->size(), 0);
  for (unsigned int i = 0; i < hits_vec[zoomlevel]->size(); i++) {
    index_mapping[i] = (*(hits_vec[zoomlevel]))[i].get_id();
    (*(hits_vec[zoomlevel]))[i].set_id(i);
  }

  ranges.clear();
  split_hits.clear();

  vote(zoomlevel);

  unsigned int n_entries = bins_vec[zoomlevel]->size();
  if (n_entries == 0) {
    return;
  }

  unsigned int n_phi = n_phi_bins[zoomlevel];
  unsigned int n_d = n_d_bins[zoomlevel];
  unsigned int n_k = n_k_bins[zoomlevel];
  unsigned int n_dzdl = n_dzdl_bins[zoomlevel];
  unsigned int n_z0 = n_z0_bins[zoomlevel];

  unsigned int count = 0;
  hits_vec[zoomlevel + 1]->clear();
  setRange((*(bins_vec[zoomlevel]))[count], zoomranges[zoomlevel],
           zoomranges[zoomlevel + 1], n_phi, n_d, n_k, n_dzdl, n_z0);
  // scan over the bins in 5-D hough space
  while (count < n_entries) {
    hits_vec[zoomlevel + 1]->push_back(
        (*(hits_vec[zoomlevel]))[(*(bins_vec[zoomlevel]))[count].entry]);
    hits_vec[zoomlevel + 1]->back().set_id(
        index_mapping[hits_vec[zoomlevel + 1]->back().get_id()]);

    count += 1;
    if ((count == n_entries) || ((*(bins_vec[zoomlevel]))[count].bin !=
                                 (*(bins_vec[zoomlevel]))[count - 1].bin)) {
      if (((*(bins_vec[zoomlevel]))[count - 1].bin != 0) &&
          (split_hits.size() == 0)) {
        for (unsigned int b = 0; b < (*(bins_vec[zoomlevel]))[count - 1].bin;
             ++b) {
          BinEntryPair5D bp;
          bp.bin = b;
          setRange(bp, zoomranges[zoomlevel], zoomranges[zoomlevel + 1], n_phi,
                   n_d, n_k, n_dzdl, n_z0);
          ranges.push_back(zoomranges[zoomlevel + 1]);
          split_hits.push_back(vector<SimpleHit3D>());
        }
        setRange((*(bins_vec[zoomlevel]))[count - 1], zoomranges[zoomlevel],
                 zoomranges[zoomlevel + 1], n_phi, n_d, n_k, n_dzdl, n_z0);
      }

      ranges.push_back(zoomranges[zoomlevel + 1]);
      split_hits.push_back(*(hits_vec[zoomlevel + 1]));
      if (count == n_entries) {
        break;
      }
      hits_vec[zoomlevel + 1]->clear();

      unsigned int bin1 = (*(bins_vec[zoomlevel]))[count - 1].bin;
      unsigned int bin2 = (*(bins_vec[zoomlevel]))[count].bin;
      for (unsigned int b = (bin1 + 1); b < bin2; ++b) {
        BinEntryPair5D bp;
        bp.bin = b;
        setRange(bp, zoomranges[zoomlevel], zoomranges[zoomlevel + 1], n_phi,
                 n_d, n_k, n_dzdl, n_z0);
        ranges.push_back(zoomranges[zoomlevel + 1]);
        split_hits.push_back(vector<SimpleHit3D>());
      }

      setRange((*(bins_vec[zoomlevel]))[count], zoomranges[zoomlevel],
               zoomranges[zoomlevel + 1], n_phi, n_d, n_k, n_dzdl, n_z0);
    }
  }

  unsigned int total_bins = n_phi * n_d * n_k * n_dzdl * n_z0;
  for (unsigned int b = (1 + (*(bins_vec[zoomlevel]))[count - 1].bin);
       b < total_bins; ++b) {
    BinEntryPair5D bp;
    bp.bin = b;
    setRange(bp, zoomranges[zoomlevel], zoomranges[zoomlevel + 1], n_phi, n_d,
             n_k, n_dzdl, n_z0);
    ranges.push_back(zoomranges[zoomlevel + 1]);
    split_hits.push_back(vector<SimpleHit3D>());
  }
}
