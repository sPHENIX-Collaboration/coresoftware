#ifndef CALORECO_PHMAKEGROUPS_H
#define CALORECO_PHMAKEGROUPS_H
// Requirements:
//
// the class type Hit needs to provide:
//
//   an operator< that sorts by ix and then by iz
//   a function is_adjacent() which returns true if the argumment hit is adjacent
//   a function ...() that returns true if the argument Hit is far enough away
//      from this one to allow breaking out of the inner loop early
//

#include <boost/bind/bind.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#pragma GCC diagnostic pop

#include <map>
#include <vector>

template <class Hit>
int PHMakeGroups(std::vector<Hit>& hits,
                 std::multimap<int, Hit>& groups)
{
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

  Graph G;

  // Process the hits in ix-then-iz order
  std::sort(hits.begin(), hits.end());

  // TODO: Since the list is sorted by channel number, it should
  // be possible to terminate the inner loop if we find the next hit
  // is more than one Z channel away (in which case the subsequent ones
  // will be adjacent neither in X nor Z).
  for (unsigned int i = 0; i < hits.size(); i++)
  {
    for (unsigned int j = i + 1; j < hits.size(); j++)
    {
      if (hits[i].is_adjacent(hits[j])) add_edge(i, j, G);
    }
    add_edge(i, i, G);
  }

  // Find the connections between the vertices of the graph (vertices are the rawhits,
  // connections are made when they are adjacent to one another)
  std::vector<int> component(num_vertices(G));
  //connected_components(G, &component[0]);
  connected_components(G, &component[0]);
  //std::cout << "Found " << num << " groups of hits" << std::endl;

  // Loop over the components(vertices) compiling a list of the unique
  // connections (ie clusters).
  std::set<int> comps;  // Number of unique components
  for (unsigned int i = 0; i < component.size(); i++)
  {
    comps.insert(component[i]);
    groups.insert(std::make_pair(component[i], hits[i]));
  }

  //       for(std::set<int>::const_iterator id=comps.begin(); id!=comps.end(); id++)
  // 	{
  // 	  std::multimap<int,SvxRawhitAdapter>::const_iterator curr, last;
  // 	  boost::tie(curr,last) = groups[iread].equal_range(*id);
  // 	  std::cout << "Group " << *id << " has " << groups[iread].count(*id) << " Rawhits:" << std::endl;
  // 	  for ( ; curr!=last; curr++)
  // 	    {
  // 	      SvxRawhitAdapter h = curr->second;
  // 	      h.hit->print();
  // 	    }
  // 	}

  return 0;
}

#endif
