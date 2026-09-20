#include "PHSiliconSeedPrunerHelper.h"

#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <gsl/gsl_randist.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace PHSiliconSeedPrunerHelper
{
  // a container for a seed's cluster keys, its index in the TrackSeedContainer, a flag for whether it has a complete MVTX triplet
  struct SeedRecord
  {
    size_t index = 0;
    std::vector<TrkrDefs::cluskey> allClusterKeys;
    std::vector<TrkrDefs::cluskey> mvtxClusterKeys;
    bool completeMvtx = false;
  };

  using ConflictGraph = std::vector<std::vector<size_t>>;

  // standard disjoint-set/union-find data structure
  class DisjointSet
  {
   public:
    explicit DisjointSet(size_t n)
      : parent_(n)
      , rank_(n, 0)
    {
      std::iota(parent_.begin(), parent_.end(), 0);
    }

    size_t Find(size_t value)
    {
      // find the root of a value, and perform path compression
      if (parent_[value] != value)
      {
        parent_[value] = Find(parent_[value]);  // the path compression step
      }
      return parent_[value];
    }

    // merge sets
    void Union(size_t lhs, size_t rhs)
    {
      size_t root_lhs = Find(lhs);
      size_t root_rhs = Find(rhs);
      if (root_lhs == root_rhs)
      {
        return;
      }

      // The lower-rank tree is attached under the higher-rank tree
      if (rank_[root_lhs] < rank_[root_rhs])
      {
        std::swap(root_lhs, root_rhs);
      }

      // combine the tree and update the rank if necessary
      parent_[root_rhs] = root_lhs;
      if (rank_[root_lhs] == rank_[root_rhs])
      {
        ++rank_[root_lhs];
      }
    }

   private:
    std::vector<size_t> parent_;
    std::vector<unsigned int> rank_;  // rank -> the depth of the tree
  };

  SeedRecord MakeSeedRecord(size_t seedIndex, const TrackSeed &seed, size_t mvtxLayerCount)
  {
    SeedRecord record;
    record.index = seedIndex;

    std::vector<unsigned int> mvtx_layer_counts(mvtxLayerCount, 0);
    for (auto iter = seed.begin_cluster_keys(); iter != seed.end_cluster_keys(); ++iter)
    {
      const TrkrDefs::cluskey key = *iter;
      record.allClusterKeys.push_back(key);

      if (TrkrDefs::getTrkrId(key) != TrkrDefs::mvtxId)
      {
        continue;
      }

      record.mvtxClusterKeys.push_back(key);
      const unsigned int layer = TrkrDefs::getLayer(key);
      if (layer < mvtxLayerCount)
      {
        ++mvtx_layer_counts[layer];
      }
    }

    record.completeMvtx =
        record.mvtxClusterKeys.size() == mvtxLayerCount &&
        mvtx_layer_counts[0] == 1 &&
        mvtx_layer_counts[1] == 1 &&
        mvtx_layer_counts[2] == 1;
    return record;
  }

  // Turn a conflict graph into a list of connected components
  // seeds that share only INTT clusters are not connected here (they are still in the same ambiguity group)
  std::vector<std::vector<size_t>> ConnectedComponents(const ConflictGraph &conflict)
  {
    std::vector<std::vector<size_t>> components;
    std::vector<bool> visited(conflict.size(), false);

    for (size_t start = 0; start < conflict.size(); ++start)
    {
      if (visited[start])
      {
        continue;
      }

      std::vector<size_t> component;
      std::vector<size_t> stack{start};
      visited[start] = true;

      while (!stack.empty())
      {
        const size_t vertex = stack.back();  // look at the top of stack
        stack.pop_back();                    // pop
        component.push_back(vertex);

        // inspect the neighbors
        for (const size_t neighbor : conflict[vertex])
        {
          if (!visited[neighbor])
          {
            visited[neighbor] = true;
            stack.push_back(neighbor);
          }
        }
      }

      components.push_back(std::move(component));
    }

    return components;
  }

  // collapse duplicate triplets
  struct TripletClass
  {
    std::vector<TrkrDefs::cluskey> mvtxKeys;  //
    std::vector<size_t> memberSeedIndices;    // global seed indices in this triplet class
  };

  // Collapse seeds that have the identical MVTX triplet into 1 class (they contribute only 1 representative)
  std::vector<TripletClass> CollapseDuplicateTriplets(  //
      const std::vector<const SeedRecord *> &completeSeeds,
      bool debugVerbosity)
  {
    if (debugVerbosity)
    {
      std::cout << __func__ << " : " << __LINE__ << " : collapsing " << completeSeeds.size() << " complete MVTX seeds into unique MVTX-triplet classes" << std::endl;
    }

    std::map<std::vector<TrkrDefs::cluskey>, size_t> classByTriplet;
    std::vector<TripletClass> classes;

    for (const SeedRecord *seed : completeSeeds)
    {
      const auto inserted = classByTriplet.emplace(seed->mvtxClusterKeys, classes.size());  // return a std::pair<iterator, bool> where the bool is true if the insertion took place
      // if not already present, create a new TripletClass
      if (inserted.second)
      {
        TripletClass cls;
        cls.mvtxKeys = seed->mvtxClusterKeys;
        classes.push_back(std::move(cls));
      }
      classes[inserted.first->second].memberSeedIndices.push_back(seed->index);
    }

    return classes;
  }

  // Build a conflict graph on the TripletClasses
  ConflictGraph BuildClassConflictGraph(const std::vector<TripletClass> &classes)
  {
    // map of cluster keys to the list of triplet classes that contain it
    std::unordered_map<TrkrDefs::cluskey, std::vector<size_t>> classesByKey;
    for (size_t c = 0; c < classes.size(); ++c)
    {
      for (const TrkrDefs::cluskey key : classes[c].mvtxKeys)
      {
        classesByKey[key].push_back(c);
      }
    }

    // build the conflict graph (nodes and edges) from the key index
    ConflictGraph conflict(classes.size());
    for (const auto &entry : classesByKey)
    {
      const std::vector<size_t> &sharing = entry.second;
      for (size_t a = 0; a < sharing.size(); ++a)
      {
        for (size_t b = a + 1; b < sharing.size(); ++b)
        {
          // make edges in both directions
          conflict[sharing[a]].push_back(sharing[b]);
          conflict[sharing[b]].push_back(sharing[a]);
        }
      }
    }

    // sort and remove duplicates from the adjacency lists
    for (std::vector<size_t> &neighbors : conflict)
    {
      std::sort(neighbors.begin(), neighbors.end());
      neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }

    return conflict;
  }

  // a container for 1 connected component of the class-conflict graph
  struct ComponentProblem
  {
    size_t m = 0;
    size_t numKeys = 0;
    std::vector<std::vector<size_t>> adj;           // local class adjacency
    std::vector<std::vector<size_t>> keyIds;        // local key ids per class (3 each)
    std::vector<std::vector<size_t>> keyToClasses;  // local key id -> local class ids
    std::vector<double> memberCount;                // class multiplicity (sampling weight)
    std::vector<size_t> globalClassId;              // local -> global class index
  };

  ComponentProblem BuildComponentProblem(const std::vector<TripletClass> &classes, const ConflictGraph &classConflict, const std::vector<size_t> &component)
  {
    ComponentProblem prob;
    prob.m = component.size();
    prob.globalClassId = component;
    prob.adj.resize(prob.m);
    prob.keyIds.resize(prob.m);
    prob.memberCount.resize(prob.m);

    std::unordered_map<size_t, size_t> localOf;  // global class index -> local class index
    for (size_t local = 0; local < component.size(); ++local)
    {
      localOf[component[local]] = local;
    }

    std::unordered_map<TrkrDefs::cluskey, size_t> keyId;  // MVTX cluster key -> local key id

    for (size_t local = 0; local < prob.m; ++local)
    {
      const size_t global = component[local];
      prob.memberCount[local] = static_cast<double>(classes[global].memberSeedIndices.size());

      for (const TrkrDefs::cluskey key : classes[global].mvtxKeys)
      {
        const auto inserted = keyId.emplace(key, keyId.size());
        prob.keyIds[local].push_back(inserted.first->second);
      }

      for (const size_t global_neighbor : classConflict[global])
      {
        const auto it = localOf.find(global_neighbor);
        if (it != localOf.end())
        {
          prob.adj[local].push_back(it->second);
        }
      }
    }

    prob.numKeys = keyId.size();
    prob.keyToClasses.assign(prob.numKeys, {});
    for (size_t local = 0; local < prob.m; ++local)
    {
      for (const size_t kid : prob.keyIds[local])
      {
        prob.keyToClasses[kid].push_back(local);
      }
    }

    return prob;
  }

  // Actual solver
  class ComponentSolver
  {
   public:
    explicit ComponentSolver(             //
        const ComponentProblem &problem,  //
        bool debugDetailed = false,       //
        size_t debugComponentIndex = 0    //
        )
      : prob_(problem)
      ,  //
      active_(problem.m, true)
      ,  //
      scratchKeyCount_(problem.numKeys, 0)
      ,  //
      covered_(problem.m, false)
      ,  //
      debugDetailed_(debugDetailed)
      ,                                          //
      debugComponentIndex_(debugComponentIndex)  //
    {
    }

    // Search with a budget (the number of nodes to explore). True if the search completes, false if the budget is hit first (then go to `BuildHeuristicSet`)
    // This gives "how big can the best non-conflicting set (i.e cardinality) of classes be?" and *one* example of a set that size "bestSet_"
    bool SolveBudgeted(size_t budget)
    {
      std::fill(active_.begin(), active_.end(), true);
      bestCardinality_ = 0;
      bestSet_.clear();
      budget_ = budget;
      nodes_ = 0;
      budgetHit_ = false;
      debugDepth_ = 0;
      std::vector<size_t> chosen;
      chosen.reserve(prob_.m);
      if (debugDetailed_)
      {
        std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " starts exact independent-set search with state budget " << budget_ << std::endl;
      }
      SearchCardinality(chosen);
      return !budgetHit_;
    }

    size_t BestCardinality() const { return bestCardinality_; }
    const std::vector<size_t> &BestSet() const { return bestSet_; }

    // If SolveBudgeted returns true, do this
    // Given that, from SolveBudgeted, we know the best cardinality (the argument "target"), this function collects sets that all have that same size, so later we can randomly pick one at random
    std::vector<std::vector<size_t>> CollectRepresentatives(size_t target, size_t maxReps, size_t budget)
    {
      std::fill(active_.begin(), active_.end(), true);
      target_ = target;
      maxReps_ = (maxReps < 1) ? 1 : maxReps;
      reps_.clear();
      budget_ = budget;
      nodes_ = 0;
      budgetHit_ = false;
      std::vector<size_t> chosen;
      chosen.reserve(target_);
      SearchRepresentatives(chosen);
      return reps_;
    }

    // Fallback when SolveBudgeted hits the budget; number of restarts is controlled by kHeuristicRestarts in the header
    std::vector<size_t> BuildHeuristicSet(size_t restarts, gsl_rng *rng) const
    {
      std::vector<bool> bestInSet(prob_.m, false);
      size_t bestSize = 0;

      // starting point 1: retain the best set found so far in the SolveBudgeted pass
      {
        std::vector<bool> inSet(prob_.m, false);
        for (const size_t v : bestSet_)
        {
          inSet[v] = true;
        }
        PolishInSet(inSet);
        const size_t sz = CountInSet(inSet);
        if (sz > bestSize)
        {
          bestSize = sz;
          bestInSet = inSet;
        }
      }

      // starting point 2: randomized searches
      std::vector<size_t> order(prob_.m);
      std::iota(order.begin(), order.end(), 0);
      for (size_t t = 0; t < restarts; ++t)  // now do the randomize search (kHeuristicRestarts times)
      {
        gsl_ran_shuffle(rng, order.data(), order.size(), sizeof(size_t));
        std::vector<bool> inSet = GreedyFromOrder(order);
        PolishInSet(inSet);
        const size_t sz = CountInSet(inSet);
        if (sz > bestSize)
        {
          bestSize = sz;
          bestInSet = inSet;
        }
      }

      std::vector<size_t> out;
      for (size_t v = 0; v < prob_.m; ++v)
      {
        if (bestInSet[v])
        {
          out.push_back(v);
        }
      }
      return out;
    }

   private:
    size_t CountActive() const
    {
      size_t count = 0;
      for (size_t v = 0; v < prob_.m; ++v)
      {
        count += active_[v];
      }
      return count;
    }

    static size_t CountInSet(const std::vector<bool> &inSet)
    {
      size_t count = 0;
      for (const bool x : inSet)
      {
        count += x;
      }
      return count;
    }

    // Used in BuildHeuristicSet, starting point 2
    // walk the list and keep every class that's still conflict-free with the classes already kept
    std::vector<bool> GreedyFromOrder(const std::vector<size_t> &order) const
    {
      std::vector<bool> inSet(prob_.m, false);
      std::vector<bool> blocked(prob_.m, false);
      for (const size_t v : order)
      {
        if (blocked[v])
        {
          continue;
        }
        inSet[v] = true;
        for (const size_t u : prob_.adj[v])
        {
          blocked[u] = true;
        }
      }
      return inSet;
    }

    // Repeatedly check whether dropping one chosen class can give two classes that are not conflicting
    // If true, drop the chosen class and add the two new classes. Repeat until no such operation exists
    void PolishInSet(std::vector<bool> &inSet) const
    {
      bool improved = true;
      const size_t passLimit = 1000;
      size_t pass = 0;
      while (improved && pass++ < passLimit)
      {
        improved = false;
        for (size_t y = 0; y < prob_.m && !improved; ++y)
        {
          if (!inSet[y])
          {
            continue;
          }

          std::vector<size_t> candidates;  // neighbors of y that are only neighbor to y (so dropping y frees them to be added)
          for (const size_t x : prob_.adj[y])
          {
            if (inSet[x])
            {
              continue;
            }
            size_t inSetNeighbors = 0;
            size_t onlyNeighbor = 0;
            for (const size_t u : prob_.adj[x])
            {
              if (inSet[u])
              {
                ++inSetNeighbors;
                onlyNeighbor = u;
              }
            }
            if (inSetNeighbors == 1 && onlyNeighbor == y)
            {
              candidates.push_back(x);
            }
          }

          // any two mutually non-adjacent candidates give a net +1
          for (size_t a = 0; a < candidates.size() && !improved; ++a)
          {
            for (size_t b = a + 1; b < candidates.size(); ++b)
            {
              bool adjacent = false;
              for (const size_t u : prob_.adj[candidates[a]])
              {
                if (u == candidates[b])
                {
                  adjacent = true;
                  break;
                }
              }
              if (!adjacent)
              {
                inSet[y] = false;
                inSet[candidates[a]] = true;
                inSet[candidates[b]] = true;
                improved = true;
                break;
              }
            }
          }
        }
      }
    }

    // Solve the clique-cover problem of a graph
    // Note that this is not an exact solution, but an approximation (clique cover does not have an exact solution in graph theory)
    size_t CliqueCoverBound()
    {
      std::fill(scratchKeyCount_.begin(), scratchKeyCount_.end(), 0);
      size_t remaining = 0;
      for (size_t v = 0; v < prob_.m; ++v)
      {
        if (!active_[v])
        {
          continue;
        }
        ++remaining;
        covered_[v] = false;
        for (const size_t kid : prob_.keyIds[v])
        {
          ++scratchKeyCount_[kid];
        }
      }
      if (remaining == 0)
      {
        return 0;
      }

      size_t cover = 0;
      while (remaining > 0)
      {
        size_t bestKey = 0;
        size_t bestCount = 0;
        for (size_t kid = 0; kid < prob_.numKeys; ++kid)
        {
          if (scratchKeyCount_[kid] > bestCount)
          {
            bestCount = scratchKeyCount_[kid];
            bestKey = kid;
          }
        }
        if (bestCount == 0)
        {
          break;  // unreachable while remaining > 0, guarded for safety
        }

        ++cover;
        for (const size_t v : prob_.keyToClasses[bestKey])
        {
          if (!active_[v] || covered_[v])
          {
            continue;
          }
          covered_[v] = true;
          --remaining;  // decrement, because it's covered
          for (const size_t kid : prob_.keyIds[v])
          {
            if (scratchKeyCount_[kid] > 0)
            {
              --scratchKeyCount_[kid];
            }
          }
        }
      }

      return cover;
    }

    // Choosing the starting point (seed) for the include-exclude search
    // The vertex (seed) with the highest degree (most conflicts) is chosen
    size_t ChooseBranchVertex() const
    {
      size_t bestVertex = 0;
      size_t bestDegree = 0;
      bool found = false;
      for (size_t v = 0; v < prob_.m; ++v)
      {
        if (!active_[v])
        {
          continue;
        }
        size_t degree = 0;
        for (const size_t neighbor : prob_.adj[v])
        {
          degree += active_[neighbor];
        }
        if (!found || degree > bestDegree)
        {
          found = true;
          bestDegree = degree;
          bestVertex = v;
        }
      }
      return bestVertex;
    }

    void FoldIsolated(std::vector<size_t> &chosen, std::vector<size_t> &folded)
    {
      for (size_t v = 0; v < prob_.m; ++v)
      {
        if (!active_[v])
        {
          continue;
        }
        bool hasActiveNeighbor = false;
        for (const size_t neighbor : prob_.adj[v])
        {
          if (active_[neighbor])
          {
            hasActiveNeighbor = true;
            break;
          }
        }
        if (!hasActiveNeighbor)
        {
          folded.push_back(v);
        }
      }
      for (const size_t v : folded)
      {
        active_[v] = false;
        chosen.push_back(v);
      }
    }

    void Unfold(std::vector<size_t> &chosen, const std::vector<size_t> &folded)
    {
      for (size_t i = folded.size(); i-- > 0;)
      {
        active_[folded[i]] = true;
        chosen.pop_back();
      }
    }

    void Deactivate(size_t vertex, std::vector<size_t> &deactivated)
    {
      active_[vertex] = false;
      deactivated.push_back(vertex);
      for (const size_t neighbor : prob_.adj[vertex])
      {
        if (active_[neighbor])
        {
          active_[neighbor] = false;
          deactivated.push_back(neighbor);
        }
      }
    }

    void Reactivate(const std::vector<size_t> &deactivated)
    {
      for (const size_t vertex : deactivated)
      {
        active_[vertex] = true;
      }
    }

    void SearchCardinality(std::vector<size_t> &chosen)
    {
      if (budgetHit_)
      {
        return;
      }
      if (++nodes_ > budget_)
      {
        budgetHit_ = true;
        if (debugDetailed_)
        {
          std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " budget hit at node " << nodes_ << " with best size " << bestCardinality_ << std::endl;
        }
        return;
      }

      const size_t activeBeforeFold = CountActive();
      if (debugDetailed_)
      {
        std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " node " << nodes_ << " depth " << debugDepth_ << " enters with chosen=" << chosen.size()
                  << ", active=" << activeBeforeFold << ", best=" << bestCardinality_ << std::endl;
      }

      std::vector<size_t> folded;
      FoldIsolated(chosen, folded);
      if (debugDetailed_ && !folded.empty())
      {
        std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " folded isolated class(es)";
        for (const size_t v : folded)
        {
          std::cout << " " << v;
        }
        std::cout << "; chosen is now " << chosen.size() << std::endl;
      }

      const size_t activeAfterFold = CountActive();
      if (activeAfterFold == 0)
      {
        if (chosen.size() > bestCardinality_)
        {
          bestCardinality_ = chosen.size();
          bestSet_ = chosen;
          if (debugDetailed_)
          {
            std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " updates best set to size " << bestCardinality_ << " at a closed leaf" << std::endl;
          }
        }
        else if (debugDetailed_)
        {
          std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " closed leaf has size " << chosen.size() << " and does not improve best " << bestCardinality_
                    << std::endl;
        }
        Unfold(chosen, folded);
        return;
      }

      const size_t cover = CliqueCoverBound();
      if (debugDetailed_)
      {
        std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " clique-cover upper bound = " << cover << "; chosen + bound = " << chosen.size() + cover
                  << std::endl;
      }
      if (chosen.size() + cover <= bestCardinality_)
      {
        if (debugDetailed_)
        {
          std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " prunes branch because " << chosen.size() + cover << " <= current best " << bestCardinality_
                    << std::endl;
        }
        Unfold(chosen, folded);
        return;
      }
      if (cover == 1)
      {
        if (chosen.size() + 1 > bestCardinality_)
        {
          size_t pick = prob_.m;
          for (size_t v = 0; v < prob_.m; ++v)
          {
            if (active_[v])
            {
              pick = v;
              break;
            }
          }
          bestCardinality_ = chosen.size() + 1;
          bestSet_ = chosen;
          if (pick < prob_.m)
          {
            bestSet_.push_back(pick);
          }
          if (debugDetailed_)
          {
            std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " single-clique shortcut picks class " << pick << " and updates best size to "
                      << bestCardinality_ << std::endl;
          }
        }
        else if (debugDetailed_)
        {
          std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " single-clique shortcut cannot improve best " << bestCardinality_ << std::endl;
        }
        Unfold(chosen, folded);
        return;
      }

      const size_t vertex = ChooseBranchVertex();
      if (debugDetailed_)
      {
        std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " branches on class " << vertex << std::endl;
      }

      std::vector<size_t> deactivated;
      Deactivate(vertex, deactivated);
      if (debugDetailed_)
      {
        std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " include class " << vertex << " and deactivate";
        for (const size_t v : deactivated)
        {
          std::cout << " " << v;
        }
        std::cout << std::endl;
      }
      chosen.push_back(vertex);
      ++debugDepth_;
      SearchCardinality(chosen);
      --debugDepth_;
      chosen.pop_back();
      Reactivate(deactivated);

      active_[vertex] = false;
      if (debugDetailed_)
      {
        std::cout << __func__ << " : " << __LINE__ << " : component " << debugComponentIndex_ << " exclude class " << vertex << " and continue residual search" << std::endl;
      }
      ++debugDepth_;
      SearchCardinality(chosen);
      --debugDepth_;
      active_[vertex] = true;

      Unfold(chosen, folded);
    }

    void Record(const std::vector<size_t> &chosen)
    {
      if (reps_.size() >= maxReps_)
      {
        return;
      }
      std::vector<size_t> sorted = chosen;
      std::sort(sorted.begin(), sorted.end());
      for (const std::vector<size_t> &existing : reps_)
      {
        if (existing == sorted)
        {
          return;
        }
      }
      reps_.push_back(std::move(sorted));
    }

    // This does slightly differently than SearchCardinality
    // SearchCardinality: does a search from scratch to discover the best possible size and track a single best set as it goes
    // SearchRepresentatives: given the size of the best set known, search and collect several different sets of that size, so that later we can randomly pick one of them
    void SearchRepresentatives(std::vector<size_t> &chosen)
    {
      if (budgetHit_ || reps_.size() >= maxReps_)
      {
        return;
      }
      if (++nodes_ > budget_)
      {
        budgetHit_ = true;
        return;
      }

      std::vector<size_t> folded;
      FoldIsolated(chosen, folded);

      if (CountActive() == 0)
      {
        if (chosen.size() == target_)
        {
          Record(chosen);
        }
        Unfold(chosen, folded);
        return;
      }

      const size_t cover = CliqueCoverBound();
      if (chosen.size() + cover < target_)  // cannot reach a maximum packing
      {
        Unfold(chosen, folded);
        return;
      }
      if (cover == 1)
      {
        for (size_t v = 0; v < prob_.m && reps_.size() < maxReps_; ++v)
        {
          if (!active_[v])
          {
            continue;
          }
          chosen.push_back(v);
          Record(chosen);
          chosen.pop_back();
        }
        Unfold(chosen, folded);
        return;
      }

      const size_t vertex = ChooseBranchVertex();

      std::vector<size_t> deactivated;
      Deactivate(vertex, deactivated);
      chosen.push_back(vertex);
      SearchRepresentatives(chosen);
      chosen.pop_back();
      Reactivate(deactivated);

      if (reps_.size() < maxReps_)
      {
        active_[vertex] = false;
        SearchRepresentatives(chosen);
        active_[vertex] = true;
      }

      Unfold(chosen, folded);
    }

    const ComponentProblem &prob_;
    std::vector<bool> active_;
    std::vector<size_t> scratchKeyCount_;
    std::vector<bool> covered_;
    size_t bestCardinality_ = 0;
    std::vector<size_t> bestSet_;
    size_t target_ = 0;
    size_t maxReps_ = 1;
    std::vector<std::vector<size_t>> reps_;
    size_t budget_ = 0;
    size_t nodes_ = 0;
    bool budgetHit_ = false;
    bool debugDetailed_ = false;
    size_t debugComponentIndex_ = 0;
    size_t debugDepth_ = 0;
  };

  std::vector<size_t> SelectLargestMvtxDisjointSubset(       //
      const std::vector<const SeedRecord *> &completeSeeds,  //
      gsl_rng *rng,                                          //
      size_t maxRepresentatives,                             //
      size_t searchBudget,                                   //
      size_t heuristicRestarts,                              //
      size_t debugMinimumGroupSize,                          //
      std::vector<size_t> &uncertifiedOut,                   //
      size_t ambiguityGroupSize,
      bool debugVerbosity,
      bool debugDetailedVerbosity)
  {
    if (completeSeeds.empty())
    {
      return {};
    }

    const bool debug = (debugVerbosity || debugDetailedVerbosity) && ambiguityGroupSize > debugMinimumGroupSize;

    const std::vector<TripletClass> classes = CollapseDuplicateTriplets(completeSeeds, debugVerbosity);
    const ConflictGraph classConflict = BuildClassConflictGraph(classes);
    if (debug)
    {
      std::cout << __func__ << " : " << __LINE__ << " : collapsed " << completeSeeds.size() << " complete MVTX seeds into " << classes.size() << " unique MVTX-triplet classes" << std::endl;
    }

    std::vector<size_t> selectedSeedIndices;

    const std::vector<std::vector<size_t>> components = ConnectedComponents(classConflict);
    if (debug)
    {
      std::cout << __func__ << " : " << __LINE__ << " : class conflict graph split into " << components.size() << " connected component(s)" << std::endl;
    }

    size_t componentIndex = 0;
    for (const std::vector<size_t> &component : components)
    {
      const ComponentProblem prob = BuildComponentProblem(classes, classConflict, component);
      ComponentSolver solver(prob, debugDetailedVerbosity && ambiguityGroupSize > debugMinimumGroupSize, componentIndex);

      if (debug)
      {
        std::cout << __func__ << " : " << __LINE__ << " : component " << componentIndex << " has " << prob.m << " class(es) and " << prob.numKeys << " MVTX key(s)" << std::endl;
      }

      const bool certified = solver.SolveBudgeted(searchBudget);
      if (debug)
      {
        std::cout << __func__ << " : " << __LINE__ << " : component " << componentIndex << " exact search " << (certified ? "closed" : "hit the state budget")
                  << "; best independent-set size = " << solver.BestCardinality() << std::endl;
      }

      std::vector<size_t> chosenLocal;
      if (certified)  // search completed, collect equally maximal sets and draw one at random
      {
        const size_t target = solver.BestCardinality();
        const std::vector<std::vector<size_t>> reps = solver.CollectRepresentatives(target, maxRepresentatives, searchBudget);
        if (debug)
        {
          std::cout << __func__ << " : " << __LINE__ << " : component " << componentIndex << " collected " << reps.size() << " representative maximum set(s) for target size " << target
                    << std::endl;
        }

        if (reps.empty())
        {
          chosenLocal = solver.BestSet();
        }
        else
        {
          size_t chosenRep = 0;
          if (reps.size() > 1)
          {
            std::vector<double> weights;
            weights.reserve(reps.size());
            for (const std::vector<size_t> &rep : reps)
            {
              double weight = 1.0;
              for (const size_t local : rep)
              {
                weight *= prob.memberCount[local];
              }
              weights.push_back(weight);
            }
            gsl_ran_discrete_t *repPick = gsl_ran_discrete_preproc(weights.size(), weights.data());
            chosenRep = gsl_ran_discrete(rng, repPick);
            gsl_ran_discrete_free(repPick);  // GNU Scientific Library function that deallocates the memory
          }
          chosenLocal = reps[chosenRep];
        }
      }
      else  // search hit budget, do heuristic search
      {
        chosenLocal = solver.BuildHeuristicSet(heuristicRestarts, rng);
        if (debug)
        {
          std::cout << __func__ << " : " << __LINE__ << " : component " << componentIndex << " heuristic polished set size = " << chosenLocal.size() << std::endl;
        }
      }

      for (const size_t local : chosenLocal)
      {
        const std::vector<size_t> &members = classes[prob.globalClassId[local]].memberSeedIndices;
        const size_t seedIndex = members[gsl_rng_uniform_int(rng, members.size())];
        selectedSeedIndices.push_back(seedIndex);
        if (!certified)
        {
          uncertifiedOut.push_back(seedIndex);
        }
      }
      if (debug)
      {
        std::cout << __func__ << " : " << __LINE__ << " : component " << componentIndex << " selected " << chosenLocal.size() << " seed(s)" << std::endl;
      }
      ++componentIndex;
    }

    return selectedSeedIndices;
  }

  std::vector<std::vector<size_t>> BuildAmbiguityGroups(const std::vector<SeedRecord> &seeds)
  {
    DisjointSet disjoint_set(seeds.size());
    std::unordered_map<TrkrDefs::cluskey, size_t> first_seed_by_cluster_key;

    for (size_t iseed = 0; iseed < seeds.size(); ++iseed)
    {
      for (const TrkrDefs::cluskey key : seeds[iseed].allClusterKeys)
      {
        const auto inserted = first_seed_by_cluster_key.emplace(key, iseed);
        if (!inserted.second)
        {
          disjoint_set.Union(iseed, inserted.first->second);
        }
      }
    }

    std::unordered_map<size_t, std::vector<size_t>> groups_by_root;
    for (size_t iseed = 0; iseed < seeds.size(); ++iseed)
    {
      groups_by_root[disjoint_set.Find(iseed)].push_back(iseed);
    }

    std::vector<std::vector<size_t>> groups;
    groups.reserve(groups_by_root.size());
    for (auto &entry : groups_by_root)
    {
      groups.push_back(entry.second);
    }

    return groups;
  }

  // randomly select seeds from a group if multiple choices are available
  std::vector<size_t> SelectRandomSeedFromGroup(  //
      const std::vector<SeedRecord> &seeds,       //
      const std::vector<size_t> &group,           //
      gsl_rng *rng                                //
  )
  {
    const size_t selected = gsl_rng_uniform_int(rng, group.size());
    return {seeds[group[selected]].index};
  }

  std::vector<size_t> SelectAmbiguityGroup(  //
      const std::vector<SeedRecord> &seeds,  //
      const std::vector<size_t> &group,      //
      gsl_rng *rng,                          //
      size_t maxRepresentatives,             //
      size_t searchBudget,                   //
      size_t heuristicRestarts,              //
      size_t debugMinimumGroupSize,          //
      std::vector<size_t> &uncertifiedOut,
      bool debugVerbosity,
      bool debugDetailedVerbosity)
  {
    const bool debug = (debugVerbosity || debugDetailedVerbosity) && group.size() > debugMinimumGroupSize;

    if (debug)
    {
      std::cout << __func__ << " : " << __LINE__ << " : ambiguity group size = " << group.size() << std::endl;
    }

    if (group.size() == 1)
    {
      return {seeds[group.front()].index};
    }

    std::vector<const SeedRecord *> complete_seeds;
    complete_seeds.reserve(group.size());

    bool all_complete = true;
    for (const size_t local_seed_index : group)
    {
      const SeedRecord &seed = seeds[local_seed_index];
      if (seed.completeMvtx)
      {
        complete_seeds.push_back(&seed);
      }
      else
      {
        all_complete = false;
      }
    }

    if (debug)
    {
      std::cout << __func__ << " : " << __LINE__ << " : complete MVTX seeds in group = " << complete_seeds.size() << "/" << group.size() << std::endl;
      // print out all seed cluster keys for debugging if debugdetail is enabled
      if (debugDetailedVerbosity)
      {
        for (const size_t local_seed_index : group)
        {
          const SeedRecord &seed = seeds[local_seed_index];
          std::cout << __func__ << " : " << __LINE__ << " : seed index " << seed.index << " completeMvtx=" << seed.completeMvtx << " cluster keys: ";
          for (const TrkrDefs::cluskey key : seed.allClusterKeys)
          {
            std::cout << key << " ";
          }
          std::cout << std::endl;
        }
      }
    }

    if (all_complete)
    {
      bool all_same_triplet = true;
      const SeedRecord &first_seed = seeds[group.front()];
      for (const size_t local_seed_index : group)
      {
        if (first_seed.mvtxClusterKeys != seeds[local_seed_index].mvtxClusterKeys)
        {
          all_same_triplet = false;
          break;
        }
      }

      if (all_same_triplet)
      {
        if (debug)
        {
          std::cout << __func__ << " : " << __LINE__ << " : all complete seeds share one MVTX triplet; selecting one random representative" << std::endl;
        }
        return SelectRandomSeedFromGroup(seeds, group, rng);
      }

      if (debug)
      {
        std::cout << __func__ << " : " << __LINE__ << " : all complete seeds have distinct MVTX triplets; selecting largest disjoint subset" << std::endl;
      }
      return SelectLargestMvtxDisjointSubset(complete_seeds, rng, maxRepresentatives, searchBudget, heuristicRestarts, debugMinimumGroupSize, uncertifiedOut, group.size(), debugVerbosity, debugDetailedVerbosity);
    }

    if (complete_seeds.empty())
    {
      if (debug)
      {
        std::cout << __func__ << " : " << __LINE__ << " : no complete MVTX seeds in group; selecting one random representative" << std::endl;
      }
      return SelectRandomSeedFromGroup(seeds, group, rng);
    }

    return SelectLargestMvtxDisjointSubset(complete_seeds, rng, maxRepresentatives, searchBudget, heuristicRestarts, debugMinimumGroupSize, uncertifiedOut, group.size(), debugVerbosity, debugDetailedVerbosity);
  }

  Result SelectPrunedSiliconSeeds(
      const std::vector<SeedRecord> &seeds,
      gsl_rng *rng,
      size_t maxRepresentatives,
      size_t searchBudget,
      size_t heuristicRestarts,
      size_t debugMinimumGroupSize,
      bool debugVerbosity,
      bool debugDetailedVerbosity)
  {
    Result result;
    if (seeds.empty())
    {
      return result;
    }

    const std::vector<std::vector<size_t>> groups = BuildAmbiguityGroups(seeds);
    for (const auto &group : groups)
    {
      const std::vector<size_t> selected_group_indices =
          SelectAmbiguityGroup(seeds, group, rng, maxRepresentatives, searchBudget, heuristicRestarts, debugMinimumGroupSize, result.uncertifiedSeedIndices, debugVerbosity, debugDetailedVerbosity);
      result.selectedSeedIndices.insert(
          result.selectedSeedIndices.end(),
          selected_group_indices.begin(),
          selected_group_indices.end());
    }

    std::sort(result.selectedSeedIndices.begin(), result.selectedSeedIndices.end());
    std::sort(result.uncertifiedSeedIndices.begin(), result.uncertifiedSeedIndices.end());
    return result;
  }

  Result SelectSeeds(
      const TrackSeedContainer &container,
      const std::vector<std::size_t> &seedIndices,
      gsl_rng *rng,
      std::size_t mvtxLayerCount,
      std::size_t maxRepresentatives,
      std::size_t searchBudget,
      std::size_t heuristicRestarts,
      std::size_t debugMinimumGroupSize,
      bool debugVerbosity,
      bool debugDetailedVerbosity)
  {
    std::vector<SeedRecord> seeds;
    seeds.reserve(seedIndices.size());
    for (const std::size_t seedIndex : seedIndices)
    {
      seeds.push_back(MakeSeedRecord(seedIndex, *container.get(seedIndex), mvtxLayerCount));
    }

    return SelectPrunedSiliconSeeds(
        seeds,
        rng,
        maxRepresentatives,
        searchBudget,
        heuristicRestarts,
        debugMinimumGroupSize,
        debugVerbosity,
        debugDetailedVerbosity);
  }
}  // namespace PHSiliconSeedPrunerHelper
