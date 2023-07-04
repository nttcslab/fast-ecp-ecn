#ifndef REL_GRAPHSIMPLIFY_HPP
#define REL_GRAPHSIMPLIFY_HPP

#include "common.hpp"
#include "graph.hpp"

#include <vector>
#include <utility>
#include <cassert>
#include <tuple>
#include <unordered_set>

using History = std::tuple<int, int, double>;

class GraphSimplify{
public:
  Graph Gnew;
  std::unordered_set<int> srcsnew;
  int tgtnew;
  std::vector<double> pinew;
  
  std::vector<int> oldtonew;
  std::vector<History> hists;
  
  GraphSimplify(){};
  
  void HHSimplify(const Graph& G, const std::vector<double>& prob, const std::unordered_set<int>& srcs, int tgt);
  void LWSimplify(const Graph& G, const std::vector<double>& prob, const std::unordered_set<int>& srcs);
};

#endif // REL_GRAPHSIMPLIFY_HPP
