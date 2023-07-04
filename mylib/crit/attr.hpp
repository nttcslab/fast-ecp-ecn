#ifndef REL_CRIT_ATTR_HPP
#define REL_CRIT_ATTR_HPP

#include "common.hpp"
#include "graph.hpp"
#include "graphsimplify.hpp"
#include "lwdp.hpp"

#include <vector>
#include <utility>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <cstdio>
#include <cstdlib>

class ATTRSolver{
public:  
  ATTRSolver(){};
  
  Real solve(const Graph& G, const std::vector<Real>& pi, const std::vector<Real>& wgt, std::vector<Real>& grad, bool reordering);
  Real solveND(const Graph& G, const std::vector<Real>& pi, const std::vector<Real>& wgt, bool reordering);
  
private:
  mutable adept::Stack stack;
};

#endif // REL_CRIT_ATTR_HPP
