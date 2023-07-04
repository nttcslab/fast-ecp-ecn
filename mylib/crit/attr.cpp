#include "common.hpp"
#include "graph.hpp"
#include "graphsimplify.hpp"
#include "../../beam_search/beam_search.hpp"
#include "lwdp.hpp"
#include "attr.hpp"

#include <vector>
#include <utility>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>

Real ATTRSolver::solve(const Graph& G, const std::vector<Real>& pi, const std::vector<Real>& wgt, std::vector<Real>& grad, bool reordering){
  std::unordered_set<int> emptysrc;
  int n = G.numV();
  int m = G.numE();
  
  std::vector<aReal> xres;
  std::vector<aReal> activepi(m);
  adept::set_values(&activepi[0], m, &pi[0]);
  std::vector<aReal> activewgt(n+1);
  adept::set_values(&activewgt[0], n+1, &wgt[0]);
  
  stack.new_recording();
  
  aReal res = 0.0;
  LWDP<aReal> LWSolver;
  LWSolver.solve(G, activepi, activewgt, emptysrc, xres, reordering);
  
  aReal wgtsum = 0.0;
  for(int i=1; i<=n; ++i){
    wgtsum += wgt[i];
  }
  wgtsum *= wgtsum;
  
  for(int i=1; i<=n; ++i){
    res += wgt[i] * (xres[i] - wgt[i]);
    wgtsum -= wgt[i] * wgt[i];
  }
  res /= wgtsum;
  
  res.set_gradient(1.0);
  stack.compute_adjoint();
  grad.resize(m);
  adept::get_gradients(&activepi[0], m, &grad[0]);
  
  return adept::value(res);
}

Real ATTRSolver::solveND(const Graph& G, const std::vector<Real>& pi, const std::vector<Real>& wgt, bool reordering){
  std::unordered_set<int> emptysrc;
  int n = G.numV();
  int m = G.numE();
  
  std::vector<Real> xres;
  LWDP<Real> LWSolver;
  
  Real res = 0.0;
  LWSolver.solve(G, pi, wgt, emptysrc, xres, reordering);
  
  Real wgtsum = 0.0;
  for(int i=1; i<=n; ++i){
    wgtsum += wgt[i];
  }
  wgtsum *= wgtsum;
  
  for(int i=1; i<=n; ++i){
    res += wgt[i] * (xres[i] - wgt[i]);
    wgtsum -= wgt[i] * wgt[i];
  }
  res /= wgtsum;
  
  return res;
}