#include "mylib/crit/common.hpp"
#include "mylib/crit/graph.hpp"
#include "mylib/crit/graphsimplify.hpp"
#include "beam_search/beam_search.hpp"

#include "tdzdd/DdSpec.hpp"
#include "tdzdd/DdStructure.hpp"
#include "tdzdd/DdEval.hpp"
#include "tdzdd/spec/FrontierBasedSearch.hpp"
#include "tdzdd/util/Graph.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

template <typename xReal>
class ProbEval : public tdzdd::DdEval<ProbEval<xReal>, xReal> {
private:
  std::vector<xReal> prob_list_;
  
public:
  ProbEval(const std::vector<xReal>& prob_list) : prob_list_(prob_list) {}
  
  void evalTerminal(xReal& p, bool one) const { p = one ? 1.0 : 0.0; }
  
  void evalNode(xReal& p, int level,
                tdzdd::DdValues<xReal, 2> const& values) const {
    xReal pc = prob_list_[prob_list_.size() - level];
    p = values.get(0) * (1 - pc) + values.get(1) * pc;
  }
};

class ALLRELSolver{
public:
  ALLRELSolver(){};
  
  Real solve(const Graph& G, const std::vector<Real>& pi, std::vector<Real>& grad, bool reordering);
  Real solveND(const Graph& G, const std::vector<Real>& pi, bool reordering);
  
private:
  mutable adept::Stack stack;
};

Real ALLRELSolver::solve(const Graph& G, const std::vector<Real>& pi, std::vector<Real>& grad, bool reordering){
  int n = G.numV();
  int m = G.numE();
  
  std::vector<aReal> activepi(m);
  adept::set_values(&activepi[0], m, &pi[0]);
  
  stack.new_recording();
  
  GraphSimplify<aReal> GSLW;
  std::unordered_set<int> emptysrc;
  GSLW.LWSimplify(G, activepi, emptysrc);
  Graph& Gn = GSLW.Gnew;
  std::vector<aReal>& pin = GSLW.pinew;
  int nn = Gn.numV();
  int mn = Gn.numE();
  
  Graph reoG;
  std::vector<aReal> reopi;
  
  if(reordering){
    std::vector<Edge> beforder;
    beforder.reserve(mn);
    for(const auto& edg : Gn.e){
      beforder.emplace_back(edg.first-1, edg.second-1);
    }
    std::vector<Edge> aftorder = ordering(nn, beforder);
    for(const auto& edg : aftorder){
      reoG.addEdge(edg.first+1, edg.second+1);
    }
    reopi.resize(mn);
    for(int i=0; i<mn; ++i){
      reopi[i] = pin[Gn.etovar(reoG.e[i].first, reoG.e[i].second)];
    }
    Gn = reoG;
    pin = reopi;
  }
  
  aReal res;
  
  tdzdd::Graph tG;
  for(const auto& edg : Gn.e){
    tG.addEdge(std::to_string(edg.first), std::to_string(edg.second));
  }
  for(int i=1; i<=nn; ++i){
    tG.setColor(std::to_string(i), 1);
  }
  tG.update();
  tdzdd::FrontierBasedSearch fbs(tG, -1, false, false);
  
  tdzdd::DdStructure<2> dd(fbs);
  dd.useMultiProcessors(false);
  res = dd.evaluate(ProbEval<aReal>(pin));
  
  for(const auto& ele : GSLW.hists){
    res *= std::get<2>(ele);
  }
  
  res.set_gradient(1.0);
  stack.compute_adjoint();
  grad.resize(m);
  adept::get_gradients(&activepi[0], m, &grad[0]);
  
  return adept::value(res);
}

Real ALLRELSolver::solveND(const Graph& G, const std::vector<Real>& pi, bool reordering){
  GraphSimplify<Real> GSLW;
  std::unordered_set<int> emptysrc;
  GSLW.LWSimplify(G, pi, emptysrc);
  Graph& Gn = GSLW.Gnew;
  std::vector<Real>& pin = GSLW.pinew;
  int nn = Gn.numV();
  int mn = Gn.numE();
  
  Graph reoG;
  std::vector<Real> reopi;
  
  if(reordering){
    std::vector<Edge> beforder;
    beforder.reserve(mn);
    for(const auto& edg : Gn.e){
      beforder.emplace_back(edg.first-1, edg.second-1);
    }
    std::vector<Edge> aftorder = ordering(nn, beforder);
    for(const auto& edg : aftorder){
      reoG.addEdge(edg.first+1, edg.second+1);
    }
    reopi.resize(mn);
    for(int i=0; i<mn; ++i){
      reopi[i] = pin[Gn.etovar(reoG.e[i].first, reoG.e[i].second)];
    }
    Gn = reoG;
    pin = reopi;
  }
  
  Real res;
  
  tdzdd::Graph tG;
  for(const auto& edg : Gn.e){
    tG.addEdge(std::to_string(edg.first), std::to_string(edg.second));
  }
  for(int i=1; i<=nn; ++i){
    tG.setColor(std::to_string(i), 1);
  }
  tG.update();
  tdzdd::FrontierBasedSearch fbs(tG, -1, false, false);
  
  tdzdd::DdStructure<2> dd(fbs);
  dd.useMultiProcessors(false);
  res = dd.evaluate(ProbEval<Real>(pin));
  
  for(const auto& ele : GSLW.hists){
    res *= std::get<2>(ele);
  }
  return res;
}

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [order_file] [mode]\n", fil);
}

int main(int argc, char **argv){
  if(argc < 5){
    fprintf(stderr, "ERROR: too few arguments.\n");
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  Graph G;
  int n, m;
  int mode = 0;
  bool reordering = false;
  std::vector<Real> pi;
  
  {
    Graph H;
    if(!H.readfromFile(argv[1])){
      fprintf(stderr, "ERROR: reading graph file %s failed.\n", argv[1]);
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
    
    n = H.numV();
    m = H.numE();
    
    std::vector<Real> prob(m);
    pi.resize(m);
    {
      FILE *fp;
      if((fp = fopen(argv[2], "r")) == NULL){
        fprintf(stderr, "ERROR: reading probability file %s failed.\n", argv[2]);
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
      }
      
      for(size_t i=0; i<m; ++i){
        fscanf(fp, "%lf", &prob[i]);
      }
      fclose(fp);
    }
    
    if(argv[3][0] == '!'){
      reordering = true;
      for(size_t i=0; i<m; ++i){
        G.addEdge(H.e[i].first, H.e[i].second);
        pi[i] = prob[i];
      }
    }else{
      if(!G.readfromFile(argv[3])){
        fprintf(stderr, "ERROR: reading order file %s failed.\n", argv[3]);
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
      }
    
      for(size_t i=0; i<m; ++i){
        pi[i] = prob[H.etovar(G.e[i].first, G.e[i].second)];
      }
    }
  }
  mode = atoi(argv[4]);
  
  Real res;
  std::vector<Real> grad;
  
  auto cstart = std::chrono::system_clock::now();
  
  ALLRELSolver ALLRELMac;
  res = ALLRELMac.solve(G, pi, grad, reordering);
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  fprintf(stderr, "ALLREL = %.15lf\n", res);
  for(size_t i=0; i<m; ++i){
    Real val;
    switch(mode){
    case 1:
      val = 1.0 - (res - grad[i] * pi[i]);
      break;
    case 2:
      val = res + grad[i] * (1.0 - pi[i]);
      break;
    case 3:
      val = pi[i] * (res + grad[i] * (1.0 - pi[i])) / res;
      break;
    default:
      val = grad[i];
      break;
    }
    printf("%.15lf\n", val);
  }
  
  fprintf(stderr, "calc time: %.6lf ms\n", ctime);
  
  return 0;
}