#include "mylib/common.hpp"
#include "mylib/graph.hpp"
#include "mylib/graphsimplify.hpp"
#include "beam_search/beam_search.hpp"

#include "tdzdd/DdSpec.hpp"
#include "tdzdd/DdStructure.hpp"
#include "tdzdd/DdEval.hpp"
#include "tdzdd/spec/FrontierBasedSearch.hpp"
#include "tdzdd/util/Graph.hpp"

#include "SAPPOROBDD/include/BDD.h"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <chrono>
#include <queue>
#include <unordered_map>
#include <unordered_set>

class ToBDD: public tdzdd::DdEval<ToBDD, BDD> {
  int const offset;
  
public:
  ToBDD(int offset = 0)
  : offset(offset) {
  }
  
  bool isThreadSafe() const {
    return false;
  }
  
  void initialize(int topLevel) const {
    while (BDD_VarUsed() < topLevel + offset) {
      BDD_NewVar();
    }
  }
  
  void evalTerminal(BDD& f, int value) const {
    f = BDD(value);
  }
  
  void evalNode(BDD& f, int level, tdzdd::DdValues<BDD, 2> const& values) const {
    f = values.get(0);
    f &= ~(BDDvar(BDD_VarOfLev(level + offset)));
    if (level + offset > 0) {
      BDD f1 = values.get(1);
      f |= (BDDvar(BDD_VarOfLev(level + offset))) & f1;
    }
  }
};

class BDDRel{
public:
  BDDRel(){};
          
  double compute(BDD f, const std::vector<double> &_wp){
    mc.clear();
    wp.reserve(_wp.size() + 1);
    std::copy(_wp.begin(), _wp.end(), wp.begin()+1);
    mc.reserve(f.Size() * 2);
    int v = f.Top();
    double res = computeInner(f);
    return res;
  }
  
private:
  std::vector<double> wp;
  std::unordered_map<bddword, double> mc;
  
  double computeInner(BDD f)
  {
    // base case: f is terminal
    if (f == BDD(0)) return 0.0;
    if (f == BDD(1)) return 1.0;
    assert(f.Top() != 0);
    
    // get ID, var and level
    BDD h = f;
    int v = h.Top();
    int lev = BDD_LevOfVar(v);
    bddword id = h.GetID();
    bool neg = false;
    // negative edge processing
    if (id & 1){
      h = ~h;
      id = h.GetID();
      neg = true;
    }
    // cache search
    auto it = mc.find(id);
    if (it != mc.end()){
      if (!neg) return it->second;
      else      return 1.0 - it->second;
    }
    
    double res = 0.0;
    // 0-edge processing
    {
      BDD h0 = h.At0(v);
      res = (1.0 - wp[v]) * computeInner(h0);
    }
    // 1-edge processing
    {
      BDD h1 = h.At1(v);
      res += wp[v] * computeInner(h1);
    }
    // cache entry
    mc.emplace(id, res);
    if (!neg) return res;
    else      return 1.0 - res;
  }
};

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [order_file] [source_file] <weight_file>\n", fil);
}

int main(int argc, char **argv){
  if(argc < 4){
    fprintf(stderr, "ERROR: too few arguments.\n");
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  Graph G;
  int n, m;
  bool reordering = false;
  std::vector<double> pi;
  std::vector<double> wgt;
  std::unordered_set<int> srcs;
  
  {
    Graph H;
    if(!H.readfromFile(argv[1])){
      fprintf(stderr, "ERROR: reading graph file %s failed.\n", argv[1]);
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
    
    n = H.numV();
    m = H.numE();
    
    std::vector<double> prob(m);
    pi.resize(m);
    wgt.resize(n+1);
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
  
  {
    FILE *fp;
    if((fp = fopen(argv[4], "r")) == NULL){
      fprintf(stderr, "ERROR: reading source file %s failed.\n", argv[4]);
      exit(EXIT_FAILURE);
    }
    int src;
    while(fscanf(fp, "%d", &src) != EOF){
      srcs.emplace(src);
    }
    fclose(fp);
  }
  
  if(argc >= 6){
    FILE *fp;
    if((fp = fopen(argv[5], "r")) == NULL){
      fprintf(stderr, "ERROR: reading weight file %s failed.\n", argv[5]);
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
    for(size_t i=1; i<=n; ++i){
      fscanf(fp, "%lf", &wgt[i]);
    }
    fclose(fp);
  }else{
    for(size_t i=1; i<=n; ++i){
      wgt[i] = 1.0;
    }
  }
  
  std::vector<double> res;
  
  BDD_Init(1ULL << 24, 1ULL << 32);
  
  auto cstart = std::chrono::system_clock::now();
  
  GraphSimplify GSLW;
  GSLW.LWSimplify(G, pi, srcs);
  Graph& Gn = GSLW.Gnew;
  std::vector<double>& pin = GSLW.pinew;
  auto& srcsn = GSLW.srcsnew;
  int nn = Gn.numV();
  int mn = Gn.numE();
  
  Graph reoG;
  std::vector<double> reopi;
  
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
  
  std::vector<double> pirev;
  pirev.reserve(pin.size());
  for(auto itr=pin.rbegin(); itr!=pin.rend(); ++itr) pirev.emplace_back(*itr);
  
  BDDRel BR;
  
  tdzdd::Graph tG;
  for(const auto& edg : Gn.e){
    tG.addEdge(std::to_string(edg.first), std::to_string(edg.second));
  }
  tG.update();
  
  std::vector<std::vector<double>> nrels(nn+1);
  for(int i=0; i<=nn; ++i) nrels[i].resize(nn+1);
  
  for(int v=1; v<=nn; ++v){
    if(srcsn.count(v)){
      for(int u=0; u<=nn; ++u) nrels[u][v] = 1.0;
      continue;
    }
    BDD BTv;
    for(const auto& src : srcsn){
      tG.clearColors();
      tG.setColor(std::to_string(src), 1);
      tG.setColor(std::to_string(v), 1);
      tG.update();
      
      tdzdd::FrontierBasedSearch fbs(tG, -1, false, false);
      
      tdzdd::DdStructure<2> dd(fbs);
      dd.useMultiProcessors(false);
      
      BDD Bsrcv = dd.evaluate(ToBDD());
      BTv |= Bsrcv;
    }
    nrels[0][v] = BR.compute(BTv, pirev);
    
    for(int u=1; u<=nn; ++u){
      if(u == v){
        nrels[u][v] = 1.0;
        continue;
      }else if(srcsn.count(u)){
        nrels[u][v] = nrels[0][v];
        continue;
      }
      tG.clearColors();
      tG.setColor(std::to_string(u), 1);
      tG.setColor(std::to_string(v), 1);
      tG.update();
      
      tdzdd::FrontierBasedSearch fbs(tG, -1, false, false);
      
      tdzdd::DdStructure<2> dd(fbs);
      dd.useMultiProcessors(false);
      dd.bddReduce();
      
      BDD Buv = dd.evaluate(ToBDD());
      Buv |= BTv;
      nrels[u][v] = BR.compute(Buv, pirev);
    }
  }
  
  std::vector<std::vector<double>> rels(n+1);
  for(int i=0; i<=n; ++i) rels[i].resize(n+1);
  
  for(int j=1; j<=n; ++j){
    double mult = 1.0;
    int v = j;
    auto hitr = GSLW.hists.begin();
    while(!GSLW.oldtonew[v]){
      if(std::get<0>(*hitr) == v){
        v = std::get<1>(*hitr);
        mult *= std::get<2>(*hitr);
      }
      if(srcs.count(v)) break;
      ++hitr;
    }
    if(srcs.count(v)) rels[0][j] = mult;
    else              rels[0][j] = mult * nrels[0][GSLW.oldtonew[v]];
  }
  
  for(int i=1; i<=n; ++i){
    for(int j=1; j<=n; ++j){
      if(i == j || srcs.count(j)){
        rels[i][j] = 1.0;
        continue;
      }
      if(srcs.count(i)){
        rels[i][j] = rels[0][j];
        continue;
      }
      double mult = 1.0;
      double pbase = 0.0;
      int u = i;
      int v = j;
      auto hitr = GSLW.hists.begin();
      while(!GSLW.oldtonew[u] || !GSLW.oldtonew[v]){
        if(std::get<0>(*hitr) == u){
          u = std::get<1>(*hitr);
          pbase += mult * (1.0 - std::get<2>(*hitr)) * rels[0][v];
          mult *= std::get<2>(*hitr);
        }else if(std::get<0>(*hitr) == v){
          v = std::get<1>(*hitr);
          mult *= std::get<2>(*hitr);
        }
        if(u == v || srcs.count(v)) break;
        ++hitr;
      }
      if(u == v || srcs.count(v)) rels[i][j] = mult + pbase;
      else                        rels[i][j] = mult * nrels[GSLW.oldtonew[u]][GSLW.oldtonew[v]] + pbase;
    }
  }
  
  res.resize(n+1);
  for(int i=1; i<=n; ++i){
    for(int j=1; j<=n; ++j){
      res[i] += wgt[j] * rels[i][j];
    }
  }
  puts("");
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  for(size_t i=1; i<=n; ++i){
    printf("%3zu : %.15lf\n", i, res[i]);
  }
  
  fprintf(stderr, "calc time: %.6lf ms\n", ctime);
  
  return 0;
}