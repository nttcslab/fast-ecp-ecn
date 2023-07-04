#include "mylib/common.hpp"
#include "mylib/graph.hpp"
#include "mylib/graphsimplify.hpp"
#include "mylib/csrel.hpp"
#include "beam_search/beam_search.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [order_file] <weight_file>\n", fil);
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
  
  if(argc >= 5){
    FILE *fp;
    if((fp = fopen(argv[4], "r")) == NULL){
      fprintf(stderr, "ERROR: reading weight file %s failed.\n", argv[4]);
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
  
  auto cstart = std::chrono::system_clock::now();
  
  GraphSimplify GSLW;
  GSLW.LWSimplify(G, pi, srcs);
  Graph& Gn = GSLW.Gnew;
  std::vector<double>& pin = GSLW.pinew;
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
  
  std::vector<std::vector<double>> nrels(nn+1);
  for(int i=1; i<=nn; ++i) nrels[i].resize(nn+1);
  
  Gn.buildFrontiers();
  
  CSREL CSSolver;
  
  for(int u=1; u<nn; ++u){
    std::vector<double> csres;
    std::unordered_set<int> tgts = {u};
    
    CSSolver.solve(Gn, pin, tgts, csres);
    
    nrels[u][u] = 1.0;
    for(int v=u+1; v<=nn; ++v){
      nrels[u][v] = nrels[v][u] = csres[v];
    }
  }
  
  std::vector<std::vector<double>> rels(n+1);
  for(int i=1; i<=n; ++i) rels[i].resize(n+1);
  
  for(int i=1; i<=n; ++i){
    rels[i][i] = 1.0;
    for(int j=i+1; j<=n; ++j){
      double mult = 1.0;
      int u = i;
      int v = j;
      auto hitr = GSLW.hists.begin();
      while(!GSLW.oldtonew[u] || !GSLW.oldtonew[v]){
        if(std::get<0>(*hitr) == u){
          u = std::get<1>(*hitr);
          mult *= std::get<2>(*hitr);
        }else if(std::get<0>(*hitr) == v){
          v = std::get<1>(*hitr);
          mult *= std::get<2>(*hitr);
        }
        if(u == v) break;
        ++hitr;
      }
      if(u == v) rels[i][j] = rels[j][i] = mult;
      else       rels[i][j] = rels[j][i] = mult * nrels[GSLW.oldtonew[u]][GSLW.oldtonew[v]];
    }
  }
  
  res.resize(n+1);
  for(int i=1; i<=n; ++i){
    for(int j=1; j<=n; ++j)res[i] += wgt[j] * rels[i][j];
  }
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  for(size_t i=1; i<=n; ++i){
    printf("%3zu : %.15lf\n", i, res[i]);
  }
  
  fprintf(stderr, "calc time: %.6lf ms\n", ctime);
  
  return 0;
}