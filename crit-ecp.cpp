#include "mylib/crit/common.hpp"
#include "mylib/crit/graph.hpp"
#include "mylib/crit/attr.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [order_file] [mode] <weight_file>\n", fil);
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
  std::vector<Real> wgt;
  
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
  mode = atoi(argv[4]);
  
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
  
  Real res;
  std::vector<Real> grad;
  
  auto cstart = std::chrono::system_clock::now();
  
  ATTRSolver ATTRMac;
  res = ATTRMac.solve(G, pi, wgt, grad, reordering);
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  fprintf(stderr, "ATTR = %.15lf\n", res);
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