#include "mylib/common.hpp"
#include "mylib/graph.hpp"
#include "mylib/graphsimplify.hpp"
#include "mylib/lwdp.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [order_file] [source_file] <weight_file>\n", fil);
}

int main(int argc, char **argv){
  if(argc < 5){
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
  
  auto cstart = std::chrono::system_clock::now();
  
  LWDP LWSolver;
  LWSolver.solve(G, pi, wgt, srcs, res, reordering);
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  for(size_t i=1; i<=n; ++i){
    printf("%3zu : %.15lf\n", i, res[i]);
  }
  
  fprintf(stderr, "calc time: %.6lf ms\n", ctime);
  fprintf(stderr, "#(states): %zu\n", LWSolver.snum);
  
  return 0;
}