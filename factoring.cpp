#include "mylib/common.hpp"
#include "mylib/graph.hpp"
#include "mylib/graphsimplify.hpp"
#include "beam_search/beam_search.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <chrono>
#include <queue>
#include <unordered_map>
#include <unordered_set>

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [source_file] <weight_file>\n", fil);
}

class FactorSolver{
public:
  std::vector<std::vector<std::pair<int, int>>> adj;
  std::vector<double> pi;
  std::unordered_set<int> u;
  int v;
  std::vector<int> states;
  double res = 0.0;
  int n;
  int m;
  
  FactorSolver(){};
  
  double solve(const Graph& G, const std::vector<double>& _pi, const std::unordered_set<int> _u, int _v);
  
private:
  void solveInner(double mult);
  bool shortest_path(std::vector<int>& rpath);
};

double FactorSolver::solve(const Graph& G, const std::vector<double>& _pi, const std::unordered_set<int> _u, int _v){
  n = G.numV();
  m = G.numE();
  adj.resize(n+1);
  for(int i=0; i<m; ++i){
    adj[G.e[i].first].emplace_back(i, G.e[i].second);
    adj[G.e[i].second].emplace_back(i, G.e[i].first);
  }
  pi.resize(m);
  std::copy(_pi.begin(), _pi.end(), pi.begin());
  u = _u;
  v = _v;
  states.resize(m);
  
  solveInner(1.0);
  
  return res;
}

void FactorSolver::solveInner(double mult){
  std::vector<int> rpath;
  if(!shortest_path(rpath)){
    return;
  }
  double padd = mult;
  std::vector<int> undecided;
  for(const auto& ele : rpath){
    if(!states[ele]){
      undecided.emplace_back(ele);
      padd *= pi[ele];
    }
  }
  res += padd;
  double mmult = mult;
  for(const auto& ele : undecided){
    states[ele] = -1;
    solveInner(mmult * (1.0 - pi[ele]));
    states[ele] = 1;
    mmult *= pi[ele];
  }
  for(const auto& ele : undecided) states[ele] = 0;
}

bool FactorSolver::shortest_path(std::vector<int>& rpath){
  std::vector<std::pair<int, int>> bef(n+1);
  std::vector<int> visited(n+1);
  std::queue<int> que;
  int finish = 0;
  que.emplace(v);
  visited[v] = 1;
  while(!que.empty()){
    int now = que.front();
    que.pop();
    for(const auto& ele : adj[now]){
      if(states[ele.first] >= 0 && !visited[ele.second]){
        visited[ele.second] = 1;
        bef[ele.second] = std::make_pair(ele.first, now);
        que.emplace(ele.second);
        if(u.count(ele.second)){
          finish = ele.second;
          goto FINISH;
        }
      }
    }
  }
  FINISH:
  if(!finish){
    return false;
  }
  std::vector<int> rrpath;
  while(finish != v){
    rrpath.emplace_back(bef[finish].first);
    finish = bef[finish].second;
  }
  size_t rpath_size = rrpath.size();
  rpath.reserve(rpath_size);
  for(size_t i=0; i<rpath_size; ++i){
    rpath.emplace_back(rrpath[rpath_size - 1 - i]);
  }
  return true;
}


int main(int argc, char **argv){
  if(argc < 4){
    fprintf(stderr, "ERROR: too few arguments.\n");
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  Graph G;
  int n, m;
  std::vector<double> pi;
  std::vector<double> wgt;
  std::unordered_set<int> srcs;
  
  {
    if(!G.readfromFile(argv[1])){
      fprintf(stderr, "ERROR: reading graph file %s failed.\n", argv[1]);
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
    
    n = G.numV();
    m = G.numE();
    
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
        fscanf(fp, "%lf", &pi[i]);
      }
      fclose(fp);
    }
  }
  
  {
    FILE *fp;
    if((fp = fopen(argv[3], "r")) == NULL){
      fprintf(stderr, "ERROR: reading source file %s failed.\n", argv[3]);
      exit(EXIT_FAILURE);
    }
    int src;
    while(fscanf(fp, "%d", &src) != EOF){
      srcs.emplace(src);
    }
    fclose(fp);
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
  auto& srcsn = GSLW.srcsnew;
  int nn = Gn.numV();
  int mn = Gn.numE();
  
  std::vector<std::vector<double>> nrels(nn+1);
  for(int i=0; i<=nn; ++i) nrels[i].resize(nn+1);
  
  for(int u=0; u<=nn; ++u){
    for(int v=1; v<=nn; ++v){
      if(u == v || srcsn.count(v)){
        nrels[u][v] = 1.0;
        continue;
      }
      if(srcsn.count(u)){
        nrels[u][v] = nrels[0][v];
        continue;
      }
      double rel = 1.0;
      std::unordered_set<int> uu(srcsn);
      if(u) uu.emplace(u);
      if(uu.empty()){
        nrels[u][v] = 0.0;
        continue;
      }
      
      GraphSimplify GS;
      GS.HHSimplify(Gn, pin, uu, v);
      Graph& Guv = GS.Gnew;
      std::vector<double>& piuv = GS.pinew;
      auto& uuv = GS.srcsnew;
      int vuv = GS.tgtnew;
      int nuv = Guv.numV();
      int muv = Guv.numE();
      
      if(!uuv.count(vuv)){
        FactorSolver FS;
        rel = FS.solve(Guv, piuv, uuv, vuv);
      }
      for(const auto& ele : GS.hists){
        rel *= std::get<2>(ele);
      }
      nrels[u][v] = rel;
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