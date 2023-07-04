#include "common.hpp"
#include "graph.hpp"

#include <cstdio>
#include <cstdlib>
#include <utility>
#include <queue>

void Graph::addEdge(int u, int v){
  if(u > v) std::swap(u, v);
  e.emplace_back(u, v);
  etopos.emplace(std::make_pair(u, v), m);
  n = std::max(n, u);
  n = std::max(n, v);
  ++m;
}

bool Graph::readfromFile(const char *filename){
  FILE *fp;
  if((fp = fopen(filename, "r")) == NULL){
    return false;
  }
  m = n = 0;
  int u, v;
  while(fscanf(fp, "%d%d", &u, &v) != EOF){
    addEdge(u, v);
  }
  fclose(fp);
  e.shrink_to_fit();
  return true;
}

void Graph::buildFrontiers(){
  fros.resize(m+1);
  mfros.resize(m);
  vpos.resize(m);
  fent.resize(m);
  flve.resize(m);
  std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> que;
  for(size_t i=0; i<n; ++i) que.push(i);
  
  std::vector<int> deg1(n+1);
  for(size_t i=0; i<m; ++i){
    ++deg1[e[i].first];
    ++deg1[e[i].second];
  }
  std::vector<int> deg2(n+1);
  fros[0].resize(n);
  for(size_t i=0; i<m; ++i){
    mfros[i].resize(n);
    std::copy(fros[i].begin(), fros[i].end(), mfros[i].begin());
    if(!deg2[e[i].first]){
      mfros[i][que.top()] = e[i].first;
      fent[i].emplace_back(que.top());
      que.pop();
    }
    if(!deg2[e[i].second]){
      mfros[i][que.top()] = e[i].second;
      fent[i].emplace_back(que.top());
      que.pop();
    }
    for(size_t j=0; j<mfros[i].size(); ++j){
      if(mfros[i][j] == e[i].first)  vpos[i].first  = j;
      if(mfros[i][j] == e[i].second) vpos[i].second = j;
    }
    ++deg2[e[i].first];
    ++deg2[e[i].second];
    
    fros[i+1].resize(n);
    std::copy(mfros[i].begin(), mfros[i].end(), fros[i+1].begin());
    --deg1[e[i].first];
    --deg1[e[i].second];
    if(!deg1[e[i].first]){
      fros[i+1][vpos[i].first] = 0;
      flve[i].emplace_back(vpos[i].first);
      que.push(vpos[i].first);
    }
    if(!deg1[e[i].second]){
      fros[i+1][vpos[i].second] = 0;
      flve[i].emplace_back(vpos[i].second);
      que.push(vpos[i].second);
    }
  }
  
  for(auto&& v : fros){
    while(!v.empty() && !v.back()) v.pop_back();
    v.shrink_to_fit();
  }
  for(auto&& v : mfros){
    while(!v.empty() && !v.back()) v.pop_back();
    v.shrink_to_fit();
  }
}

int Graph::maxFWidth() const{
  int res = 0;
  for(const auto& vec : fros){
    int tmp = 0;
    for(const auto& ele : vec){
      tmp += ele ? 1 : 0;
    }
    res = std::max(res, tmp);
  }
  return res;
}

void Graph::printFWidth() const{
  for(const auto& vec : fros){
    int tmp = 0;
    for(const auto& ele : vec){
      tmp += ele ? 1 : 0;
    }
    printf("%d ", tmp);
  }
  puts("");
}
