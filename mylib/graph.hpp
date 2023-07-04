#ifndef REL_GRAPH_HPP
#define REL_GRAPH_HPP

#include "common.hpp"

#include <unordered_map>
#include <vector>
#include <utility>
#include <cassert>

using Edge = std::pair<int, int>;

class Graph{
public:
  int n;
  int m;
  std::vector<Edge> e;
  std::unordered_map<std::pair<int, int>, int, HashPI> etopos;
  
  std::vector<std::vector<int>> fros;
  std::vector<std::vector<int>> mfros;
  std::vector<std::pair<size_t, size_t>> vpos;
  std::vector<std::vector<size_t>> fent;
  std::vector<std::vector<size_t>> flve;
  
  Graph(int _n): n(_n), m(0) {};
  Graph(): n(0), m(0) {};
  
  int numE() const{
    return m;
  }
  int numV() const{
    return n;
  }
  int maxFWidth() const;
  void printFWidth() const;
  int etovar(int u, int v) const{
    if(u > v) std::swap(u, v);
    auto key = std::make_pair(u, v);
    auto it = etopos.find(key);
    if(it == etopos.end()) return -1;
    else return it->second;
  }
  void addEdge(int u, int v);
  
  bool readfromFile(const char *filename);
  void buildFrontiers();
};

#endif // REL_GRAPH_HPP
