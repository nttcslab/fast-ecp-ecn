#ifndef REL_LWDP_HPP
#define REL_LWDP_HPP

#include "common.hpp"
#include "graph.hpp"
#include "graphsimplify.hpp"

#include <vector>
#include <utility>
#include <cassert>
#include <unordered_set>
#include <unordered_map>

using State = std::vector<int8_t>;

namespace std{
  template <>
  struct hash<State>{
    public:
    uint64_t operator()(const State& s) const{
      uint64_t h = FNV_OFFSET_BASIS_64;
      for(const auto& v : s) h = FNV_PRIME_64* h ^ v;
      return h;
    }
  };
}

class DPBlock{
public:
  State s;
  int level;
  double p = 0.0;
  std::vector<double> c;
  std::vector<std::vector<double>> q;
  std::vector<double> d;
  int8_t cnum;
  size_t lo;
  size_t hi;
  std::vector<int8_t> vlo;
  std::vector<int8_t> vhi;
  
  DPBlock() {};
  DPBlock(const State& _s, int _level, int8_t _siz): s(_s), level(_level), cnum(_siz), c(_siz), q(_siz), d(_siz), vlo(_siz), vhi(_siz) {
    for(auto&& v : q) v.resize(_siz);
  };
  DPBlock(State&& _s, int _level, int8_t _siz): s(_s), level(_level), cnum(_siz), c(_siz), q(_siz), d(_siz), vlo(_siz), vhi(_siz) {
    for(auto&& v : q) v.resize(_siz);
  };
};

class LWDP{
public:
  std::vector<DPBlock> dp;
  std::vector<std::unordered_map<State, size_t>> maps;
  std::vector<size_t> secs;
  std::vector<std::pair<int, int8_t>> vtoipos;
  size_t snum;
  
  LWDP(){};
  
  void solve(const Graph& G, const std::vector<double>& pi, const std::vector<double>& wgt, const std::unordered_set<int>& srcs, std::vector<double>& res, bool reordering);
  void innerSolve(const Graph& G, const std::vector<double>& pi, const std::vector<double>& wgt, const std::unordered_set<int>& srcs, std::vector<double>& res);
  void innerSolveEmpty(const Graph& G, const std::vector<double>& pi, const std::vector<double>& wgt, std::vector<double>& res);
  void allClear(){
    dp.clear();
    secs.clear();
    snum = 0;
  }
};

#endif // REL_LWDP_HPP
