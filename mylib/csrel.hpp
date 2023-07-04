#ifndef REL_CSREL_HPP
#define REL_CSREL_HPP

#include "common.hpp"
#include "graph.hpp"

#include <vector>
#include <utility>
#include <cassert>
#include <unordered_set>
#include <unordered_map>

using StateCS = std::pair<uint64_t, std::vector<int8_t>>;

namespace std{
  template <>
  struct hash<StateCS>{
    public:
    uint64_t operator()(const StateCS& s) const{
      uint64_t h = FNV_OFFSET_BASIS_64;
      h = FNV_PRIME_64* h ^ s.first;
      for(const auto& v : s.second) h = FNV_PRIME_64* h ^ v;
      return h;
    }
  };
}

class DPBlockCS{
public:
  StateCS s;
  int level;
  double p = 0.0;
  std::vector<double> q;
  int8_t cnum;
  size_t lo;
  size_t hi;
  std::vector<int8_t> vlo;
  std::vector<int8_t> vhi;
  
  DPBlockCS() {};
  DPBlockCS(const StateCS& _s, int _level, int8_t _siz): s(_s), level(_level), cnum(_siz), q(_siz), vlo(_siz), vhi(_siz) {};
  DPBlockCS(StateCS&& _s, int _level, int8_t _siz): s(_s), level(_level), cnum(_siz), q(_siz), vlo(_siz), vhi(_siz) {};
};

class CSREL{
public:
  std::vector<DPBlockCS> dp;
  std::vector<std::unordered_map<StateCS, size_t>> maps;
  std::vector<size_t> secs;
  size_t snum;
  
  CSREL(){};
  
  void solve(const Graph& G, const std::vector<double>& pi, const std::unordered_set<int>& srcs, std::vector<double>& res);
  void allClear(){
    dp.clear();
    secs.clear();
    if(!maps.empty()){
      for(auto&& mp : maps){
        mp.clear();
      }
    }
    maps.clear();
    snum = 0;
  }
};

#endif // REL_CSREL_HPP
