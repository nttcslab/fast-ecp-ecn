#ifndef REL_CRIT_LWDP_HPP
#define REL_CRIT_LWDP_HPP

#include "common.hpp"
#include "graph.hpp"
#include "graphsimplify.hpp"
#include "../../beam_search/beam_search.hpp"

#include <vector>
#include <utility>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <cstdio>
#include <cstdlib>

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

template <typename xReal>
class DPBlock{
public:
  State s;
  int level;
  xReal p = 0.0;
  std::vector<xReal> c;
  std::vector<std::vector<xReal>> q;
  std::vector<xReal> d;
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

template <typename xReal>
class LWDP{
public:
  std::vector<DPBlock<xReal>> dp;
  std::vector<std::unordered_map<State, size_t>> maps;
  std::vector<size_t> secs;
  std::vector<std::pair<int, int8_t>> vtoipos;
  size_t snum;
  
  LWDP(){};
  
  void solve(const Graph& G, const std::vector<xReal>& pi, const std::vector<xReal>& wgt, const std::unordered_set<int>& srcs, std::vector<xReal>& res, bool reordering);
  void innerSolve(const Graph& G, const std::vector<xReal>& pi, const std::vector<xReal>& wgt, const std::unordered_set<int>& srcs, std::vector<xReal>& res);
  void innerSolveEmpty(const Graph& G, const std::vector<xReal>& pi, const std::vector<xReal>& wgt, std::vector<xReal>& res);
  void allClear(){
    dp.clear();
    secs.clear();
    snum = 0;
  }
};

template <typename xReal>
void LWDP<xReal>::solve(const Graph& G, const std::vector<xReal>& pi, const std::vector<xReal>& wgt, const std::unordered_set<int>& srcs, std::vector<xReal>& res, bool reordering){
  allClear();
  
  GraphSimplify<xReal> GS;
  GS.LWSimplify(G, pi, srcs);
  
  Graph& Gn = GS.Gnew;
  std::vector<xReal>& pin = GS.pinew;
  int no = G.numV();
  int mo = G.numE();
  int nn = Gn.numV();
  int mn = Gn.numE();
  std::vector<xReal> innerres(nn+1);
  std::vector<xReal> innerwgt0(no+1);
  for(int i=1; i<=no; ++i) innerwgt0[i] = wgt[i];
  std::vector<xReal> innerwgt(nn+1);
  
  for(const auto& ent : GS.hists){
    innerwgt0[std::get<1>(ent)] += innerwgt0[std::get<0>(ent)] * std::get<2>(ent);
  }
  for(int i=1; i<=no; ++i){
    if(GS.oldtonew[i] != 0){
      innerwgt[GS.oldtonew[i]] = innerwgt0[i];
    }
  }
  
  Graph reoG;
  std::vector<xReal> reopi;
  
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
    fprintf(stderr, "Reordering finished.\n");
  }
  
  Gn.buildFrontiers();
  fprintf(stderr, "max_frontier_width = %d\n", Gn.maxFWidth());
  //Gn.printFWidth();
  
  if(!GS.srcsnew.empty()) innerSolve(Gn, GS.pinew, innerwgt, GS.srcsnew, innerres);
  else                    innerSolveEmpty(Gn, GS.pinew, innerwgt, innerres);
  
  res.resize(no+1);
  res[0] = innerres[0];
  for(int i=1; i<=no; ++i){
    if(GS.oldtonew[i] != 0){
      res[i] = innerres[GS.oldtonew[i]];
    }
  }
  
  std::vector<xReal> pbase(no+1);
  std::vector<int> pbase_computed(no+1);
  for(auto itr = GS.hists.rbegin(); itr != GS.hists.rend(); ++itr){
    int x = std::get<0>(*itr);
    int y = std::get<1>(*itr);
    xReal pri = std::get<2>(*itr);
    if(!pbase_computed[y]){
      if(GS.srcsnew.empty()) pbase[y] = 0.0;
      else{
        int newy = GS.oldtonew[y];
        if(GS.srcsnew.count(newy)) pbase[y] = 1.0;
        else{
          int nowi = vtoipos[newy].first;
          int8_t nowpos = vtoipos[newy].second;
          for(size_t k=secs[nowi]; k<secs[nowi+1]; ++k){
            const DPBlock<xReal>& dp_now = dp[k];
            pbase[y] += dp_now.p * dp_now.q[0][dp_now.s[nowpos]];
          }
        }
      }
      pbase_computed[y] = 1;
    }
    pbase[x] = pri * pbase[y];
    pbase_computed[x] = 1;
    
    res[x] = innerwgt0[x] * (1.0 - pri * pri - pri * (1.0 - pri) * pbase[y]) + pri * res[y] + (1.0 - pri) * res[0];
  }
}

template <typename xReal>
void LWDP<xReal>::innerSolve(const Graph& G, const std::vector<xReal>& pi, const std::vector<xReal>& wgt, const std::unordered_set<int>& srcs, std::vector<xReal>& res){
  int n = G.numV();
  int m = G.numE();
  
  maps.resize(m+1);
  secs.resize(m+2);
  
  std::vector<int> issrc(n+1);
  for(const auto& val : srcs){
    issrc[val] = 1;
  }
  
  // id=0: root node
  {
    State root;
    maps[0].emplace(root, 0);
    dp.emplace_back(root, 0, 1);
    dp[0].p = 1.0;
  }
  secs[0] = 0;
  snum = 1;
  
  for(size_t i=0; i<m; ++i){
    const auto& now_fro = G.fros[i];
    const auto& med_fro = G.mfros[i];
    const auto& next_fro = G.fros[i+1];
    const auto& now_vpos = G.vpos[i];
    const auto& now_ent = G.fent[i];
    const auto& now_lve = G.flve[i];
    size_t kk = now_fro.size();
    size_t tt = med_fro.size();
    size_t ll = next_fro.size();
    secs[i+1] = snum;
    
    for(const auto& ent : maps[i]){
      size_t now_id = ent.second;
      const State& now_state = ent.first;
      
      // generate intermediate state
      State med_state;
      med_state.resize(tt);
      int8_t cc = dp[now_id].cnum;
      int8_t cc_old = cc;
      
      std::copy(now_state.begin(), now_state.end(), med_state.begin());
      for(const auto& pos : now_ent){
        med_state[pos] = issrc[med_fro[pos]] ? 0 : cc++;
      }
      
      // lo_state processing
      {
        // generate lo_state
        State lo_state;
        lo_state.resize(ll);
        std::copy(med_state.begin(), med_state.begin() + ll, lo_state.begin());
        for(const auto& pos : now_lve){
          if(pos < ll) lo_state[pos] = -1;
        }
        std::vector<int8_t> renum(cc, -1);
        renum[0] = 0;
        int8_t cc_new = 1;
        for(auto&& val : lo_state){
          if(val < 0) continue;
          if(renum[val] < 0) renum[val] = cc_new++;
          val = renum[val];
        }
        
        // find or generate id
        size_t lo_id;
        auto it = maps[i+1].find(lo_state);
        if(it != maps[i+1].end()){ // maps[i+1] has already had entry
          lo_id = it->second;
        }else{                     // there is no entry
          maps[i+1].emplace(lo_state, snum);
          lo_id = snum++;
          dp.emplace_back(lo_state, i+1, cc_new);
        }
        dp[now_id].lo = lo_id;
        
        // vlo equals renum
        std::copy(renum.begin(), renum.begin() + cc_old, dp[now_id].vlo.begin());
      }
      
      // hi_state processing
      int8_t cat_to   = med_state[now_vpos.first];
      int8_t cat_from = med_state[now_vpos.second];
      {
        // generate hi_state
        State hi_state;
        hi_state.resize(ll);
        std::copy(med_state.begin(), med_state.begin() + ll, hi_state.begin());
        for(const auto& pos : now_lve){
          if(pos < ll) hi_state[pos] = -1;
        }
        std::vector<int8_t> renum(cc, -1);
        renum[0] = 0;
        if(cat_to == 0)   renum[cat_from] = 0;
        if(cat_from == 0) renum[cat_to]   = 0;
        int8_t cc_new = 1;
        for(auto&& val : hi_state){
          if(val < 0) continue;
          if(renum[val] < 0){
            renum[val] = cc_new++;
            if(val == cat_to)        renum[cat_from] = renum[val];
            else if(val == cat_from) renum[cat_to]   = renum[val];
          }
          val = renum[val];
        }
        
        // find or generate id
        size_t hi_id;
        auto it = maps[i+1].find(hi_state);
        if(it != maps[i+1].end()){ // maps[i+1] has already had entry
          hi_id = it->second;
        }else{                     // there is no entry
          maps[i+1].emplace(hi_state, snum);
          hi_id = snum++;
          dp.emplace_back(hi_state, i+1, cc_new);
        }
        dp[now_id].hi = hi_id;
        
        // vhi equals renum
        std::copy(renum.begin(), renum.begin() + cc_old, dp[now_id].vhi.begin());
      }
    }
    //maps[i].clear();
  }
  //maps[m].clear();
  secs[m+1] = snum;
  dp[snum-1].q[0][0] = 1.0;
  
  // dp calculation
  // degree 1 check
  std::vector<std::vector<int>> deg1(m);
  for(int i=0; i<m; ++i){
    deg1[i].resize(G.mfros[i].size());
    for(const auto& pos : G.flve[i]){
      for(const auto& pos2 : G.fent[i]){
        if(pos == pos2){
          deg1[i][pos] = 1;
        }
      }
    }
  }
  
  // top-down calc. of p
  for(size_t ind=0; ind<snum-1; ++ind){
    DPBlock<xReal>& dp_now = dp[ind];
    dp[dp_now.lo].p += (1.0 - pi[dp_now.level]) * dp_now.p;
    dp[dp_now.hi].p += pi[dp_now.level]         * dp_now.p;
  }
  
  // top-down calc. of c
  for(size_t ind=0; ind<snum-1; ++ind){
    DPBlock<xReal>& dp_now = dp[ind];
    for(int8_t b=0; b<dp_now.cnum; ++b){
      if(dp_now.vlo[b] >= 0) dp[dp_now.lo].c[dp_now.vlo[b]] += (1.0 - pi[dp_now.level]) * dp_now.c[b];
      if(dp_now.vhi[b] >= 0) dp[dp_now.hi].c[dp_now.vhi[b]] += pi[dp_now.level]         * dp_now.c[b];
    }
    for(const auto& pos : G.flve[dp_now.level]){
      int8_t b = deg1[dp_now.level][pos] ? 0 : dp_now.s[pos];
      int vid = G.mfros[dp_now.level][pos];
      if(dp_now.vlo[b] >= 0) dp[dp_now.lo].c[dp_now.vlo[b]] += wgt[vid] * (1.0 - pi[dp_now.level]) * dp_now.p;
      if(dp_now.vhi[b] >= 0) dp[dp_now.hi].c[dp_now.vhi[b]] += wgt[vid] * pi[dp_now.level]         * dp_now.p;
    }
  }
  
  // bottom-up calc. of q
  for(size_t ind=snum-2; ind>=1; --ind){
    DPBlock<xReal>& dp_now = dp[ind];
    for(int8_t a=0; a<dp_now.cnum; ++a){
      for(int8_t b=0; b<dp_now.cnum; ++b){
        if(a == b)      dp_now.q[a][b] = 1.0;
        else if(b == 0) dp_now.q[a][b] = 1.0;
        else{
          if(dp_now.vlo[b] >= 0){
            if(dp_now.vlo[a] >= 0) dp_now.q[a][b]  = (1.0 - pi[dp_now.level]) * dp[dp_now.lo].q[dp_now.vlo[a]][dp_now.vlo[b]];
            else                   dp_now.q[a][b]  = (1.0 - pi[dp_now.level]) * dp[dp_now.lo].q[0][dp_now.vlo[b]];
          }
          if(dp_now.vhi[b] >= 0){
            if(dp_now.vhi[a] >= 0) dp_now.q[a][b] += pi[dp_now.level]         * dp[dp_now.hi].q[dp_now.vhi[a]][dp_now.vhi[b]];
            else                   dp_now.q[a][b] += pi[dp_now.level]         * dp[dp_now.hi].q[0][dp_now.vhi[b]];
          }
        }
      }
    }
    if(G.flve[dp_now.level].size() == 2){
      int8_t cat1 = dp_now.s[G.flve[dp_now.level][0]];
      int8_t cat2 = dp_now.s[G.flve[dp_now.level][1]];
      if(cat1 != cat2 && dp_now.vhi[cat1] < 0 && dp_now.vhi[cat2] < 0){
        dp_now.q[cat1][cat2] = dp_now.q[cat2][cat1] = pi[dp_now.level];
      }
    }
  }
  
  // bottom-up calc. of d
  for(size_t ind=snum-2; ind>=1; --ind){
    DPBlock<xReal>& dp_now = dp[ind];
    for(int8_t b=0; b<dp_now.cnum; ++b){
      dp_now.d[b]  = (1.0 - pi[dp_now.level]) * dp[dp_now.lo].d[std::max(dp_now.vlo[b], static_cast<int8_t>(0))];
      dp_now.d[b] += pi[dp_now.level]         * dp[dp_now.hi].d[std::max(dp_now.vhi[b], static_cast<int8_t>(0))];
    }
    for(const auto& pos : G.flve[dp_now.level]){
      int8_t a = deg1[dp_now.level][pos] ? 0 : dp_now.s[pos];
      int vid = G.mfros[dp_now.level][pos];
      for(int8_t b=0; b<dp_now.cnum; ++b){
        dp_now.d[b] += wgt[vid] * dp_now.q[b][a];
      }
    }
  }
  
  // levelwise computation
  res.resize(n+1);
  std::vector<int> computed(n+1);
  vtoipos.resize(n+1);
  computed[0] = 1;
  for(const auto& u : srcs) computed[u] = 1;
  for(size_t i=0; i<m; ++i){
    std::vector<size_t> compute_pos;
    const auto& now_fro = G.fros[i];
    for(size_t pos=0; pos<now_fro.size(); ++pos){
      int now_v = now_fro[pos];
      if(!computed[now_v]){
        compute_pos.emplace_back(pos);
        vtoipos[now_v] = std::make_pair(i, pos);
        computed[now_v] = 1;
        res[now_v] = 0.0;
      }
    }
    for(size_t k=secs[i]; k<secs[i+1]; ++k){
      const DPBlock<xReal>& dp_now = dp[k];
      for(size_t pos : compute_pos){
        int now_v = now_fro[pos];
        int8_t cat = dp_now.s[pos];
        for(int8_t b=0; b<dp_now.cnum; ++b){
          res[now_v] += dp_now.c[b] * dp_now.q[cat][b];
        }
        res[now_v] += dp_now.p * dp_now.d[cat];
      }
    }
  }
  res[0] = dp[snum-1].c[0];
  for(const auto& u : srcs) res[u] = res[0];
}

template <typename xReal>
void LWDP<xReal>::innerSolveEmpty(const Graph& G, const std::vector<xReal>& pi, const std::vector<xReal>& wgt, std::vector<xReal>& res){
  int n = G.numV();
  int m = G.numE();
  
  maps.resize(m+1);
  secs.resize(m+2);
  
  // id=0: root node
  {
    State root;
    maps[0].emplace(root, 0);
    dp.emplace_back(root, 0, 0);
    dp[0].p = 1.0;
  }
  secs[0] = 0;
  snum = 1;
  
  for(size_t i=0; i<m; ++i){
    const auto& now_fro = G.fros[i];
    const auto& med_fro = G.mfros[i];
    const auto& next_fro = G.fros[i+1];
    const auto& now_vpos = G.vpos[i];
    const auto& now_ent = G.fent[i];
    const auto& now_lve = G.flve[i];
    size_t kk = now_fro.size();
    size_t tt = med_fro.size();
    size_t ll = next_fro.size();
    secs[i+1] = snum;
    
    for(const auto& ent : maps[i]){
      size_t now_id = ent.second;
      const State& now_state = ent.first;
      
      // generate intermediate state
      State med_state;
      med_state.resize(tt);
      int8_t cc = dp[now_id].cnum;
      int8_t cc_old = cc;
      
      std::copy(now_state.begin(), now_state.end(), med_state.begin());
      for(const auto& pos : now_ent){
        med_state[pos] = cc++;
      }
      
      // lo_state processing
      {
        // generate lo_state
        State lo_state;
        lo_state.resize(ll);
        std::copy(med_state.begin(), med_state.begin() + ll, lo_state.begin());
        for(const auto& pos : now_lve){
          if(pos < ll) lo_state[pos] = -1;
        }
        std::vector<int8_t> renum(cc, -1);
        int8_t cc_new = 0;
        for(auto&& val : lo_state){
          if(val < 0) continue;
          if(renum[val] < 0) renum[val] = cc_new++;
          val = renum[val];
        }
        
        // find or generate id
        size_t lo_id;
        auto it = maps[i+1].find(lo_state);
        if(it != maps[i+1].end()){ // maps[i+1] has already had entry
          lo_id = it->second;
        }else{                     // there is no entry
          maps[i+1].emplace(lo_state, snum);
          lo_id = snum++;
          dp.emplace_back(lo_state, i+1, cc_new);
        }
        dp[now_id].lo = lo_id;
        
        // vlo equals renum
        std::copy(renum.begin(), renum.begin() + cc_old, dp[now_id].vlo.begin());
      }
      
      // hi_state processing
      int8_t cat_to   = med_state[now_vpos.first];
      int8_t cat_from = med_state[now_vpos.second];
      {
        // generate hi_state
        State hi_state;
        hi_state.resize(ll);
        std::copy(med_state.begin(), med_state.begin() + ll, hi_state.begin());
        for(const auto& pos : now_lve){
          if(pos < ll) hi_state[pos] = -1;
        }
        std::vector<int8_t> renum(cc, -1);
        int8_t cc_new = 0;
        for(auto&& val : hi_state){
          if(val < 0) continue;
          if(renum[val] < 0){
            renum[val] = cc_new++;
            if(val == cat_to)        renum[cat_from] = renum[val];
            else if(val == cat_from) renum[cat_to]   = renum[val];
          }
          val = renum[val];
        }
        
        // find or generate id
        size_t hi_id;
        auto it = maps[i+1].find(hi_state);
        if(it != maps[i+1].end()){ // maps[i+1] has already had entry
          hi_id = it->second;
        }else{                     // there is no entry
          maps[i+1].emplace(hi_state, snum);
          hi_id = snum++;
          dp.emplace_back(hi_state, i+1, cc_new);
        }
        dp[now_id].hi = hi_id;
        
        // vhi equals renum
        std::copy(renum.begin(), renum.begin() + cc_old, dp[now_id].vhi.begin());
      }
    }
    maps[i].clear();
  }
  maps[m].clear();
  secs[m+1] = snum;
  
  // dp calculation
  // top-down calc. of p
  for(size_t ind=0; ind<snum-1; ++ind){
    DPBlock<xReal>& dp_now = dp[ind];
    dp[dp_now.lo].p += (1.0 - pi[dp_now.level]) * dp_now.p;
    dp[dp_now.hi].p += pi[dp_now.level]         * dp_now.p;
  }
  
  // top-down calc. of c
  for(size_t ind=0; ind<snum-1; ++ind){
    DPBlock<xReal>& dp_now = dp[ind];
    for(int8_t b=0; b<dp_now.cnum; ++b){
      if(dp_now.vlo[b] >= 0) dp[dp_now.lo].c[dp_now.vlo[b]] += (1.0 - pi[dp_now.level]) * dp_now.c[b];
      if(dp_now.vhi[b] >= 0) dp[dp_now.hi].c[dp_now.vhi[b]] += pi[dp_now.level]         * dp_now.c[b];
    }
    for(const auto& pos : G.flve[dp_now.level]){
      int8_t b = dp_now.s[pos];
      int vid = G.mfros[dp_now.level][pos];
      if(dp_now.vlo[b] >= 0) dp[dp_now.lo].c[dp_now.vlo[b]] += wgt[vid] * (1.0 - pi[dp_now.level]) * dp_now.p;
      if(dp_now.vhi[b] >= 0) dp[dp_now.hi].c[dp_now.vhi[b]] += wgt[vid] * pi[dp_now.level]         * dp_now.p;
    }
  }
  
  // bottom-up calc. of q
  for(size_t ind=snum-2; ind>=1; --ind){
    DPBlock<xReal>& dp_now = dp[ind];
    for(int8_t a=0; a<dp_now.cnum; ++a){
      for(int8_t b=0; b<a; ++b){
        dp_now.q[a][b] = dp_now.q[b][a];
      }
      dp_now.q[a][a] = 1.0;
      for(int8_t b=a+1; b<dp_now.cnum; ++b){
        if(dp_now.vlo[a] >= 0 && dp_now.vlo[b] >= 0) dp_now.q[a][b]  = (1.0 - pi[dp_now.level]) * dp[dp_now.lo].q[dp_now.vlo[a]][dp_now.vlo[b]];
        if(dp_now.vhi[a] >= 0 && dp_now.vhi[b] >= 0) dp_now.q[a][b] += pi[dp_now.level]         * dp[dp_now.hi].q[dp_now.vhi[a]][dp_now.vhi[b]];
      }
    }
    if(G.flve[dp_now.level].size() == 2){
      int8_t cat1 = dp_now.s[G.flve[dp_now.level][0]];
      int8_t cat2 = dp_now.s[G.flve[dp_now.level][1]];
      if(cat1 != cat2 && dp_now.vhi[cat1] < 0 && dp_now.vhi[cat2] < 0){
        dp_now.q[cat1][cat2] = dp_now.q[cat2][cat1] = pi[dp_now.level];
      }
    }
  }
  
  // bottom-up calc. of d
  for(size_t ind=snum-2; ind>=1; --ind){
    DPBlock<xReal>& dp_now = dp[ind];
    for(int8_t b=0; b<dp_now.cnum; ++b){
      if(dp_now.vlo[b] >= 0) dp_now.d[b]  = (1.0 - pi[dp_now.level]) * dp[dp_now.lo].d[dp_now.vlo[b]];
      if(dp_now.vhi[b] >= 0) dp_now.d[b] += pi[dp_now.level]         * dp[dp_now.hi].d[dp_now.vhi[b]];
    }
    for(const auto& pos : G.flve[dp_now.level]){
      int8_t a = dp_now.s[pos];
      int vid = G.mfros[dp_now.level][pos];
      for(int8_t b=0; b<dp_now.cnum; ++b){
        dp_now.d[b] += wgt[vid] * dp_now.q[a][b];
      }
    }
  }
  
  // levelwise calculation
  res.resize(n+1);
  res[0] = 0.0;
  std::vector<int> computed(n+1);
  computed[0] = 1;
  for(size_t i=0; i<m; ++i){
    std::vector<size_t> compute_pos;
    const auto& now_fro = G.fros[i];
    for(size_t pos=0; pos<now_fro.size(); ++pos){
      int now_v = now_fro[pos];
      if(!computed[now_v]){
        compute_pos.emplace_back(pos);
        computed[now_v] = 1;
        res[now_v] = 0.0;
      }
    }
    for(size_t k=secs[i]; k<secs[i+1]; ++k){
      const DPBlock<xReal>& dp_now = dp[k];
      for(size_t pos : compute_pos){
        int now_v = now_fro[pos];
        int8_t cat = dp_now.s[pos];
        for(int8_t b=0; b<dp_now.cnum; ++b){
          res[now_v] += dp_now.c[b] * dp_now.q[b][cat];
        }
        res[now_v] += dp_now.p * dp_now.d[cat];
      }
    }
  }
}

#endif // REL_CRIT_LWDP_HPP
