#include "common.hpp"
#include "graph.hpp"
#include "graphsimplify.hpp"
#include "csrel.hpp"

#include <vector>
#include <utility>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>

void CSREL::solve(const Graph& G, const std::vector<double>& pi, const std::unordered_set<int>& srcs, std::vector<double>& res){
  allClear();
  
  int n = G.numV();
  int m = G.numE();
  
  maps.resize(m+1);
  secs.resize(m+2);
  
  size_t src_final = 0;
  for(size_t i=0; i<m; ++i){
    for(const auto& pos : G.fent[i]){
      if(srcs.count(G.mfros[i][pos])) src_final = i;
    }
  }
  
  // id=0: dummy node (sink node)
  dp.emplace_back(StateCS(), m, 2);
  dp[0].q[0] = 0.0;
  dp[0].q[1] = 1.0;
  
  // id=1: root node
  {
    StateCS root;
    maps[0].emplace(root, 1);
    dp.emplace_back(root, 0, 0);
    dp[1].p = 1.0;
  }
  secs[0] = 1;
  
  snum = 2;
  
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
    maps[i+1].reserve(maps[i].size() * 2);
    
    // build ast_ent and ast_lve
    std::vector<size_t> ast_ent, ast_lve;
    for(const auto& pos : now_ent){
      if(srcs.count(med_fro[pos])) ast_ent.emplace_back(pos);
    }
    for(const auto& pos : now_lve){
      if(srcs.count(med_fro[pos])) ast_ent.emplace_back(pos);
    }
    
    for(const auto& ent : maps[i]){
      size_t now_id = ent.second;
      const StateCS& now_state = ent.first;
      const auto& now_comp = now_state.second;
      
      // generate intermediate state
      StateCS med_state;
      auto& med_comp = med_state.second;
      med_comp.resize(tt);
      int8_t cc = dp[now_id].cnum;
      int8_t cc_old = cc;
      
      med_state.first = now_state.first;
      std::copy(now_comp.begin(), now_comp.end(), med_comp.begin());
      for(const auto& pos : now_ent){
        med_comp[pos] = cc++;
      }
      // introduce asterisk when src enters the (intermediate-)frontier
      for(const auto& pos : ast_ent){
        med_state.first |= 1ULL << med_comp[pos];
      }
      
      // lo_state processing
      {
        // generate lo_state
        StateCS lo_state;
        auto& lo_comp = lo_state.second;
        lo_comp.resize(ll);
        std::copy(med_comp.begin(), med_comp.begin() + ll, lo_comp.begin());
        for(const auto& pos : now_lve){
          if(pos < ll) lo_comp[pos] = -1;
        }
        std::vector<int8_t> renum(cc, -1);
        int8_t cc_new = 0;
        for(auto&& val : lo_comp){
          if(val < 0) continue;
          if(renum[val] < 0) renum[val] = cc_new++;
          val = renum[val];
        }
        // asterisk renumbering
        bool prune = false;
        lo_state.first = 0ULL;
        {
          uint64_t astset = med_state.first;
          while(astset){
            uint64_t ast_c2ton = astset & (-astset);
            astset ^= ast_c2ton;
            int8_t ast_c = renum[log2ton(ast_c2ton)];
            if(ast_c < 0){  // asterisk component leaves without being connected to other components
              prune = true; // if so, prune this state
              break;
            }
            lo_state.first |= 1ULL << ast_c;
          }
        }
        
        // find or generate id
        size_t lo_id = 0;
        if(prune){ // prune lo_state if prune = true
          dp[now_id].lo = lo_id;
          if(i >= src_final && is2ton(med_state.first)) dp[now_id].vlo[log2ton(med_state.first)] = 1;
        }else{     // not pruned
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
      }
      
      // hi_state processing
      int8_t cat_to   = med_comp[now_vpos.first];
      int8_t cat_from = med_comp[now_vpos.second];
      {
        // generate hi_state
        StateCS hi_state;
        auto& hi_comp = hi_state.second;
        hi_comp.resize(ll);
        std::copy(med_comp.begin(), med_comp.begin() + ll, hi_comp.begin());
        for(const auto& pos : now_lve){
          if(pos < ll) hi_comp[pos] = -1;
        }
        std::vector<int8_t> renum(cc, -1);
        int8_t cc_new = 0;
        for(auto&& val : hi_comp){
          if(val < 0) continue;
          if(renum[val] < 0){
            renum[val] = cc_new++;
            if(val == cat_to)        renum[cat_from] = renum[val];
            else if(val == cat_from) renum[cat_to]   = renum[val];
          }
          val = renum[val];
        }
        
        // asterisk renumbering
        bool prune = false;
        hi_state.first = 0ULL;
        {
          uint64_t astset = med_state.first;
          while(astset){
            uint64_t ast_c2ton = astset & (-astset);
            astset ^= ast_c2ton;
            int8_t ast_c = renum[log2ton(ast_c2ton)];
            if(ast_c < 0){  // asterisk component leaves without being connected to other components
              prune = true; // if so, prune this state
              break;
            }
            hi_state.first |= 1ULL << ast_c;
          }
        }
        
        // find or generate id
        size_t hi_id = 0;
        if(prune){ // prune hi_state if prune = true
          dp[now_id].hi = hi_id;
          if(i >= src_final && is2ton(med_state.first)) dp[now_id].vhi[log2ton(med_state.first)] = 1;
          // process corner case
          if(i >= src_final && now_lve.size() == 2 && cat_from != cat_to){
            if     (med_state.first == 1ULL << cat_from) dp[now_id].vhi[cat_to]   = 1;
            else if(med_state.first == 1ULL << cat_to)   dp[now_id].vhi[cat_from] = 1;
            else if(med_state.first == ((1ULL << cat_from) | (1ULL << cat_to))){
              dp[now_id].vhi[cat_to]   = 1;
              dp[now_id].vhi[cat_from] = 1;
            }
          }
        }else{     // not pruned
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
    }
    maps[i].clear();
  }
  maps[m].clear();
  secs[m+1] = snum;
  
  // dp calculation
  // top-down calc. of p
  for(size_t ind=1; ind<snum; ++ind){
    DPBlockCS& dp_now = dp[ind];
    dp[dp_now.lo].p += (1.0 - pi[dp_now.level]) * dp_now.p;
    dp[dp_now.hi].p += pi[dp_now.level]         * dp_now.p;
  }
  
  // bottom-up calc. of q
  for(size_t ind=snum-1; ind>=2; --ind){
    DPBlockCS& dp_now = dp[ind];
    for(size_t c=0; c<dp_now.cnum; ++c){
      dp_now.q[c]  = dp_now.vlo[c] >= 0 ? (1.0 - pi[dp_now.level]) * dp[dp_now.lo].q[dp_now.vlo[c]] : 0.0;
      dp_now.q[c] += dp_now.vhi[c] >= 0 ? pi[dp_now.level]         * dp[dp_now.hi].q[dp_now.vhi[c]] : 0.0;
    }
  }
  
  // levelwise calculation
  res.resize(n+1);
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
      const DPBlockCS& dp_now = dp[k];
      for(size_t pos : compute_pos){
        int now_v = now_fro[pos];
        int8_t cat = dp_now.s.second[pos];
        res[now_v] += dp_now.p * dp_now.q[cat];
      }
    }
  }
}
