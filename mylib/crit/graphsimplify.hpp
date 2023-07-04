#ifndef REL_CRIT_GRAPHSIMPLIFY_HPP
#define REL_CRIT_GRAPHSIMPLIFY_HPP

#include "common.hpp"
#include "graph.hpp"

#include <vector>
#include <utility>
#include <cassert>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

template <typename xReal>
using History = std::tuple<int, int, xReal>;

template <typename xReal>
class GraphSimplify{
public:
  Graph Gnew;
  std::unordered_set<int> srcsnew;
  int tgtnew;
  std::vector<xReal> pinew;
  
  std::vector<int> oldtonew;
  std::vector<History<xReal>> hists;
  
  GraphSimplify(){};
  
  void HHSimplify(const Graph& G, const std::vector<xReal>& prob, const std::unordered_set<int>& srcs, int tgt);
  void LWSimplify(const Graph& G, const std::vector<xReal>& prob, const std::unordered_set<int>& srcs);
};

template <typename xReal>
void GraphSimplify<xReal>::HHSimplify(const Graph& G, const std::vector<xReal>& prob, const std::unordered_set<int>& srcs, int tgt){
  int mold = G.m;
  int nold = G.n;
  std::vector<std::unordered_map<int, int>> adjs(nold+1);
  std::vector<xReal> tmpi(prob);
  bool updated = false;
  auto srcsnew0 = srcs;
  tgtnew = tgt;
  oldtonew.resize(nold+1);
  
  // build adjacency list
  for(int i=0; i<mold; ++i){
    adjs[G.e[i].first].emplace(G.e[i].second, i);
    adjs[G.e[i].second].emplace(G.e[i].first, i);
  }
  do{
    updated = false;
    for(int k=1; k<=nold; ++k){
      if(adjs[k].size() == 1){
        const auto kitr = adjs[k].begin();
        int kk = kitr->first;
        const auto kkitr = adjs[kk].find(k);
        const auto srcsitr = srcsnew0.find(k);
        if(srcsitr == srcsnew0.cend() && tgtnew != k){
          adjs[k].erase(kitr);
          adjs[kk].erase(kkitr);
          updated = true;
        }else{
          bool innerupdated = false;
          if(tgtnew == k){
            tgtnew = kk;
            hists.emplace_back(k, kk, tmpi[kitr->second]);
            innerupdated = true;
          }
          if(srcsitr != srcsnew0.cend() && srcsnew0.size() == 1){
            srcsnew0.erase(srcsitr);
            srcsnew0.emplace(kk);
            hists.emplace_back(k, kk, tmpi[kitr->second]);
            innerupdated = true;
          }
          if(innerupdated){
            adjs[k].erase(kitr);
            adjs[kk].erase(kkitr);
          }
          updated = updated || innerupdated;
        }
      }
      else if(adjs[k].size() == 2 && srcsnew0.count(k) == 0 && tgtnew != k){
        const auto kitr1 = adjs[k].begin();
        const auto kitr2 = std::next(kitr1);
        int newid = std::min(kitr1->second, kitr2->second);
        xReal newprob = tmpi[kitr1->second] * tmpi[kitr2->second];
        int kk1 = kitr1->first;
        const auto kk1itr = adjs[kk1].find(k);
        int kk2 = kitr2->first;
        const auto kk2itr = adjs[kk2].find(k);
        adjs[k].erase(kitr1);
        adjs[k].erase(kitr2);
        adjs[kk1].erase(kk1itr);
        adjs[kk2].erase(kk2itr);
        auto kk1itrnew = adjs[kk1].find(kk2);
        auto kk2itrnew = adjs[kk2].find(kk1);
        if(kk1itrnew != adjs[kk1].end()){
          newid = std::min(newid, kk1itrnew->second);
          newprob = tmpi[kk1itrnew->second] + newprob - tmpi[kk1itrnew->second] * newprob;
          kk1itrnew->second = newid;
          kk2itrnew->second = newid;
        }else{
          adjs[kk1].emplace(kk2, newid);
          adjs[kk2].emplace(kk1, newid);
        }
        tmpi[newid] = newprob;
        updated = true;
      }
      if(srcsnew0.count(tgtnew) == 1){
        updated = false;
        break;
      }
    }
  }while(updated);
  
  int nnew = 0;
  for(int k=1; k<=nold; ++k){
    if(adjs[k].size() >= 1){
      oldtonew[k] = ++nnew;
    }
  }
  std::vector<Edge> newe(mold);
  for(int k=1; k<=nold; ++k){
    for(const auto& elem : adjs[k]){
      if(k < elem.first) newe[elem.second] = std::make_pair(oldtonew[k], oldtonew[elem.first]);
    }
  }
  for(int i=0; i<mold; ++i){
    if(newe[i].first != 0 && newe[i].second != 0){
      Gnew.addEdge(newe[i].first, newe[i].second);
      pinew.emplace_back(tmpi[i]);
    }
  }
  for(const auto& ele : srcsnew0){
    srcsnew.emplace(oldtonew[ele]);
  }
  tgtnew = oldtonew[tgtnew];
}

template <typename xReal>
void GraphSimplify<xReal>::LWSimplify(const Graph& G, const std::vector<xReal>& prob, const std::unordered_set<int>& srcs){
  int mold = G.m;
  int nold = G.n;
  std::vector<std::unordered_map<int, int>> adjs(nold+1);
  bool updated = false;
  oldtonew.resize(nold+1);
  
  // build adjacency list
  for(int i=0; i<mold; ++i){
    adjs[G.e[i].first].emplace(G.e[i].second, i);
    adjs[G.e[i].second].emplace(G.e[i].first, i);
  }
  do{
    updated = false;
    for(int k=1; k<=nold; ++k){
      if(adjs[k].size() == 1){
        const auto kitr = adjs[k].begin();
        int kk = kitr->first;
        const auto kkitr = adjs[kk].find(k);
        const auto srcsitr = srcs.find(k);
        if(srcsitr == srcs.cend()){
          adjs[k].erase(kitr);
          adjs[kk].erase(kkitr);
          hists.emplace_back(k, kk, prob[kitr->second]);
          updated = true;
        }
      }
    }
  }while(updated);
  
  int nnew = 0;
  for(int k=1; k<=nold; ++k){
    if(adjs[k].size() >= 1){
      oldtonew[k] = ++nnew;
    }
  }
  std::vector<Edge> newe(mold);
  for(int k=1; k<=nold; ++k){
    for(const auto& elem : adjs[k]){
      if(k < elem.first) newe[elem.second] = std::make_pair(oldtonew[k], oldtonew[elem.first]);
    }
  }
  for(int i=0; i<mold; ++i){
    if(newe[i].first != 0 && newe[i].second != 0){
      Gnew.addEdge(newe[i].first, newe[i].second);
      pinew.emplace_back(prob[i]);
    }
  }
  for(const auto& ele : srcs){
    srcsnew.emplace(oldtonew[ele]);
  }
}

#endif // REL_CRIT_GRAPHSIMPLIFY_HPP
