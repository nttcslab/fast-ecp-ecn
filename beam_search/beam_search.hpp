#ifndef BSEARCH_BEAM_SEARCH_HPP
#define BSEARCH_BEAM_SEARCH_HPP

#include <vector>
#include <utility>

typedef std::pair<int,int> Bedge;

std::vector<Bedge> ordering(int n, const std::vector<Bedge> &g);

#endif // BSEARCH_BEAM_SEARCH_HPP