#pragma once
#include <utility>
#include <memory>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "_main.hxx"
#include "vertices.hxx"
#include "edges.hxx"
#include "csr.hxx"
#include "rak.hxx"

using std::unique_ptr;
using std::tuple;
using std::vector;
using std::swap;




// RAK-MOVE-ITERATION
// ------------------

/**
 * Move each vertex to its best community.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param vcom community each vertex belongs to (updated)
 * @param x original graph
 * @param vdom community each vertex belonged to
 * @returns number of changed vertices
 */
template <bool STRICT=false, class G, class K, class V, class FA, class FP>
K rakMoveIterationOmp(vector<unique_ptr<vector<K>>>& vcs, vector<unique_ptr<vector<V>>>& vcout, vector<K>& vcom, const G& x, FA fa, FP fp) {
  K a = K();
  K S = x.span();
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (K u=0; u<S; ++u) {
    int t = omp_get_thread_num();
    if (!x.hasVertex(u)) continue;
    if (!fa(u)) continue;
    K d = vcom[u];
    rakClearScan(*vcs[t], *vcout[t]);
    rakScanCommunities(*vcs[t], *vcout[t], x, u, vcom);
    auto [c, w] = rakChooseCommunity<STRICT>(x, u, vcom, *vcs[t], *vcout[t]);
    if (c && c!=d) { vcom[u] = c; ++a; fp(u); }
  }
  return a;
}




// RAK-OMP
// -------

template <bool STRICT=false, class G, class K, class FA, class FP>
RakResult<K> rakOmp(const G& x, const vector<K>* q, const RakOptions& o, FA fa, FP fp) {
  using V = typename G::edge_value_type;
  int l = 0;
  int T = omp_get_num_threads();
  K S = x.span();
  K N = x.order();
  vector<K> vcom(S);
  vector<unique_ptr<vector<K>>> vcs(T);
  vector<unique_ptr<vector<V>>> vcout(T);
  for (int t=0; t<T; ++t) {
    vcs[t]   = new vector<K>();
    vcout[t] = new vector<V>(S);
  }
  float t = measureDuration([&]() {
    rakInitialize(vcom, x);
    for (l=0; l<o.maxIterations;) {
      K n = rakMoveIterationOmp<STRICT>(vcs, vcout, vcom, x, fa, fp); ++l;
      PRINTFD("rakOmp(): l=%d, n=%d, N=%d, n/N=%f\n", l, n, N, float(n)/N);
      if (float(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  return {vcom, l, t};
}
template <bool STRICT=false, class G, class K, class FA>
inline RakResult<K> rakOmp(const G& x, const vector<K>* q, const RakOptions& o, FA fa) {
  auto fp = [](auto u) {};
  return rakOmp<STRICT>(x, q, o, fa, fp);
}
template <bool STRICT=false, class G, class K>
inline RakResult<K> rakOmp(const G& x, const vector<K>* q, const RakOptions& o) {
  auto fa = [](auto u) { return true; };
  return rakOmp<STRICT>(x, q, o, fa);
}




// RAK-OMP-STATIC
// --------------

template <bool STRICT=false, class G, class K>
inline RakResult<K> rakOmpStatic(const G& x, const vector<K>* q=nullptr, const RakOptions& o={}) {
  return rakOmp<STRICT>(x, q, o);
}




// RAK-OMP-DYNAMIC-DELTA-SCREENING
// -------------------------------

template <bool STRICT=false, class G, class K, class V>
inline RakResult<K> rakOmpDynamicDeltaScreening(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const RakOptions& o={}) {
  K S = x.span();
  const vector<K>& vcom = *q;
  auto vaff = rakAffectedVerticesDeltaScreening<STRICT>(x, deletions, insertions, vcom);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return rakOmp<STRICT>(x, q, o, fa);
}




// RAK-OMP-DYNAMIC-FRONTIER
// ------------------------

template <bool STRICT=false, class G, class K, class V>
inline RakResult<K> rakOmpDynamicFrontier(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const RakOptions& o={}) {
  K S = x.span();
  const vector<K>& vcom = *q;
  auto vaff = rakAffectedVerticesFrontier(x, deletions, insertions, vcom);
  auto fa = [&](auto u) { return vaff[u]==true; };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vaff[v] = true; }); };
  return rakOmp<STRICT>(x, q, o, fa, fp);
}
