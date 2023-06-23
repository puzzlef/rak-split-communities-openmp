#pragma once
#include <utility>
#include <vector>
#include "_main.hxx"
#include "update.hxx"

#ifdef OPENMP
#include <omp.h>
#endif

using std::pair;
using std::tuple;
using std::vector;
using std::make_pair;
using std::swap;
using std::move;
using std::get;




// RAK OPTIONS
// -----------

struct RakOptions {
  int    repeat;
  double tolerance;
  int    maxIterations;

  RakOptions(int repeat=1, double tolerance=0.05, int maxIterations=20) :
  repeat(repeat), tolerance(tolerance), maxIterations(maxIterations) {}
};

// Weight to be using in hashtable.
#define RAK_WEIGHT_TYPE double




// RAK RESULT
// ----------

template <class K>
struct RakResult {
  vector<K> membership;
  int   iterations;
  float time;
  float preprocessingTime;

  RakResult(vector<K>&& membership, int iterations=0, float time=0, float preprocessingTime=0) :
  membership(membership), iterations(iterations), time(time), preprocessingTime(preprocessingTime) {}

  RakResult(vector<K>& membership, int iterations=0, float time=0, float preprocessingTime=0) :
  membership(move(membership)), iterations(iterations), time(time), preprocessingTime(preprocessingTime) {}
};




// RAK HASHTABLES
// --------------

/**
 * Allocate a number of hashtables.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param S size of each hashtable
 */
template <class K, class W>
inline void rakAllocateHashtables(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, size_t S) {
  size_t N = vcs.size();
  for (size_t i=0; i<N; ++i) {
    vcs[i]   = new vector<K>();
    vcout[i] = new vector<W>(S);
  }
}


/**
 * Free a number of hashtables.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 */
template <class K, class W>
inline void rakFreeHashtables(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout) {
  size_t N = vcs.size();
  for (size_t i=0; i<N; ++i) {
    delete vcs[i];
    delete vcout[i];
  }
}




// RAK INITIALIZE
// --------------

/**
 * Initialize communities such that each vertex is its own community.
 * @param vcom community each vertex belongs to (updated)
 * @param x original graph
 */
template <class G, class K>
inline void rakInitialize(vector<K>& vcom, const G& x) {
  x.forEachVertexKey([&](auto u) { vcom[u] = u; });
}

#ifdef OPENMP
template <class G, class K>
inline void rakInitializeOmp(vector<K>& vcom, const G& x) {
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    vcom[u] = u;
  }
}
#endif


/**
 * Initialize communities from given initial communities.
 * @param vcom community each vertex belongs to (updated)
 * @param x original graph
 * @param q initial community each vertex belongs to
 */
template <class G, class K>
inline void rakInitializeFrom(vector<K>& vcom, const G& x, const vector<K>& q) {
  x.forEachVertexKey([&](auto u) { vcom[u] = q[u]; });
}

#ifdef OPENMP
template <class G, class K>
inline void rakInitializeFromOmp(vector<K>& vcom, const G& x, const vector<K>& q) {
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    vcom[u] = q[u];
  }
}
#endif




// RAK CHOOSE COMMUNITY
// --------------------

/**
 * Scan an edge community connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param u given vertex
 * @param v outgoing edge vertex
 * @param w outgoing edge weight
 * @param vcom community each vertex belongs to
 */
template <bool SELF=false, class K, class V, class W>
inline void rakScanCommunity(vector<K>& vcs, vector<W>& vcout, K u, K v, V w, const vector<K>& vcom) {
  if (!SELF && u==v) return;
  K c = vcom[v];
  if (!vcout[c]) vcs.push_back(c);
  vcout[c] += w;
}


/**
 * Scan communities connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 */
template <bool SELF=false, class G, class K, class W>
inline void rakScanCommunities(vector<K>& vcs, vector<W>& vcout, const G& x, K u, const vector<K>& vcom) {
  x.forEachEdge(u, [&](auto v, auto w) { rakScanCommunity<SELF>(vcs, vcout, u, v, w, vcom); });
}


/**
 * Clear communities scan data.
 * @param vcs total edge weight from vertex u to community C (updated)
 * @param vcout communities vertex u is linked to (updated)
 */
template <class K, class W>
inline void rakClearScan(vector<K>& vcs, vector<W>& vcout) {
  for (K c : vcs)
    vcout[c] = W();
  vcs.clear();
}


/**
 * Choose connected community with most weight.
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 * @param vcs communities vertex u is linked to
 * @param vcout total edge weight from vertex u to community C
 * @returns [best community, best edge weight to community]
 */
template <bool STRICT=false, class G, class K, class W>
inline pair<K, W> rakChooseCommunity(const G& x, K u, const vector<K>& vcom, const vector<K>& vcs, const vector<W>& vcout) {
  K d = vcom[u];
  K cmax = K();
  W wmax = W();
  for (K c : vcs) {
    // Do some basic randomization if multiple labels have max weight.
    if (vcout[c]>wmax || (!STRICT && vcout[c]==wmax && (c & 2))) { cmax = c; wmax = vcout[c]; }
  }
  return make_pair(cmax, wmax);
}




// RAK AFFECTED VERTICES DELTA-SCREENING
// -------------------------------------
// Using delta-screening approach.
// - All edge batches are undirected, and sorted by source vertex-id.
// - For edge additions across communities with source vertex `i` and highest modularity changing edge vertex `j*`,
//   `i`'s neighbors and `j*`'s community is marked as affected.
// - For edge deletions within the same community `i` and `j`,
//   `i`'s neighbors and `j`'s community is marked as affected.

/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param x original graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community each vertex belongs to
 * @returns flags for each vertex marking whether it is affected
 */
template <bool STRICT=false, class B, class G, class K, class V, class W>
inline void rakAffectedVerticesDeltaScreening(vector<K>& vcs, vector<W>& vcout, vector<B>& vertices, vector<B>& neighbors, vector<B>& communities, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  fillValueU(vertices,    B());
  fillValueU(neighbors,   B());
  fillValueU(communities, B());
  for (const auto& [u, v] : deletions) {
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
    neighbors[u] = 1;
    communities[vcom[v]] = 1;
  }
  for (size_t i=0; i<insertions.size();) {
    K u = get<0>(insertions[i]);
    rakClearScan(vcs, vcout);
    for (; i<insertions.size() && get<0>(insertions[i])==u; ++i) {
      K v = get<1>(insertions[i]);
      V w = get<2>(insertions[i]);
      if (vcom[u] == vcom[v]) continue;
      rakScanCommunity(vcs, vcout, u, v, w, vcom);
    }
    auto [c, w] = rakChooseCommunity<STRICT>(x, u, vcom, vcs, vcout);
    if (w<=0 || c==vcom[u]) continue;
    vertices[u]  = 1;
    neighbors[u] = 1;
    communities[c] = 1;
  }
  x.forEachVertexKey([&](auto u) {
    if (neighbors[u]) x.forEachEdgeKey(u, [&](auto v) { vertices[v] = 1; });
    if (communities[vcom[u]]) vertices[u] = 1;
  });
}


#ifdef OPENMP
template <bool STRICT=false, class B, class G, class K, class V, class W>
inline void rakAffectedVerticesDeltaScreeningOmp(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, vector<B>& vertices, vector<B>& neighbors, vector<B>& communities, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  size_t S = x.span();
  size_t D = deletions.size();
  size_t I = insertions.size();
  fillValueOmpU(vertices,    B());
  fillValueOmpU(neighbors,   B());
  fillValueOmpU(communities, B());
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<D; ++i) {
    K u = get<0>(deletions[i]);
    K v = get<1>(deletions[i]);
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
    neighbors[u] = 1;
    communities[vcom[v]] = 1;
  }
  #pragma omp parallel
  {
    int T = omp_get_num_threads();
    int t = omp_get_thread_num();
    K  u0 = I>0? get<0>(insertions[0]) : 0;
    for (size_t i=0, n=0; i<I;) {
      K u = get<0>(insertions[i]);
      if (u!=u0) { ++n; u0 = u; }
      if (n % T != t) { ++i; continue; }
      rakClearScan(*vcs[t], *vcout[t]);
      for (; i<I && get<0>(insertions[i])==u; ++i) {
        K v = get<1>(insertions[i]);
        V w = get<2>(insertions[i]);
        if (vcom[u] == vcom[v]) continue;
        rakScanCommunity(*vcs[t], *vcout[t], u, v, w, vcom);
      }
      auto [c, w] = rakChooseCommunity<STRICT>(x, u, vcom, *vcs[t], *vcout[t]);
      if (w<=0 || c==vcom[u]) continue;
      vertices[u]  = 1;
      neighbors[u] = 1;
      communities[c] = 1;
    }
  }
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    if (neighbors[u]) x.forEachEdgeKey(u, [&](auto v) { vertices[v] = 1; });
    if (communities[vcom[u]]) vertices[u] = 1;
  }
}
#endif




// RAK AFFECTED VERTICES FRONTIER
// ------------------------------
// Using frontier based approach.
// - All source and destination vertices are marked as affected for insertions and deletions.
// - For edge additions across communities with source vertex `i` and destination vertex `j`,
//   `i` is marked as affected.
// - For edge deletions within the same community `i` and `j`,
//   `i` is marked as affected.
// - Vertices whose communities change in local-moving phase have their neighbors marked as affected.

/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param x original graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community each vertex belongs to
 * @returns flags for each vertex marking whether it is affected
 */
template <class B, class G, class K, class V>
inline void rakAffectedVerticesFrontier(vector<B>& vertices, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  fillValueU(vertices, B());
  for (const auto& [u, v] : deletions) {
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
  }
  for (const auto& [u, v, w] : insertions) {
    if (vcom[u] == vcom[v]) continue;
    vertices[u]  = 1;
  }
}


#ifdef OPENMP
template <class B, class G, class K, class V>
inline void rakAffectedVerticesFrontierOmp(vector<B>& vertices, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  fillValueOmpU(vertices, B());
  size_t D = deletions.size();
  size_t I = insertions.size();
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<D; ++i) {
    K u = get<0>(deletions[i]);
    K v = get<1>(deletions[i]);
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
  }
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<I; ++i) {
    K u = get<0>(insertions[i]);
    K v = get<1>(insertions[i]);
    if (vcom[u] == vcom[v]) continue;
    vertices[u]  = 1;
  }
}
#endif




// RAK MOVE ITERATION
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
template <bool STRICT=false, class G, class K, class W, class FA, class FP>
inline size_t rakMoveIteration(vector<K>& vcs, vector<W>& vcout, vector<K>& vcom, const G& x, FA fa, FP fp) {
  size_t a = 0;
  x.forEachVertexKey([&](auto u) {
    if (!fa(u)) return;
    K d = vcom[u];
    rakClearScan(vcs, vcout);
    rakScanCommunities(vcs, vcout, x, u, vcom);
    auto [c, w] = rakChooseCommunity<STRICT>(x, u, vcom, vcs, vcout);
    if (c && c!=d) { vcom[u] = c; ++a; fp(u); }
  });
  return a;
}

#ifdef OPENMP
template <bool STRICT=false, class G, class K, class W, class FA, class FP>
inline size_t rakMoveIterationOmp(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, vector<K>& vcom, const G& x, FA fa, FP fp) {
  size_t a = K();
  size_t S = x.span();
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
#endif




// RAK
// ---

template <bool STRICT=false, class G, class K, class FM, class FA, class FP>
RakResult<K> rakSeq(const G& x, const vector<K>* q, const RakOptions& o, FM fm, FA fa, FP fp) {
  using V = typename G::edge_value_type;
  using W = RAK_WEIGHT_TYPE;
  int l = 0;
  size_t S = x.span();
  size_t N = x.order();
  vector<K> vcom(S), vcs;
  vector<W> vcout(S);
  float tm = 0;
  float t  = measureDuration([&]() {
    tm += measureDuration(fm);
    if (q) rakInitializeFrom(vcom, x, *q);
    else   rakInitialize(vcom, x);
    for (l=0; l<o.maxIterations;) {
      size_t n = rakMoveIteration<STRICT>(vcs, vcout, vcom, x, fa, fp); ++l;
      if (double(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  return {vcom, l, t, tm/o.repeat};
}


#ifdef OPENMP
template <bool STRICT=false, class G, class K, class FM, class FA, class FP>
RakResult<K> rakOmp(const G& x, const vector<K>* q, const RakOptions& o, FM fm, FA fa, FP fp) {
  using V = typename G::edge_value_type;
  using W = RAK_WEIGHT_TYPE;
  int l = 0;
  int T = omp_get_max_threads();
  size_t S = x.span();
  size_t N = x.order();
  vector<K> vcom(S);
  vector<vector<K>*> vcs(T);
  vector<vector<W>*> vcout(T);
  rakAllocateHashtables(vcs, vcout, S);
  float tm = 0;
  float t  = measureDuration([&]() {
    tm += measureDuration(fm);
    if (q) rakInitializeFromOmp(vcom, x, *q);
    else   rakInitializeOmp(vcom, x);
    for (l=0; l<o.maxIterations;) {
      size_t n = rakMoveIterationOmp<STRICT>(vcs, vcout, vcom, x, fa, fp); ++l;
      if (double(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  rakFreeHashtables(vcs, vcout);
  return {vcom, l, t, tm/o.repeat};
}
#endif




// RAK STATIC
// ----------

template <bool STRICT=false, class G, class K>
inline RakResult<K> rakStaticSeq(const G& x, const vector<K>* q=nullptr, const RakOptions& o={}) {
  auto fm = []() {};
  auto fa = [](auto u) { return u; };
  auto fp = [](auto u) {};
  return rakSeq<STRICT>(x, q, o, fm, fa, fp);
}

#ifdef OPENMP
template <bool STRICT=false, class G, class K>
inline RakResult<K> rakStaticOmp(const G& x, const vector<K>* q=nullptr, const RakOptions& o={}) {
  auto fm = []() {};
  auto fa = [](auto u) { return u; };
  auto fp = [](auto u) {};
  return rakOmp<STRICT>(x, q, o, fm, fa, fp);
}
#endif




// RAK DYNAMIC DELTA-SCREENING
// ---------------------------

template <bool STRICT=false, class FLAG=char, class G, class K, class V>
inline RakResult<K> rakDynamicDeltaScreeningSeq(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const RakOptions& o={}) {
  using  W = RAK_WEIGHT_TYPE;
  using  B = FLAG;
  size_t S = x.span();
  const vector<K>& vcom = *q;
  vector<K> vcs; vector<W> vcout(S);
  vector<B> vertices(S), neighbors(S), communities(S);
  auto fm = [&]() { rakAffectedVerticesDeltaScreening<STRICT>(vcs, vcout, vertices, neighbors, communities, x, deletions, insertions, vcom); vcs.clear(); vcout.clear(); };
  auto fa = [&](auto u) { return vertices[u]==B(1); };
  auto fp = [](auto u) {};
  return rakSeq<STRICT>(x, q, o, fm, fa, fp);
}


#ifdef OPENMP
template <bool STRICT=false, class FLAG=char, class G, class K, class V>
inline RakResult<K> rakDynamicDeltaScreeningOmp(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const RakOptions& o={}) {
  using  W = RAK_WEIGHT_TYPE;
  using  B = FLAG;
  size_t S = x.span();
  int    T = omp_get_max_threads();
  const vector<K>& vcom = *q;
  vector<B> vertices(S), neighbors(S), communities(S);
  vector<vector<K>*> vcs(T);
  vector<vector<W>*> vcout(T);
  rakAllocateHashtables(vcs, vcout, S);
  auto fm = [&]() { rakAffectedVerticesDeltaScreeningOmp<STRICT>(vcs, vcout, vertices, neighbors, communities, x, deletions, insertions, vcom); rakFreeHashtables(vcs, vcout); };
  auto fa = [&](auto u) { return vertices[u]==B(1); };
  auto fp = [](auto u) {};
  return rakOmp<STRICT>(x, q, o, fm, fa, fp);
}
#endif




// RAK DYNAMIC FRONTIER
// --------------------

template <bool STRICT=false, class FLAG=char, class G, class K, class V>
inline RakResult<K> rakDynamicFrontierSeq(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const RakOptions& o={}) {
  using  B = FLAG;
  size_t S = x.span();
  const vector<K>& vcom = *q;
  vector<B> vertices(S);
  auto fm = [&]() { rakAffectedVerticesFrontier(vertices, x, deletions, insertions, vcom); };
  auto fa = [&](auto u) { return vertices[u]==B(1); };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vertices[v] = 1; }); };
  return rakSeq<STRICT>(x, q, o, fm, fa, fp);
}


#ifdef OPENMP
template <bool STRICT=false, class FLAG=char, class G, class K, class V>
inline RakResult<K> rakDynamicFrontierOmp(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const RakOptions& o={}) {
  using  B = FLAG;
  size_t S = x.span();
  const vector<K>& vcom = *q;
  vector<B> vertices(S);
  auto fm = [&]() { rakAffectedVerticesFrontier(vertices, x, deletions, insertions, vcom); };
  auto fa = [&](auto u) { return vertices[u]==B(1); };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vertices[v] = 1; }); };
  return rakOmp<STRICT>(x, q, o, fm, fa, fp);
}
#endif
